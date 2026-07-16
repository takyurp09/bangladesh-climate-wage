## 12_asymmetry_test.R
## Test asymmetric yield pass-through: attached (3-meal) vs casual (0-meal)
## across positive vs negative first-stage yield shocks
## Output: output/stage2/asymmetry/

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(fixest)
  library(ggplot2)
  library(kableExtra)
  library(patchwork)
})

ROOT <- here::here()
OUT  <- file.path(ROOT, "output/stage2/asymmetry")
for (d in c("", "tables", "figures", "models", "summary"))
  dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)

cat("=== YIELD PASS-THROUGH ASYMMETRY TEST ===\n")
cat("Output directory:", OUT, "\n\n")

## ── Helpers ───────────────────────────────────────────────────────────────── ##
fmt_p <- function(p) {
  vapply(p, function(x) {
    if (is.na(x)) return("NA")
    if (x < 0.001) return("<0.001")
    sprintf("%.3f", x)
  }, character(1))
}

sig_flag <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 0.01, "significant ***",
                ifelse(p < 0.05, "significant **",
                       ifelse(p < 0.10, "significant *",
                              "imprecise (not significant)"))))
}

extract_coef <- function(mod, param = "diff_log_yield_hat") {
  ct <- summary(mod)$coeftable
  if (!param %in% rownames(ct)) return(NULL)
  data.frame(
    coef  = ct[param, "Estimate"],
    se    = ct[param, "Std. Error"],
    p     = ct[param, "Pr(>|t|)"],
    ci_lo = ct[param, "Estimate"] - 1.96 * ct[param, "Std. Error"],
    ci_hi = ct[param, "Estimate"] + 1.96 * ct[param, "Std. Error"],
    nobs  = nobs(mod),
    stringsAsFactors = FALSE
  )
}

run_split_reg <- function(df, meal, shock_direction, sample_label, fe_rhs) {
  shock_val <- if (shock_direction == "positive") 1L else 0L
  sub <- df %>%
    filter(meal_type == meal, positive_shock == shock_val)
  if (nrow(sub) < 30) {
    warning(sprintf("Small sample: %s / %s / %s (N=%d)",
                    sample_label, meal, shock_direction, nrow(sub)))
  }
  fml <- as.formula(paste0("diff_real_wage ~ diff_log_yield_hat | ", fe_rhs))
  mod <- feols(fml, data = sub, cluster = ~District)
  est <- extract_coef(mod)
  if (is.null(est)) return(NULL)
  est %>%
    mutate(
      sample       = sample_label,
      meal_type    = meal,
      meal_label   = ifelse(meal == "Three", "3 meals (attached)", "0 meals (casual)"),
      shock_sign   = shock_direction,
      shock_label  = ifelse(shock_direction == "positive", "Positive shock", "Negative shock"),
      cell_id      = paste(meal_label, shock_label, sep = " | "),
      sd_yield_hat = sd(sub$diff_log_yield_hat, na.rm = TRUE),
      mean_wage    = mean(sub$real_wage, na.rm = TRUE),
      n_districts  = length(unique(sub$District)),
      sig_status   = sig_flag(p),
      stringsAsFactors = FALSE
    )
}

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 1 — Define yield shocks from first-stage residuals                    ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 1 — Construct yield shock variables\n")

yield_shocks <- read_csv(
  file.path(ROOT, "output/stage1/fitted/yield_hat_2017_2023.csv"),
  show_col_types = FALSE
) %>%
  rename(District = district, growing_season = season) %>%
  mutate(
    yield_shock    = residual,
    shock_abs      = abs(residual),
    positive_shock = as.integer(residual > 0),
    shock_sign     = ifelse(residual > 0, "positive", "negative")
  ) %>%
  distinct(District, growing_season, year, .keep_all = TRUE)

cat(sprintf("  District-season-year cells: N = %d\n", nrow(yield_shocks)))
cat(sprintf("  Positive shocks: %d (%.1f%%)\n",
            sum(yield_shocks$positive_shock),
            100 * mean(yield_shocks$positive_shock)))
cat(sprintf("  Negative shocks: %d (%.1f%%)\n",
            sum(yield_shocks$positive_shock == 0),
            100 * mean(yield_shocks$positive_shock == 0)))
cat(sprintf("  Mean shock: %.4f | SD: %.4f | Min: %.4f | Max: %.4f\n",
            mean(yield_shocks$yield_shock), sd(yield_shocks$yield_shock),
            min(yield_shocks$yield_shock), max(yield_shocks$yield_shock)))

## Balance by season (district-season-year level)
shock_by_season <- yield_shocks %>%
  group_by(growing_season, shock_sign) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(growing_season) %>%
  mutate(share = n_cells / sum(n_cells)) %>%
  ungroup()

cat("\n  Shock distribution by season (district-season-year cells):\n")
for (s in unique(shock_by_season$growing_season)) {
  sub <- shock_by_season %>% filter(growing_season == s)
  pos <- sub$share[sub$shock_sign == "positive"]
  if (length(pos) == 0) pos <- 0
  flag <- if (pos < 0.25 || pos > 0.75) " ** SEVERELY IMBALANCED **" else ""
  cat(sprintf("    %s: %.1f%% positive, %.1f%% negative%s\n",
              s, 100 * pos, 100 * (1 - pos), flag))
}

## Merge shocks into wage panel (0-meal and 3-meal only)
df_wage <- read_csv(
  file.path(ROOT, "data/Regression_data/df_2_merged_v2.csv"),
  show_col_types = FALSE
) %>%
  filter(meal_type %in% c("None", "Three"),
         !is.na(diff_real_wage), !is.na(diff_log_yield_hat)) %>%
  left_join(
    yield_shocks %>% select(District, growing_season, year,
                            yield_shock, shock_abs, positive_shock, shock_sign),
    by = c("District", "growing_season", "year")
  ) %>%
  filter(!is.na(yield_shock)) %>%
  mutate(
    meal3 = as.integer(meal_type == "Three"),
    meal_label = ifelse(meal_type == "Three", "3 meals (attached)", "0 meals (casual)")
  )

cat(sprintf("\n  Wage observations (0- and 3-meal): N = %d\n", nrow(df_wage)))
cat(sprintf("  Positive shock obs: %d (%.1f%%)\n",
            sum(df_wage$positive_shock), 100 * mean(df_wage$positive_shock)))
cat(sprintf("  Negative shock obs: %d (%.1f%%)\n\n",
            sum(df_wage$positive_shock == 0), 100 * mean(df_wage$positive_shock == 0)))

## Balance by season × meal (wage obs level)
shock_balance <- df_wage %>%
  group_by(growing_season, meal_label, shock_sign) %>%
  summarise(n_obs = n(), .groups = "drop") %>%
  group_by(growing_season, meal_label) %>%
  mutate(share = n_obs / sum(n_obs)) %>%
  ungroup() %>%
  mutate(
    imbalance_flag = ifelse(share < 0.20 | share > 0.80,
                            "SEVERELY IMBALANCED", "OK")
  )

imbalanced <- shock_balance %>% filter(imbalance_flag == "SEVERELY IMBALANCED")
if (nrow(imbalanced) > 0) {
  cat("  ** FLAG: Severely imbalanced shock split in subgroups:\n")
  for (i in seq_len(nrow(imbalanced))) {
    r <- imbalanced[i, ]
    cat(sprintf("    %s | %s | %s: %.1f%% of obs\n",
                r$growing_season, r$meal_label, r$shock_sign, 100 * r$share))
  }
} else {
  cat("  Shock split reasonably balanced across season × meal subgroups.\n")
}

write_csv(yield_shocks,    file.path(OUT, "tables", "yield_shocks_dsy.csv"))
write_csv(shock_by_season, file.path(OUT, "tables", "shock_distribution_by_season.csv"))
write_csv(shock_balance,   file.path(OUT, "tables", "shock_balance_by_season_meal.csv"))

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 2 — Split-sample regressions (pooled, all seasons)                    ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("\nSTEP 2 — Split-sample regressions (pooled seasons)\n")
cat("  Spec: diff_real_wage ~ diff_log_yield_hat | year + District^growing_season\n")

FE_POOLED <- "year + District^growing_season"
FE_BORO   <- "year + District"

split_pooled <- bind_rows(
  run_split_reg(df_wage, "Three", "positive", "Pooled", FE_POOLED),
  run_split_reg(df_wage, "Three", "negative", "Pooled", FE_POOLED),
  run_split_reg(df_wage, "None",  "positive", "Pooled", FE_POOLED),
  run_split_reg(df_wage, "None",  "negative", "Pooled", FE_POOLED)
)

for (i in seq_len(nrow(split_pooled))) {
  r <- split_pooled[i, ]
  cat(sprintf("  %s: coef = %.1f  SE = %.1f  p = %s  [%s]  N = %d\n",
              r$cell_id, r$coef, r$se, fmt_p(r$p), r$sig_status, r$nobs))
}

## Key comparison: 3-meal differential by shock direction
three_pos <- split_pooled %>% filter(meal_type == "Three", shock_sign == "positive")
three_neg <- split_pooled %>% filter(meal_type == "Three", shock_sign == "negative")
none_pos  <- split_pooled %>% filter(meal_type == "None",  shock_sign == "positive")
none_neg  <- split_pooled %>% filter(meal_type == "None",  shock_sign == "negative")

diff_attached <- three_neg$coef - three_pos$coef
diff_casual   <- none_neg$coef  - none_pos$coef
differential_neg <- three_neg$coef - none_neg$coef
differential_pos <- three_pos$coef - none_pos$coef

cat(sprintf("\n  3-meal coef (neg shock) - 3-meal coef (pos shock): %.1f\n", diff_attached))
cat(sprintf("  0-meal coef (neg shock) - 0-meal coef (pos shock): %.1f\n", diff_casual))
cat(sprintf("  Attached/casual differential, NEG shock: %.1f\n", differential_neg))
cat(sprintf("  Attached/casual differential, POS shock: %.1f\n", differential_pos))
if (abs(differential_neg) > abs(differential_pos) * 1.5) {
  cat("  >> Asymmetry pattern: differential LARGER in negative shock subsample.\n")
} else if (abs(differential_pos - differential_neg) < 50) {
  cat("  >> Differential similar across shock directions — symmetric pattern.\n")
}
cat("\n")

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 3 — Formal asymmetry test (triple interaction)                        ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 3 — Formal asymmetry test (interaction model)\n")
cat("  Spec: wage ~ yield_hat + meal3 + pos_shock + yield×meal3 + yield×pos_shock",
    "       + yield×meal3×pos_shock | FE\n")

m_asym <- feols(
  diff_real_wage ~ diff_log_yield_hat + meal3 + positive_shock +
    diff_log_yield_hat:meal3 + diff_log_yield_hat:positive_shock +
    diff_log_yield_hat:meal3:positive_shock |
    year + District^growing_season,
  data = df_wage,
  cluster = ~District
)

ct <- summary(m_asym)$coeftable
triple_nm <- "diff_log_yield_hat:meal3:positive_shock"
meal_int_nm <- "diff_log_yield_hat:meal3"

asym_results <- data.frame(
  term = c(
    "diff_log_yield_hat (0-meal, neg shock)",
    "diff_log_yield_hat:meal3 (attached/casual diff, neg shock)",
    "diff_log_yield_hat:positive_shock (0-meal, pos vs neg)",
    "diff_log_yield_hat:meal3:positive_shock (TRIPLE: asymmetry test)"
  ),
  coef = ct[c("diff_log_yield_hat", meal_int_nm,
              "diff_log_yield_hat:positive_shock", triple_nm), "Estimate"],
  se   = ct[c("diff_log_yield_hat", meal_int_nm,
              "diff_log_yield_hat:positive_shock", triple_nm), "Std. Error"],
  p    = ct[c("diff_log_yield_hat", meal_int_nm,
              "diff_log_yield_hat:positive_shock", triple_nm), "Pr(>|t|)"],
  stringsAsFactors = FALSE
) %>%
  mutate(sig_status = sig_flag(p))

joint_meal <- wald(m_asym, c(meal_int_nm, triple_nm))

cat(sprintf("  Attached/casual diff (neg shock): %.1f (SE %.1f) p = %s [%s]\n",
            asym_results$coef[2], asym_results$se[2],
            fmt_p(asym_results$p[2]), asym_results$sig_status[2]))
cat(sprintf("  TRIPLE interaction (asymmetry):   %.1f (SE %.1f) p = %s [%s]\n",
            asym_results$coef[4], asym_results$se[4],
            fmt_p(asym_results$p[4]), asym_results$sig_status[4]))
cat(sprintf("  Joint F (meal3 interactions): F(%d,%d) = %.3f  p = %s\n\n",
            joint_meal$df1, joint_meal$df2, joint_meal$stat, fmt_p(joint_meal$p)))

if (asym_results$coef[4] > 0 && asym_results$p[4] < 0.10) {
  cat("  >> Positive significant triple interaction: differential SHRINKS in good years",
      "(one-sided risk transfer).\n")
} else if (abs(asym_results$coef[4]) < asym_results$se[4]) {
  cat("  >> Triple interaction near zero: consistent with symmetric pass-through.\n")
}

write_csv(asym_results, file.path(OUT, "tables", "asymmetry_interaction_model.csv"))
write_csv(data.frame(
  test = "Joint F: meal3 interactions",
  F_stat = joint_meal$stat, df1 = joint_meal$df1,
  df2 = joint_meal$df2, p_value = joint_meal$p
), file.path(OUT, "tables", "asymmetry_joint_f_test.csv"))

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 4 — Boro season split-sample analysis                                 ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 4 — Boro season split-sample regressions\n")

df_boro <- df_wage %>% filter(growing_season == "Boro")

split_boro <- bind_rows(
  run_split_reg(df_boro, "Three", "positive", "Boro", FE_BORO),
  run_split_reg(df_boro, "Three", "negative", "Boro", FE_BORO),
  run_split_reg(df_boro, "None",  "positive", "Boro", FE_BORO),
  run_split_reg(df_boro, "None",  "negative", "Boro", FE_BORO)
)

for (i in seq_len(nrow(split_boro))) {
  r <- split_boro[i, ]
  cat(sprintf("  %s: coef = %.1f  SE = %.1f  p = %s  [%s]  N = %d\n",
              r$cell_id, r$coef, r$se, fmt_p(r$p), r$sig_status, r$nobs))
}

boro_none_pos <- split_boro %>% filter(meal_type == "None", shock_sign == "positive")
boro_none_neg <- split_boro %>% filter(meal_type == "None", shock_sign == "negative")
boro_diff_none <- boro_none_pos$coef - boro_none_neg$coef

cat(sprintf("\n  Boro 0-meal: pos shock coef (%.1f) - neg shock coef (%.1f) = %.1f\n",
            boro_none_pos$coef, boro_none_neg$coef, boro_diff_none))
if (boro_none_pos$coef > 0 && boro_none_neg$coef < 0) {
  cat("  >> Boro 0-meal: gains in good years, loses in bad — pure demand story.\n")
} else if (abs(boro_diff_none) < 100) {
  cat("  >> Boro 0-meal: similar pass-through in both shock directions.\n")
}
cat("\n")

split_all <- bind_rows(split_pooled, split_boro)
write_csv(split_all, file.path(OUT, "tables", "split_sample_coefficients.csv"))

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 5 — Visualizations                                                    ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 5 — Visualizations\n")

## 5a. 2×2 coefficient panel (pooled split-sample)
plot_split <- split_pooled %>%
  mutate(
    x_pos = ifelse(shock_sign == "positive", 1, 0),
    facet_meal = meal_label,
    label_txt = sprintf("%.0f\n(p=%s)", coef, fmt_p(p))
  )

p_2x2 <- ggplot(plot_split, aes(x = shock_label, y = coef, color = shock_sign)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.15, linewidth = 0.8) +
  geom_point(size = 4) +
  geom_text(aes(label = label_txt), vjust = -1.3, size = 3, color = "grey25") +
  facet_wrap(~ facet_meal, nrow = 1) +
  scale_color_manual(values = c("positive" = "#009E73", "negative" = "#D55E00"),
                     labels = c("Negative shock", "Positive shock"), guide = "none") +
  labs(
    title = "Pass-Through by Meal Type and Yield Shock Direction",
    subtitle = "Split-sample regressions (pooled seasons) | 95% CI | Clustered SE",
    x = "First-stage yield shock direction",
    y = "Pass-through coefficient (BDT/day)"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold")) +
  coord_cartesian(ylim = range(c(plot_split$ci_lo, plot_split$ci_hi, 0)) * c(1.25, 1.25))

ggsave(file.path(OUT, "figures", "asymmetry_2x2_split_sample.png"),
       p_2x2, width = 9, height = 5, dpi = 300)
ggsave(file.path(OUT, "figures", "asymmetry_2x2_split_sample.pdf"),
       p_2x2, width = 9, height = 5)

## 5b. Side-by-side pass-through by meal, shock direction
plot_side <- split_pooled %>%
  mutate(
    group = interaction(meal_label, shock_label, sep = "\n"),
    group = factor(group, levels = c(
      "0 meals (casual)\nNegative shock",
      "0 meals (casual)\nPositive shock",
      "3 meals (attached)\nNegative shock",
      "3 meals (attached)\nPositive shock"
    ))
  )

p_side <- ggplot(plot_side, aes(x = group, y = coef, fill = shock_sign)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2, linewidth = 0.7) +
  geom_col(width = 0.55, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.0f (p=%s)", coef, fmt_p(p))),
            vjust = -0.8, size = 2.8, color = "grey20") +
  scale_fill_manual(values = c("negative" = "#D55E00", "positive" = "#009E73"),
                    name = "Shock") +
  labs(
    title = "Yield Pass-Through Asymmetry: Attached vs Casual Workers",
    subtitle = "Implied pass-through in positive vs negative shock subsamples",
    x = NULL, y = "Pass-through coefficient (BDT/day)"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(size = 9), panel.grid.minor = element_blank()) +
  coord_cartesian(ylim = range(c(plot_side$ci_lo, plot_side$ci_hi, 0)) * c(1.3, 1.3))

ggsave(file.path(OUT, "figures", "asymmetry_side_by_side.png"),
       p_side, width = 9, height = 5.5, dpi = 300)

## 5c. Boro 2×2 panel
plot_boro <- split_boro %>%
  mutate(label_txt = sprintf("%.0f\n(p=%s)", coef, fmt_p(p)))

p_boro <- ggplot(plot_boro, aes(x = shock_label, y = coef, color = shock_sign)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.15, linewidth = 0.8) +
  geom_point(size = 4) +
  geom_text(aes(label = label_txt), vjust = -1.3, size = 3, color = "grey25") +
  facet_wrap(~ meal_label, nrow = 1) +
  scale_color_manual(values = c("positive" = "#009E73", "negative" = "#D55E00"), guide = "none") +
  labs(
    title = "Boro Season: Pass-Through by Meal Type and Shock Direction",
    subtitle = "FE: year + District | 95% CI",
    x = "Shock direction", y = "Pass-through coefficient (BDT/day)"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))

ggsave(file.path(OUT, "figures", "asymmetry_boro_2x2.png"),
       p_boro, width = 9, height = 5, dpi = 300)

## 5d. Yield shock distribution
p_density <- ggplot(yield_shocks, aes(x = yield_shock, fill = growing_season)) +
  geom_histogram(bins = 40, alpha = 0.65, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
  facet_wrap(~ growing_season, nrow = 1) +
  scale_fill_manual(values = c(Boro = "#0072B2", Aus = "#E69F00", Aman = "#CC79A7"),
                    name = "Season") +
  labs(
    title = "Distribution of First-Stage Yield Shocks",
    subtitle = "Residual = actual Δlog(yield) − predicted Δlog(yield) | district-season-year cells",
    x = "Yield shock (first-stage residual)", y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave(file.path(OUT, "figures", "yield_shock_distribution.png"),
       p_density, width = 10, height = 5, dpi = 300)

p_shock_bar <- yield_shocks %>%
  count(growing_season, shock_sign) %>%
  group_by(growing_season) %>%
  mutate(share = n / sum(n)) %>%
  ggplot(aes(x = growing_season, y = share, fill = shock_sign)) +
  geom_col(position = "stack", width = 0.6) +
  scale_fill_manual(values = c("negative" = "#D55E00", "positive" = "#009E73"),
                    labels = c("Negative", "Positive"), name = "Shock") +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  labs(
    title = "Share of Positive vs Negative Yield Shocks by Season",
    x = "Season", y = "Share of district-season-year cells"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(OUT, "figures", "yield_shock_balance_by_season.png"),
       p_shock_bar, width = 7, height = 4.5, dpi = 300)

cat("  Saved asymmetry_2x2_split_sample.png, asymmetry_side_by_side.png,\n")
cat("       asymmetry_boro_2x2.png, yield_shock_distribution.png\n\n")

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 6 — Economic interpretation table                                       ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 6 — Economic interpretation table\n")

econ_tbl <- split_all %>%
  mutate(
    wage_effect_1sd = coef * sd_yield_hat,
    pct_mean_wage   = 100 * wage_effect_1sd / mean_wage,
    sig_status      = sig_flag(p)
  ) %>%
  select(
    sample, cell_id, meal_type, meal_label, shock_sign,
    coef, se, p, sig_status, nobs, n_districts,
    sd_yield_hat, mean_wage,
    wage_effect_1sd, pct_mean_wage
  ) %>%
  arrange(sample, meal_type, shock_sign)

for (i in seq_len(nrow(econ_tbl))) {
  r <- econ_tbl[i, ]
  cat(sprintf("  %s [%s]: 1-SD shock → %.1f BDT/day (%.1f%% of mean wage) [%s]\n",
              r$cell_id, r$sample, r$wage_effect_1sd, r$pct_mean_wage, r$sig_status))
}

write_csv(econ_tbl, file.path(OUT, "tables", "asymmetry_economic_interpretation.csv"))

econ_html <- econ_tbl %>%
  mutate(
    coef_se = sprintf("%.1f (%.1f)", coef, se),
    p_fmt   = fmt_p(p),
    effect  = sprintf("%.1f BDT (%.1f%% mean wage)", wage_effect_1sd, pct_mean_wage)
  ) %>%
  select(Sample = sample, Cell = cell_id, `Coef (SE)` = coef_se,
         `p-value` = p_fmt, Status = sig_status, N = nobs, `1-SD effect` = effect) %>%
  kbl(caption = "Asymmetry Test: Economic Magnitudes of Pass-Through",
      format = "html", booktabs = TRUE) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = TRUE, font_size = 12) %>%
  footnote(
    general = paste(
      "1-SD effect = coefficient × SD(Δ fitted log yield) within each subsample.",
      "Pct mean wage = 1-SD effect / mean daily wage × 100.",
      "Shock split based on sign of first-stage yield residual.",
      "Pooled FE: year + District×Season. Boro FE: year + District."
    ),
    general_title = "Note: ",
    footnote_as_chunk = TRUE
  )

writeLines(as.character(econ_html),
           file.path(OUT, "tables", "asymmetry_economic_interpretation.html"))

## ── Summary report ──────────────────────────────────────────────────────────── ##
interpretation <- if (asym_results$p[4] < 0.10 && asym_results$coef[4] > 0) {
  "One-sided risk transfer: attached/casual differential shrinks significantly in positive shock years."
} else if (asym_results$p[4] < 0.10 && asym_results$coef[4] < 0) {
  "Asymmetric pass-through: differential expands in positive shock years."
} else if (abs(differential_neg) > abs(differential_pos) * 1.25) {
  "Directionally suggestive asymmetry (larger differential in bad years) but triple interaction imprecise."
} else {
  "No strong evidence of asymmetry — differential similar across shock directions."
}

report <- c(
  "# Yield Pass-Through Asymmetry Test",
  sprintf("Generated: %s", Sys.time()),
  "",
  "## Step 1 — Shock distribution",
  sprintf("- District-season-year cells: %d", nrow(yield_shocks)),
  sprintf("- Positive shocks: %.1f%%", 100 * mean(yield_shocks$positive_shock)),
  sprintf("- Negative shocks: %.1f%%", 100 * mean(1 - yield_shocks$positive_shock)),
  if (nrow(imbalanced) > 0) "- **FLAG**: Severely imbalanced subgroups — see shock_balance_by_season_meal.csv" else "- Shock split reasonably balanced.",
  "",
  "## Step 2 — Split-sample (pooled)",
  capture.output(print(split_pooled %>% select(cell_id, coef, se, p, sig_status, nobs))),
  sprintf("- Attached/casual differential, NEG shock: %.1f", differential_neg),
  sprintf("- Attached/casual differential, POS shock: %.1f", differential_pos),
  "",
  "## Step 3 — Formal asymmetry test",
  capture.output(print(asym_results)),
  sprintf("- Joint F (meal3 interactions): p = %s", fmt_p(joint_meal$p)),
  "",
  "## Step 4 — Boro split-sample",
  capture.output(print(split_boro %>% select(cell_id, coef, se, p, sig_status, nobs))),
  sprintf("- Boro 0-meal pos vs neg shock difference: %.1f", boro_diff_none),
  "",
  "## Interpretation",
  interpretation,
  "",
  "## Outputs (output/stage2/asymmetry/)",
  "- tables/split_sample_coefficients.csv",
  "- tables/asymmetry_interaction_model.csv",
  "- tables/asymmetry_economic_interpretation.csv/.html",
  "- figures/asymmetry_2x2_split_sample.png",
  "- figures/asymmetry_side_by_side.png",
  "- figures/asymmetry_boro_2x2.png",
  "- figures/yield_shock_distribution.png"
)

writeLines(report, file.path(OUT, "summary", "asymmetry_test_report.md"))

save(m_asym, split_all, asym_results, econ_tbl, yield_shocks,
     file = file.path(OUT, "models", "asymmetry_models.RData"))

cat("\nSaved asymmetry_test_report.md\n")
cat("=== ASYMMETRY TEST COMPLETE ===\n")
