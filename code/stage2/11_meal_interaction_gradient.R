## 11_meal_interaction_gradient.R
## Investigate pass-through gradient via interaction model vs separate regressions
## Compares full sample and season-specific (Boro, Aus, Aman) estimates
## Output: output/stage2/meal_monotonicity/

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
OUT  <- file.path(ROOT, "output/stage2/meal_monotonicity")
for (d in c("", "tables", "figures", "models", "summary"))
  dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)

MEAL_LEVELS <- c("None", "One", "Two", "Three")
MEAL_LABELS <- c("0 meals", "1 meal", "2 meals", "3 meals")
names(MEAL_LABELS) <- MEAL_LEVELS
SEASONS     <- c("Boro", "Aus", "Aman")
SAMPLES     <- c("Pooled", SEASONS)

cat("=== MEAL INTERACTION GRADIENT ANALYSIS ===\n")
cat("Output directory:", OUT, "\n\n")

## ── Helpers ───────────────────────────────────────────────────────────────── ##
fmt_p <- function(p) {
  vapply(p, function(x) {
    if (is.na(x)) return("NA")
    if (x < 0.001) return("<0.001")
    sprintf("%.3f", x)
  }, character(1))
}

sig_flag <- function(p, alpha = 0.10) {
  ifelse(is.na(p), "NA",
         ifelse(p < 0.01, "significant ***",
                ifelse(p < 0.05, "significant **",
                       ifelse(p < 0.10, "significant *",
                              "imprecise (not significant)"))))
}

coef_cell <- function(est, se, p = NULL) {
  if (is.null(p)) sprintf("%.1f (%.1f)", est, se)
  else sprintf("%.1f (%.1f) %s", est, se, sig_flag(p))
}

## Implied total pass-through per meal from interaction model (delta-method SE)
implied_totals <- function(mod, meal_levels = MEAL_LEVELS) {
  ct  <- summary(mod)$coeftable
  bnm <- "diff_log_yield_hat"
  if (!bnm %in% rownames(ct)) stop("Baseline yield term not found.")
  V   <- vcov(mod)
  b0  <- ct[bnm, "Estimate"]
  df2 <- mod$nobs - length(coef(mod))

  bind_rows(lapply(meal_levels, function(m) {
    if (m == "None") {
      est <- b0
      se  <- sqrt(V[bnm, bnm])
    } else {
      int_nm <- paste0(bnm, ":meal_type", m)
      est <- b0 + coef(mod)[int_nm]
      se  <- sqrt(V[bnm, bnm] + V[int_nm, int_nm] + 2 * V[bnm, int_nm])
    }
    data.frame(
      meal_type  = m,
      meal_label = MEAL_LABELS[m],
      meal_num   = which(meal_levels == m) - 1L,
      coef       = est,
      se         = se,
      p          = 2 * pt(-abs(est / se), df = df2),
      ci_lo      = est - 1.96 * se,
      ci_hi      = est + 1.96 * se,
      stringsAsFactors = FALSE
    )
  }))
}

## Separate regression per meal category
separate_by_meal <- function(df, fe_rhs) {
  fml <- as.formula(paste0("diff_real_wage ~ diff_log_yield_hat | ", fe_rhs))
  bind_rows(lapply(MEAL_LEVELS, function(m) {
    sub <- df %>% filter(meal_type == m)
    mod <- feols(fml, data = sub, cluster = ~District)
    ct <- summary(mod)$coeftable
    nm <- "diff_log_yield_hat"
    data.frame(
      meal_type  = m,
      meal_label = MEAL_LABELS[m],
      meal_num   = which(MEAL_LEVELS == m) - 1L,
      coef       = ct[nm, "Estimate"],
      se         = ct[nm, "Std. Error"],
      p          = ct[nm, "Pr(>|t|)"],
      ci_lo      = ct[nm, "Estimate"] - 1.96 * ct[nm, "Std. Error"],
      ci_hi      = ct[nm, "Estimate"] + 1.96 * ct[nm, "Std. Error"],
      n_obs      = nobs(mod),
      stringsAsFactors = FALSE
    )
  }))
}

## Fit interaction model and extract all Step-1/2 outputs
fit_interaction_block <- function(df, sample_label, fe_rhs) {
  fml <- as.formula(sprintf(
    "diff_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal_type | %s",
    fe_rhs
  ))
  mod <- feols(fml, data = df, cluster = ~District)

  ct <- summary(mod)$coeftable
  bnm <- "diff_log_yield_hat"
  int_nms <- paste0(bnm, ":meal_type", MEAL_LEVELS[-1])

  joint <- wald(mod, int_nms)

  block <- data.frame(
    sample          = sample_label,
    n_obs           = nobs(mod),
    n_districts     = length(unique(df$District)),
    baseline_coef   = ct[bnm, "Estimate"],
    baseline_se     = ct[bnm, "Std. Error"],
    baseline_p      = ct[bnm, "Pr(>|t|)"],
    int_one_coef    = ct[int_nms[1], "Estimate"],
    int_one_se      = ct[int_nms[1], "Std. Error"],
    int_one_p       = ct[int_nms[1], "Pr(>|t|)"],
    int_two_coef    = ct[int_nms[2], "Estimate"],
    int_two_se      = ct[int_nms[2], "Std. Error"],
    int_two_p       = ct[int_nms[2], "Pr(>|t|)"],
    int_three_coef  = ct[int_nms[3], "Estimate"],
    int_three_se    = ct[int_nms[3], "Std. Error"],
    int_three_p     = ct[int_nms[3], "Pr(>|t|)"],
    joint_F         = joint$stat,
    joint_p         = joint$p,
    joint_df1       = joint$df1,
    joint_df2       = joint$df2,
    stringsAsFactors = FALSE
  )

  list(
    model   = mod,
    block   = block,
    implied = implied_totals(mod) %>% mutate(sample = sample_label),
    separate = separate_by_meal(df, fe_rhs) %>% mutate(sample = sample_label)
  )
}

## Gradient plot: interaction implied totals with optional separate overlay
plot_gradient <- function(implied_df, separate_df, title, subtitle, filename) {
  imp <- implied_df %>%
    mutate(
      meal_label = factor(meal_label, levels = MEAL_LABELS),
      sig_label  = sig_flag(p),
      label_txt  = sprintf("%.0f\n(p=%s)", coef, fmt_p(p))
    )
  sep <- separate_df %>%
    mutate(meal_label = factor(meal_label, levels = MEAL_LABELS))

  meal_cols <- c("0 meals" = "#E69F00", "1 meal" = "#56B4E9",
                 "2 meals" = "#009E73", "3 meals" = "#D55E00")

  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
    geom_errorbar(data = imp,
                  aes(x = meal_label, ymin = ci_lo, ymax = ci_hi),
                  width = 0.15, linewidth = 0.9, color = "#0072B2") +
    geom_point(data = imp,
               aes(x = meal_label, y = coef, color = meal_label),
               size = 4, shape = 16) +
    geom_errorbar(data = sep,
                  aes(x = meal_label, ymin = ci_lo, ymax = ci_hi),
                  width = 0.25, linewidth = 0.6, color = "grey55", linetype = "dotted") +
    geom_point(data = sep,
               aes(x = meal_label, y = coef),
               size = 3, shape = 17, color = "grey40") +
    geom_text(data = imp,
              aes(x = meal_label, y = coef, label = label_txt),
              vjust = -1.5, size = 2.8, color = "grey20", lineheight = 0.85) +
    scale_color_manual(values = meal_cols, guide = "none") +
    scale_x_discrete(labels = MEAL_LABELS) +
    labs(
      title    = title,
      subtitle = subtitle,
      x        = "Meal provision",
      y        = "Pass-through coefficient (BDT/day per unit change in fitted log yield)",
      caption  = "Filled circles + solid bars = interaction model (implied total). Triangles + dotted bars = separate regressions."
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold")
    ) +
    coord_cartesian(ylim = range(c(imp$ci_lo, imp$ci_hi, sep$ci_lo, sep$ci_hi, 0)) * c(1.2, 1.2))

  ggsave(file.path(OUT, "figures", paste0(filename, ".png")),
         p, width = 8, height = 5.5, dpi = 300)
  ggsave(file.path(OUT, "figures", paste0(filename, ".pdf")),
         p, width = 8, height = 5.5)
  p
}

## ════════════════════════════════════════════════════════════════════════════ ##
## Load panel (same sample as monotonicity / main analysis)                    ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("Loading panel...\n")

df_all <- read_csv(
  file.path(ROOT, "data/Regression_data/df_2_merged_v2.csv"),
  show_col_types = FALSE
) %>%
  mutate(meal_type = relevel(factor(meal_type), ref = "None")) %>%
  filter(!is.na(diff_real_wage), !is.na(diff_log_yield_hat))

cat(sprintf("  Full panel: N = %d | districts = %d | years = %s\n\n",
            nrow(df_all), length(unique(df_all$District)),
            paste(sort(unique(df_all$year)), collapse = ", ")))

## FE structure: pooled uses District×Season (main spec); within-season uses District
FE_POOLED <- "year + District^growing_season"
FE_SEASON <- "year + District"

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 1 — Interaction model: full sample                                     ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 1 — Interaction model (full sample, pooled seasons)\n")
cat(sprintf("  Spec: diff_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal_type | %s\n",
            FE_POOLED))

res_pooled <- fit_interaction_block(df_all, "Pooled", FE_POOLED)
b <- res_pooled$block

cat(sprintf("  N = %d | districts = %d\n", b$n_obs, b$n_districts))
cat(sprintf("  Baseline (0-meal): coef = %.2f  SE = %.2f  p = %s  [%s]\n",
            b$baseline_coef, b$baseline_se, fmt_p(b$baseline_p), sig_flag(b$baseline_p)))
cat(sprintf("  Interaction 1-meal:  coef = %.2f  SE = %.2f  p = %s  [%s]\n",
            b$int_one_coef, b$int_one_se, fmt_p(b$int_one_p), sig_flag(b$int_one_p)))
cat(sprintf("  Interaction 2-meal:  coef = %.2f  SE = %.2f  p = %s  [%s]\n",
            b$int_two_coef, b$int_two_se, fmt_p(b$int_two_p), sig_flag(b$int_two_p)))
cat(sprintf("  Interaction 3-meal:  coef = %.2f  SE = %.2f  p = %s  [%s]\n",
            b$int_three_coef, b$int_three_se, fmt_p(b$int_three_p), sig_flag(b$int_three_p)))
cat(sprintf("  Joint F-test (all interactions = 0): F(%d,%d) = %.3f  p = %s\n\n",
            b$joint_df1, b$joint_df2, b$joint_F, fmt_p(b$joint_p)))

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 2 — Season-specific interaction models                                 ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 2 — Season-specific interaction models\n")

res_season <- lapply(SEASONS, function(s) {
  df_s <- df_all %>% filter(growing_season == s)
  cat(sprintf("\n  --- %s ---\n", s))
  cat(sprintf("  N = %d | districts = %d\n", nrow(df_s), length(unique(df_s$District))))
  out <- fit_interaction_block(df_s, s, FE_SEASON)
  bk <- out$block
  cat(sprintf("  Baseline (0-meal): coef = %.2f  SE = %.2f  p = %s  [%s]\n",
              bk$baseline_coef, bk$baseline_se, fmt_p(bk$baseline_p), sig_flag(bk$baseline_p)))
  if (s == "Boro" && bk$baseline_coef > 0)
    cat("  >> Boro 0-meal coefficient is POSITIVE — consistent with casual-labor demand channel.\n")
  cat(sprintf("  Interaction 1-meal:  coef = %.2f  p = %s  [%s]\n",
              bk$int_one_coef, fmt_p(bk$int_one_p), sig_flag(bk$int_one_p)))
  cat(sprintf("  Interaction 2-meal:  coef = %.2f  p = %s  [%s]\n",
              bk$int_two_coef, fmt_p(bk$int_two_p), sig_flag(bk$int_two_p)))
  cat(sprintf("  Interaction 3-meal:  coef = %.2f  p = %s  [%s]\n",
              bk$int_three_coef, fmt_p(bk$int_three_p), sig_flag(bk$int_three_p)))
  cat(sprintf("  Joint F-test: F(%d,%d) = %.3f  p = %s\n",
              bk$joint_df1, bk$joint_df2, bk$joint_F, fmt_p(bk$joint_p)))
  out
})
names(res_season) <- SEASONS

interaction_blocks <- bind_rows(
  res_pooled$block,
  bind_rows(lapply(res_season, `[[`, "block"))
)
cat("\n")

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 3 — Implied pass-through plots                                         ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 3 — Implied pass-through coefficient plots\n")

all_implied <- bind_rows(
  res_pooled$implied,
  bind_rows(lapply(res_season, `[[`, "implied"))
)
all_separate <- bind_rows(
  res_pooled$separate,
  bind_rows(lapply(res_season, `[[`, "separate"))
)

plot_gradient(
  res_pooled$implied, res_pooled$separate,
  title    = "Implied Pass-Through Gradient — Full Sample (Pooled Seasons)",
  subtitle = "Interaction model implied totals vs. separate regressions | 95% CI",
  filename = "implied_gradient_pooled"
)

for (s in SEASONS) {
  plot_gradient(
    res_season[[s]]$implied, res_season[[s]]$separate,
    title    = sprintf("Implied Pass-Through Gradient — %s Season", s),
    subtitle = sprintf("%s | FE: year + District | 95%% CI", s),
    filename = sprintf("implied_gradient_%s", tolower(s))
  )
}

## Combined faceted plot
facet_df <- all_implied %>%
  mutate(
    sample = factor(sample, levels = SAMPLES),
    meal_label = factor(meal_label, levels = MEAL_LABELS),
    label_txt = sprintf("%.0f\n(p=%s)", coef, fmt_p(p))
  )

p_facet <- ggplot(facet_df, aes(x = meal_label, y = coef, color = meal_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.12, linewidth = 0.7) +
  geom_point(size = 2.8) +
  geom_text(aes(label = label_txt), vjust = -1.2, size = 2.2, color = "grey25") +
  facet_wrap(~ sample, nrow = 1) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#D55E00"), guide = "none") +
  labs(
    title = "Interaction Model: Implied Pass-Through by Meal Category",
    subtitle = "Full sample and season-specific estimates | 95% CI",
    x = "Meal provision", y = "Implied pass-through coefficient"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))

ggsave(file.path(OUT, "figures", "implied_gradient_all_samples.png"),
       p_facet, width = 14, height = 5, dpi = 300)
cat("  Saved implied_gradient_*.png and implied_gradient_all_samples.png\n\n")

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 4 — 0-meal pass-through by season                                      ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 4 — 0-meal pass-through by season (separate regressions)\n")

none_by_season <- bind_rows(lapply(SEASONS, function(s) {
  sub <- df_all %>% filter(meal_type == "None", growing_season == s)
  mod <- feols(
    diff_real_wage ~ diff_log_yield_hat | year + District,
    data = sub, cluster = ~District
  )
  ct <- summary(mod)$coeftable
  nm <- "diff_log_yield_hat"
  data.frame(
    season     = s,
    coef       = ct[nm, "Estimate"],
    se         = ct[nm, "Std. Error"],
    p          = ct[nm, "Pr(>|t|)"],
    ci_lo      = ct[nm, "Estimate"] - 1.96 * ct[nm, "Std. Error"],
    ci_hi      = ct[nm, "Estimate"] + 1.96 * ct[nm, "Std. Error"],
    n_obs      = nobs(mod),
    n_districts = length(unique(sub$District)),
    flag       = sig_flag(ct[nm, "Pr(>|t|)"]),
    stringsAsFactors = FALSE
  )
}))

for (i in seq_len(nrow(none_by_season))) {
  r <- none_by_season[i, ]
  note <- if (r$coef > 0) " [POSITIVE]" else ""
  cat(sprintf("  %s: coef = %.2f  SE = %.2f  p = %s  [%s]%s\n",
              r$season, r$coef, r$se, fmt_p(r$p), r$flag, note))
}

pos_seasons <- none_by_season$season[none_by_season$coef > 0]
if (length(pos_seasons) == 1 && pos_seasons == "Boro") {
  cat("\n  >> Positive 0-meal coefficient is concentrated in Boro — demand-channel story.\n")
} else if (length(pos_seasons) == length(SEASONS)) {
  cat("\n  ** FLAG: Positive 0-meal coefficient in ALL seasons — investigate further.\n")
} else if (length(pos_seasons) == 0) {
  cat("\n  >> No season shows a positive 0-meal coefficient in separate regressions.\n")
} else {
  cat(sprintf("\n  >> Positive 0-meal coefficient in: %s\n", paste(pos_seasons, collapse = ", ")))
}
cat("\n")

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 5 — Gradient comparison table                                          ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 5 — Gradient comparison table (separate vs interaction)\n")

comparison_long <- bind_rows(
  all_separate %>%
    transmute(sample, meal_type, meal_label, meal_num,
              approach = "Separate regression",
              coef, se, p, ci_lo, ci_hi),
  all_implied %>%
    transmute(sample, meal_type, meal_label, meal_num,
              approach = "Interaction model (implied)",
              coef, se, p, ci_lo, ci_hi)
) %>%
  left_join(
    all_separate %>% select(sample, meal_type, se_sep = se),
    by = c("sample", "meal_type")
  ) %>%
  mutate(
    sig_status = sig_flag(p),
    se_ratio   = ifelse(approach == "Interaction model (implied)",
                        se_sep / se, NA_real_),
    coef_se    = sprintf("%.1f (%.1f)", coef, se),
    cell       = sprintf("%.1f (%.1f) [%s]", coef, se, sig_flag(p))
  )

## Wide table for export
comparison_wide <- comparison_long %>%
  select(sample, meal_label, approach, coef, se, p, sig_status, se_ratio) %>%
  pivot_wider(
    id_cols = c(sample, meal_label),
    names_from = approach,
    values_from = c(coef, se, p, sig_status, se_ratio),
    names_glue = "{approach}_{.value}"
  ) %>%
  mutate(
    se_gain_pct = (`Separate regression_se` - `Interaction model (implied)_se`) /
      `Separate regression_se` * 100,
    gradient_sep = `Separate regression_coef` >= lag(`Separate regression_coef`),
    gradient_int = `Interaction model (implied)_coef` >= lag(`Interaction model (implied)_coef`)
  ) %>%
  group_by(sample) %>%
  mutate(
    monotone_sep = all(gradient_sep, na.rm = TRUE),
    monotone_int = all(gradient_int, na.rm = TRUE)
  ) %>%
  ungroup()

## Display table (coef + SE cells)
display_tbl <- comparison_long %>%
  select(sample, meal_label, approach, cell, se_ratio) %>%
  pivot_wider(
    id_cols = c(sample, meal_label),
    names_from = approach,
    values_from = c(cell, se_ratio),
    names_glue = "{approach}_{.value}"
  ) %>%
  arrange(match(sample, SAMPLES), meal_label)

write_csv(comparison_wide, file.path(OUT, "tables", "gradient_comparison.csv"))
write_csv(display_tbl,   file.path(OUT, "tables", "gradient_comparison_display.csv"))
write_csv(interaction_blocks, file.path(OUT, "tables", "interaction_model_results.csv"))
write_csv(all_implied,  file.path(OUT, "tables", "implied_passthrough_by_sample.csv"))
write_csv(none_by_season, file.path(OUT, "tables", "none_meal_by_season.csv"))

tbl_html <- display_tbl %>%
  rename(
    `Meal category` = meal_label,
    `Separate regression` = `Separate regression_cell`,
    `Interaction (implied)` = `Interaction model (implied)_cell`,
    `SE ratio (sep/int)` = `Interaction model (implied)_se_ratio`
  ) %>%
  select(Sample = sample, `Meal category`, `Separate regression`,
         `Interaction (implied)`, `SE ratio (sep/int)`) %>%
  kbl(caption = "Pass-Through Gradient: Separate Regressions vs Interaction Model",
      format = "html", booktabs = TRUE, escape = FALSE) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = TRUE, font_size = 12) %>%
  pack_rows(index = setNames(table(display_tbl$sample), names(table(display_tbl$sample)))) %>%
  footnote(
    general = paste(
      "Coefficients are pass-through of change in fitted log yield on change in real daily wage (BDT/day).",
      "Bracketed labels: significant * p<0.10, ** p<0.05, *** p<0.01; otherwise imprecise.",
      "SE ratio > 1 means interaction model has smaller SE (more efficient).",
      "Pooled FE: year + District×Season. Season FE: year + District."
    ),
    general_title = "Note: ",
    footnote_as_chunk = TRUE
  )

writeLines(as.character(tbl_html), file.path(OUT, "tables", "gradient_comparison.html"))
cat("  Saved gradient_comparison.csv/.html\n\n")

## ── Efficiency summary ─────────────────────────────────────────────────────── ##
eff_summary <- comparison_long %>%
  filter(approach == "Interaction model (implied)") %>%
  group_by(sample) %>%
  summarise(
    mean_se_ratio = mean(se_ratio, na.rm = TRUE),
    pct_tighter   = mean(se_ratio > 1, na.rm = TRUE) * 100,
    .groups = "drop"
  )

cat("SE efficiency (interaction vs separate, mean ratio > 1 = tighter):\n")
for (i in seq_len(nrow(eff_summary))) {
  r <- eff_summary[i, ]
  cat(sprintf("  %s: mean SE ratio = %.2f | interaction tighter in %.0f%% of meal categories\n",
              r$sample, r$mean_se_ratio, r$pct_tighter))
}

## ── Key findings summary ───────────────────────────────────────────────────── ##
boro_int <- res_season$Boro$implied
boro_joint_p <- res_season$Boro$block$joint_p
pooled_joint_p <- res_pooled$block$joint_p

gradient_verdict <- function(implied_df, joint_p) {
  mono <- all(diff(implied_df$coef) <= 0)
  any_sig <- any(implied_df$p < 0.10)
  three_sig <- implied_df$p[implied_df$meal_type == "Three"] < 0.10
  list(
    monotone = mono,
    any_sig  = any_sig,
    three_sig = three_sig,
    joint_sig = joint_p < 0.10,
    verdict = if (mono && (any_sig || joint_p < 0.10)) {
      "gradient present and partially significant"
    } else if (mono) {
      "directionally consistent but imprecise"
    } else {
      "non-monotonic or flat — rely on binary 3-vs-0 contrast"
    }
  )
}

verdicts <- list(
  Pooled = gradient_verdict(res_pooled$implied, pooled_joint_p),
  Boro   = gradient_verdict(res_season$Boro$implied, boro_joint_p),
  Aus    = gradient_verdict(res_season$Aus$implied, res_season$Aus$block$joint_p),
  Aman   = gradient_verdict(res_season$Aman$implied, res_season$Aman$block$joint_p)
)

## ── Save models & report ────────────────────────────────────────────────────── ##
save(
  res_pooled, res_season,
  interaction_blocks, all_implied, all_separate,
  none_by_season, comparison_wide, verdicts,
  file = file.path(OUT, "models", "interaction_gradient_models.RData")
)

report_lines <- c(
  "# Interaction Model Gradient Analysis",
  sprintf("Generated: %s", Sys.time()),
  "",
  "## Step 1 — Full sample interaction model",
  sprintf("- Baseline (0-meal): %.1f (SE %.1f), p = %s [%s]",
          b$baseline_coef, b$baseline_se, fmt_p(b$baseline_p), sig_flag(b$baseline_p)),
  sprintf("- Interaction 1-meal: %.1f (SE %.1f), p = %s [%s]",
          b$int_one_coef, b$int_one_se, fmt_p(b$int_one_p), sig_flag(b$int_one_p)),
  sprintf("- Interaction 2-meal: %.1f (SE %.1f), p = %s [%s]",
          b$int_two_coef, b$int_two_se, fmt_p(b$int_two_p), sig_flag(b$int_two_p)),
  sprintf("- Interaction 3-meal: %.1f (SE %.1f), p = %s [%s]",
          b$int_three_coef, b$int_three_se, fmt_p(b$int_three_p), sig_flag(b$int_three_p)),
  sprintf("- Joint F-test (all interactions = 0): F = %.3f, p = %s [%s]",
          b$joint_F, fmt_p(b$joint_p),
          if (b$joint_p < 0.10) "significant" else "not significant"),
  "",
  "## Step 2 — Season-specific interaction models",
  capture.output(print(interaction_blocks %>%
    select(sample, baseline_coef, baseline_p, int_three_coef, int_three_p, joint_F, joint_p))),
  "",
  "## Step 3 — Implied pass-through (interaction model)",
  capture.output(print(all_implied %>%
    select(sample, meal_label, coef, se, p) %>%
    mutate(sig = sig_flag(p)))),
  "",
  "## Step 4 — 0-meal coefficient by season",
  capture.output(print(none_by_season)),
  if (length(pos_seasons) == 1 && pos_seasons == "Boro")
    "- **Finding**: Positive 0-meal coefficient concentrated in Boro (demand channel)."
  else if (length(pos_seasons) == length(SEASONS))
    "- **FLAG**: Positive 0-meal coefficient in all seasons."
  else "",
  "",
  "## Verdicts",
  sprintf("- Pooled: %s", verdicts$Pooled$verdict),
  sprintf("- Boro: %s", verdicts$Boro$verdict),
  sprintf("- Aus: %s", verdicts$Aus$verdict),
  sprintf("- Aman: %s", verdicts$Aman$verdict),
  "",
  "## SE efficiency (interaction vs separate)",
  capture.output(print(eff_summary)),
  "",
  "## Outputs",
  "- tables/gradient_comparison.csv/.html",
  "- tables/interaction_model_results.csv",
  "- tables/implied_passthrough_by_sample.csv",
  "- tables/none_meal_by_season.csv",
  "- figures/implied_gradient_pooled.png",
  "- figures/implied_gradient_boro.png",
  "- figures/implied_gradient_aus.png",
  "- figures/implied_gradient_aman.png",
  "- figures/implied_gradient_all_samples.png"
)

writeLines(report_lines, file.path(OUT, "summary", "interaction_gradient_report.md"))
cat("\nSaved interaction_gradient_report.md\n")
cat("=== INTERACTION GRADIENT ANALYSIS COMPLETE ===\n")
