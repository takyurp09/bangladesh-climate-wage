## 13_extensive_margin.R
## Extensive margin: does share of attached (3-meal) workers fall in bad yield years?
## Output: output/stage2/extensive_margin/

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
OUT  <- file.path(ROOT, "output/stage2/extensive_margin")
for (d in c("", "tables", "figures", "models", "summary"))
  dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)

SEASONS <- c("Boro", "Aus", "Aman")

cat("=== EXTENSIVE MARGIN: ATTACHED LABOR SHARE ===\n")
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

extract_coef <- function(mod, param) {
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

run_ext_reg <- function(df, outcome, regressor, sample_label, fe_rhs) {
  fml <- as.formula(sprintf("%s ~ %s | %s", outcome, regressor, fe_rhs))
  mod <- tryCatch(
    feols(fml, data = df, cluster = ~District),
    error = function(e) { cat("  Model failed:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(mod)) return(NULL)
  est <- extract_coef(mod, regressor)
  est$outcome <- outcome
  est$sample  <- sample_label
  est$sig_status <- sig_flag(est$p)
  est
}

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 1 — Data structure audit                                              ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 1 — Data structure audit\n\n")

wage_raw <- read_csv(
  file.path(ROOT, "data/Regression_data/wage_by_growing_season.csv"),
  show_col_types = FALSE
) %>%
  filter(meal_type != "No_info", !is.na(real_wage), year >= 2017, year <= 2023)

merged <- read_csv(
  file.path(ROOT, "data/Regression_data/df_2_merged_v2.csv"),
  show_col_types = FALSE
)

## Raw wage: rows per district-season-year
dsy_counts_raw <- wage_raw %>%
  count(District, year, growing_season, name = "n_monthly_rows")

## Raw: rows per district-season-year-meal
dsy_meal_counts <- wage_raw %>%
  count(District, year, growing_season, meal_type, name = "n_monthly_rows")

## Merged panel: rows per district-season-year
dsy_counts_merged <- merged %>%
  count(District, year, growing_season, name = "n_panel_rows")

## Meal types present per dsy cell (raw)
meal_types_per_cell <- wage_raw %>%
  group_by(District, year, growing_season) %>%
  summarise(
    n_meal_types_present = n_distinct(meal_type),
    meal_types           = paste(sort(unique(meal_type)), collapse = ", "),
    .groups = "drop"
  )

structure_report <- c(
  "# Extensive Margin — Data Structure Audit",
  sprintf("Generated: %s", Sys.time()),
  "",
  "## Raw wage data (wage_by_growing_season.csv, 2017–2023)",
  sprintf("- Total rows: %d", nrow(wage_raw)),
  sprintf("- Unit: monthly wage survey records"),
  sprintf("- Columns include: District, year, growing_season, gender, meal_type, real_wage"),
  sprintf("- NO worker-count or respondent-weight variable available"),
  "",
  "### Rows per district-season-year (raw monthly records)",
  capture.output(print(summary(dsy_counts_raw$n_monthly_rows))),
  sprintf("- Distribution: %s",
          paste(names(table(dsy_counts_raw$n_monthly_rows)), table(dsy_counts_raw$n_monthly_rows),
                sep = "=", collapse = ", ")),
  "",
  "### Rows per district-season-year-meal (raw)",
  capture.output(print(summary(dsy_meal_counts$n_monthly_rows))),
  "",
  "### Meal categories present per district-season-year cell",
  capture.output(print(table(meal_types_per_cell$n_meal_types_present))),
  sprintf("- Cells with all 4 meal types: %d (%.1f%%)",
          sum(meal_types_per_cell$n_meal_types_present == 4),
          100 * mean(meal_types_per_cell$n_meal_types_present == 4)),
  sprintf("- Cells with only 1 meal type: %d (%.1f%%)",
          sum(meal_types_per_cell$n_meal_types_present == 1),
          100 * mean(meal_types_per_cell$n_meal_types_present == 1)),
  "",
  "## Merged analysis panel (df_2_merged_v2.csv)",
  sprintf("- Total rows: %d", nrow(merged)),
  sprintf("- Unit: one row per district × season × year × gender × meal_type"),
  sprintf("- After first-differencing for wage regressions; used here for presence checks"),
  "",
  "### Rows per district-season-year (merged panel)",
  capture.output(print(summary(dsy_counts_merged$n_panel_rows))),
  sprintf("- Max 8 rows = 4 meal types × 2 genders (when fully observed)"),
  "",
  "## Implications for share construction",
  "- Primary share_attached: fraction of monthly wage records in each district-season-year",
  "  cell that belong to the 3-meal (Three) category. This weights by underlying survey",
  "  records (months) rather than treating each gender-meal row equally.",
  "- Binary presence_attached: indicator that 3-meal has ≥1 non-missing wage record in cell.",
  "- LIMITATION: shares reflect survey record composition, not actual workforce headcount.",
  "- LIMITATION: 704 district-season-year cells have only 1 meal type in raw data — shares",
  "  in sparse cells are mechanically 0 or 1 and should be interpreted cautiously.",
  "- No continuous worker-count weighting available; monthly record counts are the finest",
  "  granularity offered by the data."
)

writeLines(structure_report, file.path(OUT, "summary", "data_structure_audit.md"))
write_csv(dsy_counts_raw,     file.path(OUT, "tables", "dsy_row_counts_raw.csv"))
write_csv(meal_types_per_cell, file.path(OUT, "tables", "meal_types_per_dsy_cell.csv"))

cat("  Raw wage: median", median(dsy_counts_raw$n_monthly_rows),
    "monthly rows per district-season-year\n")
cat("  Merged panel: median", median(dsy_counts_merged$n_panel_rows),
    "rows per district-season-year (gender × meal)\n")
cat("  Cells with all 4 meal types:", sum(meal_types_per_cell$n_meal_types_present == 4),
    "/", nrow(meal_types_per_cell), "\n")
cat("  Saved data_structure_audit.md\n\n")

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 2 — Construct extensive margin variables at district-season-year      ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 2 — Construct share_attached, share_casual, attached_to_casual_ratio\n")

## Monthly-record shares (primary)
meal_counts <- wage_raw %>%
  group_by(District, year, growing_season, meal_type) %>%
  summarise(n_monthly = n(), .groups = "drop") %>%
  group_by(District, year, growing_season) %>%
  mutate(
    share = n_monthly / sum(n_monthly),
    total_monthly = sum(n_monthly)
  ) %>%
  ungroup()

dsy_panel <- meal_counts %>%
  select(District, year, growing_season, total_monthly) %>%
  distinct() %>%
  left_join(
    meal_counts %>%
      filter(meal_type == "Three") %>%
      select(District, year, growing_season,
             share_attached = share, n_attached = n_monthly),
    by = c("District", "year", "growing_season")
  ) %>%
  left_join(
    meal_counts %>%
      filter(meal_type == "None") %>%
      select(District, year, growing_season,
             share_casual = share, n_casual = n_monthly),
    by = c("District", "year", "growing_season")
  ) %>%
  mutate(
    share_attached = replace_na(share_attached, 0),
    share_casual   = replace_na(share_casual, 0),
    n_attached     = replace_na(n_attached, 0),
    n_casual       = replace_na(n_casual, 0),
    presence_attached = as.integer(n_attached > 0),
    presence_casual   = as.integer(n_casual > 0),
    attached_to_casual_ratio = ifelse(n_casual > 0, n_attached / n_casual, NA_real_)
  )

## Binary presence from merged panel (alternative / robustness)
presence_merged <- merged %>%
  group_by(District, year, growing_season) %>%
  summarise(
    share_attached_rows = mean(meal_type == "Three"),
    share_casual_rows   = mean(meal_type == "None"),
    presence_attached_rows = as.integer(any(meal_type == "Three")),
    .groups = "drop"
  )

dsy_panel <- dsy_panel %>%
  left_join(presence_merged, by = c("District", "year", "growing_season"))

## Merge yield variables
yield_hat <- read_csv(
  file.path(ROOT, "output/stage1/fitted/yield_hat_2017_2023.csv"),
  show_col_types = FALSE
) %>%
  rename(District = district, growing_season = season) %>%
  distinct(District, growing_season, year, .keep_all = TRUE) %>%
  arrange(District, growing_season, year) %>%
  group_by(District, growing_season) %>%
  mutate(diff_log_yield_hat = yield_hat - lag(yield_hat)) %>%
  ungroup() %>%
  mutate(
    yield_shock    = residual,
    positive_shock = as.integer(residual > 0)
  ) %>%
  select(District, growing_season, year,
         yield_hat, diff_log_yield_hat, yield_shock, positive_shock, residual)

dsy_panel <- dsy_panel %>%
  left_join(yield_hat, by = c("District", "year", "growing_season")) %>%
  filter(!is.na(diff_log_yield_hat))

## Mean 3-meal wage for Step 6 selection-bias check
mean_three_wage <- wage_raw %>%
  filter(meal_type == "Three") %>%
  group_by(District, year, growing_season) %>%
  summarise(mean_three_wage = mean(real_wage, na.rm = TRUE), .groups = "drop")

dsy_panel <- dsy_panel %>%
  left_join(mean_three_wage, by = c("District", "year", "growing_season"))

## Shock severity quartiles within season (among negative shocks only)
neg_shock_labels <- dsy_panel %>%
  filter(yield_shock < 0) %>%
  group_by(growing_season) %>%
  mutate(
    neg_q = ntile(yield_shock, 4),
    shock_severity = case_when(
      neg_q == 4 ~ "mild_negative",
      neg_q == 1 ~ "severe_negative",
      TRUE       ~ "moderate_negative"
    )
  ) %>%
  ungroup() %>%
  select(District, year, growing_season, shock_severity, neg_q)

dsy_panel <- dsy_panel %>%
  left_join(neg_shock_labels, by = c("District", "year", "growing_season"))

## Summary statistics
summarise_var <- function(df, var, label) {
  data.frame(
    variable = label,
    sample   = "All",
    mean     = mean(df[[var]], na.rm = TRUE),
    sd       = sd(df[[var]], na.rm = TRUE),
    min      = min(df[[var]], na.rm = TRUE),
    max      = max(df[[var]], na.rm = TRUE),
    n        = sum(!is.na(df[[var]])),
    stringsAsFactors = FALSE
  )
}

sum_stats_all <- bind_rows(
  summarise_var(dsy_panel, "share_attached", "share_attached"),
  summarise_var(dsy_panel, "share_casual", "share_casual"),
  summarise_var(dsy_panel, "attached_to_casual_ratio", "attached_to_casual_ratio"),
  summarise_var(dsy_panel, "presence_attached", "presence_attached (binary)")
)

sum_stats_season <- bind_rows(lapply(SEASONS, function(s) {
  sub <- dsy_panel %>% filter(growing_season == s)
  bind_rows(
    summarise_var(sub, "share_attached", "share_attached"),
    summarise_var(sub, "share_casual", "share_casual"),
    summarise_var(sub, "attached_to_casual_ratio", "attached_to_casual_ratio")
  ) %>% mutate(sample = s)
}))

cat("  share_attached: mean =", round(mean(dsy_panel$share_attached), 3),
    "SD =", round(sd(dsy_panel$share_attached), 3), "\n")
cat("  share_casual:   mean =", round(mean(dsy_panel$share_casual), 3),
    "SD =", round(sd(dsy_panel$share_casual), 3), "\n")
cat("  N district-season-year cells:", nrow(dsy_panel), "\n\n")

write_csv(dsy_panel, file.path(OUT, "tables", "extensive_margin_panel.csv"))
write_csv(sum_stats_all,   file.path(OUT, "tables", "share_summary_stats_all.csv"))
write_csv(sum_stats_season, file.path(OUT, "tables", "share_summary_stats_by_season.csv"))

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 3 — Main extensive margin regressions                                 ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 3 — Main extensive margin regressions\n")
cat("  Spec: outcome ~ diff_log_yield_hat | year + District^growing_season\n")
cat("  (diff_log_yield_hat = change in first-stage predicted yield)\n\n")

FE_POOLED <- "year + District^growing_season"
FE_SEASON <- "year + District"
REG       <- "diff_log_yield_hat"
OUTCOMES  <- c("share_attached", "share_casual", "attached_to_casual_ratio")

dsy_ratio <- dsy_panel %>% filter(n_casual > 0, !is.na(attached_to_casual_ratio))

main_results <- bind_rows(
  lapply(c("share_attached", "share_casual"), function(y) {
    run_ext_reg(dsy_panel, y, REG, "Pooled", FE_POOLED)
  }),
  run_ext_reg(dsy_ratio, "attached_to_casual_ratio", REG, "Pooled", FE_POOLED),
  unlist(lapply(SEASONS, function(s) {
    sub <- dsy_panel %>% filter(growing_season == s)
    sub_r <- dsy_ratio %>% filter(growing_season == s)
    c(
      list(run_ext_reg(sub, "share_attached", REG, s, FE_SEASON)),
      list(run_ext_reg(sub, "share_casual",   REG, s, FE_SEASON)),
      list(run_ext_reg(sub_r, "attached_to_casual_ratio", REG, s, FE_SEASON))
    )
  }), recursive = FALSE)
)

for (i in seq_len(nrow(main_results))) {
  r <- main_results[i, ]
  cat(sprintf("  [%s] %s: coef = %.4f  SE = %.4f  p = %s  [%s]\n",
              r$sample, r$outcome, r$coef, r$se, fmt_p(r$p), r$sig_status))
}
cat("\n")

## Also run with yield_shock (residual) as regressor for reference
shock_results <- bind_rows(
  lapply(c("share_attached", "share_casual"), function(y) {
    run_ext_reg(dsy_panel, y, "yield_shock", "Pooled", FE_POOLED)
  }),
  run_ext_reg(dsy_ratio, "attached_to_casual_ratio", "yield_shock", "Pooled", FE_POOLED)
)

write_csv(main_results,  file.path(OUT, "tables", "extensive_margin_main.csv"))
write_csv(shock_results, file.path(OUT, "tables", "extensive_margin_yield_shock.csv"))

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 4 — Asymmetric extensive margin (positive vs negative shocks)         ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 4 — Asymmetric extensive margin (split by shock sign)\n")

run_asym_ext <- function(df, outcome, shock_dir, sample_label, fe_rhs) {
  shock_val <- if (shock_dir == "positive") 1L else 0L
  sub <- df %>% filter(positive_shock == shock_val)
  run_ext_reg(sub, outcome, REG, paste(sample_label, shock_dir, sep = "_"), fe_rhs)
}

asym_results <- bind_rows(
  lapply(c("positive", "negative"), function(sh) {
    bind_rows(
      run_asym_ext(dsy_panel, "share_attached", sh, "Pooled", FE_POOLED),
      run_asym_ext(dsy_panel, "share_casual",   sh, "Pooled", FE_POOLED),
      run_asym_ext(dsy_ratio, "attached_to_casual_ratio", sh, "Pooled", FE_POOLED)
    )
  })
)

for (i in seq_len(nrow(asym_results))) {
  r <- asym_results[i, ]
  cat(sprintf("  [%s] %s: coef = %.4f  p = %s  [%s]\n",
              r$sample, r$outcome, r$coef, fmt_p(r$p), r$sig_status))
}

neg_att <- asym_results %>% filter(grepl("negative", sample), outcome == "share_attached")
pos_att <- asym_results %>% filter(grepl("positive", sample), outcome == "share_attached")
if (neg_att$coef < 0 && neg_att$p < 0.10 && abs(pos_att$coef) < abs(neg_att$coef)) {
  cat("\n  >> Hysteresis pattern: share_attached falls in bad years, not restored in good years.\n")
} else if (neg_att$p >= 0.10 && pos_att$p >= 0.10) {
  cat("\n  >> No significant composition shift in either shock direction — durable relationships.\n")
}
cat("\n")

write_csv(asym_results, file.path(OUT, "tables", "extensive_margin_asymmetry.csv"))

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 5 — Interaction with shock severity (quartiles)                     ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 5 — Shock severity quartiles (negative shocks only)\n")

severity_results <- bind_rows(
  lapply(SEASONS, function(s) {
    bind_rows(
      run_ext_reg(
        dsy_panel %>% filter(growing_season == s, shock_severity == "mild_negative"),
        "share_attached", REG, paste(s, "mild_neg"), FE_SEASON
      ),
      run_ext_reg(
        dsy_panel %>% filter(growing_season == s, shock_severity == "severe_negative"),
        "share_attached", REG, paste(s, "severe_neg"), FE_SEASON
      )
    )
  }),
  run_ext_reg(
    dsy_panel %>% filter(shock_severity == "mild_negative"),
    "share_attached", REG, "Pooled_mild_neg", FE_POOLED
  ),
  run_ext_reg(
    dsy_panel %>% filter(shock_severity == "severe_negative"),
    "share_attached", REG, "Pooled_severe_neg", FE_POOLED
  )
)

for (i in seq_len(nrow(severity_results))) {
  r <- severity_results[i, ]
  cat(sprintf("  [%s]: coef = %.4f  p = %s  [%s]  N = %d\n",
              r$sample, r$coef, fmt_p(r$p), r$sig_status, r$nobs))
}
cat("\n")

write_csv(severity_results, file.path(OUT, "tables", "extensive_margin_severity.csv"))

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 6 — Selection bias implication for wage estimates                     ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 6 — Selection bias check (mean 3-meal wage)\n")

df_sel <- dsy_panel %>% filter(!is.na(mean_three_wage))

m_wage_base <- feols(
  mean_three_wage ~ diff_log_yield_hat | year + District^growing_season,
  data = df_sel, cluster = ~District
)
m_wage_ctrl <- feols(
  mean_three_wage ~ diff_log_yield_hat + share_attached | year + District^growing_season,
  data = df_sel, cluster = ~District
)

sel_compare <- data.frame(
  specification = c("Without share_attached", "With share_attached"),
  yield_coef    = c(coef(m_wage_base)["diff_log_yield_hat"],
                    coef(m_wage_ctrl)["diff_log_yield_hat"]),
  yield_se      = c(se(m_wage_base)["diff_log_yield_hat"],
                    se(m_wage_ctrl)["diff_log_yield_hat"]),
  yield_p       = c(pvalue(m_wage_base)["diff_log_yield_hat"],
                    pvalue(m_wage_ctrl)["diff_log_yield_hat"]),
  share_coef    = c(NA, coef(m_wage_ctrl)["share_attached"]),
  share_se      = c(NA, se(m_wage_ctrl)["share_attached"]),
  share_p       = c(NA, pvalue(m_wage_ctrl)["share_attached"]),
  nobs          = c(nobs(m_wage_base), nobs(m_wage_ctrl))
)

pct_change <- 100 * (sel_compare$yield_coef[2] - sel_compare$yield_coef[1]) /
  abs(sel_compare$yield_coef[1])

cat(sprintf("  Yield coef without share_attached: %.1f (SE %.1f, p = %s)\n",
            sel_compare$yield_coef[1], sel_compare$yield_se[1],
            fmt_p(sel_compare$yield_p[1])))
cat(sprintf("  Yield coef with share_attached:    %.1f (SE %.1f, p = %s)\n",
            sel_compare$yield_coef[2], sel_compare$yield_se[2],
            fmt_p(sel_compare$yield_p[2])))
cat(sprintf("  Change in yield coef: %.1f%%\n", pct_change))
cat(sprintf("  share_attached coef: %.1f (SE %.1f, p = %s)\n",
            sel_compare$share_coef[2], sel_compare$share_se[2],
            fmt_p(sel_compare$share_p[2])))

if (abs(pct_change) > 20) {
  cat("  >> Substantial change (>20%) — evidence of selection bias in wage estimates.\n")
} else {
  cat("  >> Modest change — limited evidence of selection bias from composition shifts.\n")
}
cat("\n")

write_csv(sel_compare, file.path(OUT, "tables", "selection_bias_wage_check.csv"))

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 7 — Visualizations                                                    ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 7 — Visualizations\n")

## 7a. Scatter: yield shock vs share_attached by season
p_scatter <- dsy_panel %>%
  ggplot(aes(x = yield_shock, y = share_attached)) +
  geom_point(alpha = 0.35, size = 1.5, color = "#0072B2") +
  geom_smooth(method = "lm", se = TRUE, color = "#D55E00", linewidth = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ growing_season, nrow = 1) +
  labs(
    title = "Extensive Margin: Yield Shock vs Share of Attached (3-Meal) Workers",
    subtitle = "District-season-year cells | Fitted OLS line per season",
    x = "First-stage yield shock (residual)",
    y = "share_attached (monthly record fraction)"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))

ggsave(file.path(OUT, "figures", "scatter_shock_vs_share_attached.png"),
       p_scatter, width = 11, height = 4.5, dpi = 300)

## 7b. Coefficient plot: extensive margin across seasons and shock directions
coef_plot_df <- bind_rows(
  main_results %>% filter(outcome == "share_attached") %>%
    mutate(shock_direction = "All shocks", type = "Main"),
  asym_results %>% filter(outcome == "share_attached") %>%
    mutate(
      shock_direction = ifelse(grepl("positive", sample), "Positive shock", "Negative shock"),
      sample = "Pooled",
      type = "Split sample"
    )
) %>%
  mutate(
    label = paste(sample, shock_direction, sep = " | "),
    label = factor(label, levels = rev(unique(label)))
  )

p_coef <- ggplot(coef_plot_df, aes(x = coef, y = label, color = type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.2, linewidth = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Main" = "#0072B2", "Split sample" = "#D55E00")) +
  labs(
    title = "Extensive Margin Coefficients: share_attached on Predicted Yield",
    subtitle = "95% CI | Outcome: share of 3-meal monthly records",
    x = "Coefficient on diff_log_yield_hat", y = NULL, color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(OUT, "figures", "coefplot_extensive_margin.png"),
       p_coef, width = 8, height = 4.5, dpi = 300)

## Season-specific main coefficients
coef_season <- main_results %>%
  filter(outcome == "share_attached", sample %in% SEASONS) %>%
  mutate(sample = factor(sample, levels = SEASONS))

p_season_coef <- ggplot(coef_season, aes(x = sample, y = coef, color = sample)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.15, linewidth = 0.8) +
  geom_point(size = 4) +
  geom_text(aes(label = sprintf("%.3f\n(p=%s)", coef, fmt_p(p))),
            vjust = -1.2, size = 3, color = "grey25") +
  scale_color_manual(values = c(Boro = "#0072B2", Aus = "#E69F00", Aman = "#CC79A7"),
                     guide = "none") +
  labs(
    title = "Extensive Margin by Season",
    subtitle = "share_attached ~ diff_log_yield_hat | year + District",
    x = "Season", y = "Coefficient"
  ) +
  theme_minimal(base_size = 12) +
  coord_cartesian(ylim = range(c(coef_season$ci_lo, coef_season$ci_hi, 0)) * c(1.3, 1.3))

ggsave(file.path(OUT, "figures", "coefplot_extensive_margin_by_season.png"),
       p_season_coef, width = 7, height = 5, dpi = 300)

## 7c. Time series of share_attached over 2017-2023
ts_df <- dsy_panel %>%
  group_by(year) %>%
  summarise(
    mean_share_attached = mean(share_attached, na.rm = TRUE),
    median_share_attached = median(share_attached, na.rm = TRUE),
    .groups = "drop"
  )

p_ts <- ggplot(ts_df, aes(x = year)) +
  geom_line(aes(y = mean_share_attached), color = "#0072B2", linewidth = 1) +
  geom_point(aes(y = mean_share_attached), color = "#0072B2", size = 2.5) +
  geom_line(aes(y = median_share_attached), color = "#D55E00", linewidth = 0.8, linetype = "dashed") +
  scale_x_continuous(breaks = 2017:2023) +
  labs(
    title = "Prevalence of Attached (3-Meal) Labor Over Time",
    subtitle = "Mean and median share_attached across district-season-year cells",
    x = "Year", y = "share_attached",
    caption = "Blue = mean; dashed orange = median"
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(OUT, "figures", "timeseries_share_attached.png"),
       p_ts, width = 8, height = 4.5, dpi = 300)

## District-level spaghetti (subset for readability)
set.seed(42)
sample_districts <- sample(unique(dsy_panel$District), min(15, length(unique(dsy_panel$District))))
p_spaghetti <- dsy_panel %>%
  filter(District %in% sample_districts) %>%
  ggplot(aes(x = year, y = share_attached, group = interaction(District, growing_season), color = growing_season)) +
  geom_line(alpha = 0.5, linewidth = 0.5) +
  geom_point(size = 0.8, alpha = 0.6) +
  scale_color_manual(values = c(Boro = "#0072B2", Aus = "#E69F00", Aman = "#CC79A7"),
                     name = "Season") +
  labs(
    title = "share_attached Over Time (15 Sample Districts)",
    subtitle = "Each line = one district-season cell",
    x = "Year", y = "share_attached"
  ) +
  theme_minimal(base_size = 11)

ggsave(file.path(OUT, "figures", "timeseries_share_attached_districts.png"),
       p_spaghetti, width = 9, height = 5, dpi = 300)

cat("  Saved scatter, coefficient plots, and time series figures\n\n")

## ── HTML summary table ──────────────────────────────────────────────────────── ##
main_html <- main_results %>%
  mutate(
    coef_se = sprintf("%.4f (%.4f)", coef, se),
    p_fmt   = fmt_p(p)
  ) %>%
  select(Sample = sample, Outcome = outcome, `Coef (SE)` = coef_se,
         `p-value` = p_fmt, Status = sig_status, N = nobs) %>%
  kbl(caption = "Extensive Margin Regressions", format = "html", booktabs = TRUE) %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)

writeLines(as.character(main_html), file.path(OUT, "tables", "extensive_margin_main.html"))

## ── Final report ───────────────────────────────────────────────────────────── ##
share_att_pooled <- main_results %>%
  filter(sample == "Pooled", outcome == "share_attached")

interpretation <- if (share_att_pooled$coef < 0 && share_att_pooled$p < 0.10) {
  "Significant extensive margin: share_attached falls when predicted yield falls. Wage estimates may understate true risk exposure due to worker dropout."
} else if (share_att_pooled$coef < 0 && share_att_pooled$p < 0.20) {
  "Directionally consistent (share falls in bad years) but imprecise. Wage channel may be primary adjustment mechanism."
} else {
  "No evidence of extensive margin adjustment — attached relationships appear durable. Wage channel is the primary (or sole) adjustment margin."
}

report <- c(
  structure_report,
  "",
  "---",
  "",
  "# Extensive Margin Analysis Results",
  sprintf("Generated: %s", Sys.time()),
  "",
  "## Step 2 — Summary statistics (share_attached)",
  capture.output(print(sum_stats_all)),
  capture.output(print(sum_stats_season)),
  "",
  "## Step 3 — Main regressions",
  capture.output(print(main_results)),
  "",
  "## Step 4 — Asymmetry",
  capture.output(print(asym_results)),
  "",
  "## Step 5 — Severity quartiles",
  capture.output(print(severity_results)),
  "",
  "## Step 6 — Selection bias",
  capture.output(print(sel_compare)),
  sprintf("- Percent change in yield coef when controlling for share_attached: %.1f%%", pct_change),
  "",
  "## Interpretation",
  interpretation,
  "",
  "## Outputs",
  "- summary/data_structure_audit.md",
  "- tables/extensive_margin_panel.csv",
  "- tables/extensive_margin_main.csv/.html",
  "- tables/extensive_margin_asymmetry.csv",
  "- tables/extensive_margin_severity.csv",
  "- tables/selection_bias_wage_check.csv",
  "- figures/scatter_shock_vs_share_attached.png",
  "- figures/coefplot_extensive_margin.png",
  "- figures/timeseries_share_attached.png"
)

writeLines(report, file.path(OUT, "summary", "extensive_margin_report.md"))

save(dsy_panel, main_results, asym_results, severity_results, sel_compare,
     m_wage_base, m_wage_ctrl,
     file = file.path(OUT, "models", "extensive_margin_models.RData"))

cat("Saved extensive_margin_report.md\n")
cat("=== EXTENSIVE MARGIN ANALYSIS COMPLETE ===\n")
