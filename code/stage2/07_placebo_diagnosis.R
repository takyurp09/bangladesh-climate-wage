## stage2_placebo_diagnosis.R
## Diagnose why ROB3 placebo (lagged yield) is significant (p=0.048)
## Five tests: wage residual AR(1), yield AR(1), placebo+controls,
##             Driscoll-Kraay SE, lagged-wage control
## Output: output/stage2/summary/PLACEBO_DIAGNOSIS.md

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(fixest)
})

ROOT    <- here::here()
out_sum <- file.path(ROOT, "output/stage2/summary")
dir.create(out_sum, recursive = TRUE, showWarnings = FALSE)

cat("=== PLACEBO DIAGNOSIS ===\n")

## ── Load & prepare ────────────────────────────────────────────────────────── ##
df <- read_csv(file.path(ROOT, "data/Regression_data/df_2_merged_v2.csv"),
               show_col_types = FALSE) %>%
  mutate(
    gender    = relevel(factor(gender),    ref = "Female"),
    meal_type = relevel(factor(meal_type), ref = "None")
  ) %>%
  arrange(District, growing_season, gender, meal_type, year)

## Build lagged variables within District × growing_season × gender × meal_type
df <- df %>%
  group_by(District, growing_season, gender, meal_type) %>%
  mutate(
    lag_diff_real_wage      = lag(diff_real_wage),
    lag_diff_log_yield_hat  = lag(diff_log_yield_hat)   # lagged instrument (placebo)
  ) %>%
  ungroup()

## Also need lagged yield at District×season level (for yield AR(1) test)
df_yield_level <- df %>%
  distinct(District, growing_season, year, diff_log_yield_hat) %>%
  arrange(District, growing_season, year) %>%
  group_by(District, growing_season) %>%
  mutate(lag_yield = lag(diff_log_yield_hat)) %>%
  ungroup()

cat(sprintf("Full panel: N=%d\n", nrow(df)))

## ─────────────────────────────────────────────────────────────────────────── ##
## TEST 1: AR(1) in wage residuals                                             ##
## ─────────────────────────────────────────────────────────────────────────── ##
cat("\n--- TEST 1: AR(1) in wage residuals ---\n")

m3_main <- feols(diff_real_wage ~ diff_log_yield_hat + gender + meal_type |
                   year + District^growing_season,
                 data = df, cluster = ~District)

## Attach residuals
df_resid <- df %>%
  filter(!is.na(diff_real_wage) & !is.na(diff_log_yield_hat)) %>%
  mutate(resid_m3 = residuals(m3_main))

## Lag residuals within panel group
df_resid <- df_resid %>%
  arrange(District, growing_season, gender, meal_type, year) %>%
  group_by(District, growing_season, gender, meal_type) %>%
  mutate(lag_resid_m3 = lag(resid_m3)) %>%
  ungroup()

ar1_wage <- feols(resid_m3 ~ lag_resid_m3 |
                    District^growing_season + gender + meal_type,
                  data = df_resid, cluster = ~District)

t1_coef <- coef(ar1_wage)["lag_resid_m3"]
t1_se   <- se(ar1_wage)["lag_resid_m3"]
t1_p    <- pvalue(ar1_wage)["lag_resid_m3"]

cat(sprintf("AR(1) wage residuals: coef=%.4f SE=%.4f p=%.4f\n", t1_coef, t1_se, t1_p))
cat(sprintf("→ %s\n",
    if (!is.na(t1_p) && t1_p < 0.05)
      "SIGNIFICANT — serial correlation in wage residuals present"
    else
      "Not significant at 5%"))

## ─────────────────────────────────────────────────────────────────────────── ##
## TEST 2: AR(1) in diff_log_yield_hat                                         ##
## ─────────────────────────────────────────────────────────────────────────── ##
cat("\n--- TEST 2: AR(1) in yield_hat ---\n")

df_yield_ar <- df_yield_level %>%
  filter(!is.na(diff_log_yield_hat) & !is.na(lag_yield))

ar1_yield <- feols(diff_log_yield_hat ~ lag_yield |
                     District^growing_season,
                   data = df_yield_ar, cluster = ~District)

t2_coef <- coef(ar1_yield)["lag_yield"]
t2_se   <- se(ar1_yield)["lag_yield"]
t2_p    <- pvalue(ar1_yield)["lag_yield"]

cat(sprintf("AR(1) yield_hat: coef=%.4f SE=%.4f p=%.4f\n", t2_coef, t2_se, t2_p))
cat(sprintf("→ %s\n",
    if (!is.na(t2_p) && t2_p < 0.05)
      "SIGNIFICANT — yield instrument is autocorrelated"
    else
      "Not significant at 5%"))

## ─────────────────────────────────────────────────────────────────────────── ##
## TEST 3: Placebo (lagged yield) with and without controls                    ##
## ─────────────────────────────────────────────────────────────────────────── ##
cat("\n--- TEST 3: Placebo p with/without controls ---\n")

df_lag <- df %>% filter(!is.na(lag_diff_log_yield_hat))

## ROB3 baseline (M3 spec, lagged yield)
rob3_base <- feols(diff_real_wage ~ lag_diff_log_yield_hat + gender + meal_type |
                     year + District^growing_season,
                   data = df_lag, cluster = ~District)

t3a_coef <- coef(rob3_base)["lag_diff_log_yield_hat"]
t3a_se   <- se(rob3_base)["lag_diff_log_yield_hat"]
t3a_p    <- pvalue(rob3_base)["lag_diff_log_yield_hat"]

cat(sprintf("ROB3 (M3, no controls):  coef=%.4f SE=%.4f p=%.4f  N=%d\n",
            t3a_coef, t3a_se, t3a_p, nobs(rob3_base)))

## ROB3 + M3_extended controls
df_lag_ext <- df_lag %>%
  filter(!is.na(share_Boro) & !is.na(share_Aus) &
           !is.na(pop_density) & !is.na(ratio_double_cropped))

rob3_ext <- feols(diff_real_wage ~ lag_diff_log_yield_hat + gender + meal_type +
                    share_Boro + share_Aus + pop_density + ratio_double_cropped |
                    year + District^growing_season,
                  data = df_lag_ext, cluster = ~District)

t3b_coef <- coef(rob3_ext)["lag_diff_log_yield_hat"]
t3b_se   <- se(rob3_ext)["lag_diff_log_yield_hat"]
t3b_p    <- pvalue(rob3_ext)["lag_diff_log_yield_hat"]

cat(sprintf("ROB3 (M3_extended):      coef=%.4f SE=%.4f p=%.4f  N=%d\n",
            t3b_coef, t3b_se, t3b_p, nobs(rob3_ext)))
cat(sprintf("→ Placebo %s with controls\n",
    if (!is.na(t3b_p) && t3b_p >= 0.05) "LOSES significance" else "remains significant"))

## ─────────────────────────────────────────────────────────────────────────── ##
## TEST 4: Driscoll-Kraay (spatial HAC) SE on M3                              ##
## fixest ≥0.11: summary(model, vcov = "DK") or vcov = "NW"                  ##
## ─────────────────────────────────────────────────────────────────────────── ##
cat("\n--- TEST 4: Driscoll-Kraay HAC SE on M3 ---\n")

# Set panel id for DK SE
m3_dk <- tryCatch({
  setFixest_estimation(panel.id = ~District + year)
  feols(diff_real_wage ~ diff_log_yield_hat + gender + meal_type |
          year + District^growing_season,
        data = df, vcov = "DK")
}, error = function(e) {
  cat("DK with panel.id failed, trying NW...\n")
  tryCatch(
    feols(diff_real_wage ~ diff_log_yield_hat + gender + meal_type |
            year + District^growing_season,
          data = df, vcov = "NW"),
    error = function(e2) {
      cat("NW also failed:", conditionMessage(e2), "\n"); NULL
    }
  )
})

if (!is.null(m3_dk)) {
  t4_coef <- coef(m3_dk)["diff_log_yield_hat"]
  t4_se   <- se(m3_dk)["diff_log_yield_hat"]
  t4_p    <- pvalue(m3_dk)["diff_log_yield_hat"]
  t4_label <- "Driscoll-Kraay/NW"
} else {
  ## Fallback: standard errors clustered two-way (District + year)
  cat("Falling back to two-way cluster SE (District + year)...\n")
  m3_2w  <- feols(diff_real_wage ~ diff_log_yield_hat + gender + meal_type |
                    year + District^growing_season,
                  data = df, cluster = ~District + year)
  t4_coef <- coef(m3_2w)["diff_log_yield_hat"]
  t4_se   <- se(m3_2w)["diff_log_yield_hat"]
  t4_p    <- pvalue(m3_2w)["diff_log_yield_hat"]
  t4_label <- "Two-way cluster (District+year)"
}

cat(sprintf("M3 clustered (District):     coef=%.4f SE=%.4f p=%.4f\n",
            coef(m3_main)["diff_log_yield_hat"],
            se(m3_main)["diff_log_yield_hat"],
            pvalue(m3_main)["diff_log_yield_hat"]))
cat(sprintf("M3 %s: coef=%.4f SE=%.4f p=%.4f\n",
            t4_label, t4_coef, t4_se, t4_p))
cat(sprintf("→ SE change: %.1f → %.1f (%+.1f%%)\n",
            se(m3_main)["diff_log_yield_hat"], t4_se,
            (t4_se / se(m3_main)["diff_log_yield_hat"] - 1) * 100))

## ─────────────────────────────────────────────────────────────────────────── ##
## TEST 5: M3 + lagged wage as control; placebo with lagged wage               ##
## ─────────────────────────────────────────────────────────────────────────── ##
cat("\n--- TEST 5: Add lag(diff_real_wage) as control ---\n")

df_lagw <- df %>% filter(!is.na(lag_diff_real_wage))

## M3 + lagged wage
m3_lagw <- feols(diff_real_wage ~ diff_log_yield_hat + gender + meal_type +
                   lag_diff_real_wage |
                   year + District^growing_season,
                 data = df_lagw, cluster = ~District)

t5_coef <- coef(m3_lagw)["diff_log_yield_hat"]
t5_se   <- se(m3_lagw)["diff_log_yield_hat"]
t5_p    <- pvalue(m3_lagw)["diff_log_yield_hat"]
t5_lagw_coef <- coef(m3_lagw)["lag_diff_real_wage"]
t5_lagw_p    <- pvalue(m3_lagw)["lag_diff_real_wage"]

cat(sprintf("M3 + lag wage: yield coef=%.4f SE=%.4f p=%.4f  N=%d\n",
            t5_coef, t5_se, t5_p, nobs(m3_lagw)))
cat(sprintf("  lag_diff_real_wage: coef=%.4f p=%.4f → %s\n",
            t5_lagw_coef, t5_lagw_p,
            if (t5_lagw_p < 0.05) "SIGNIFICANT (serial corr. confirmed)" else "Not significant"))

## Placebo with lagged wage control
df_lagw_lag <- df_lagw %>% filter(!is.na(lag_diff_log_yield_hat))

rob3_lagw <- feols(diff_real_wage ~ lag_diff_log_yield_hat + gender + meal_type +
                     lag_diff_real_wage |
                     year + District^growing_season,
                   data = df_lagw_lag, cluster = ~District)

t5b_coef <- coef(rob3_lagw)["lag_diff_log_yield_hat"]
t5b_se   <- se(rob3_lagw)["lag_diff_log_yield_hat"]
t5b_p    <- pvalue(rob3_lagw)["lag_diff_log_yield_hat"]

cat(sprintf("Placebo + lag wage: coef=%.4f SE=%.4f p=%.4f  N=%d\n",
            t5b_coef, t5b_se, t5b_p, nobs(rob3_lagw)))
cat(sprintf("→ Placebo %s after adding lagged wage control\n",
    if (!is.na(t5b_p) && t5b_p >= 0.05) "LOSES significance" else "remains significant"))

## ─────────────────────────────────────────────────────────────────────────── ##
## Diagnosis summary                                                            ##
## ─────────────────────────────────────────────────────────────────────────── ##
cat("\n=== DIAGNOSIS SUMMARY ===\n")

diagnose_serial <- !is.na(t1_p) && t1_p < 0.05
diagnose_yield_ar <- !is.na(t2_p) && t2_p < 0.05
placebo_survives_controls <- !is.na(t3b_p) && t3b_p < 0.05
placebo_survives_lagw     <- !is.na(t5b_p) && t5b_p < 0.05

if (diagnose_serial && diagnose_yield_ar) {
  diagnosis <- "BOTH wage residuals and yield_hat are serially correlated. The placebo is likely spurious due to temporal autocorrelation propagating through the FD regression. Recommended fix: include lagged wage as control OR use Driscoll-Kraay SE."
} else if (diagnose_serial) {
  diagnosis <- "Wage residuals are serially correlated but yield_hat is not. The placebo captures persistence in wages not absorbed by FD. Recommended fix: include lagged wage as control."
} else if (diagnose_yield_ar) {
  diagnosis <- "Yield_hat is autocorrelated but wage residuals are not. The 'placebo' lagged yield is correlated with current yield changes (ARMA structure). This is a misspecification of the placebo test, not an IV validity problem. Recommended fix: report this finding and note the placebo test is biased."
} else {
  diagnosis <- "No strong serial correlation detected. Placebo significance may be due to small-sample bias (G=63 clusters) or insufficient within-variation. Recommended fix: verify cluster bootstrap p-value on placebo."
}

cat(paste0(strwrap(diagnosis, 80), collapse="\n"), "\n")

## ─────────────────────────────────────────────────────────────────────────── ##
## Save report                                                                  ##
## ─────────────────────────────────────────────────────────────────────────── ##
report <- c(
  "# Placebo Diagnosis Report",
  sprintf("Generated: %s", Sys.time()),
  "",
  "## Context",
  "ROB3 in prior robustness run: lagged yield_hat as instrument → p=0.048 (significant).",
  "This script diagnoses whether the significance is due to serial correlation,",
  "autocorrelated yield, or another structural issue.",
  "",
  "---",
  "",
  "## TEST 1 — Serial correlation in M3 wage residuals",
  sprintf("Regress resid_t on resid_(t-1) within District×season×gender×meal_type, FE: District×season + gender + meal_type, cluster=~District"),
  sprintf("| Statistic | Value |"),
  sprintf("|-----------|-------|"),
  sprintf("| AR(1) coefficient | %.4f |", t1_coef),
  sprintf("| SE (clustered) | %.4f |", t1_se),
  sprintf("| p-value | %.4f |", t1_p),
  sprintf("| N obs | %d |", nobs(ar1_wage)),
  sprintf("| **Verdict** | **%s** |",
          if (diagnose_serial) "Significant — serial correlation present" else "Not significant"),
  "",
  "---",
  "",
  "## TEST 2 — Serial correlation in ΔlogYieldHat",
  sprintf("Regress Δyield_t on Δyield_(t-1) within District×season, FE: District×season, cluster=~District"),
  sprintf("| Statistic | Value |"),
  sprintf("|-----------|-------|"),
  sprintf("| AR(1) coefficient | %.4f |", t2_coef),
  sprintf("| SE (clustered) | %.4f |", t2_se),
  sprintf("| p-value | %.4f |", t2_p),
  sprintf("| N obs | %d |", nobs(ar1_yield)),
  sprintf("| **Verdict** | **%s** |",
          if (diagnose_yield_ar) "Significant — yield instrument autocorrelated" else "Not significant"),
  "",
  "---",
  "",
  "## TEST 3 — Placebo p-value with/without controls",
  "Spec: diff_real_wage ~ lag(diff_log_yield_hat) + gender + meal_type [+ controls]",
  sprintf("| Spec | Coef | SE | p-value | N |"),
  sprintf("|------|------|----|---------|---|"),
  sprintf("| M3 (no controls) | %.4f | %.4f | %.4f | %d |",
          t3a_coef, t3a_se, t3a_p, nobs(rob3_base)),
  sprintf("| M3_extended (+controls) | %.4f | %.4f | %.4f | %d |",
          t3b_coef, t3b_se, t3b_p, nobs(rob3_ext)),
  sprintf("| **Verdict** | **Placebo %s with controls** | | | |",
          if (!placebo_survives_controls) "LOSES significance" else "remains significant"),
  "",
  "---",
  "",
  sprintf("## TEST 4 — %s SE on M3", t4_label),
  sprintf("| SE type | Coef | SE | p-value |"),
  sprintf("|---------|------|----|---------|"),
  sprintf("| Clustered by District | %.4f | %.4f | %.4f |",
          coef(m3_main)["diff_log_yield_hat"],
          se(m3_main)["diff_log_yield_hat"],
          pvalue(m3_main)["diff_log_yield_hat"]),
  sprintf("| %s | %.4f | %.4f | %.4f |", t4_label, t4_coef, t4_se, t4_p),
  sprintf("| SE inflation | | %+.1f%% | |",
          (t4_se / se(m3_main)["diff_log_yield_hat"] - 1) * 100),
  sprintf("| **Verdict** | Inference %s | | |",
          if (!is.na(t4_p) && t4_p < 0.05) "CHANGES (now significant)" else "unchanged (still insignificant)"),
  "",
  "---",
  "",
  "## TEST 5 — Add lag(diff_real_wage) as control",
  sprintf("| Spec | Yield coef | SE | p-value | Lag wage p | N |"),
  sprintf("|------|------------|----|---------|-----------|----|"),
  sprintf("| M3 (baseline) | %.4f | %.4f | %.4f | — | %d |",
          coef(m3_main)["diff_log_yield_hat"],
          se(m3_main)["diff_log_yield_hat"],
          pvalue(m3_main)["diff_log_yield_hat"],
          nobs(m3_main)),
  sprintf("| M3 + lag wage | %.4f | %.4f | %.4f | %.4f | %d |",
          t5_coef, t5_se, t5_p, t5_lagw_p, nobs(m3_lagw)),
  sprintf("| Placebo + lag wage | %.4f | %.4f | %.4f | — | %d |",
          t5b_coef, t5b_se, t5b_p, nobs(rob3_lagw)),
  sprintf("| **Verdict** | Placebo %s after adding lag wage | | | | |",
          if (!placebo_survives_lagw) "LOSES significance" else "remains significant"),
  "",
  "---",
  "",
  "## Overall Diagnosis",
  "",
  paste0(strwrap(diagnosis, 100), collapse = "\n"),
  "",
  "## Recommended Fix",
  "",
  if (!placebo_survives_lagw) {
    c("**Add lag(diff_real_wage) as control** to main spec.",
      "- Absorbs wage persistence that the FD transformation does not fully remove.",
      "- Placebo becomes insignificant once this control is included.",
      "- This is the minimal fix with the largest diagnostic payoff.")
  } else if (diagnose_yield_ar) {
    c("**Clarify placebo test interpretation.**",
      "- The placebo test is invalid because ΔlogYieldHat is itself autocorrelated.",
      "- A 'lagged' instrument then predicts current wages mechanically via yield AR structure.",
      "- This is a test design flaw, not an IV validity failure.",
      "- Recommended: report that placebo test is biased due to yield autocorrelation;",
      "  use Driscoll-Kraay SE as the primary robustness check instead.")
  } else {
    c("**Use Driscoll-Kraay SE** as a robustness check.",
      "- The covariance structure is more complex than district clustering alone.",
      "  A spatial-HAC approach is more conservative and appropriate.")
  }
)

writeLines(report, file.path(out_sum, "PLACEBO_DIAGNOSIS.md"))
cat("\nSaved PLACEBO_DIAGNOSIS.md\n")
cat("=== PLACEBO DIAGNOSIS COMPLETE ===\n")
