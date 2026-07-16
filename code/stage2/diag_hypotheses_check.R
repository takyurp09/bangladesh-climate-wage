## diag_hypotheses_check.R
## Diagnostic only — do NOT modify any paper files.
## Verifies M5 total-effect SEs and p-values via delta method on VCV.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(fixest)
})

ROOT <- here::here()
cat("Root:", ROOT, "\n")

## ── Load data ────────────────────────────────────────────────────────────── ##
df <- read_csv(file.path(ROOT, "data/Regression_data/df_2_merged_v2.csv"),
               show_col_types = FALSE) %>%
  mutate(
    gender    = relevel(factor(gender),    ref = "Female"),
    meal_type = relevel(factor(meal_type), ref = "None")
  )
cat(sprintf("N = %d | districts = %d\n", nrow(df), length(unique(df$District))))

## ── Delta-method helper ───────────────────────────────────────────────────── ##
# Computes total effect = b1 + b2 (e.g., None baseline + interaction)
# SE via delta method: sqrt(var(b1) + var(b2) + 2*cov(b1,b2))
# p-value: two-sided t with G-1 = 62 df (63 districts)
total_effect <- function(model, nm1, nm2, df_resid = 62) {
  b  <- coef(model)
  V  <- vcov(model)
  est <- b[nm1] + b[nm2]
  var_e <- V[nm1, nm1] + V[nm2, nm2] + 2 * V[nm1, nm2]
  se_e  <- sqrt(var_e)
  t_e   <- est / se_e
  p_e   <- 2 * pt(-abs(t_e), df = df_resid)
  list(coef = est, se = se_e, t = t_e, p = p_e)
}

## ── Fit M5 ───────────────────────────────────────────────────────────────── ##
m5 <- feols(diff_real_wage ~ diff_log_yield_hat +
              diff_log_yield_hat:meal_type + gender |
              year + District^growing_season,
            data = df, cluster = ~District)

cat("\n--- M5 raw coeftable ---\n")
print(summary(m5)$coeftable)
cat(sprintf("\nNumber of clusters (G): %d\n", length(unique(df$District))))
cat(sprintf("Residual df used for t-test (G-1): %d\n", length(unique(df$District)) - 1))

## ── Three-meal total effect ───────────────────────────────────────────────── ##
cat("\n=== THREE-MEAL TOTAL EFFECT (delta method) ===\n")
r_three <- total_effect(m5, "diff_log_yield_hat", "diff_log_yield_hat:meal_typeThree")
cat(sprintf("  Coef (None + Three_interaction): %.4f\n", r_three$coef))
cat(sprintf("  SE (delta method):               %.4f\n", r_three$se))
cat(sprintf("  t-stat:                          %.4f\n", r_three$t))
cat(sprintf("  p-value (two-sided, df=62):      %.4f\n", r_three$p))

## ── One-meal total effect ─────────────────────────────────────────────────── ##
cat("\n=== ONE-MEAL TOTAL EFFECT (delta method) ===\n")
r_one <- total_effect(m5, "diff_log_yield_hat", "diff_log_yield_hat:meal_typeOne")
cat(sprintf("  Coef: %.4f | SE: %.4f | t: %.4f | p: %.4f\n",
            r_one$coef, r_one$se, r_one$t, r_one$p))

## ── Female-only M5 ───────────────────────────────────────────────────────── ##
cat("\n=== FEMALE-ONLY M5 ===\n")
df_f <- df %>% filter(gender == "Female")
m5_f <- feols(diff_real_wage ~ diff_log_yield_hat +
                diff_log_yield_hat:meal_type |
                year + District^growing_season,
              data = df_f, cluster = ~District)
G_f <- length(unique(df_f$District))
cat(sprintf("Female N=%d, G=%d\n", nrow(df_f), G_f))
cat("Female M5 coeftable:\n")
print(summary(m5_f)$coeftable)
r_f <- total_effect(m5_f, "diff_log_yield_hat", "diff_log_yield_hat:meal_typeThree",
                    df_resid = G_f - 1)
cat(sprintf("  Total: coef=%.4f | SE=%.4f | t=%.4f | p=%.4f\n",
            r_f$coef, r_f$se, r_f$t, r_f$p))

## ── Male-only M5 ─────────────────────────────────────────────────────────── ##
cat("\n=== MALE-ONLY M5 ===\n")
df_m <- df %>% filter(gender == "Male")
m5_m <- feols(diff_real_wage ~ diff_log_yield_hat +
                diff_log_yield_hat:meal_type |
                year + District^growing_season,
              data = df_m, cluster = ~District)
G_m <- length(unique(df_m$District))
cat(sprintf("Male N=%d, G=%d\n", nrow(df_m), G_m))
cat("Male M5 coeftable:\n")
print(summary(m5_m)$coeftable)
r_m <- total_effect(m5_m, "diff_log_yield_hat", "diff_log_yield_hat:meal_typeThree",
                    df_resid = G_m - 1)
cat(sprintf("  Total: coef=%.4f | SE=%.4f | t=%.4f | p=%.4f\n",
            r_m$coef, r_m$se, r_m$t, r_m$p))

## ── DIAGNOSTIC REPORT ────────────────────────────────────────────────────── ##
cat("\n\nDIAGNOSTIC REPORT — M5 Total Effect Inference\n")
cat("==============================================\n\n")

cat("Three-meal:\n")
cat(sprintf("  Currently reported: coef=-103.7, SE=95.0, p=0.022\n"))
cat(sprintf("  delta-method result: coef=%.4f, SE=%.4f, t=%.4f, p=%.4f\n",
            r_three$coef, r_three$se, r_three$t, r_three$p))
three_match <- abs(r_three$p - 0.022) < 0.005
cat(sprintf("  MATCH: %s\n\n", ifelse(three_match, "YES", "NO")))

cat("One-meal:\n")
cat(sprintf("  Currently reported: coef=-55.5, SE=102.4, p=0.094\n"))
cat(sprintf("  delta-method result: coef=%.4f, SE=%.4f, t=%.4f, p=%.4f\n",
            r_one$coef, r_one$se, r_one$t, r_one$p))
one_match <- abs(r_one$p - 0.094) < 0.010
cat(sprintf("  MATCH: %s\n\n", ifelse(one_match, "YES", "NO")))

cat("Female three-meal:\n")
cat(sprintf("  Currently reported: coef=-682, SE=145, p=0.020\n"))
cat(sprintf("  delta-method result: coef=%.4f, SE=%.4f, t=%.4f, p=%.4f\n",
            r_f$coef, r_f$se, r_f$t, r_f$p))
f_match <- abs(r_f$p - 0.020) < 0.010
cat(sprintf("  MATCH: %s\n\n", ifelse(f_match, "YES", "NO")))

cat("Male three-meal:\n")
cat(sprintf("  Currently reported: coef=+314, SE=89, p=0.122\n"))
cat(sprintf("  delta-method result: coef=%.4f, SE=%.4f, t=%.4f, p=%.4f\n",
            r_m$coef, r_m$se, r_m$t, r_m$p))
m_match <- abs(r_m$p - 0.122) < 0.015
cat(sprintf("  MATCH: %s\n\n", ifelse(m_match, "YES", "NO")))

all_match <- three_match && one_match && f_match && m_match
cat("VERDICT:\n")
if (all_match) {
  cat("  Inference is correct. p-values are valid.\n")
} else {
  flags <- c()
  if (!three_match) flags <- c(flags, "Three-meal")
  if (!one_match)   flags <- c(flags, "One-meal")
  if (!f_match)     flags <- c(flags, "Female three-meal")
  if (!m_match)     flags <- c(flags, "Male three-meal")
  cat(sprintf("  INFERENCE ERROR DETECTED. Paper must be corrected before submission.\n"))
  cat(sprintf("  Affected results: %s\n", paste(flags, collapse = ", ")))
}

## ── SE convention note ────────────────────────────────────────────────────── ##
ct <- summary(m5)$coeftable
cat("\n=== SE CONVENTION NOTE ===\n")
cat(sprintf("Interaction term (Three):  coef=%.4f, SE=%.4f, p=%.4f  [currently used as 'total SE' in paper]\n",
            ct["diff_log_yield_hat:meal_typeThree","Estimate"],
            ct["diff_log_yield_hat:meal_typeThree","Std. Error"],
            ct["diff_log_yield_hat:meal_typeThree","Pr(>|t|)"]))
cat(sprintf("Baseline (None):           coef=%.4f, SE=%.4f\n",
            ct["diff_log_yield_hat","Estimate"],
            ct["diff_log_yield_hat","Std. Error"]))
cat(sprintf("Cov(None, Three_int):      %.4f\n",
            vcov(m5)["diff_log_yield_hat","diff_log_yield_hat:meal_typeThree"]))
cat(sprintf("Naive SE (no cov term):    %.4f\n",
            sqrt(ct["diff_log_yield_hat:meal_typeThree","Std. Error"]^2 +
                   ct["diff_log_yield_hat","Std. Error"]^2)))
cat(sprintf("Delta-method SE (with cov): %.4f\n", r_three$se))

cat("\n=== DONE ===\n")


cat("\n--- M5 raw coeftable ---\n")
print(summary(m5)$coeftable)

## ── Check coefficient names ───────────────────────────────────────────────── ##
cat("\nCoefficient names in M5:\n")
print(names(coef(m5)))

## ── Three-meal total effect via hypotheses() ─────────────────────────────── ##
cat("\n=== THREE-MEAL TOTAL EFFECT (hypotheses) ===\n")
h_three <- hypotheses(m5,
  "diff_log_yield_hat + diff_log_yield_hat:meal_typeThree = 0")
print(h_three)

## ── One-meal total effect via hypotheses() ───────────────────────────────── ##
cat("\n=== ONE-MEAL TOTAL EFFECT (hypotheses) ===\n")
h_one <- hypotheses(m5,
  "diff_log_yield_hat + diff_log_yield_hat:meal_typeOne = 0")
print(h_one)

## ── Female-only M5 ───────────────────────────────────────────────────────── ##
cat("\n=== FEMALE-ONLY M5 ===\n")
df_f <- df %>% filter(gender == "Female")
m5_f <- feols(diff_real_wage ~ diff_log_yield_hat +
                diff_log_yield_hat:meal_type |
                year + District^growing_season,
              data = df_f, cluster = ~District)
cat("Female M5 coeftable:\n")
print(summary(m5_f)$coeftable)
h_f <- hypotheses(m5_f,
  "diff_log_yield_hat + diff_log_yield_hat:meal_typeThree = 0")
print(h_f)

## ── Male-only M5 ─────────────────────────────────────────────────────────── ##
cat("\n=== MALE-ONLY M5 ===\n")
df_m <- df %>% filter(gender == "Male")
m5_m <- feols(diff_real_wage ~ diff_log_yield_hat +
                diff_log_yield_hat:meal_type |
                year + District^growing_season,
              data = df_m, cluster = ~District)
cat("Male M5 coeftable:\n")
print(summary(m5_m)$coeftable)
h_m <- hypotheses(m5_m,
  "diff_log_yield_hat + diff_log_yield_hat:meal_typeThree = 0")
print(h_m)

## ── Extract scalars helper ────────────────────────────────────────────────── ##
# hypotheses() returns a data.frame; pull coef, se, tstat, pval
get_hyp <- function(h) {
  # column names vary slightly by fixest version
  nm <- colnames(h)
  coef_c <- h[[nm[grepl("Estim|coef", nm, ignore.case=TRUE)][1]]]
  se_c   <- h[[nm[grepl("Std|se", nm, ignore.case=TRUE)][1]]]
  t_c    <- h[[nm[grepl("^t$|t.stat|statistic", nm, ignore.case=TRUE)][1]]]
  p_c    <- h[[nm[grepl("Pr|p.val|pval", nm, ignore.case=TRUE)][1]]]
  list(coef=coef_c, se=se_c, t=t_c, p=p_c)
}

r_three <- get_hyp(h_three)
r_one   <- get_hyp(h_one)
r_f     <- get_hyp(h_f)
r_m     <- get_hyp(h_m)

## ── DIAGNOSTIC REPORT ────────────────────────────────────────────────────── ##
cat("\n\nDIAGNOSTIC REPORT — M5 Total Effect Inference\n")
cat("==============================================\n\n")

cat("Three-meal:\n")
cat(sprintf("  Currently reported: coef=-103.7, SE=95.0, p=0.022\n"))
cat(sprintf("  hypotheses() result: coef=%.4f, SE=%.4f, t=%.4f, p=%.4f\n",
            r_three$coef, r_three$se, r_three$t, r_three$p))
three_match <- abs(r_three$p - 0.022) < 0.005
cat(sprintf("  MATCH: %s\n\n", ifelse(three_match, "YES", "NO")))

cat("One-meal:\n")
cat(sprintf("  Currently reported: coef=-55.5, SE=102.4, p=0.094\n"))
cat(sprintf("  hypotheses() result: coef=%.4f, SE=%.4f, t=%.4f, p=%.4f\n",
            r_one$coef, r_one$se, r_one$t, r_one$p))
one_match <- abs(r_one$p - 0.094) < 0.010
cat(sprintf("  MATCH: %s\n\n", ifelse(one_match, "YES", "NO")))

cat("Female three-meal:\n")
cat(sprintf("  Currently reported: coef=-682, SE=145, p=0.020\n"))
cat(sprintf("  hypotheses() result: coef=%.4f, SE=%.4f, t=%.4f, p=%.4f\n",
            r_f$coef, r_f$se, r_f$t, r_f$p))
f_match <- abs(r_f$p - 0.020) < 0.010
cat(sprintf("  MATCH: %s\n\n", ifelse(f_match, "YES", "NO")))

cat("Male three-meal:\n")
cat(sprintf("  Currently reported: coef=+314, SE=89, p=0.122\n"))
cat(sprintf("  hypotheses() result: coef=%.4f, SE=%.4f, t=%.4f, p=%.4f\n",
            r_m$coef, r_m$se, r_m$t, r_m$p))
m_match <- abs(r_m$p - 0.122) < 0.015
cat(sprintf("  MATCH: %s\n\n", ifelse(m_match, "YES", "NO")))

all_match <- three_match && one_match && f_match && m_match
cat("VERDICT:\n")
if (all_match) {
  cat("  Inference is correct. p-values are valid.\n")
} else {
  flags <- c()
  if (!three_match) flags <- c(flags, "Three-meal")
  if (!one_match)   flags <- c(flags, "One-meal")
  if (!f_match)     flags <- c(flags, "Female three-meal")
  if (!m_match)     flags <- c(flags, "Male three-meal")
  cat(sprintf("  INFERENCE ERROR DETECTED. Paper must be corrected before submission.\n"))
  cat(sprintf("  Affected results: %s\n", paste(flags, collapse=", ")))
}

cat("\n=== NOTE ON SE CONVENTION ===\n")
ct <- summary(m5)$coeftable
cat(sprintf("Interaction coef (Three): %.4f, SE: %.4f, p: %.4f\n",
            ct["diff_log_yield_hat:meal_typeThree","Estimate"],
            ct["diff_log_yield_hat:meal_typeThree","Std. Error"],
            ct["diff_log_yield_hat:meal_typeThree","Pr(>|t|)"]))
cat(sprintf("Baseline (None):          %.4f, SE: %.4f\n",
            ct["diff_log_yield_hat","Estimate"],
            ct["diff_log_yield_hat","Std. Error"]))
cat(sprintf("Total naive SE (sqrt sum of squares, ignoring cov): %.4f\n",
            sqrt(ct["diff_log_yield_hat:meal_typeThree","Std. Error"]^2 +
                   ct["diff_log_yield_hat","Std. Error"]^2)))
cat(sprintf("hypotheses() SE (delta method, with covariance):    %.4f\n",
            r_three$se))

cat("\n=== DONE ===\n")
