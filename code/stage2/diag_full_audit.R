## diag_full_audit.R
## Read-only diagnostic audit of ALL M5 specs.
## For each spec: interaction coef/SE/p vs delta-method total effect coef/SE/p.
## DO NOT modify any paper files.

suppressPackageStartupMessages({
  library(here); library(dplyr); library(readr); library(fixest)
})
ROOT <- here::here()

## ── Data ─────────────────────────────────────────────────────────────────── ##
df <- read_csv(file.path(ROOT, "data/Regression_data/df_2_merged_v2.csv"),
               show_col_types = FALSE) %>%
  mutate(
    meal_type      = relevel(factor(meal_type), ref = "None"),
    gender         = relevel(factor(gender),    ref = "Female"),
    growing_season = factor(growing_season, levels = c("Boro", "Aus", "Aman"))
  )
cat(sprintf("Main data: N=%d, G=%d\n\n", nrow(df), n_distinct(df$District)))

## ── Delta-method helper ───────────────────────────────────────────────────── ##
## Returns total=base+interaction; SE via full delta method including cov term
## p-value: two-sided t with G-1 df (clusters - 1)
dm_total <- function(model, base_nm, int_nm, G = NULL) {
  if (is.null(model)) return(list(coef=NA, se=NA, t=NA, p=NA))
  cf <- coef(model); V <- vcov(model)
  if (!base_nm %in% names(cf) || !int_nm %in% names(cf))
    return(list(coef=NA, se=NA, t=NA, p=NA))
  est <- cf[base_nm] + cf[int_nm]
  var_e <- V[base_nm,base_nm] + V[int_nm,int_nm] + 2*V[base_nm,int_nm]
  se_e <- sqrt(var_e)
  t_e  <- est / se_e
  if (is.null(G)) G <- length(unique(model$fixef_id[[1]]))
  p_e  <- 2 * pt(-abs(t_e), df = G - 1)
  list(coef = as.numeric(est), se = as.numeric(se_e),
       t = as.numeric(t_e), p = as.numeric(p_e))
}

## Interaction-only extractor (what paper currently uses as "total SE/p")
int_only <- function(model, base_nm, int_nm) {
  if (is.null(model)) return(list(coef=NA, int_coef=NA, se=NA, p=NA, n=NA))
  cf <- coef(model); sv <- se(model); pv <- pvalue(model)
  if (!base_nm %in% names(cf) || !int_nm %in% names(cf))
    return(list(coef=NA, int_coef=NA, se=NA, p=NA, n=NA))
  list(
    coef     = as.numeric(cf[base_nm] + cf[int_nm]),  # total point estimate
    int_coef = as.numeric(cf[int_nm]),                # interaction only
    se       = as.numeric(sv[int_nm]),                # SE of interaction
    p        = as.numeric(pv[int_nm]),                # p of interaction
    n        = nobs(model)
  )
}

## Safe feols
sf <- function(fml, data, clust = ~District) {
  tryCatch(
    feols(fml, data=data, cluster=clust, warn=FALSE, notes=FALSE),
    error = function(e) { cat("  FAILED:", conditionMessage(e), "\n"); NULL }
  )
}

## Audit one spec: print row
audit_row <- function(label, model, base_nm = "diff_log_yield_hat",
                      int_nm = "diff_log_yield_hat:meal_typeThree",
                      G = NULL, reported_coef, reported_se, reported_p,
                      meal = "THREE") {
  cat(sprintf("\nROW: %s [%s]\n", label, meal))
  if (is.null(model)) {
    cat("  Model failed — cannot audit\n")
    return(invisible(NULL))
  }
  io  <- int_only(model, base_nm, int_nm)
  dm  <- dm_total(model, base_nm, int_nm, G)
  # What does the paper say?
  cat(sprintf("  Interaction coef: %.4f, SE=%.4f, p=%.4f\n",
              io$int_coef, io$se, io$p))
  cat(sprintf("  Total effect (delta-method): coef=%.4f, SE=%.4f, p=%.4f\n",
              dm$coef, dm$se, dm$p))
  cat(sprintf("  Currently reported in paper: coef=%.4f, SE=%.4f, p=%.4f\n",
              reported_coef, reported_se, reported_p))
  # Determine what the paper is using
  se_matches_int   <- abs(io$se   - reported_se)  < 0.1
  p_matches_int    <- abs(io$p    - reported_p)   < 0.005
  se_matches_delta <- abs(dm$se   - reported_se)  < 1.0
  p_matches_delta  <- abs(dm$p    - reported_p)   < 0.01
  if (p_matches_int && se_matches_int) {
    cat(sprintf("  SE source (paper): INTERACTION\n"))
    cat(sprintf("  MATCH on p-value: YES (interaction p)\n"))
  } else if (p_matches_delta) {
    cat(sprintf("  SE source (paper): TOTAL (delta-method)\n"))
    cat(sprintf("  MATCH on p-value: YES (delta-method p)\n"))
  } else {
    cat(sprintf("  SE source (paper): UNKNOWN — neither interaction nor delta-method matches\n"))
    cat(sprintf("  MATCH on p-value: NO\n"))
  }
  invisible(list(label=label, meal=meal, io=io, dm=dm,
                 reported_coef=reported_coef, reported_se=reported_se,
                 reported_p=reported_p,
                 p_matches_int=p_matches_int, p_matches_delta=p_matches_delta))
}

## ═══════════════════════════════════════════════════════════════════════════ ##
## FIT ALL SPECS                                                                ##
## ═══════════════════════════════════════════════════════════════════════════ ##

## R1: Add controls
census <- read_csv(file.path(ROOT, "data/agricultural_census_2019.csv"),
                   show_col_types = FALSE) %>%
  mutate(irrigation_share = Net_Irrigated_Area / Net_Cultivated_Area,
         avg_holdings     = Number_of_Holdings / Net_Cultivated_Area * 1000,
         crop_intensity   = Intensity_of_Cropping) %>%
  select(District, irrigation_share, avg_holdings, crop_intensity)
df_r1 <- df %>% left_join(census, by="District") %>%
  filter(!is.na(irrigation_share), !is.na(avg_holdings), !is.na(crop_intensity))

## R3: Levels
df_r3 <- df %>% arrange(District, growing_season, year) %>%
  group_by(District, growing_season) %>%
  mutate(log_yield_hat = cumsum(diff_log_yield_hat)) %>% ungroup()

## R5: Drop top-5 variance
top5 <- df %>% group_by(District) %>%
  summarise(v=var(diff_log_yield_hat, na.rm=TRUE)) %>%
  arrange(desc(v)) %>% slice_head(n=5) %>% pull(District)
cat("Top-5 variance districts:", paste(top5, collapse=", "), "\n")
df_r5 <- df %>% filter(!District %in% top5)

## R6
df_r6a <- df %>% filter(year <= 2020)
df_r6b <- df %>% filter(year >= 2021)

## R10: log wage (diff_log_real_wage column)
cat("Checking log wage column: ")
cat(ifelse("diff_log_real_wage" %in% names(df), "diff_log_real_wage OK\n", "MISSING\n"))

## R11: nominal wage
wage_raw <- read_csv(file.path(ROOT, "data/Regression_data/wage_by_growing_season.csv"),
                     show_col_types = FALSE) %>% filter(!is.na(wage))
wage_agg <- wage_raw %>%
  group_by(District, growing_season, gender, meal_type, year) %>%
  summarise(nominal_wage = mean(wage, na.rm=TRUE), .groups="drop")
wage_fd <- wage_agg %>%
  arrange(District, growing_season, gender, meal_type, year) %>%
  group_by(District, growing_season, gender, meal_type) %>%
  mutate(diff_nominal_wage = nominal_wage - lag(nominal_wage)) %>%
  ungroup() %>% filter(!is.na(diff_nominal_wage))
df_r11 <- df %>% select(District, growing_season, gender, meal_type, year,
                         diff_log_yield_hat) %>%
  inner_join(wage_fd %>% select(District, growing_season, gender, meal_type,
                                 year, diff_nominal_wage),
             by=c("District","growing_season","gender","meal_type","year")) %>%
  mutate(meal_type = relevel(factor(meal_type), ref="None"),
         gender    = relevel(factor(gender),    ref="Female"))

## Gender subsamples
df_f <- df %>% filter(gender == "Female")
df_m <- df %>% filter(gender == "Male")
df_boro <- df %>% filter(growing_season == "Boro")

## Fit all models
## IMPORTANT: paper's M5 uses +: form (no meal_type main effects).
## stage2_robustness_m5.R converts * to +: before running.
## rob_r1_r3_r6_r11.R uses * form directly (slight inconsistency — noted in report).
cat("\n--- Fitting all models ---\n")
## +: form (matches paper's 01_main_regressions.R and stage2_robustness_m5.R)
fml_base <- diff_real_wage  ~ diff_log_yield_hat + diff_log_yield_hat:meal_type + gender | year + District^growing_season
fml_r2   <- diff_real_wage  ~ diff_log_yield_hat + diff_log_yield_hat:meal_type + gender | year + District
fml_log  <- diff_log_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal_type + gender | year + District^growing_season
fml_boro <- diff_real_wage  ~ diff_log_yield_hat + diff_log_yield_hat:meal_type + gender | year + District
fml_gen  <- diff_real_wage  ~ diff_log_yield_hat + diff_log_yield_hat:meal_type | year + District^growing_season
## * form (matches rob_r1_r3_r6_r11.R — note inconsistency with baseline)
fml_r1s  <- diff_real_wage ~ diff_log_yield_hat * meal_type + gender | year + District^growing_season
fml_r3s  <- real_wage ~ log_yield_hat * meal_type + gender | year + District^growing_season
fml_r11s <- diff_nominal_wage ~ diff_log_yield_hat * meal_type + gender | year + District^growing_season

m_base  <- sf(fml_base,  df);      cat("  Baseline (+: form): OK\n")
m_r1    <- sf(fml_r1s,   df_r1);   cat("  R1 (* form):        OK\n")
m_r2    <- sf(fml_r2,    df);      cat("  R2 (+: form):       OK\n")
m_r3    <- sf(fml_r3s,   df_r3);   cat("  R3 (* form):        OK\n")
m_r4    <- sf(fml_boro,  df_boro); cat("  R4 (+: form):       OK\n")
m_r5    <- sf(fml_base,  df_r5);   cat("  R5 (+: form):       OK\n")
m_r6a   <- sf(fml_r1s,   df_r6a);  cat("  R6a (* form):       OK\n")
m_r6b   <- sf(fml_r1s,   df_r6b);  cat("  R6b (* form):       OK\n")
m_r10   <- sf(fml_log,   df);      cat("  R10 (+: form):      OK\n")
m_r11   <- sf(fml_r11s,  df_r11);  cat("  R11 (* form):       OK\n")
m_f     <- sf(fml_gen,   df_f);    cat("  Female (+: form):   OK\n")
m_m     <- sf(fml_gen,   df_m);    cat("  Male (+: form):     OK\n")

## DK SE and two-way cluster (R8/R9): refit with different VCV
## R8: Driscoll-Kraay lag=2
m_r8  <- tryCatch(
  feols(fml_base, data=df, se="driscoll_kraay", dof=dof(adj=FALSE), notes=FALSE,
        panel.id = ~District+year),
  error = function(e) {
    cat("  R8 DK feols failed:", conditionMessage(e), "\n")
    # Try alternative
    tryCatch(
      feols(fml_base, data=df, vcov=vcov_DK(lag=2), notes=FALSE),
      error = function(e2) { cat("  R8 alt also failed:", conditionMessage(e2), "\n"); NULL }
    )
  }
)

## R9: Two-way cluster
m_r9  <- tryCatch(
  feols(fml_base, data=df, cluster=~District+year, warn=FALSE, notes=FALSE),
  error = function(e) { cat("  R9 two-way cluster failed:", conditionMessage(e), "\n"); NULL }
)

cat("\n")

## ═══════════════════════════════════════════════════════════════════════════ ##
## AUDIT TABLE                                                                  ##
## ═══════════════════════════════════════════════════════════════════════════ ##
cat("\n\n========================================================\n")
cat("FULL AUDIT: Interaction vs. Delta-Method Total Effect\n")
cat("========================================================\n")
cat("Currently reported = what appears in paper's rob_m5_threemeal.tex\n")

G_all <- n_distinct(df$District)       # 63

results <- list()

## ── BASELINE ──
results[[1]] <- audit_row("Baseline M5", m_base,
  "diff_log_yield_hat", "diff_log_yield_hat:meal_typeThree", G=G_all,
  reported_coef=-103.69, reported_se=94.99, reported_p=0.022)

## ── R1 ──
G_r1 <- n_distinct(df_r1$District)
results[[2]] <- audit_row("R1: Add controls", m_r1,
  "diff_log_yield_hat", "diff_log_yield_hat:meal_typeThree", G=G_r1,
  reported_coef=-283.23, reported_se=104.25, reported_p=0.045)

## ── R2 ──
results[[3]] <- audit_row("R2: District FE only", m_r2,
  "diff_log_yield_hat", "diff_log_yield_hat:meal_typeThree", G=G_all,
  reported_coef=-84.12, reported_se=95.29, reported_p=0.021)

## ── R3 (levels) ──
results[[4]] <- audit_row("R3: Levels", m_r3,
  "log_yield_hat", "log_yield_hat:meal_typeThree", G=G_all,
  reported_coef=29.51, reported_se=20.95, reported_p=0.103)

## ── R4 (Boro only) ──
G_boro <- n_distinct(df_boro$District)
results[[5]] <- audit_row("R4: Boro only", m_r4,
  "diff_log_yield_hat", "diff_log_yield_hat:meal_typeThree", G=G_boro,
  reported_coef=-760.02, reported_se=137.71, reported_p=0.060)

## ── R5 (drop top-5) ──
G_r5 <- n_distinct(df_r5$District)
results[[6]] <- audit_row("R5: Drop top-5 var", m_r5,
  "diff_log_yield_hat", "diff_log_yield_hat:meal_typeThree", G=G_r5,
  reported_coef=80.61, reported_se=97.28, reported_p=0.056)

## ── R6a ──
G_r6a <- n_distinct(df_r6a$District)
results[[7]] <- audit_row("R6a: Early 2017-2020", m_r6a,
  "diff_log_yield_hat", "diff_log_yield_hat:meal_typeThree", G=G_r6a,
  reported_coef=-532.51, reported_se=137.37, reported_p=0.005)

## ── R6b ──
G_r6b <- n_distinct(df_r6b$District)
results[[8]] <- audit_row("R6b: Late 2021-2023", m_r6b,
  "diff_log_yield_hat", "diff_log_yield_hat:meal_typeThree", G=G_r6b,
  reported_coef=239.15, reported_se=127.13, reported_p=0.305)

## ── R8: DK SE ──
cat("\nROW: R8: DK SE (lag=2) [THREE]\n")
if (!is.null(m_r8)) {
  io8 <- int_only(m_r8, "diff_log_yield_hat", "diff_log_yield_hat:meal_typeThree")
  cat(sprintf("  Interaction coef: %.4f, SE=%.4f, p=%.4f\n",
              io8$int_coef, io8$se, io8$p))
  cat(sprintf("  Currently reported: coef=-103.69, SE=95.75, p=0.020\n"))
  cat(sprintf("  MATCH: %s\n", ifelse(abs(io8$p-0.020)<0.005, "YES","NO")))
  cat("  NOTE: delta-method NA for DK VCV (not applicable for non-cluster inference)\n")
  results[[9]] <- list(label="R8: DK SE", io=io8, dm=NULL,
                        reported_p=0.020, p_matches_int=abs(io8$p-0.020)<0.005)
} else {
  cat("  Model not available\n")
  results[[9]] <- NULL
}

## ── R9: Two-way cluster ──
cat("\nROW: R9: Two-way cluster [THREE]\n")
if (!is.null(m_r9)) {
  io9 <- int_only(m_r9, "diff_log_yield_hat", "diff_log_yield_hat:meal_typeThree")
  cat(sprintf("  Interaction coef: %.4f, SE=%.4f, p=%.4f\n",
              io9$int_coef, io9$se, io9$p))
  cat(sprintf("  Currently reported: coef=-103.69, SE=101.17, p=0.028\n"))
  cat(sprintf("  MATCH: %s\n", ifelse(abs(io9$p-0.028)<0.005, "YES","NO")))
  cat("  NOTE: delta-method not computed for two-way cluster\n")
  results[[10]] <- list(label="R9: Two-way cluster", io=io9, dm=NULL,
                         reported_p=0.028, p_matches_int=abs(io9$p-0.028)<0.005)
} else {
  cat("  Model not available\n")
  results[[10]] <- NULL
}

## ── R10: Log wage ──
results[[11]] <- audit_row("R10: Log wage", m_r10,
  "diff_log_yield_hat", "diff_log_yield_hat:meal_typeThree", G=G_all,
  reported_coef=-0.56, reported_se=0.25, reported_p=0.031)

## ── R11: Nominal wage ──
G_r11 <- n_distinct(df_r11$District)
results[[12]] <- audit_row("R11: Nominal wage", m_r11,
  "diff_log_yield_hat", "diff_log_yield_hat:meal_typeThree", G=G_r11,
  reported_coef=-78.81, reported_se=89.75, reported_p=0.043)

## ── Female ──
G_f <- n_distinct(df_f$District)
results[[13]] <- audit_row("Female only", m_f,
  "diff_log_yield_hat", "diff_log_yield_hat:meal_typeThree", G=G_f,
  reported_coef=-682.03, reported_se=144.75, reported_p=0.019)

## ── Male ──
G_m <- n_distinct(df_m$District)
results[[14]] <- audit_row("Male only", m_m,
  "diff_log_yield_hat", "diff_log_yield_hat:meal_typeThree", G=G_m,
  reported_coef=313.51, reported_se=89.31, reported_p=0.122)

## ═══════════════════════════════════════════════════════════════════════════ ##
## ONE-MEAL CHECK: Baseline, R2, R4, R10                                        ##
## ═══════════════════════════════════════════════════════════════════════════ ##
cat("\n\n=== ONE-MEAL TOTAL EFFECT AUDIT (selected specs) ===\n")
one_audit <- function(label, model, G, reported_coef, reported_se, reported_p,
                       base_nm="diff_log_yield_hat") {
  int_nm <- paste0(base_nm, ":meal_typeOne")
  cat(sprintf("\nROW: %s [ONE]\n", label))
  if (is.null(model)) { cat("  Model failed\n"); return(invisible(NULL)) }
  io <- int_only(model, base_nm, int_nm)
  dm <- dm_total(model, base_nm, int_nm, G)
  cat(sprintf("  Interaction coef: %.4f, SE=%.4f, p=%.4f\n",
              io$int_coef, io$se, io$p))
  cat(sprintf("  Total effect (delta-method): coef=%.4f, SE=%.4f, p=%.4f\n",
              dm$coef, dm$se, dm$p))
  cat(sprintf("  Currently reported: coef=%.4f, SE=%.4f, p=%.4f\n",
              reported_coef, reported_se, reported_p))
  p_match_int <- abs(io$p - reported_p) < 0.005
  cat(sprintf("  SE source (paper): %s\n",
              ifelse(p_match_int, "INTERACTION", "OTHER/TOTAL")))
  cat(sprintf("  MATCH on p-value: %s\n", ifelse(p_match_int, "YES","NO")))
}
one_audit("Baseline M5", m_base, G_all,  -55.47, 102.44, 0.094)
one_audit("R2: District FE", m_r2,  G_all,  -35.17, 103.37, 0.091)
one_audit("R4: Boro only",   m_r4,  G_boro, -855.62, 134.39, 0.009)
one_audit("R10: Log wage",   m_r10, G_all,  -0.42, 0.25, 0.108)

## ═══════════════════════════════════════════════════════════════════════════ ##
## R5 and R6 INTERACTION vs TOTAL CHECK                                         ##
## ═══════════════════════════════════════════════════════════════════════════ ##
cat("\n\n=== STEP 6: R5 and R6 — INTERACTION vs TOTAL sign ===\n")

check_sign <- function(label, model, base_nm="diff_log_yield_hat",
                        int_nm="diff_log_yield_hat:meal_typeThree", G=63,
                        baseline_m=m_base) {
  if (is.null(model)) { cat(sprintf("  %s: model failed\n", label)); return(invisible(NULL)) }
  cf <- coef(model)
  io <- int_only(model, base_nm, int_nm)
  dm <- dm_total(model, base_nm, int_nm, G)
  base_coef <- cf[base_nm]
  int_coef  <- cf[int_nm]
  total     <- base_coef + int_coef
  base_int  <- coef(baseline_m)[int_nm]
  cat(sprintf("\n%s:\n", label))
  cat(sprintf("  None baseline:     %.4f\n", base_coef))
  cat(sprintf("  Three interaction: %.4f  SE=%.4f  p=%.4f\n", int_coef, io$se, io$p))
  cat(sprintf("  Three total:       %.4f  SE=%.4f  p=%.4f\n", dm$coef, dm$se, dm$p))
  cat(sprintf("  INTERACTION sign:  %s | TOTAL sign: %s\n",
              ifelse(int_coef < 0, "NEGATIVE", "POSITIVE"),
              ifelse(total < 0, "NEGATIVE", "POSITIVE")))
  cat(sprintf("  Does INTERACTION reverse vs baseline? baseline int=%.4f, this int=%.4f -> %s\n",
              base_int, int_coef,
              ifelse(sign(int_coef) != sign(base_int),
                     "YES — SIGN REVERSAL", "NO — same sign")))
}

check_sign("R5: Drop top-5 var",    m_r5,  G=G_r5)
## R6a/R6b use * form — baseline reference should also use * form for fair comparison
m_base_star <- sf(fml_r1s, df)  # * form baseline for R6 comparison
check_sign("R6a: Early 2017-2020",  m_r6a, G=G_r6a, baseline_m=m_base_star)
check_sign("R6b: Late 2021-2023",   m_r6b, G=G_r6b, baseline_m=m_base_star)

## ═══════════════════════════════════════════════════════════════════════════ ##
## SUMMARY                                                                      ##
## ═══════════════════════════════════════════════════════════════════════════ ##
cat("\n\nSUMMARY\n")
cat("=======\n\n")

## Collect interaction p-values and delta-method p-values from results
valid_results <- Filter(Negate(is.null), results)
int_sig_labs  <- c()
dm_sig_labs   <- c()
mismatch_labs <- c()
int_coefs     <- c()

for (r in valid_results) {
  if (is.null(r)) next
  io_p <- if (!is.null(r$io)) r$io$p else NA
  dm_p <- if (!is.null(r$dm)) r$dm$p else NA
  int_c <- if (!is.null(r$io)) r$io$int_coef else NA
  if (!is.na(io_p) && io_p < 0.10) int_sig_labs <- c(int_sig_labs, r$label)
  if (!is.na(dm_p) && dm_p < 0.10) dm_sig_labs  <- c(dm_sig_labs,  r$label)
  if (!is.na(io_p) && !is.na(dm_p) && abs(io_p - dm_p) > 0.05)
    mismatch_labs <- c(mismatch_labs, r$label)
  if (!is.na(int_c)) int_coefs <- c(int_coefs, setNames(int_c, r$label))
}

cat("Results where INTERACTION p-value is valid (differential claim, p < 0.10):\n")
if (length(int_sig_labs)==0) cat("  None\n") else
  for (l in int_sig_labs) cat(sprintf("  - %s\n", l))

cat("\nResults where TOTAL EFFECT p-value is valid (delta-method, p < 0.10):\n")
if (length(dm_sig_labs)==0) cat("  None\n") else
  for (l in dm_sig_labs) cat(sprintf("  - %s\n", l))

cat("\nResults where paper reports interaction p as total p (MISMATCH > 0.05):\n")
if (length(mismatch_labs)==0) cat("  None\n") else
  for (l in mismatch_labs) cat(sprintf("  - %s\n", l))

cat("\nINTERACTION STABILITY (Three-meal interaction, not total):\n")
n_sig_int <- sum(!is.na(int_coefs) & sapply(valid_results, function(r) {
  if (is.null(r) || is.null(r$io)) return(FALSE)
  !is.na(r$io$p) && r$io$p < 0.10
}))
n_total <- length(int_coefs)
cat(sprintf("  Significant at p < 0.10: %d of %d specs\n", n_sig_int, n_total))
cat(sprintf("  Range of interaction coefs: [%.2f, %.2f]\n",
            min(int_coefs, na.rm=TRUE), max(int_coefs, na.rm=TRUE)))
cat(sprintf("  Baseline interaction coef: %.4f\n",
            coef(m_base)["diff_log_yield_hat:meal_typeThree"]))
n_same_sign <- sum(int_coefs < 0, na.rm=TRUE)
cat(sprintf("  Same sign as baseline (negative): %d of %d\n", n_same_sign, n_total))
cat(sprintf("  INTERACTION: %s\n",
            ifelse(n_sig_int >= n_total * 0.7 &&
                     n_same_sign >= n_total * 0.7, "STABLE", "UNSTABLE")))

cat("\nAudit complete. Awaiting rewrite decision.\n")
