## 17_identification_diagnostics.R
## Identification diagnostics: first-stage strength, reduced form heterogeneity,
## RF interaction tests, and joint bootstrap for generated-regressor inference
## Output: output/stage2/identification_diagnostics/

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(fixest)
  library(ggplot2)
  library(kableExtra)
})

ROOT <- here::here()
OUT  <- file.path(ROOT, "output/stage2/identification_diagnostics")
for (d in c("tables", "figures", "logs", "models", "summary"))
  dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(OUT, "logs/identification_diagnostics_log.txt")
log_con  <- file(log_path, open = "wt")
on.exit(close(log_con), add = TRUE)

CLIM_VARS <- c("diff_gdd_10_30", "diff_edd_30", "diff_precip", "diff_precip_sq")
FE_YIELD  <- "year + district^season"
FE_WAGE   <- "year + District^growing_season"
STOCK_YOGO_10PCT <- 16.85  ## 1 endogenous regressor, 4 instruments (approx.)

log_msg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  cat(msg, "\n", file = log_con)
}

fmt_p <- function(p) {
  vapply(p, function(x) {
    if (is.na(x)) return("NA")
    if (x < 0.001) return("<0.001")
    sprintf("%.3f", x)
  }, character(1))
}

sig_flag <- function(p) {
  ifelse(is.na(p), "imprecise",
         ifelse(p < 0.05, "significant",
                ifelse(p < 0.10, "marginally significant", "imprecise")))
}

extract_climate_coefs <- function(mod, vars = CLIM_VARS) {
  ct <- summary(mod)$coeftable
  out <- as.data.frame(ct[vars, , drop = FALSE])
  out$variable <- rownames(out)
  rownames(out) <- NULL
  names(out) <- c("estimate", "std_error", "t_value", "p_value", "variable")
  out[, c("variable", "estimate", "std_error", "p_value")]
}

first_stage_diag <- function(data, sample_label) {
  fml <- as.formula(paste(
    "diff_log_yield ~", paste(CLIM_VARS, collapse = " + "), "|", FE_YIELD
  ))
  mod <- feols(fml, data = data, cluster = ~district)
  w   <- wald(mod, CLIM_VARS)
  r2  <- as.numeric(fitstat(mod, "r2"))
  wr2 <- as.numeric(fitstat(mod, "wr2"))
  coefs <- extract_climate_coefs(mod)
  coefs$sample <- sample_label
  list(
    sample = sample_label,
    nobs = nobs(mod),
    coefs = coefs,
    partial_f = w$stat,
    partial_f_p = w$p,
    r_squared = r2,
    within_r_squared = wr2,
    clears_f10 = w$stat >= 10,
    clears_stock_yogo = w$stat >= STOCK_YOGO_10PCT,
    model = mod
  )
}

rf_by_meal <- function(data, meal_level, meal_label) {
  sub <- data %>% filter(meal_type == meal_level)
  fml <- as.formula(paste(
    "diff_real_wage ~", paste(CLIM_VARS, collapse = " + "), "|", FE_WAGE
  ))
  mod <- feols(fml, data = sub, cluster = ~District)
  ct  <- summary(mod)$coeftable
  sd_gdd <- sd(sub$diff_gdd_10_30, na.rm = TRUE)
  sd_edd <- sd(sub$diff_edd_30, na.rm = TRUE)
  data.frame(
    meal_type = meal_level,
    meal_label = meal_label,
    n_obs = nobs(mod),
    gdd_coef = ct["diff_gdd_10_30", "Estimate"],
    gdd_se   = ct["diff_gdd_10_30", "Std. Error"],
    gdd_p    = ct["diff_gdd_10_30", "Pr(>|t|)"],
    edd_coef = ct["diff_edd_30", "Estimate"],
    edd_se   = ct["diff_edd_30", "Std. Error"],
    edd_p    = ct["diff_edd_30", "Pr(>|t|)"],
    precip_coef = ct["diff_precip", "Estimate"],
    precip_se   = ct["diff_precip", "Std. Error"],
    precip_p    = ct["diff_precip", "Pr(>|t|)"],
    precip_sq_coef = ct["diff_precip_sq", "Estimate"],
    precip_sq_se   = ct["diff_precip_sq", "Std. Error"],
    precip_sq_p    = ct["diff_precip_sq", "Pr(>|t|)"],
    gdd_effect_1sd = ct["diff_gdd_10_30", "Estimate"] * sd_gdd,
    edd_effect_1sd = ct["diff_edd_30", "Estimate"] * sd_edd,
    stringsAsFactors = FALSE
  )
}

run_rf_interaction <- function(data, sample_label) {
  int_vars <- paste0(CLIM_VARS, ":meal3")
  fml <- as.formula(paste(
    "diff_real_wage ~", paste(c(CLIM_VARS, int_vars), collapse = " + "),
    "|", FE_WAGE
  ))
  mod <- feols(fml, data = data, cluster = ~District)
  ct  <- summary(mod)$coeftable
  w_joint <- wald(mod, int_vars)
  rows <- lapply(int_vars, function(v) {
    data.frame(
      sample = sample_label,
      term = v,
      estimate = ct[v, "Estimate"],
      std_error = ct[v, "Std. Error"],
      p_value = ct[v, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
  rows$joint_f <- w_joint$stat
  rows$joint_f_p <- w_joint$p
  list(model = mod, interactions = rows, joint = data.frame(
    sample = sample_label,
    joint_f = w_joint$stat,
    joint_f_p = w_joint$p,
    nobs = nobs(mod)
  ))
}

extract_m5_diff <- function(mod) {
  ct <- summary(mod)$coeftable
  int <- "diff_log_yield_hat:meal_typeThree"
  data.frame(
    differential = ct[int, "Estimate"],
    se = ct[int, "Std. Error"],
    p = ct[int, "Pr(>|t|)"],
    nobs = nobs(mod)
  )
}

log_msg("=== IDENTIFICATION DIAGNOSTICS ===")
log_msg("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

## ── Build yield (first stage) panel ─────────────────────────────────────────
df_yield <- read_csv(
  file.path(ROOT, "data/Regression_data/bangladesh_rice_regression_panel.csv"),
  show_col_types = FALSE
) %>%
  arrange(district, season, year) %>%
  group_by(district, season) %>%
  mutate(
    log_yield      = log(yield_per_ha),
    diff_log_yield = log_yield - lag(log_yield),
    diff_gdd_10_30 = gdd_10_30 - lag(gdd_10_30)
  ) %>%
  ungroup() %>%
  mutate(
    diff_precip    = diff_pr1,
    diff_precip_sq = diff_pr1^2
  ) %>%
  filter(
    !is.na(diff_log_yield),
    !is.na(diff_gdd_10_30),
    !is.na(diff_edd_30),
    !is.na(diff_precip),
    !is.na(diff_precip_sq)
  )

## ── Build wage + climate reduced-form panel ───────────────────────────────────
df_wage <- read_csv(
  file.path(ROOT, "data/Regression_data/df_2_merged_v2.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    meal_type = relevel(factor(meal_type), ref = "None"),
    gender    = relevel(factor(gender),    ref = "Female")
  ) %>%
  filter(!is.na(diff_real_wage))

climate_merge <- df_yield %>%
  rename(District = district, growing_season = season) %>%
  group_by(District, growing_season, year) %>%
  summarise(across(all_of(CLIM_VARS), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

df_rf <- df_wage %>%
  left_join(climate_merge, by = c("District", "growing_season", "year")) %>%
  filter(if_all(all_of(CLIM_VARS), ~ !is.na(.x)))

log_msg("Yield FS panel N = ", nrow(df_yield))
log_msg("Wage RF panel N = ", nrow(df_rf))

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 1 — First stage F-statistics by season                               ##
## ════════════════════════════════════════════════════════════════════════════ ##
log_msg("\nSTEP 1 — First stage by season")

fs_samples <- list(
  Pooled = df_yield,
  Boro   = df_yield %>% filter(season == "Boro"),
  Aus    = df_yield %>% filter(season == "Aus"),
  Aman   = df_yield %>% filter(season == "Aman")
)

fs_results <- lapply(names(fs_samples), function(nm) {
  first_stage_diag(fs_samples[[nm]], nm)
})
names(fs_results) <- names(fs_samples)

fs_coef_tbl <- bind_rows(lapply(fs_results, `[[`, "coefs"))
fs_summary <- bind_rows(lapply(fs_results, function(x) {
  data.frame(
    sample = x$sample,
    partial_f = x$partial_f,
    partial_f_p = x$partial_f_p,
    r_squared = x$r_squared,
    within_r_squared = x$within_r_squared,
    nobs = x$nobs,
    clears_f10 = x$clears_f10,
    clears_stock_yogo_10pct = x$clears_stock_yogo,
    stringsAsFactors = FALSE
  )
}))

write_csv(fs_coef_tbl, file.path(OUT, "tables/first_stage_coefs_by_season.csv"))
write_csv(fs_summary, file.path(OUT, "tables/first_stage_summary_by_season.csv"))

for (x in fs_results) {
  log_msg(sprintf("  %s: partial F = %.2f (p=%s) | R2=%.3f | within R2=%.3f | N=%d | F>=10: %s",
                  x$sample, x$partial_f, fmt_p(x$partial_f_p),
                  x$r_squared, x$within_r_squared, x$nobs, x$clears_f10))
  for (i in seq_len(nrow(x$coefs))) {
    r <- x$coefs[i, ]
    log_msg(sprintf("    %s: %.4f (SE %.4f, p=%s)", r$variable, r$estimate,
                    r$std_error, fmt_p(r$p_value)))
  }
}

boro_f <- fs_summary$partial_f[fs_summary$sample == "Boro"]
aus_f  <- fs_summary$partial_f[fs_summary$sample == "Aus"]
aman_f <- fs_summary$partial_f[fs_summary$sample == "Aman"]
q1_answer <- boro_f > max(aus_f, aman_f)

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 2 — Reduced form by meal category                                    ##
## ════════════════════════════════════════════════════════════════════════════ ##
log_msg("\nSTEP 2 — Reduced form by meal category")

meal_levels <- c("None", "One", "Two", "Three")
meal_labels <- c("0-meal", "1-meal", "2-meal", "3-meal")

rf_meal <- bind_rows(Map(rf_by_meal, meal_levels, meal_labels, MoreArgs = list(data = df_rf)))
rf_meal <- rf_meal %>%
  mutate(
    gdd_sig = sig_flag(gdd_p),
    edd_sig = sig_flag(edd_p)
  )

write_csv(rf_meal, file.path(OUT, "tables/reduced_form_by_meal.csv"))

for (i in seq_len(nrow(rf_meal))) {
  r <- rf_meal[i, ]
  log_msg(sprintf("  %s (N=%d): GDD=%.2f (p=%s, 1SD=%.1f BDT) | EDD=%.2f (p=%s, 1SD=%.1f BDT)",
                  r$meal_label, r$n_obs, r$gdd_coef, fmt_p(r$gdd_p), r$gdd_effect_1sd,
                  r$edd_coef, fmt_p(r$edd_p), r$edd_effect_1sd))
}

gdd_0 <- rf_meal$gdd_coef[rf_meal$meal_type == "None"]
gdd_3 <- rf_meal$gdd_coef[rf_meal$meal_type == "Three"]
edd_0 <- rf_meal$edd_coef[rf_meal$meal_type == "None"]
edd_3 <- rf_meal$edd_coef[rf_meal$meal_type == "Three"]

q2_pattern_gdd <- gdd_3 < 0 && abs(gdd_0) < abs(gdd_3)
q2_pattern_edd <- edd_3 < 0 && abs(edd_0) < abs(edd_3)

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 3 — Pooled RF differential (3-meal vs 0-meal)                        ##
## ════════════════════════════════════════════════════════════════════════════ ##
log_msg("\nSTEP 3 — Reduced form interaction (3-meal vs 0-meal, pooled)")

df_30 <- df_rf %>%
  filter(meal_type %in% c("None", "Three")) %>%
  mutate(meal3 = as.integer(meal_type == "Three"))

rf_int_pooled <- run_rf_interaction(df_30, "Pooled")
write_csv(rf_int_pooled$interactions, file.path(OUT, "tables/rf_interaction_pooled.csv"))
write_csv(rf_int_pooled$joint, file.path(OUT, "tables/rf_interaction_pooled_joint.csv"))

for (i in seq_len(nrow(rf_int_pooled$interactions))) {
  r <- rf_int_pooled$interactions[i, ]
  log_msg(sprintf("  %s: %.2f (SE %.2f, p=%s)", r$term, r$estimate, r$std_error, fmt_p(r$p_value)))
}
log_msg(sprintf("  Joint F-test (all interactions): F=%.2f, p=%s",
                rf_int_pooled$joint$joint_f, fmt_p(rf_int_pooled$joint$joint_f_p)))

gdd_int <- rf_int_pooled$interactions %>%
  filter(term == "diff_gdd_10_30:meal3")
edd_int <- rf_int_pooled$interactions %>%
  filter(term == "diff_edd_30:meal3")
## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 4 — Boro-specific RF differential                                    ##
## ════════════════════════════════════════════════════════════════════════════ ##
log_msg("\nSTEP 4 — Reduced form interaction (Boro only)")

df_30_boro <- df_30 %>% filter(growing_season == "Boro")
rf_int_boro <- run_rf_interaction(df_30_boro, "Boro")
write_csv(rf_int_boro$interactions, file.path(OUT, "tables/rf_interaction_boro.csv"))
write_csv(rf_int_boro$joint, file.path(OUT, "tables/rf_interaction_boro_joint.csv"))

for (i in seq_len(nrow(rf_int_boro$interactions))) {
  r <- rf_int_boro$interactions[i, ]
  log_msg(sprintf("  %s: %.2f (SE %.2f, p=%s)", r$term, r$estimate, r$std_error, fmt_p(r$p_value)))
}
log_msg(sprintf("  Joint F-test: F=%.2f, p=%s",
                rf_int_boro$joint$joint_f, fmt_p(rf_int_boro$joint$joint_f_p)))

q3_sig_pooled <- (gdd_int$p_value < 0.05) || (edd_int$p_value < 0.05) ||
  rf_int_pooled$joint$joint_f_p < 0.05
q3_sig_boro <- rf_int_boro$joint$joint_f_p < 0.05

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 5 — Joint cluster bootstrap (FS + SS)                                ##
## NOTE: FD bootstrap below is legacy. Canonical levels bootstrap:            ##
##       Rscript code/paper/levels_bootstrap.R                                ##
## ════════════════════════════════════════════════════════════════════════════ ##
log_msg("\nSTEP 5 — Joint cluster bootstrap (FD legacy; see levels_bootstrap.R)")

B_BOOT <- 9999
set.seed(42)

districts_yield <- unique(df_yield$district)
## Use merged panel with yield_hat for non-bootstrap benchmark
df_ts <- read_csv(file.path(ROOT, "data/Regression_data/df_2_merged_v2.csv"),
                  show_col_types = FALSE) %>%
  mutate(
    meal_type = relevel(factor(meal_type), ref = "None"),
    gender    = relevel(factor(gender), ref = "Female")
  ) %>%
  filter(!is.na(diff_log_yield_hat), !is.na(diff_real_wage))

m5_orig <- feols(
  diff_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal_type + gender |
    year + District^growing_season,
  data = df_ts,
  cluster = ~District
)
orig_diff <- extract_m5_diff(m5_orig)

## Verify bootstrap estimand: with None baseline, interaction coef = attached-casual differential
cf_orig <- coef(m5_orig)
base_none <- unname(cf_orig["diff_log_yield_hat"])
int_three <- unname(cf_orig["diff_log_yield_hat:meal_typeThree"])
total_three <- base_none + int_three
differential_alt <- total_three - base_none
log_msg(sprintf(
  "Bootstrap estimand verification: base (0-meal)=%.4f; interaction (3-meal)=%.4f; total (3-meal)=%.4f; differential (3-0)=%.4f; extract_m5_diff=%.4f",
  base_none, int_three, total_three, differential_alt, orig_diff$differential
))
stopifnot(
  abs(int_three - differential_alt) < 1e-8,
  abs(orig_diff$differential - int_three) < 1e-8
)

boot_diffs <- numeric(B_BOOT)
boot_ok    <- 0L

fs_fml <- as.formula(paste(
  "diff_log_yield ~", paste(CLIM_VARS, collapse = " + "), "|", FE_YIELD
))
ss_fml <- diff_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal_type + gender |
  year + District^growing_season

df_ts_base <- df_ts %>% select(-any_of("diff_log_yield_hat"))

for (b in seq_len(B_BOOT)) {
  samp <- sample(districts_yield, length(districts_yield), replace = TRUE)

  boot_yield <- bind_rows(lapply(seq_along(samp), function(i) {
    df_yield %>%
      filter(district == samp[i]) %>%
      mutate(boot_draw = i)
  }))

  m_fs <- tryCatch(
    feols(fs_fml, data = boot_yield, cluster = ~district, warn = FALSE, notes = FALSE),
    error = function(e) NULL
  )
  if (is.null(m_fs)) next

  boot_yield$yield_hat <- predict(m_fs, newdata = boot_yield)
  boot_hat <- boot_yield %>%
    filter(year >= 2017) %>%
    rename(District = district, growing_season = season) %>%
    arrange(boot_draw, District, growing_season, year) %>%
    group_by(boot_draw, District, growing_season) %>%
    mutate(diff_log_yield_hat = yield_hat - lag(yield_hat)) %>%
    ungroup() %>%
    filter(!is.na(diff_log_yield_hat)) %>%
    select(boot_draw, District, growing_season, year, diff_log_yield_hat)

  boot_wage <- bind_rows(lapply(seq_along(samp), function(i) {
    df_ts_base %>%
      filter(District == samp[i]) %>%
      mutate(boot_draw = i)
  })) %>%
    inner_join(boot_hat, by = c("boot_draw", "District", "growing_season", "year"))

  if (nrow(boot_wage) < 500) next

  m_ss <- tryCatch(
    feols(ss_fml, data = boot_wage, cluster = ~District, warn = FALSE, notes = FALSE),
    error = function(e) NULL
  )
  if (is.null(m_ss)) next

  cf <- coef(m_ss)
  if (!"diff_log_yield_hat:meal_typeThree" %in% names(cf)) next

  boot_diffs[b] <- cf["diff_log_yield_hat:meal_typeThree"]
  boot_ok <- boot_ok + 1L

  if (b %% 1000 == 0) log_msg("  Bootstrap draw ", b, "/", B_BOOT)
}

boot_valid <- boot_diffs[is.finite(boot_diffs)]
boot_se  <- sd(boot_valid)
boot_ci  <- quantile(boot_valid, c(0.025, 0.975), na.rm = TRUE)
boot_p   <- 2 * min(mean(boot_valid <= 0), mean(boot_valid >= 0))

boot_tbl <- data.frame(
  estimate = orig_diff$differential,
  orig_se = orig_diff$se,
  orig_p = orig_diff$p,
  boot_se = boot_se,
  boot_ci_lo = boot_ci[1],
  boot_ci_hi = boot_ci[2],
  boot_p = boot_p,
  n_boot_success = length(boot_valid),
  n_boot_requested = B_BOOT,
  se_ratio_boot_to_orig = boot_se / orig_diff$se,
  stringsAsFactors = FALSE
)
write_csv(boot_tbl, file.path(OUT, "tables/bootstrap_differential_fd.csv"))

log_msg(sprintf("  Original differential: %.1f (SE %.1f, p=%s)",
                orig_diff$differential, orig_diff$se, fmt_p(orig_diff$p)))
log_msg(sprintf("  Bootstrap SE: %.1f | 95%% CI [%.1f, %.1f] | boot p=%s | successful draws=%d",
                boot_se, boot_ci[1], boot_ci[2], fmt_p(boot_p), length(boot_valid)))
log_msg(sprintf("  SE ratio (bootstrap / original): %.2f", boot_tbl$se_ratio_boot_to_orig))

q4_se_similar <- boot_tbl$se_ratio_boot_to_orig < 1.25

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 6 — Summary diagnostic table                                         ##
## ════════════════════════════════════════════════════════════════════════════ ##
log_msg("\nSTEP 6 — Summary diagnostic table")

sec1 <- fs_summary %>%
  transmute(
    section = "First stage F-statistics",
    row = sample,
    partial_f = round(partial_f, 2),
    n = nobs,
    clears_f10 = ifelse(clears_f10, "Yes", "No"),
    extra = paste0("within R2=", round(within_r_squared, 3))
  )

sec2 <- rf_meal %>%
  transmute(
    section = "Reduced form by meal",
    row = meal_label,
    gdd_coef = round(gdd_coef, 2),
    gdd_p = fmt_p(gdd_p),
    edd_coef = round(edd_coef, 2),
    edd_p = fmt_p(edd_p),
    extra = paste0("GDD 1SD effect=", round(gdd_effect_1sd, 1), " BDT")
  )

sec3 <- rf_int_pooled$interactions %>%
  filter(term %in% c("diff_gdd_10_30:meal3", "diff_edd_30:meal3")) %>%
  transmute(
    section = "RF interaction (pooled)",
    row = term,
    coef = round(estimate, 2),
    se = round(std_error, 2),
    p = fmt_p(p_value),
    extra = paste0("Joint F p=", fmt_p(unique(joint_f_p)))
  )

sec4 <- boot_tbl %>%
  transmute(
    section = "Bootstrap comparison",
    row = "Attached/casual differential",
    orig_se = round(orig_se, 1),
    boot_se = round(boot_se, 1),
    orig_p = fmt_p(orig_p),
    boot_p = fmt_p(boot_p),
    extra = paste0("95% CI [", round(boot_ci[1], 1), ", ", round(boot_ci[2], 1), "]")
  )

write_csv(sec1, file.path(OUT, "tables/summary_section1_first_stage.csv"))
write_csv(sec2, file.path(OUT, "tables/summary_section2_rf_by_meal.csv"))
write_csv(sec3, file.path(OUT, "tables/summary_section3_rf_interaction.csv"))
write_csv(sec4, file.path(OUT, "tables/summary_section4_bootstrap.csv"))

## Diagnostic Q&A
diag_answers <- data.frame(
  question = c(
    "Q1: Is first stage stronger in Boro than other seasons?",
    "Q2: Does RF show heterogeneous meal effects despite pooled insignificance?",
    "Q3: Is RF interaction (GDD×meal3 or EDD×meal3) significant?",
    "Q4: Is bootstrap SE similar to generated-regressor SE?"
  ),
  answer = c(
    ifelse(q1_answer,
           sprintf("YES — Boro F=%.1f vs Aus=%.1f, Aman=%.1f", boro_f, aus_f, aman_f),
           sprintf("NO — Boro F=%.1f vs Aus=%.1f, Aman=%.1f", boro_f, aus_f, aman_f)),
    ifelse(q2_pattern_gdd || q2_pattern_edd,
           "YES — 3-meal climate coefficients more negative than 0-meal (heterogeneous RF)",
           "MIXED/NO — pattern not clear in all climate variables"),
    ifelse(q3_sig_pooled || q3_sig_boro,
           sprintf("PARTIAL — pooled joint p=%s; Boro joint p=%s; GDD×meal3 p=%s; EDD×meal3 p=%s",
                   fmt_p(rf_int_pooled$joint$joint_f_p),
                   fmt_p(rf_int_boro$joint$joint_f_p),
                   fmt_p(gdd_int$p_value), fmt_p(edd_int$p_value)),
           sprintf("NO at 5%% — pooled joint p=%s; GDD×meal3 p=%s; EDD×meal3 p=%s",
                   fmt_p(rf_int_pooled$joint$joint_f_p),
                   fmt_p(gdd_int$p_value), fmt_p(edd_int$p_value))),
    ifelse(q4_se_similar,
           sprintf("YES — bootstrap SE %.1f vs original %.1f (ratio %.2f)",
                   boot_se, orig_diff$se, boot_tbl$se_ratio_boot_to_orig),
           sprintf("NO — bootstrap SE %.1f > original %.1f (ratio %.2f); inference may be overstated",
                   boot_se, orig_diff$se, boot_tbl$se_ratio_boot_to_orig))
  ),
  stringsAsFactors = FALSE
)
write_csv(diag_answers, file.path(OUT, "tables/diagnostic_questions.csv"))

strongest <- if (q3_sig_boro) {
  "Boro reduced-form interaction (joint p<0.05) — season-specific RF without yield intermediate"
} else if (q4_se_similar && orig_diff$p < 0.05) {
  "Two-stage M5 with joint bootstrap validation — generated-regressor inference approximately correct"
} else if (q2_pattern_gdd || q2_pattern_edd) {
  "Heterogeneous reduced form by meal category — pooled insignificance is cancellation artifact"
} else if (q1_answer) {
  "Boro-specific analyses — relatively strongest first stage (F=3.0) though below F=10 threshold"
} else {
  "Mixed — identification relies on meal heterogeneity and Boro subsample with weak pooled FS"
}

log_msg("\n=== DIAGNOSTIC ANSWERS ===")
for (i in seq_len(nrow(diag_answers))) {
  log_msg(diag_answers$question[i])
  log_msg("  ", diag_answers$answer[i])
}
log_msg("\nStrongest identification approach: ", strongest)

summary_html <- bind_rows(
  fs_summary %>% mutate(table = "First stage") %>%
    select(table, sample, partial_f, nobs, clears_f10),
  rf_meal %>% mutate(table = "RF by meal") %>%
    select(table, meal_label, gdd_coef, gdd_p, edd_coef, edd_p),
  rf_int_pooled$interactions %>% mutate(table = "RF interaction") %>%
    select(table, term, estimate, std_error, p_value),
  boot_tbl %>% mutate(table = "Bootstrap")
) %>%
  kbl(caption = "Identification diagnostics summary") %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover")) %>%
  footnote(general = paste(
    "First stage: Δlog(yield) on Δclimate | year + district×season FE.",
    "Reduced form: Δwage on Δclimate | year + District×growing_season FE.",
    "Partial F tests joint significance of all four climate instruments.",
    "Conventional weak-IV threshold: F > 10; Stock-Yogo 10% (1 endog, 4 instruments) ≈ 16.85.",
    "Bootstrap: B=500 district cluster resamples with joint FS+SS.",
    "Strongest approach:", strongest
  ))
writeLines(summary_html, file.path(OUT, "tables/identification_diagnostics_summary.html"))

writeLines(c(
  "IDENTIFICATION DIAGNOSTICS SUMMARY",
  "==================================",
  diag_answers$question,
  diag_answers$answer,
  "",
  paste("Strongest identification approach:", strongest)
), file.path(OUT, "summary/diagnostic_conclusions.txt"))

log_msg("\nDone. Outputs: ", OUT)
