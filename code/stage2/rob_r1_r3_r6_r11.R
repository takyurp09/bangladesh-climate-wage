## rob_r1_r3_r6_r11.R
## Missing robustness tests: R1 (controls), R3 (levels), R6 (split), R11 (nominal)
## ALL feols calls use the +: form (no meal_type main effects), matching M5 baseline.
## Differential effect  = interaction coef (β_Three), SE and p from interaction term.
## Total effect         = β_None + β_Three; SE via delta method (includes cov term).
## SE and p from interaction term; cluster = District.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(fixest)
  library(knitr)
  library(kableExtra)
})

ROOT    <- here::here()
out_tbl <- file.path(ROOT, "output/stage2/tables")
dir.create(out_tbl, recursive = TRUE, showWarnings = FALSE)

cat("=== Loading main wage panel ===\n")
df <- read_csv(file.path(ROOT, "data/Regression_data/df_2_merged_v2.csv"),
               show_col_types = FALSE) %>%
  mutate(
    meal_type = relevel(factor(meal_type), ref = "None"),
    gender    = relevel(factor(gender),    ref = "Female"),
    growing_season = factor(growing_season, levels = c("Boro", "Aus", "Aman"))
  )
cat(sprintf("  N = %d, years = %s\n", nrow(df),
            paste(sort(unique(df$year)), collapse=", ")))

## ── Delta-method total effect helper ─────────────────────────────────────── ##
## Returns total = β_base + β_int, SE via full delta method (includes cov term)
## p from two-sided t with df = G-1 (clusters - 1)
dm_total_effect <- function(model, base_nm, int_nm, G = NULL) {
  null <- list(dm_coef=NA_real_, dm_se=NA_real_, dm_p=NA_real_)
  if (is.null(model)) return(null)
  cf <- coef(model); V <- vcov(model)
  if (!base_nm %in% names(cf) || !int_nm %in% names(cf)) return(null)
  est  <- cf[base_nm] + cf[int_nm]
  var_e <- V[base_nm,base_nm] + V[int_nm,int_nm] + 2*V[base_nm,int_nm]
  se_e  <- sqrt(max(var_e, 0))
  if (is.null(G)) G <- length(unique(model$fixef_id[[1]]))
  t_e  <- est / se_e
  p_e  <- 2 * pt(-abs(t_e), df = G - 1L)
  list(dm_coef = as.numeric(est), dm_se = as.numeric(se_e), dm_p = as.numeric(p_e))
}

## ── Helper: extract Three-meal differential AND total from M5 ────────────── ##
## Uses +: form (interaction only, no meal_type main effects)
extract_three <- function(m, base = "diff_log_yield_hat", G = NULL) {
  null <- list(coef=NA_real_, se=NA_real_, pval=NA_real_, n=NA_integer_,
               dm_coef=NA_real_, dm_se=NA_real_, dm_p=NA_real_)
  if (is.null(m)) return(null)
  cf <- coef(m); sv <- se(m); pv <- pvalue(m)
  inter_name <- paste0(base, ":meal_typeThree")
  if (!inter_name %in% names(cf)) {
    if (!base %in% names(cf)) return(null)
    return(c(list(coef=cf[[base]], se=sv[[base]], pval=pv[[base]], n=nobs(m)),
             list(dm_coef=cf[[base]], dm_se=sv[[base]], dm_p=pv[[base]])))
  }
  if (!base %in% names(cf)) return(null)
  dm <- dm_total_effect(m, base, inter_name, G)
  c(
    list(coef = cf[base] + cf[inter_name],  ## total point estimate
         se   = sv[inter_name],             ## differential SE
         pval = pv[inter_name],             ## differential p
         n    = nobs(m)),
    dm
  )
}

## ── Helper: save tex + html (kableExtra) ─────────────────────────────────── ##
save_tex <- function(df_tbl, stem, title, note_text) {
  kt_tex <- kable(df_tbl, format = "latex", booktabs = TRUE,
                  caption = title, align = c("l", rep("r", ncol(df_tbl)-1))) %>%
    kable_styling(latex_options = c("hold_position", "scale_down")) %>%
    footnote(general = note_text, escape = FALSE,
             general_title = "\\\\textit{Note:} ",
             footnote_as_chunk = TRUE)
  kableExtra::save_kable(kt_tex, file.path(out_tbl, paste0(stem, ".tex")))
  cat(sprintf("  Saved %s.tex\n", stem))
}

## ════════════════════════════════════════════════════════════════════════════ ##
## R1: Add controls (irrigation share, land holdings, crop intensity)          ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("\n=== R1: Add controls ===\n")

census <- read_csv(file.path(ROOT, "data/agricultural_census_2019.csv"),
                   show_col_types = FALSE) %>%
  mutate(
    irrigation_share = Net_Irrigated_Area / Net_Cultivated_Area,
    avg_holdings     = Number_of_Holdings / Net_Cultivated_Area * 1000,  ## per 1000 ha
    crop_intensity   = Intensity_of_Cropping
  ) %>%
  select(District, irrigation_share, avg_holdings, crop_intensity)

## Fix known district name mismatches (standardise to wage panel names)
## Check overlap
df_districts  <- unique(df$District)
cen_districts <- unique(census$District)
unmatched <- setdiff(cen_districts, df_districts)
if (length(unmatched) > 0)
  cat(sprintf("  Census districts not in wage panel: %s\n",
              paste(unmatched, collapse=", ")))

df_r1 <- df %>%
  left_join(census, by = "District")

n_missing_irr <- sum(is.na(df_r1$irrigation_share))
cat(sprintf("  After merge: N=%d, missing irrigation_share=%d\n",
            nrow(df_r1), n_missing_irr))

## Drop rows with missing controls
df_r1_clean <- df_r1 %>%
  filter(!is.na(irrigation_share) & !is.na(avg_holdings) & !is.na(crop_intensity))
cat(sprintf("  After dropping missing: N=%d\n", nrow(df_r1_clean)))

m_r1 <- feols(diff_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal_type + gender +
                irrigation_share + avg_holdings + crop_intensity |
                year + District^growing_season,
              data = df_r1_clean, cluster = ~District,
              warn = FALSE, notes = FALSE)

G_r1 <- n_distinct(df_r1_clean$District)
te_r1 <- extract_three(m_r1, G = G_r1)
cat(sprintf("  R1 Three-meal: diff=%.4f  SE=%.4f  p=%.4f  |  total=%.4f  dm_SE=%.4f  dm_p=%.4f  N=%d\n",
            te_r1$coef - coef(m_r1)["diff_log_yield_hat"],
            te_r1$se, te_r1$pval,
            te_r1$dm_coef, te_r1$dm_se, te_r1$dm_p,
            te_r1$n))

## Build display table
stars_fn <- function(p) ifelse(p<0.01,"***",ifelse(p<0.05,"**",ifelse(p<0.10,"*","")))

r1_tbl <- data.frame(
  Specification = "R1: Add controls",
  `Three-meal coef`  = sprintf("%.2f%s", te_r1$coef, stars_fn(te_r1$pval)),
  `Three-meal SE`    = sprintf("%.2f", te_r1$se),
  `Three-meal p`     = sprintf("%.3f", te_r1$pval),
  N = format(te_r1$n, big.mark=","),
  FE = "yr+Dist^season",
  Cluster = "District",
  check.names = FALSE, stringsAsFactors = FALSE
)

save_tex(r1_tbl, "rob_r1_controls",
         "ROB-R1: Three-Meal Pass-Through with Additional Controls",
         paste("M5 specification with three census-based controls: irrigation share",
               "(Net irrigated/net cultivated area, 2019 agricultural census),",
               "land holdings per 1{,}000 ha, and cropping intensity.",
               "Total effect = baseline + interaction; SE and $p$ from interaction term.",
               "FE: year + District$\\\\times$growing season.",
               "$^{*}$p$<$0.1, $^{**}$p$<$0.05, $^{***}$p$<$0.01."))

## ════════════════════════════════════════════════════════════════════════════ ##
## R3: Levels specification                                                     ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("\n=== R3: Levels specification ===\n")

## Reconstruct log_yield_hat in levels from cumsum of diff_log_yield_hat
## within each (District, growing_season) group, sorted by year.
## With year + District^Season FE, the baseline level is absorbed.
df_r3 <- df %>%
  arrange(District, growing_season, year) %>%
  group_by(District, growing_season) %>%
  mutate(
    log_yield_hat = cumsum(diff_log_yield_hat)   ## running sum, starts at 0 in 2017
  ) %>%
  ungroup()

## NOTE: R3 uses * form (not +: form). In a LEVELS regression, meal_type main
## effects must be included — first-differencing does not remove them.
## The +: form without meal_type intercepts gives biased results in levels
## because District^Season FE does NOT absorb meal_type level wage differences.
m_r3 <- feols(real_wage ~ log_yield_hat * meal_type + gender |
                year + District^growing_season,
              data = df_r3, cluster = ~District,
              warn = FALSE, notes = FALSE)

G_r3 <- n_distinct(df_r3$District)
te_r3 <- extract_three(m_r3, base = "log_yield_hat", G = G_r3)
cat(sprintf("  R3 Three-meal: diff=%.4f  SE=%.4f  p=%.4f  |  total=%.4f  dm_SE=%.4f  dm_p=%.4f  N=%d\n",
            te_r3$coef - coef(m_r3)["log_yield_hat"],  ## interaction coef
            te_r3$se, te_r3$pval,
            te_r3$dm_coef, te_r3$dm_se, te_r3$dm_p,
            te_r3$n))

r3_tbl <- data.frame(
  Specification = "R3: Levels (not FD)",
  `Three-meal coef`  = sprintf("%.2f%s", te_r3$coef, stars_fn(te_r3$pval)),
  `Three-meal SE`    = sprintf("%.2f", te_r3$se),
  `Three-meal p`     = sprintf("%.3f", te_r3$pval),
  N = format(te_r3$n, big.mark=","),
  FE = "yr+Dist^season",
  Cluster = "District",
  check.names = FALSE, stringsAsFactors = FALSE
)

save_tex(r3_tbl, "rob_r3_levels",
         "ROB-R3: Levels Specification (Wage in Levels, Yield Hat Cumulative)",
         paste("Levels specification: outcome is real wage (BDT/day, not first-differenced).",
               "Yield hat is constructed as the running sum of",
               "$\\\\Delta\\\\log(\\\\widehat{\\\\text{yield}})$ within each",
               "District$\\\\times$Season group (baseline year 2017 = 0).",
               "FE: year + District$\\\\times$growing season absorb level differences.",
               "Baseline FD specification uses first differences of the same variables.",
               "$^{*}$p$<$0.1, $^{**}$p$<$0.05, $^{***}$p$<$0.01."))

## ════════════════════════════════════════════════════════════════════════════ ##
## R6: Early vs late period split                                               ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("\n=== R6: Early/late split ===\n")

df_early <- df %>% filter(year <= 2020)
df_late  <- df %>% filter(year >= 2021)
cat(sprintf("  Early: N=%d, years=%s\n", nrow(df_early),
            paste(sort(unique(df_early$year)), collapse=",")))
cat(sprintf("  Late:  N=%d, years=%s\n", nrow(df_late),
            paste(sort(unique(df_late$year)), collapse=",")))

m_r6_early <- tryCatch(
  feols(diff_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal_type + gender |
          year + District^growing_season,
        data = df_early, cluster = ~District, warn = FALSE, notes = FALSE),
  error = function(e) { cat("  R6 early error:", conditionMessage(e), "\n"); NULL }
)

m_r6_late <- tryCatch(
  feols(diff_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal_type + gender |
          year + District^growing_season,
        data = df_late, cluster = ~District, warn = FALSE, notes = FALSE),
  error = function(e) { cat("  R6 late error:", conditionMessage(e), "\n"); NULL }
)

G_r6a <- n_distinct(df_early$District)
G_r6b <- n_distinct(df_late$District)
te_r6e <- if (!is.null(m_r6_early)) extract_three(m_r6_early, G = G_r6a) else
  list(coef=NA, se=NA, pval=NA, n=NA, dm_coef=NA, dm_se=NA, dm_p=NA)
te_r6l <- if (!is.null(m_r6_late)) extract_three(m_r6_late, G = G_r6b) else
  list(coef=NA, se=NA, pval=NA, n=NA, dm_coef=NA, dm_se=NA, dm_p=NA)

cat(sprintf("  R6 Early Three-meal: diff=%.4f  SE=%.4f  p=%.4f  |  total=%.4f  dm_SE=%.4f  dm_p=%.4f  N=%s\n",
            te_r6e$coef - coef(m_r6_early)["diff_log_yield_hat"],
            te_r6e$se, te_r6e$pval,
            te_r6e$dm_coef, te_r6e$dm_se, te_r6e$dm_p,
            ifelse(is.na(te_r6e$n), "NA", as.character(te_r6e$n))))
cat(sprintf("  R6 Late  Three-meal: diff=%.4f  SE=%.4f  p=%.4f  |  total=%.4f  dm_SE=%.4f  dm_p=%.4f  N=%s\n",
            te_r6l$coef - coef(m_r6_late)["diff_log_yield_hat"],
            te_r6l$se, te_r6l$pval,
            te_r6l$dm_coef, te_r6l$dm_se, te_r6l$dm_p,
            ifelse(is.na(te_r6l$n), "NA", as.character(te_r6l$n))))

r6_tbl <- data.frame(
  Panel      = c("Panel A: 2017-2020", "Panel B: 2021-2023"),
  Specification = c("R6a: Early period", "R6b: Late period"),
  `Three-meal coef` = c(
    sprintf("%.2f%s", te_r6e$coef, stars_fn(te_r6e$pval)),
    sprintf("%.2f%s", te_r6l$coef, stars_fn(te_r6l$pval))
  ),
  `Three-meal SE`   = c(sprintf("%.2f", te_r6e$se),  sprintf("%.2f", te_r6l$se)),
  `Three-meal p`    = c(sprintf("%.3f", te_r6e$pval), sprintf("%.3f", te_r6l$pval)),
  N = c(format(te_r6e$n, big.mark=","), format(te_r6l$n, big.mark=",")),
  FE = c("yr+Dist^season", "yr+Dist^season"),
  check.names = FALSE, stringsAsFactors = FALSE
)

save_tex(r6_tbl, "rob_r6_split",
         "ROB-R6: Early vs.~Late Period Split (M5 Three-Meal Pass-Through)",
         paste("M5 specification estimated on subsamples.",
               "Panel A: 2017--2020 (4 years).",
               "Panel B: 2021--2023 (3 years).",
               "FE: year + District$\\\\times$growing season.",
               "$^{*}$p$<$0.1, $^{**}$p$<$0.05, $^{***}$p$<$0.01."))

## ════════════════════════════════════════════════════════════════════════════ ##
## R11: Nominal wage outcome                                                    ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("\n=== R11: Nominal wage ===\n")

## Load raw wage data to get nominal wages
## wage_by_growing_season.csv has: District, growing_season, gender, meal_type,
## year, wage (nominal), cpi, real_wage
wage_raw <- read_csv(file.path(ROOT, "data/Regression_data/wage_by_growing_season.csv"),
                     show_col_types = FALSE) %>%
  filter(!is.na(wage)) %>%
  mutate(
    meal_type      = relevel(factor(meal_type), ref = "None"),
    gender         = relevel(factor(gender),    ref = "Female"),
    growing_season = factor(growing_season, levels = c("Boro", "Aus", "Aman"))
  )

cat(sprintf("  wage_raw: N=%d, cols=%s\n", nrow(wage_raw),
            paste(names(wage_raw), collapse=",")))

## Aggregate nominal wage to district-season-gender-meal-year level
## (matching the structure of df_2_merged_v2)
wage_agg <- wage_raw %>%
  group_by(District, growing_season, gender, meal_type, year) %>%
  summarise(nominal_wage = mean(wage, na.rm = TRUE), .groups = "drop")

## Compute first differences within (District, growing_season, gender, meal_type)
wage_fd <- wage_agg %>%
  arrange(District, growing_season, gender, meal_type, year) %>%
  group_by(District, growing_season, gender, meal_type) %>%
  mutate(diff_nominal_wage = nominal_wage - lag(nominal_wage)) %>%
  ungroup() %>%
  filter(!is.na(diff_nominal_wage))

cat(sprintf("  FD nominal wage: N=%d\n", nrow(wage_fd)))

## Merge with main panel to get diff_log_yield_hat
df_r11 <- df %>%
  select(District, growing_season, gender, meal_type, year, diff_log_yield_hat) %>%
  inner_join(wage_fd %>% select(District, growing_season, gender, meal_type, year,
                                 diff_nominal_wage),
             by = c("District", "growing_season", "gender", "meal_type", "year")) %>%
  mutate(
    meal_type      = relevel(factor(meal_type), ref = "None"),
    gender         = relevel(factor(gender),    ref = "Female"),
    growing_season = factor(growing_season, levels = c("Boro", "Aus", "Aman"))
  )

cat(sprintf("  After merge: N=%d\n", nrow(df_r11)))

m_r11 <- tryCatch(
  feols(diff_nominal_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal_type + gender |
          year + District^growing_season,
        data = df_r11, cluster = ~District, warn = FALSE, notes = FALSE),
  error = function(e) { cat("  R11 error:", conditionMessage(e), "\n"); NULL }
)

if (!is.null(m_r11)) {
  G_r11 <- n_distinct(df_r11$District)
  te_r11 <- extract_three(m_r11, G = G_r11)
  cat(sprintf("  R11 Three-meal: diff=%.4f  SE=%.4f  p=%.4f  |  total=%.4f  dm_SE=%.4f  dm_p=%.4f  N=%d\n",
              te_r11$coef - coef(m_r11)["diff_log_yield_hat"],
              te_r11$se, te_r11$pval,
              te_r11$dm_coef, te_r11$dm_se, te_r11$dm_p,
              te_r11$n))
} else {
  te_r11 <- list(coef=NA, se=NA, pval=NA, n=NA, dm_coef=NA, dm_se=NA, dm_p=NA)
  cat("  R11: model failed\n")
}

r11_tbl <- data.frame(
  Specification = "R11: Nominal wage",
  `Three-meal coef`  = sprintf("%.2f%s", te_r11$coef, stars_fn(te_r11$pval)),
  `Three-meal SE`    = sprintf("%.2f", te_r11$se),
  `Three-meal p`     = sprintf("%.3f", te_r11$pval),
  N = ifelse(is.na(te_r11$n), "NA", format(te_r11$n, big.mark=",")),
  FE = "yr+Dist^season",
  Cluster = "District",
  check.names = FALSE, stringsAsFactors = FALSE
)

save_tex(r11_tbl, "rob_r11_nominal",
         "ROB-R11: Nominal Wage Outcome (Not CPI-Deflated)",
         paste("M5 specification with nominal wage (BDT/day, not CPI-deflated)",
               "as outcome, in first differences.",
               "Nominal wages are monthly averages from the BBS Agricultural Wage Survey,",
               "aggregated to district-season-gender-meal-year cells.",
               "FE: year + District$\\\\times$growing season.",
               "$^{*}$p$<$0.1, $^{**}$p$<$0.05, $^{***}$p$<$0.01."))

## ════════════════════════════════════════════════════════════════════════════ ##
## Summary                                                                      ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("\n=== SUMMARY ===\n")
fmt <- function(x, digits=3) ifelse(is.na(x), "NA", sprintf(paste0("%.",digits,"f"), x))
cat(sprintf("R1 (controls):  coef=%-10s  SE=%-8s  p=%s  N=%s\n",
    fmt(te_r1$coef,2), fmt(te_r1$se,2), fmt(te_r1$pval), te_r1$n))
cat(sprintf("R3 (levels):    coef=%-10s  SE=%-8s  p=%s  N=%s\n",
    fmt(te_r3$coef,2), fmt(te_r3$se,2), fmt(te_r3$pval), te_r3$n))
cat(sprintf("R6a (early):    coef=%-10s  SE=%-8s  p=%s  N=%s\n",
    fmt(te_r6e$coef,2), fmt(te_r6e$se,2), fmt(te_r6e$pval), te_r6e$n))
cat(sprintf("R6b (late):     coef=%-10s  SE=%-8s  p=%s  N=%s\n",
    fmt(te_r6l$coef,2), fmt(te_r6l$se,2), fmt(te_r6l$pval), te_r6l$n))
cat(sprintf("R11 (nominal):  coef=%-10s  SE=%-8s  p=%s  N=%s\n",
    fmt(te_r11$coef,2), fmt(te_r11$se,2), fmt(te_r11$pval), te_r11$n))

cat("\nAll tables saved to output/stage2/tables/\n")
