## 15_district_heterogeneity_sdg.R
## District heterogeneity in attached/casual wage pass-through differential
## World Bank Bangladesh Spatial Database (Zila-SDG-v2.xlsx)
## Output: output/stage2/district_heterogeneity/

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(readxl)
  library(purrr)
  library(fixest)
  library(ggplot2)
  library(kableExtra)
})

ROOT <- here::here()
OUT  <- file.path(ROOT, "output/stage2/district_heterogeneity")
WB   <- file.path(ROOT, "data/Zila-SDG-v2.xlsx")

for (d in c("tables", "figures", "logs"))
  dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(OUT, "logs/heterogeneity_log.txt")
log_con  <- file(log_path, open = "wt")
on.exit(close(log_con), add = TRUE)

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

source(file.path(ROOT, "code/utils/district_names.R"))

extract_wb_indicator <- function(sheet, indicator, var_name) {
  read_excel(WB, sheet = sheet) %>%
    filter(Indicator == indicator) %>%
    transmute(
      ZilaCode = as.integer(ZilaCode),
      Zila,
      !!var_name := as.numeric(Estimate)
    )
}

extract_meal_pass <- function(mod) {
  ct  <- summary(mod)$coeftable
  V   <- vcov(mod)
  bnm <- "diff_log_yield_hat"
  int <- "diff_log_yield_hat:meal_typeThree"
  df2 <- mod$nobs - length(coef(mod))

  b0 <- ct[bnm, "Estimate"]
  diff_est <- ct[int, "Estimate"]
  diff_se  <- ct[int, "Std. Error"]
  diff_p   <- ct[int, "Pr(>|t|)"]

  list(
    casual_coef     = b0,
    casual_se       = ct[bnm, "Std. Error"],
    casual_p        = ct[bnm, "Pr(>|t|)"],
    attached_coef   = b0 + diff_est,
    attached_se     = sqrt(V[bnm, bnm] + V[int, int] + 2 * V[bnm, int]),
    differential    = diff_est,
    differential_se = diff_se,
    differential_p  = diff_p,
    nobs            = nobs(mod)
  )
}

wald_linear <- function(mod, terms, weights) {
  b <- coef(mod); V <- vcov(mod)
  idx <- match(terms, names(b))
  L <- numeric(length(b)); L[idx] <- weights
  est <- sum(L * b)
  se  <- sqrt(as.numeric(t(L) %*% V %*% L))
  df2 <- mod$nobs - length(b)
  list(estimate = est, se = se, p = 2 * pt(-abs(est / se), df = df2))
}

run_split_heterogeneity <- function(df, split_var, split_label_high, split_label_low) {
  base_rhs <- "diff_log_yield_hat + diff_log_yield_hat:meal_type + gender"
  triple_nm <- paste0("diff_log_yield_hat:meal3:", split_var)

  m_high <- feols(
    as.formula(paste("diff_real_wage ~", base_rhs, "| year + District^growing_season")),
    data = df %>% filter(.data[[split_var]] == 1),
    cluster = ~District
  )
  m_low <- feols(
    as.formula(paste("diff_real_wage ~", base_rhs, "| year + District^growing_season")),
    data = df %>% filter(.data[[split_var]] == 0),
    cluster = ~District
  )
  m_pool <- feols(
    as.formula(paste(
      "diff_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal3 +",
      triple_nm, "+ gender | year + District^growing_season"
    )),
    data = df,
    cluster = ~District
  )

  r_high <- extract_meal_pass(m_high)
  r_low  <- extract_meal_pass(m_low)
  diff_est <- r_high$differential - r_low$differential

  diff_test <- if (triple_nm %in% names(coef(m_pool))) {
    ct <- summary(m_pool)$coeftable
    list(
      estimate = diff_est,
      se       = ct[triple_nm, "Std. Error"],
      p        = ct[triple_nm, "Pr(>|t|)"]
    )
  } else {
    list(
      estimate = diff_est,
      se       = sqrt(r_high$differential_se^2 + r_low$differential_se^2),
      p        = 2 * pt(-abs(diff_est / sqrt(r_high$differential_se^2 +
                                               r_low$differential_se^2)),
                        df = min(r_high$nobs, r_low$nobs))
    )
  }

  data.frame(
    split_var           = split_var,
    group_high_label    = split_label_high,
    group_low_label     = split_label_low,
    high_differential   = r_high$differential,
    high_se             = r_high$differential_se,
    high_p              = r_high$differential_p,
    high_n              = r_high$nobs,
    low_differential    = r_low$differential,
    low_se              = r_low$differential_se,
    low_p               = r_low$differential_p,
    low_n               = r_low$nobs,
    diff_high_minus_low = diff_test$estimate,
    diff_se             = diff_test$se,
    diff_p              = diff_test$p,
    stringsAsFactors    = FALSE
  )
}

log_msg("=== DISTRICT HETEROGENEITY: IMPLICIT CONTRACT ENVIRONMENT ===")
log_msg("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 1 — Extract World Bank district variables                               ##
## ════════════════════════════════════════════════════════════════════════════ ##
log_msg("\nSTEP 1 — Extract district-level heterogeneity variables")

poverty   <- extract_wb_indicator("SDG1",  "Poverty rate (Percent)",
                                  "poverty_rate")
agri_emp  <- extract_wb_indicator("SDG8",  "Primary employment: Agriculture (Percent)",
                                  "agri_employment_share")
low_edu   <- extract_wb_indicator("SDG4",  "Adults with incomplete primary (Percent)",
                                  "low_education_share")
dist_dh   <- extract_wb_indicator("Geography", "Distance to Dhaka (km)",
                                  "dist_dhaka")
mobile    <- extract_wb_indicator("SDG9",  "Owns mobile phone (Percent)",
                                  "mobile_ownership")
firms     <- extract_wb_indicator("Geography", "Density of Firms (Number per km2)",
                                  "firm_density")

wb_district <- poverty %>%
  full_join(select(agri_emp, -Zila), by = "ZilaCode") %>%
  full_join(select(low_edu, -Zila), by = "ZilaCode") %>%
  full_join(select(dist_dh, -Zila), by = "ZilaCode") %>%
  full_join(select(mobile, -Zila), by = "ZilaCode") %>%
  full_join(select(firms, -Zila), by = "ZilaCode")

bbs_districts <- read_csv(
  file.path(ROOT, "data/Regression_data/wage_by_growing_season.csv"),
  show_col_types = FALSE
) %>% pull(District) %>% unique() %>% sort()

wb_district <- wb_district %>%
  mutate(District = map_wb_district(Zila, bbs_districts))

unmatched_wb  <- wb_district %>% filter(is.na(District) | !District %in% bbs_districts)
unmatched_bbs <- setdiff(bbs_districts, wb_district$District)

log_msg("  WB districts extracted: ", nrow(wb_district))
if (nrow(unmatched_wb) > 0)
  log_msg("  FLAG: Unmatched WB districts: ", paste(unmatched_wb$Zila, collapse = ", "))
if (length(unmatched_bbs) > 0)
  log_msg("  FLAG: Unmatched BBS districts: ", paste(unmatched_bbs, collapse = ", "))

write_csv(wb_district, file.path(OUT, "tables/wb_district_variables.csv"))

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 2 — Binary splits and composite index                                 ##
## ════════════════════════════════════════════════════════════════════════════ ##
log_msg("\nSTEP 2 — Construct heterogeneity split variables")

het_vars <- wb_district %>%
  filter(District %in% bbs_districts)

med <- het_vars %>%
  summarise(
    poverty_rate          = median(poverty_rate, na.rm = TRUE),
    agri_employment_share = median(agri_employment_share, na.rm = TRUE),
    low_education_share   = median(low_education_share, na.rm = TRUE),
    dist_dhaka            = median(dist_dhaka, na.rm = TRUE),
    mobile_ownership      = median(mobile_ownership, na.rm = TRUE),
    firm_density          = median(firm_density, na.rm = TRUE)
  )

het_vars <- het_vars %>%
  mutate(
    high_poverty  = as.integer(poverty_rate >= med$poverty_rate),
    high_agri     = as.integer(agri_employment_share >= med$agri_employment_share),
    high_low_edu  = as.integer(low_education_share >= med$low_education_share),
    high_remote   = as.integer(dist_dhaka >= med$dist_dhaka),
    low_mobile    = as.integer(mobile_ownership < med$mobile_ownership),
    low_firms     = as.integer(firm_density < med$firm_density),
    implicit_contract_environment = high_poverty + high_agri + high_low_edu +
      high_remote + low_mobile + low_firms
  )

write_csv(het_vars, file.path(OUT, "tables/district_heterogeneity_splits.csv"))
log_msg("  Medians — poverty: ", round(med$poverty_rate, 2),
        " | agri emp: ", round(med$agri_employment_share, 2),
        " | low edu: ", round(med$low_education_share, 2),
        " | dist Dhaka: ", round(med$dist_dhaka, 1), " km",
        " | mobile: ", round(med$mobile_ownership, 2),
        " | firm density: ", round(med$firm_density, 4))
log_msg("  Composite index range: ", min(het_vars$implicit_contract_environment),
        "–", max(het_vars$implicit_contract_environment))

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 3 — Main heterogeneity regressions (six binary splits)                  ##
## ════════════════════════════════════════════════════════════════════════════ ##
log_msg("\nSTEP 3 — Binary split heterogeneity regressions")

df <- read_csv(
  file.path(ROOT, "data/Regression_data/df_2_merged_v2.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    meal_type = relevel(factor(meal_type), ref = "None"),
    gender    = relevel(factor(gender),    ref = "Female")
  ) %>%
  left_join(
    het_vars %>%
      select(District, high_poverty, high_agri, high_low_edu, high_remote,
             low_mobile, low_firms, implicit_contract_environment),
    by = "District"
  ) %>%
  mutate(meal3 = as.integer(meal_type == "Three")) %>%
  filter(!is.na(diff_log_yield_hat))

split_specs <- tibble(
  split_var         = c("high_poverty", "high_agri", "high_low_edu",
                      "high_remote", "low_mobile", "low_firms"),
  variable_label    = c("Poverty rate", "Agricultural employment share",
                        "Low education share", "Distance to Dhaka",
                        "Low mobile ownership", "Low firm density"),
  group_high_label  = c("High poverty", "High agri employment",
                        "High low-education", "Remote (far from Dhaka)",
                        "Low mobile ownership", "Low firm density"),
  group_low_label   = c("Low poverty", "Low agri employment",
                        "Low education", "Near Dhaka",
                        "High mobile ownership", "High firm density"),
  theory_direction  = "Larger |differential| in constrained group (more negative)"
)

split_results <- pmap_dfr(
  split_specs,
  function(split_var, variable_label, group_high_label, group_low_label, theory_direction) {
    res <- run_split_heterogeneity(df, split_var, group_high_label, group_low_label)
    res$variable_label <- variable_label
    res$theory_direction <- theory_direction
    res
  }
)

split_results <- split_results %>%
  mutate(
    direction_consistent = high_differential < low_differential,
    high_sig             = sig_flag(high_p),
    low_sig              = sig_flag(low_p),
    diff_sig             = sig_flag(diff_p),
    interpretation = ifelse(
      direction_consistent,
      ifelse(diff_p < 0.05,
             "Consistent with theory (significant)",
             "Directionally consistent with theory (imprecise)"),
      "Not consistent with theory"
    )
  )

write_csv(split_results, file.path(OUT, "tables/split_heterogeneity_results.csv"))

for (i in seq_len(nrow(split_results))) {
  r <- split_results[i, ]
  log_msg(sprintf("  %s:", r$variable_label))
  log_msg(sprintf("    High: %.1f (SE %.1f, p=%s) | Low: %.1f (SE %.1f, p=%s) | Diff: %.1f (p=%s) [%s]",
                  r$high_differential, r$high_se, fmt_p(r$high_p),
                  r$low_differential, r$low_se, fmt_p(r$low_p),
                  r$diff_high_minus_low, fmt_p(r$diff_p), r$interpretation))
}

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 4 — Composite index analysis                                          ##
## ════════════════════════════════════════════════════════════════════════════ ##
log_msg("\nSTEP 4 — Composite index analysis")

df <- df %>%
  mutate(
    ic_group = case_when(
      implicit_contract_environment <= 2 ~ "Low (0-2)",
      implicit_contract_environment <= 4 ~ "Medium (3-4)",
      TRUE ~ "High (5-6)"
    )
  )

composite_specs <- c("Low (0-2)", "Medium (3-4)", "High (5-6)")
composite_results <- lapply(composite_specs, function(g) {
  sub <- df %>% filter(ic_group == g)
  m <- feols(
    diff_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal_type + gender |
      year + District^growing_season,
    data = sub, cluster = ~District
  )
  r <- extract_meal_pass(m)
  data.frame(
    ic_group = g,
    n_districts = n_distinct(sub$District),
    differential = r$differential,
    differential_se = r$differential_se,
    differential_p = r$differential_p,
    nobs = r$nobs,
    sig_status = sig_flag(r$differential_p)
  )
}) %>% bind_rows()

write_csv(composite_results, file.path(OUT, "tables/composite_index_groups.csv"))

for (i in seq_len(nrow(composite_results))) {
  r <- composite_results[i, ]
  log_msg(sprintf("  %s (n=%d dist): differential = %.1f (SE %.1f, p=%s) [%s]",
                  r$ic_group, r$n_districts, r$differential, r$differential_se,
                  fmt_p(r$differential_p), r$sig_status))
}

m_ic_cont <- feols(
  diff_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal3 +
    diff_log_yield_hat:meal3:implicit_contract_environment + gender |
    year + District^growing_season,
  data = df,
  cluster = ~District
)

triple_ic <- "diff_log_yield_hat:meal3:implicit_contract_environment"
ic_triple <- summary(m_ic_cont)$coeftable[triple_ic, ]

log_msg(sprintf("  Continuous triple interaction: %.1f (SE %.1f, p=%s)",
                ic_triple["Estimate"], ic_triple["Std. Error"],
                fmt_p(ic_triple["Pr(>|t|)"])))
log_msg("  (Negative coef => differential more negative as composite index rises)")

ic_continuous <- data.frame(
  term = triple_ic,
  coefficient = ic_triple["Estimate"],
  std_error = ic_triple["Std. Error"],
  p_value = ic_triple["Pr(>|t|)"],
  nobs = nobs(m_ic_cont)
)
write_csv(ic_continuous, file.path(OUT, "tables/composite_index_continuous.csv"))

monotone <- all(diff(composite_results$differential) < 0)
log_msg("  Monotonic increase in |differential| from low to high index: ",
        ifelse(monotone, "YES", "NO"))

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 5 — Summary table and coefficient plot                                ##
## ════════════════════════════════════════════════════════════════════════════ ##
log_msg("\nSTEP 5 — Summary outputs")

summary_tbl <- split_results %>%
  transmute(
    variable = variable_label,
    theoretical_direction = theory_direction,
    high_group_diff = round(high_differential, 1),
    high_se = round(high_se, 1),
    high_p = fmt_p(high_p),
    low_group_diff = round(low_differential, 1),
    low_se = round(low_se, 1),
    low_p = fmt_p(low_p),
    difference_high_minus_low = round(diff_high_minus_low, 1),
    diff_se = round(diff_se, 1),
    diff_p = fmt_p(diff_p),
    direction_consistent_with_theory = ifelse(direction_consistent, "Yes", "No"),
    notes = interpretation
  )

composite_row <- composite_results %>%
  transmute(
    variable = paste0("Composite index: ", ic_group),
    theoretical_direction = "Larger |differential| in high-index districts",
    high_group_diff = round(differential, 1),
    high_se = round(differential_se, 1),
    high_p = fmt_p(differential_p),
    low_group_diff = NA_real_,
    low_se = NA_real_,
    low_p = NA_character_,
    difference_high_minus_low = NA_real_,
    diff_se = NA_real_,
    diff_p = NA_character_,
    direction_consistent_with_theory = NA_character_,
    notes = paste0(sig_status, "; ", n_districts, " districts")
  )

summary_tbl <- bind_rows(summary_tbl, composite_row)
write_csv(summary_tbl, file.path(OUT, "tables/summary_heterogeneity.csv"))

summary_html <- summary_tbl %>%
  kbl(caption = "District heterogeneity in attached/casual wage pass-through differential") %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover")) %>%
  footnote(
    general = paste(
      "Cross-sectional splits using time-invariant district characteristics (HIES 2016, MICS 2019, Geography 2013)",
      "applied to the 2017-2023 BBS wage panel. Interpret as descriptive heterogeneity, not causal identification.",
      "Baseline spec: diff_real_wage ~ diff_log_yield_hat x meal_type + gender | year + District x growing_season;",
      "clustered SE at district level."
    )
  )
writeLines(summary_html, file.path(OUT, "tables/summary_heterogeneity.html"))

plot_df <- split_results %>%
  select(variable_label, group_high_label, group_low_label,
         high_differential, high_se, low_differential, low_se) %>%
  pivot_longer(
    cols = c(high_differential, low_differential),
    names_to = "which",
    values_to = "differential"
  ) %>%
  mutate(
    se = ifelse(which == "high_differential", high_se, low_se),
    group = ifelse(which == "high_differential", group_high_label, group_low_label),
    variable_label = factor(variable_label, levels = split_specs$variable_label)
  )

p_coef <- ggplot(plot_df, aes(x = differential, y = group, color = group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = differential - 1.96 * se, xmax = differential + 1.96 * se),
                 height = 0.2) +
  facet_wrap(~variable_label, scales = "free_y", ncol = 2) +
  labs(
    title = "Attached/casual pass-through differential by district characteristic",
    subtitle = "High-constraint vs low-constraint groups (median split); 95% CI",
    x = "Attached/casual differential (BDT per unit log yield)",
    y = NULL,
    caption = paste(
      "Descriptive heterogeneity: time-invariant district splits on 2017-2023 panel.",
      "Theory predicts more negative differential in constrained groups."
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

ggsave(file.path(OUT, "figures/coefplot_split_heterogeneity.png"),
       p_coef, width = 11, height = 9, dpi = 150)

p_composite <- ggplot(composite_results,
                      aes(x = factor(ic_group, levels = composite_specs),
                          y = differential)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3, color = "#0072B2") +
  geom_errorbar(aes(ymin = differential - 1.96 * differential_se,
                    ymax = differential + 1.96 * differential_se),
                width = 0.15, color = "#0072B2") +
  labs(
    title = "Pass-through differential by implicit contract environment index",
    x = "Composite index group (sum of 6 binary indicators)",
    y = "Attached/casual differential (BDT)"
  ) +
  theme_minimal(base_size = 11)

ggsave(file.path(OUT, "figures/coefplot_composite_index.png"),
       p_composite, width = 8, height = 5, dpi = 150)

n_consistent <- sum(split_results$direction_consistent)
n_sig_consistent <- sum(split_results$direction_consistent & split_results$diff_p < 0.05)

log_msg("\n=== SUMMARY ===")
log_msg("  Directionally consistent splits: ", n_consistent, "/6")
log_msg("  Statistically significant group differences: ", n_sig_consistent, "/6")
log_msg("  Composite index monotone (more negative at higher index): ", monotone)
log_msg("  NOTE: Cross-sectional splits with time-invariant district characteristics;",
        " descriptive heterogeneity only.")
log_msg("Done. Outputs: ", OUT)
