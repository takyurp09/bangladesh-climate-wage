## 10_meal_monotonicity.R
## Internal validation: monotonic yield pass-through across meal categories (0–3 meals)
## Steps: (1) first stage yield_hat, (2) separate second stages by meal,
##        (3) coefficient plot, (4) formal monotonicity tests,
##        (5) wage-variance ordering, (6) summary table
## Output: output/stage2/meal_monotonicity/

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
OUT  <- file.path(ROOT, "output/stage2/meal_monotonicity")
for (d in c("", "tables", "figures", "models", "summary"))
  dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)

MEAL_LEVELS <- c("None", "One", "Two", "Three")
MEAL_LABELS <- c("0 meals", "1 meal", "2 meals", "3 meals")
names(MEAL_LABELS) <- MEAL_LEVELS

cat("=== MEAL MONOTONICITY VALIDATION ===\n")
cat("Output directory:", OUT, "\n\n")

## ── Helpers ───────────────────────────────────────────────────────────────── ##
fmt_p <- function(p) {
  vapply(p, function(x) {
    if (is.na(x)) return("NA")
    if (x < 0.001) return("<0.001")
    sprintf("%.3f", x)
  }, character(1))
}

extract_yield_coef <- function(mod, param = "diff_log_yield_hat") {
  ct <- summary(mod)$coeftable
  if (!param %in% rownames(ct)) stop("Parameter not found: ", param)
  data.frame(
    coef  = ct[param, "Estimate"],
    se    = ct[param, "Std. Error"],
    t     = ct[param, "t value"],
    p     = ct[param, "Pr(>|t|)"],
    ci_lo = ct[param, "Estimate"] - 1.96 * ct[param, "Std. Error"],
    ci_hi = ct[param, "Estimate"] + 1.96 * ct[param, "Std. Error"],
    nobs  = nobs(mod),
    stringsAsFactors = FALSE
  )
}

one_sided_p <- function(est, se, df, direction = "less") {
  t_stat <- est / se
  if (direction == "less") pt(t_stat, df = df, lower.tail = TRUE)
  else pt(t_stat, df = df, lower.tail = FALSE)
}

## Brown-Forsythe Levene test (median-centered), no extra dependencies
wald_linear <- function(mod, terms, weights, df2 = NULL) {
  b <- coef(mod)
  V <- vcov(mod)
  idx <- match(terms, names(b))
  if (any(is.na(idx))) stop("Term not found in model: ", paste(terms[is.na(idx)], collapse = ", "))
  L <- numeric(length(b))
  L[idx] <- weights
  est <- sum(L * b)
  var <- as.numeric(t(L) %*% V %*% L)
  se  <- sqrt(var)
  if (is.null(df2)) df2 <- mod$nobs - length(b)
  t_stat <- est / se
  list(
    estimate = est,
    se       = se,
    stat     = t_stat^2,
    p_value  = 2 * pt(-abs(t_stat), df = df2)
  )
}

levene_test <- function(y, group) {
  d <- data.frame(y = y, g = factor(group))
  d <- d[complete.cases(d), , drop = FALSE]
  d$center <- ave(d$y, d$g, FUN = median)
  d$z <- abs(d$y - d$center)
  fit <- lm(z ~ g, data = d)
  a <- anova(fit)
  data.frame(
    F_stat  = a$`F value`[1],
    df1     = a$Df[1],
    df2     = a$Df[2],
    p_value = a$`Pr(>F)`[1],
    stringsAsFactors = FALSE
  )
}

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 1 — First stage: predict yield from climate (same spec as main)        ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 1 — First stage (climate → yield)\n")

df_yield_raw <- read_csv(
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
  )

df_yield_est <- df_yield_raw %>%
  filter(
    !is.na(diff_log_yield),
    !is.na(diff_gdd_10_30),
    !is.na(diff_edd_30),
    !is.na(diff_precip),
    !is.na(diff_precip_sq)
  )

m_first <- feols(
  diff_log_yield ~ diff_gdd_10_30 + diff_edd_30 + diff_precip + diff_precip_sq |
    year + district^season,
  data = df_yield_est,
  cluster = ~district
)

cat(sprintf("  First-stage N = %d | districts = %d\n",
            nobs(m_first), length(unique(df_yield_est$district))))

df_yield_est$yield_hat <- fitted(m_first)

df_yield_hat <- df_yield_est %>%
  filter(year >= 2017) %>%
  rename(District = district, growing_season = season) %>%
  arrange(District, growing_season, year) %>%
  group_by(District, growing_season) %>%
  mutate(diff_log_yield_hat = yield_hat - lag(yield_hat)) %>%
  ungroup() %>%
  filter(!is.na(diff_log_yield_hat)) %>%
  select(District, growing_season, year, yield_hat, diff_log_yield_hat)

cat(sprintf("  Fitted yield_hat rows (2017+, after FD): %d\n", nrow(df_yield_hat)))

## ── Build Stage 2 panel (mirrors 00_dataprep.R) ───────────────────────────── ##
wage_raw <- read_csv(
  file.path(ROOT, "data/Regression_data/wage_by_growing_season.csv"),
  show_col_types = FALSE
) %>%
  filter(meal_type != "No_info", !is.na(real_wage)) %>%
  group_by(District, year, growing_season, gender, meal_type) %>%
  summarise(real_wage = mean(real_wage, na.rm = TRUE), .groups = "drop") %>%
  arrange(District, growing_season, gender, meal_type, year) %>%
  group_by(District, growing_season, gender, meal_type) %>%
  mutate(
    diff_real_wage = real_wage - lag(real_wage)
  ) %>%
  ungroup() %>%
  filter(!is.na(diff_real_wage)) %>%
  mutate(
    meal_type = relevel(factor(meal_type), ref = "None"),
    meal_num  = as.integer(factor(meal_type, levels = MEAL_LEVELS)) - 1L
  )

df <- wage_raw %>%
  left_join(df_yield_hat, by = c("District", "growing_season", "year")) %>%
  filter(!is.na(diff_log_yield_hat))

cat(sprintf("  Stage 2 panel: N = %d | districts = %d | years = %s\n\n",
            nrow(df), length(unique(df$District)),
            paste(sort(unique(df$year)), collapse = ", ")))

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 2 — Four separate second-stage regressions (one per meal category)     ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 2 — Separate pass-through regressions by meal category\n")
cat("  Spec: diff_real_wage ~ diff_log_yield_hat | year + District^growing_season\n")

sep_models <- lapply(MEAL_LEVELS, function(m) {
  sub <- df %>% filter(meal_type == m)
  feols(
    diff_real_wage ~ diff_log_yield_hat | year + District^growing_season,
    data = sub,
    cluster = ~District
  )
})
names(sep_models) <- MEAL_LEVELS

sep_results <- bind_rows(lapply(MEAL_LEVELS, function(m) {
  sub  <- df %>% filter(meal_type == m)
  est  <- extract_yield_coef(sep_models[[m]])
  data.frame(
    meal_type      = m,
    meal_label     = MEAL_LABELS[m],
    meal_num       = which(MEAL_LEVELS == m) - 1L,
    n_obs          = nrow(sub),
    n_districts    = length(unique(sub$District)),
    mean_wage      = mean(sub$real_wage, na.rm = TRUE),
    pass_coef      = est$coef,
    pass_se        = est$se,
    pass_p         = est$p,
    pass_ci_lo     = est$ci_lo,
    pass_ci_hi     = est$ci_hi,
    stringsAsFactors = FALSE
  )
}))

for (i in seq_len(nrow(sep_results))) {
  r <- sep_results[i, ]
  cat(sprintf("  %s: coef=%.2f SE=%.2f p=%s N=%d districts=%d\n",
              r$meal_label, r$pass_coef, r$pass_se, fmt_p(r$pass_p),
              r$n_obs, r$n_districts))
}

is_monotone_sep <- all(diff(sep_results$pass_coef) <= 0)
cat(sprintf("\n  Separate-regression ordering β(0) ≥ β(1) ≥ β(2) ≥ β(3): %s\n\n",
            if (is_monotone_sep) "HOLDS" else "** VIOLATED — FLAG **"))

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 3 — Coefficient plot                                                   ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 3 — Coefficient plot\n")

plot_df <- sep_results %>%
  mutate(
    meal_label = factor(meal_label, levels = MEAL_LABELS),
    label_txt  = sprintf("%.1f\n(p=%s)", pass_coef, fmt_p(pass_p))
  )

meal_cols <- c("0 meals" = "#E69F00", "1 meal" = "#56B4E9",
               "2 meals" = "#009E73", "3 meals" = "#D55E00")

p_coef <- ggplot(plot_df, aes(x = meal_label, y = pass_coef, color = meal_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.5) +
  geom_errorbar(aes(ymin = pass_ci_lo, ymax = pass_ci_hi),
                width = 0.12, linewidth = 0.9, color = "grey30") +
  geom_point(size = 4) +
  geom_text(aes(label = label_txt), vjust = -1.4, size = 3.2, color = "grey20",
            lineheight = 0.9) +
  scale_color_manual(values = meal_cols, guide = "none") +
  scale_x_discrete(labels = MEAL_LABELS) +
  labs(
    title    = "Yield Pass-Through by Meal Provision",
    subtitle = "Separate regressions per meal category | 95% CI | SE clustered by district",
    x        = "Meal provision (contract attachment proxy)",
    y        = "Pass-through coefficient on change in fitted log yield",
    caption  = if (!is_monotone_sep) {
      "Note: coefficient ordering is non-monotonic — interpret with caution."
    } else {
      "Expected pattern: coefficients become more negative as meal provision increases."
    }
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  coord_cartesian(ylim = range(c(plot_df$pass_ci_lo, plot_df$pass_ci_hi, 0)) * c(1.15, 1.15))

ggsave(file.path(OUT, "figures", "passthrough_by_meal.png"),
       p_coef, width = 8, height = 5.5, dpi = 300)
ggsave(file.path(OUT, "figures", "passthrough_by_meal.pdf"),
       p_coef, width = 8, height = 5.5)
cat("  Saved passthrough_by_meal.png/.pdf\n\n")

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 4 — Formal monotonicity tests                                          ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 4 — Formal monotonicity tests\n")

## 4a. Joint interaction model (0-meal baseline) — mirrors M5 without gender
m_joint <- feols(
  diff_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal_type |
    year + District^growing_season,
  data = df,
  cluster = ~District
)

m_joint_ct <- summary(m_joint)$coeftable
base_coef  <- m_joint_ct["diff_log_yield_hat", "Estimate"]

## Wald tests vs 0-meal baseline (interaction coef = 0)
wald_vs_base <- lapply(MEAL_LEVELS[-1], function(m) {
  int_nm <- paste0("diff_log_yield_hat:meal_type", m)
  w <- wald(m_joint, int_nm)
  data.frame(
    comparison = sprintf("%s-meal vs 0-meal", which(MEAL_LEVELS == m) - 1L),
    hypothesis = sprintf("β(%s) = β(0)", m),
    constraint = int_nm,
    stat       = w$stat,
    p_value    = w$p,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

## 4b. Pairwise adjacent tests (linear combinations via vcov)
adj_pairs <- list(
  list(label = "1 vs 0",
       terms = "diff_log_yield_hat:meal_typeOne",
       weights = 1),
  list(label = "2 vs 1",
       terms = c("diff_log_yield_hat:meal_typeTwo", "diff_log_yield_hat:meal_typeOne"),
       weights = c(1, -1)),
  list(label = "3 vs 2",
       terms = c("diff_log_yield_hat:meal_typeThree", "diff_log_yield_hat:meal_typeTwo"),
       weights = c(1, -1))
)

## Total coefficients with proper delta-method SE via vcov
V <- vcov(m_joint)
df_resid <- m_joint$nobs - length(coef(m_joint))

total_joint <- lapply(MEAL_LEVELS, function(m) {
  if (m == "None") {
    est <- base_coef; se <- sqrt(V["diff_log_yield_hat", "diff_log_yield_hat"])
  } else {
    int_nm <- paste0("diff_log_yield_hat:meal_type", m)
    est <- base_coef + coef(m_joint)[int_nm]
    se  <- sqrt(V["diff_log_yield_hat", "diff_log_yield_hat"] +
                  V[int_nm, int_nm] +
                  2 * V["diff_log_yield_hat", int_nm])
  }
  data.frame(
    meal_type = m, total_coef = est, total_se = se,
    total_p = 2 * pt(-abs(est / se), df = m_joint$nobs - length(coef(m_joint)))
  )
}) %>% bind_rows()

wald_adjacent <- bind_rows(lapply(seq_along(adj_pairs), function(i) {
  ap <- adj_pairs[[i]]
  w  <- wald_linear(m_joint, ap$terms, ap$weights, df2 = df_resid)
  lo <- total_joint$total_coef[i]
  hi <- total_joint$total_coef[i + 1]
  data.frame(
    comparison       = ap$label,
    coef_lower_meal  = lo,
    coef_higher_meal = hi,
    diff_observed    = hi - lo,
    stat             = w$stat,
    p_value          = w$p_value,
    stringsAsFactors = FALSE
  )
}))

## 4c. One-sided monotonicity constraints: β(k+1) ≤ β(k)
##     Equivalently: interaction_One ≤ 0; int_Two - int_One ≤ 0; int_Three - int_Two ≤ 0
mono_constraints <- list(
  list(label = "β(1) ≤ β(0)", expr = "diff_log_yield_hat:meal_typeOne",
       est = coef(m_joint)["diff_log_yield_hat:meal_typeOne"],
       se  = se(m_joint)["diff_log_yield_hat:meal_typeOne"]),
  list(label = "β(2) ≤ β(1)",
       expr = "diff_log_yield_hat:meal_typeTwo - diff_log_yield_hat:meal_typeOne",
       est = coef(m_joint)["diff_log_yield_hat:meal_typeTwo"] -
         coef(m_joint)["diff_log_yield_hat:meal_typeOne"],
       se  = {
         v <- V["diff_log_yield_hat:meal_typeTwo", "diff_log_yield_hat:meal_typeTwo"] +
           V["diff_log_yield_hat:meal_typeOne", "diff_log_yield_hat:meal_typeOne"] -
           2 * V["diff_log_yield_hat:meal_typeTwo", "diff_log_yield_hat:meal_typeOne"]
         sqrt(v)
       }),
  list(label = "β(3) ≤ β(2)",
       expr = "diff_log_yield_hat:meal_typeThree - diff_log_yield_hat:meal_typeTwo",
       est = coef(m_joint)["diff_log_yield_hat:meal_typeThree"] -
         coef(m_joint)["diff_log_yield_hat:meal_typeTwo"],
       se  = {
         v <- V["diff_log_yield_hat:meal_typeThree", "diff_log_yield_hat:meal_typeThree"] +
           V["diff_log_yield_hat:meal_typeTwo", "diff_log_yield_hat:meal_typeTwo"] -
           2 * V["diff_log_yield_hat:meal_typeThree", "diff_log_yield_hat:meal_typeTwo"]
         sqrt(v)
       })
)

mono_one_sided <- bind_rows(lapply(mono_constraints, function(c) {
  data.frame(
    constraint    = c$label,
    expression    = c$expr,
    estimate      = c$est,
    se            = c$se,
    t_stat        = c$est / c$se,
    p_one_sided   = one_sided_p(c$est, c$se, df_resid, direction = "less"),
    satisfied     = c$est <= 0,
    stringsAsFactors = FALSE
  )
}))

## Joint monotonicity: all three one-sided constraints hold
## Report intersection-union style: max t-stat across constraints
mono_joint_stat <- max(mono_one_sided$t_stat, na.rm = TRUE)
mono_joint_p    <- min(mono_one_sided$p_one_sided, na.rm = TRUE)
mono_holm_p     <- p.adjust(mono_one_sided$p_one_sided, method = "holm")

mono_one_sided$p_holm <- mono_holm_p
mono_all_satisfied    <- all(mono_one_sided$satisfied)

cat("  4a. Wald tests vs 0-meal baseline:\n")
for (i in seq_len(nrow(wald_vs_base)))
  cat(sprintf("      %s: stat=%.3f p=%s\n",
              wald_vs_base$comparison[i], wald_vs_base$stat[i],
              fmt_p(wald_vs_base$p_value[i])))

cat("  4b. Pairwise adjacent (two-sided wald):\n")
for (i in seq_len(nrow(wald_adjacent)))
  cat(sprintf("      %s: Δ=%.2f stat=%.3f p=%s\n",
              wald_adjacent$comparison[i], wald_adjacent$diff_observed[i],
              wald_adjacent$stat[i], fmt_p(wald_adjacent$p_value[i])))

cat("  4c. One-sided monotonicity constraints:\n")
for (i in seq_len(nrow(mono_one_sided)))
  cat(sprintf("      %s: est=%.2f p(1-sided)=%s Holm=%s %s\n",
              mono_one_sided$constraint[i], mono_one_sided$estimate[i],
              fmt_p(mono_one_sided$p_one_sided[i]),
              fmt_p(mono_one_sided$p_holm[i]),
              if (mono_one_sided$satisfied[i]) "[OK]" else "[VIOLATED]"))

cat(sprintf("\n  Joint monotonicity (all constraints): %s | min p(1-sided)=%s\n\n",
            if (mono_all_satisfied) "SATISFIED" else "** VIOLATED **",
            fmt_p(mono_joint_p)))

## Save test outputs
write_csv(wald_vs_base,    file.path(OUT, "tables", "wald_vs_baseline.csv"))
write_csv(wald_adjacent,   file.path(OUT, "tables", "wald_adjacent_pairs.csv"))
write_csv(mono_one_sided,  file.path(OUT, "tables", "monotonicity_one_sided.csv"))
write_csv(total_joint,     file.path(OUT, "tables", "joint_model_total_coefs.csv"))

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 5 — Wage variance ordering                                             ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 5 — Within-district wage SD by meal category\n")

wage_levels <- read_csv(
  file.path(ROOT, "data/Regression_data/wage_by_growing_season.csv"),
  show_col_types = FALSE
) %>%
  filter(
    meal_type != "No_info",
    !is.na(real_wage),
    year >= 2017, year <= 2023
  ) %>%
  mutate(meal_type = factor(meal_type, levels = MEAL_LEVELS))

## District-year-meal mean (pool gender & season for overall SD)
wage_dy_meal <- wage_levels %>%
  group_by(District, year, meal_type) %>%
  summarise(real_wage = mean(real_wage, na.rm = TRUE), .groups = "drop")

within_sd_overall <- wage_dy_meal %>%
  group_by(District, meal_type) %>%
  summarise(
    within_sd = sd(real_wage, na.rm = TRUE),
    n_years   = n(),
    .groups = "drop"
  ) %>%
  filter(n_years >= 2)

var_summary_overall <- within_sd_overall %>%
  group_by(meal_type) %>%
  summarise(
    mean_within_sd = mean(within_sd, na.rm = TRUE),
    median_within_sd = median(within_sd, na.rm = TRUE),
    n_districts    = n(),
    .groups = "drop"
  ) %>%
  mutate(meal_label = MEAL_LABELS[as.character(meal_type)])

levene_overall <- levene_test(wage_dy_meal$real_wage, wage_dy_meal$meal_type)

## By season: district-season-meal SD across years
wage_dsy_meal <- wage_levels %>%
  group_by(District, growing_season, year, meal_type) %>%
  summarise(real_wage = mean(real_wage, na.rm = TRUE), .groups = "drop")

within_sd_season <- wage_dsy_meal %>%
  group_by(District, growing_season, meal_type) %>%
  summarise(
    within_sd = sd(real_wage, na.rm = TRUE),
    n_years   = n(),
    .groups = "drop"
  ) %>%
  filter(n_years >= 2)

var_summary_season <- within_sd_season %>%
  group_by(growing_season, meal_type) %>%
  summarise(
    mean_within_sd = mean(within_sd, na.rm = TRUE),
    n_units        = n(),
    .groups = "drop"
  ) %>%
  mutate(meal_label = MEAL_LABELS[as.character(meal_type)])

levene_by_season <- lapply(c("Boro", "Aus", "Aman"), function(s) {
  sub <- wage_dsy_meal %>% filter(growing_season == s)
  lt  <- levene_test(sub$real_wage, sub$meal_type)
  data.frame(season = s, lt, stringsAsFactors = FALSE)
}) %>% bind_rows()

is_var_monotone <- all(diff(var_summary_overall$mean_within_sd[order(
  match(var_summary_overall$meal_type, MEAL_LEVELS))]) < 0)

cat("  Mean within-district wage SD (2017–2023, pooled seasons):\n")
for (i in seq_len(nrow(var_summary_overall))) {
  r <- var_summary_overall[i, ]
  cat(sprintf("    %s: SD=%.2f (n=%d districts)\n",
              r$meal_label, r$mean_within_sd, r$n_districts))
}
cat(sprintf("  Expected SD(3) < SD(2) < SD(1) < SD(0): %s\n",
            if (is_var_monotone) "HOLDS" else "VIOLATED"))
cat(sprintf("  Levene test (overall): F=%.3f p=%s\n",
            levene_overall$F_stat[1], fmt_p(levene_overall$p_value[1])))

cat("  Levene by season:\n")
for (i in seq_len(nrow(levene_by_season)))
  cat(sprintf("    %s: F=%.3f p=%s\n",
              levene_by_season$season[i], levene_by_season$F_stat[i],
              fmt_p(levene_by_season$p_value[i])))
cat("\n")

write_csv(var_summary_overall, file.path(OUT, "tables", "wage_sd_overall.csv"))
write_csv(var_summary_season,  file.path(OUT, "tables", "wage_sd_by_season.csv"))
write_csv(levene_by_season,    file.path(OUT, "tables", "levene_by_season.csv"))

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 6 — Summary table (proxy validation)                                   ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 6 — Summary table\n")

summary_tbl <- sep_results %>%
  left_join(
    var_summary_overall %>%
      select(meal_type, within_sd = mean_within_sd),
    by = "meal_type"
  ) %>%
  transmute(
    meal_category     = meal_label,
    n_observations    = n_obs,
    n_districts       = n_districts,
    mean_wage         = round(mean_wage, 2),
    within_district_sd = round(within_sd, 2),
    passthrough_coef  = round(pass_coef, 2),
    passthrough_se    = round(pass_se, 2),
    p_value           = pass_p,
    ci_95_lo          = round(pass_ci_lo, 2),
    ci_95_hi          = round(pass_ci_hi, 2),
    monotone_order    = if (is_monotone_sep) "Yes" else "No"
  )

write_csv(summary_tbl, file.path(OUT, "tables", "proxy_validation_summary.csv"))

## HTML table
tbl_html <- summary_tbl %>%
  mutate(
    p_value  = sapply(p_value, fmt_p),
    ci_95    = sprintf("[%.2f, %.2f]", ci_95_lo, ci_95_hi)
  ) %>%
  select(-ci_95_lo, -ci_95_hi) %>%
  kbl(caption = "Meal Provision Proxy Validation: Pass-Through and Wage Dispersion",
      col.names = c("Meal category", "N obs.", "N districts", "Mean wage",
                    "Within-district SD", "Pass-through β",
                    "SE", "p-value", "95% CI", "Monotone order"),
      digits = 2, format = "html", booktabs = TRUE) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE, font_size = 13) %>%
  footnote(
    general = paste(
      "Pass-through: separate regressions of Δ real wage on Δ fitted log yield",
      "with year and District×Season FE; SE clustered by district.",
      "Within-district SD: SD of mean district-year wages (2017–2023) within each district.",
      if (!is_monotone_sep) "FLAG: pass-through coefficients are not monotonically ordered."
      else "Coefficients decrease with meal provision, consistent with risk-absorption gradient."
    ),
    general_title = "Note: ",
    footnote_as_chunk = TRUE
  )

writeLines(as.character(tbl_html), file.path(OUT, "tables", "proxy_validation_summary.html"))
cat("  Saved proxy_validation_summary.csv/.html\n\n")

## ── Save models & narrative summary ───────────────────────────────────────── ##
save(m_first, m_joint, sep_models,
     sep_results, wald_vs_base, wald_adjacent, mono_one_sided,
     var_summary_overall, levene_overall, levene_by_season, summary_tbl,
     file = file.path(OUT, "models", "meal_monotonicity_models.RData"))

summary_lines <- c(
  "# Meal Monotonicity Validation",
  sprintf("Generated: %s", Sys.time()),
  "",
  "## Pass-through gradient (separate regressions)",
  sprintf("- Monotone ordering β(0) ≥ β(1) ≥ β(2) ≥ β(3): **%s**",
          if (is_monotone_sep) "HOLDS" else "VIOLATED"),
  capture.output(print(summary_tbl %>% select(meal_category, passthrough_coef, p_value, monotone_order))),
  "",
  "## Formal tests (interaction model)",
  "- Wald vs 0-meal baseline: see tables/wald_vs_baseline.csv",
  "- Adjacent pairwise tests: see tables/wald_adjacent_pairs.csv",
  sprintf("- One-sided monotonicity constraints all satisfied: **%s** (min p = %s)",
          if (mono_all_satisfied) "Yes" else "No", fmt_p(mono_joint_p)),
  "",
  "## Wage variance ordering",
  sprintf("- Expected SD(3) < SD(2) < SD(1) < SD(0): **%s**",
          if (is_var_monotone) "HOLDS" else "VIOLATED"),
  sprintf("- Levene test (overall): F=%.3f, p=%s",
          levene_overall$F_stat[1], fmt_p(levene_overall$p_value[1])),
  "",
  "## Outputs",
  "- figures/passthrough_by_meal.png",
  "- tables/proxy_validation_summary.csv/.html",
  "- tables/wald_vs_baseline.csv",
  "- tables/wald_adjacent_pairs.csv",
  "- tables/monotonicity_one_sided.csv",
  "- tables/wage_sd_overall.csv",
  "- tables/wage_sd_by_season.csv",
  "- tables/levene_by_season.csv"
)

writeLines(summary_lines, file.path(OUT, "summary", "meal_monotonicity_report.md"))
cat("Saved meal_monotonicity_report.md\n")
cat("=== MEAL MONOTONICITY VALIDATION COMPLETE ===\n")
