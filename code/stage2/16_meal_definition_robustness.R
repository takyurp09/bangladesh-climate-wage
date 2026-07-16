## 16_meal_definition_robustness.R
## Robustness: attached/casual differential under alternative meal definitions
## Specs A–D + monotonicity gradient + formal 3-meal vs 2-meal test
## Output: output/stage2/meal_definition_robustness/

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(fixest)
  library(ggplot2)
  library(kableExtra)
})

ROOT <- here::here()
OUT  <- file.path(ROOT, "output/stage2/meal_definition_robustness")
for (d in c("tables", "figures", "logs"))
  dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(OUT, "logs/meal_robustness_log.txt")
log_con  <- file(log_path, open = "wt")
on.exit(close(log_con), add = TRUE)

PAPER_DIFF <- -222.6
PAPER_SE   <-  95.0

MEAL_LEVELS <- c("None", "One", "Two", "Three")
MEAL_LABELS <- c("0 meals", "1 meal", "2 meals", "3 meals")
names(MEAL_LABELS) <- MEAL_LEVELS

FE <- "year + District^growing_season"

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
         ifelse(p < 0.01, "significant (***)",
                ifelse(p < 0.05, "significant (**)",
                       ifelse(p < 0.10, "marginally significant (*)",
                              "imprecise (not significant)"))))
}

direction_flag <- function(diff) {
  ifelse(is.na(diff), "NA",
         ifelse(diff < 0, "directionally consistent",
                "REVERSED SIGN — investigate"))
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

run_m5_spec_a <- function(df) {
  mod <- feols(
    diff_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal_type + gender |
      year + District^growing_season,
    data = df,
    cluster = ~District
  )
  ct <- summary(mod)$coeftable
  V  <- vcov(mod)
  df2 <- mod$nobs - length(coef(mod))
  bnm <- "diff_log_yield_hat"
  int <- "diff_log_yield_hat:meal_typeThree"

  b0  <- ct[bnm, "Estimate"]
  s0  <- ct[bnm, "Std. Error"]
  p0  <- ct[bnm, "Pr(>|t|)"]
  b1  <- ct[int, "Estimate"]
  s1  <- ct[int, "Std. Error"]
  p1  <- ct[int, "Pr(>|t|)"]
  b_att <- b0 + b1
  s_att <- sqrt(V[bnm, bnm] + V[int, int] + 2 * V[bnm, int])
  p_att <- 2 * pt(-abs(b_att / s_att), df = df2)

  sd_yield   <- sd(df$diff_log_yield_hat, na.rm = TRUE)
  mean_wage  <- mean(df$real_wage, na.rm = TRUE)
  effect_1sd <- b1 * sd_yield

  data.frame(
    spec_id = "A", spec_label = "Baseline (3 vs 0)",
    attached_def = "3-meal only", casual_def = "0-meal only",
    n_attached = sum(df$meal_type == "Three"),
    n_casual   = sum(df$meal_type == "None"),
    n_total = nobs(mod),
    casual_coef = b0, casual_se = s0, casual_p = p0,
    attached_coef = b_att, attached_se = s_att, attached_p = p_att,
    differential = b1, differential_se = s1, differential_p = p1,
    significant = sig_flag(p1), direction = direction_flag(b1),
    sd_yield_hat = sd_yield, mean_daily_wage = mean_wage,
    effect_1sd_bdt = effect_1sd,
    pct_mean_wage = 100 * abs(effect_1sd) / mean_wage,
    stringsAsFactors = FALSE
  )
}

run_binary_spec <- function(df, spec_id, spec_label,
                            attached_def, casual_def,
                            attached_meals, casual_meals) {
  sub <- df %>%
    filter(meal_type %in% c(attached_meals, casual_meals)) %>%
    mutate(attached = as.integer(meal_type %in% attached_meals))

  n_att <- sum(sub$attached == 1)
  n_cas <- sum(sub$attached == 0)

  mod <- feols(
    diff_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:attached + gender |
      year + District^growing_season,
    data = sub,
    cluster = ~District
  )

  ct <- summary(mod)$coeftable
  V  <- vcov(mod)
  df2 <- mod$nobs - length(coef(mod))

  b0  <- ct["diff_log_yield_hat", "Estimate"]
  s0  <- ct["diff_log_yield_hat", "Std. Error"]
  p0  <- ct["diff_log_yield_hat", "Pr(>|t|)"]
  int <- "diff_log_yield_hat:attached"
  b1  <- ct[int, "Estimate"]
  s1  <- ct[int, "Std. Error"]
  p1  <- ct[int, "Pr(>|t|)"]

  b_att <- b0 + b1
  s_att <- sqrt(V["diff_log_yield_hat", "diff_log_yield_hat"] +
                  V[int, int] + 2 * V["diff_log_yield_hat", int])
  p_att <- 2 * pt(-abs(b_att / s_att), df = df2)

  sd_yield   <- sd(sub$diff_log_yield_hat, na.rm = TRUE)
  mean_wage  <- mean(sub$real_wage, na.rm = TRUE)
  effect_1sd <- b1 * sd_yield
  pct_wage   <- 100 * abs(effect_1sd) / mean_wage

  data.frame(
    spec_id            = spec_id,
    spec_label         = spec_label,
    attached_def       = attached_def,
    casual_def         = casual_def,
    n_attached         = n_att,
    n_casual           = n_cas,
    n_total            = nobs(mod),
    casual_coef        = b0,
    casual_se          = s0,
    casual_p           = p0,
    attached_coef      = b_att,
    attached_se        = s_att,
    attached_p         = p_att,
    differential       = b1,
    differential_se    = s1,
    differential_p     = p1,
    significant        = sig_flag(p1),
    direction          = direction_flag(b1),
    sd_yield_hat       = sd_yield,
    mean_daily_wage    = mean_wage,
    effect_1sd_bdt     = effect_1sd,
    pct_mean_wage      = pct_wage,
    stringsAsFactors   = FALSE
  )
}

log_msg("=== MEAL DEFINITION ROBUSTNESS (SPECS A–D) ===")
log_msg("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

df <- read_csv(
  file.path(ROOT, "data/Regression_data/df_2_merged_v2.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    meal_type = relevel(factor(meal_type), ref = "None"),
    gender    = relevel(factor(gender),    ref = "Female")
  ) %>%
  filter(!is.na(diff_log_yield_hat), !is.na(diff_real_wage))

log_msg("Panel N = ", nrow(df), " | districts = ", length(unique(df$District)))

## ── Four specifications ───────────────────────────────────────────────────────
spec_defs <- tribble(
  ~spec_id, ~spec_label, ~attached_def, ~casual_def,
           ~attached_meals, ~casual_meals,
  "A", "Baseline (3 vs 0)",
    "3-meal only", "0-meal only",
    list("Three"), list("None"),
  "B", "Broader casual (3 vs 0+1)",
    "3-meal only", "0-meal + 1-meal",
    list("Three"), list("None", "One"),
  "C", "Broader attached (2+3 vs 0)",
    "2-meal + 3-meal", "0-meal only",
    list("Two", "Three"), list("None"),
  "D", "Broadest split (2+3 vs 0+1)",
    "2-meal + 3-meal", "0-meal + 1-meal",
    list("Two", "Three"), list("None", "One")
)

spec_a <- run_m5_spec_a(df)

spec_bcd <- pmap_dfr(
  spec_defs %>% filter(spec_id != "A"),
  function(spec_id, spec_label, attached_def, casual_def,
           attached_meals, casual_meals) {
    run_binary_spec(df, spec_id, spec_label, attached_def, casual_def,
                    attached_meals, casual_meals)
  }
)

spec_results <- bind_rows(spec_a, spec_bcd) %>%
  arrange(match(spec_id, c("A", "B", "C", "D")))

write_csv(spec_results, file.path(OUT, "tables/four_specification_results.csv"))

## Verify Spec A against standard M5 interaction
m_m5 <- feols(
  diff_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal_type + gender |
    year + District^growing_season,
  data = df,
  cluster = ~District
)
m5_diff <- coef(m_m5)["diff_log_yield_hat:meal_typeThree"]
m5_se   <- se(m_m5)["diff_log_yield_hat:meal_typeThree"]

log_msg("\n--- Specification results ---")
for (i in seq_len(nrow(spec_results))) {
  r <- spec_results[i, ]
  log_msg(sprintf("  Spec %s (%s):", r$spec_id, r$spec_label))
  log_msg(sprintf("    Attached: %s | Casual: %s", r$attached_def, r$casual_def))
  log_msg(sprintf("    N_att = %d | N_cas = %d", r$n_attached, r$n_casual))
  log_msg(sprintf("    Casual coef:   %.1f (SE %.1f, p=%s)", r$casual_coef, r$casual_se, fmt_p(r$casual_p)))
  log_msg(sprintf("    Attached coef: %.1f (SE %.1f, p=%s)", r$attached_coef, r$attached_se, fmt_p(r$attached_p)))
  log_msg(sprintf("    Differential:  %.1f (SE %.1f, p=%s) [%s | %s]",
                  r$differential, r$differential_se, fmt_p(r$differential_p),
                  r$significant, r$direction))
  log_msg(sprintf("    1-SD shock effect: %.1f BDT/day (%.1f%% mean wage)",
                  r$effect_1sd_bdt, r$pct_mean_wage))
}

log_msg(sprintf("\n  M5 full-model 3-meal interaction (replication check): %.1f (SE %.1f)",
                m5_diff, m5_se))
log_msg(sprintf("  Spec A binary (3 vs 0 only): %.1f (SE %.1f)",
                spec_results$differential[1], spec_results$differential_se[1]))
spec_a_diff <- spec_results$differential[spec_results$spec_id == "A"]
if (abs(spec_a_diff - PAPER_DIFF) > 5) {
  log_msg("  FLAG: Spec A differs from paper benchmark by >5 BDT")
} else {
  log_msg("  Spec A matches paper benchmark (", PAPER_DIFF, " BDT)")
}

reversed <- spec_results %>% filter(differential > 0)
if (nrow(reversed) > 0)
  log_msg("  FLAG: Reversed-sign specifications: ",
          paste(reversed$spec_id, collapse = ", "))

## ── Publication-ready table (rows × columns) ────────────────────────────────
pub_rows <- c(
  "Attached definition",
  "Casual definition",
  "N (attached)",
  "N (casual)",
  "Attached coefficient (SE)",
  "Casual coefficient (SE)",
  "Differential (SE)",
  "p-value (differential)",
  "Significant?",
  "Economic magnitude (BDT per 1 SD shock)"
)

pub_tbl <- spec_results %>%
  select(spec_id, attached_def, casual_def, n_attached, n_casual,
         attached_coef, attached_se, casual_coef, casual_se,
         differential, differential_se, differential_p,
         significant, effect_1sd_bdt) %>%
  mutate(
    attached_coef_se = sprintf("%.1f (%.1f)", attached_coef, attached_se),
    casual_coef_se     = sprintf("%.1f (%.1f)", casual_coef, casual_se),
    differential_se_fmt = sprintf("%.1f (%.1f)", differential, differential_se),
    p_fmt = fmt_p(differential_p),
    econ_mag = sprintf("%.1f", effect_1sd_bdt)
  )

pub_wide <- data.frame(row = pub_rows, stringsAsFactors = FALSE)
for (sid in spec_results$spec_id) {
  r <- pub_tbl %>% filter(spec_id == sid)
  pub_wide[[paste0("Spec_", sid)]] <- c(
    r$attached_def, r$casual_def,
    as.character(r$n_attached), as.character(r$n_casual),
    r$attached_coef_se, r$casual_coef_se,
    r$differential_se_fmt, r$p_fmt, r$significant, r$econ_mag
  )
}

write_csv(pub_wide, file.path(OUT, "tables/four_specification_table.csv"))

footnote <- paste(
  "Theory predicts attenuation under broader definitions because misclassification",
  "biases the differential toward zero. A smaller but still significant or",
  "directionally consistent differential under Specs B–D confirms the baseline",
  "(Spec A) rather than undermining it. Spec A uses institutionally unambiguous",
  "endpoints (Bardhan 1979; Binswanger and Rosenzweig 1984).",
  "Regression: diff_real_wage ~ diff_log_yield_hat x attached + gender |",
  "year + District x growing_season; SE clustered by district."
)

pub_html <- pub_wide %>%
  kbl(caption = "Robustness: attached/casual differential under alternative meal definitions") %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover")) %>%
  footnote(general = footnote)
writeLines(pub_html, file.path(OUT, "tables/four_specification_table.html"))

## ── Monotonicity gradient table ───────────────────────────────────────────────
log_msg("\n--- Monotonicity gradient ---")

sep_results <- bind_rows(lapply(MEAL_LEVELS, function(m) {
  sub <- df %>% filter(meal_type == m)
  mod <- feols(
    diff_real_wage ~ diff_log_yield_hat + gender | year + District^growing_season,
    data = sub, cluster = ~District
  )
  ct <- summary(mod)$coeftable
  data.frame(
    source       = "Separate regression",
    meal_type    = m,
    meal_label   = MEAL_LABELS[m],
    coefficient  = ct["diff_log_yield_hat", "Estimate"],
    std_error    = ct["diff_log_yield_hat", "Std. Error"],
    p_value      = ct["diff_log_yield_hat", "Pr(>|t|)"],
    n_obs        = nobs(mod),
    stringsAsFactors = FALSE
  )
}))

m_interact <- feols(
  diff_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal_type + gender |
    year + District^growing_season,
  data = df,
  cluster = ~District
)

V_int <- vcov(m_interact)
b_base <- coef(m_interact)["diff_log_yield_hat"]

int_results <- bind_rows(
  data.frame(
    source = "Interaction model (implied total)",
    meal_type = "None", meal_label = "0 meals (base)",
    coefficient = b_base,
    std_error = sqrt(V_int["diff_log_yield_hat", "diff_log_yield_hat"]),
    p_value = summary(m_interact)$coeftable["diff_log_yield_hat", "Pr(>|t|)"],
    n_obs = nobs(m_interact)
  ),
  lapply(c("One", "Two", "Three"), function(m) {
    int_nm <- paste0("diff_log_yield_hat:meal_type", m)
    est <- b_base + coef(m_interact)[int_nm]
    se  <- sqrt(V_int["diff_log_yield_hat", "diff_log_yield_hat"] +
                  V_int[int_nm, int_nm] +
                  2 * V_int["diff_log_yield_hat", int_nm])
    data.frame(
      source = "Interaction model (implied total)",
      meal_type = m,
      meal_label = paste0(MEAL_LABELS[m], " (base + interaction)"),
      coefficient = est,
      std_error = se,
      p_value = 2 * pt(-abs(est / se), df = m_interact$nobs - length(coef(m_interact))),
      n_obs = nobs(m_interact)
    )
  }) %>% bind_rows()
)

int_only <- lapply(c("One", "Two", "Three"), function(m) {
  int_nm <- paste0("diff_log_yield_hat:meal_type", m)
  ct <- summary(m_interact)$coeftable
  data.frame(
    source = "Interaction model (incremental)",
    meal_type = m,
    meal_label = paste0(MEAL_LABELS[m], " interaction"),
    coefficient = ct[int_nm, "Estimate"],
    std_error = ct[int_nm, "Std. Error"],
    p_value = ct[int_nm, "Pr(>|t|)"],
    n_obs = nobs(m_interact)
  )
}) %>% bind_rows()

mono_tbl <- bind_rows(sep_results, int_only, int_results) %>%
  mutate(
    coef_se = sprintf("%.1f (%.1f)", coefficient, std_error),
    p_fmt = fmt_p(p_value),
    sig = sig_flag(p_value),
    direction = direction_flag(coefficient)
  )

write_csv(mono_tbl, file.path(OUT, "tables/monotonicity_gradient.csv"))

mono_html <- mono_tbl %>%
  select(source, meal_label, coef_se, p_fmt, n_obs, sig, direction) %>%
  kbl(col.names = c("Source", "Meal category", "Coefficient (SE)", "p-value",
                    "N", "Significance", "Direction"),
      caption = "Monotonicity gradient: pass-through by meal category") %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))
writeLines(mono_html, file.path(OUT, "tables/monotonicity_gradient.html"))

for (i in seq_len(nrow(sep_results))) {
  r <- sep_results[i, ]
  log_msg(sprintf("  Separate %s: %.1f (SE %.1f, p=%s)",
                  r$meal_label, r$coefficient, r$std_error, fmt_p(r$p_value)))
}

## ── Formal test: 3-meal vs 2-meal interaction ───────────────────────────────
log_msg("\n--- Formal test: 3-meal vs 2-meal interaction ---")

test_3v2 <- wald_linear(
  m_interact,
  c("diff_log_yield_hat:meal_typeThree", "diff_log_yield_hat:meal_typeTwo"),
  c(1, -1)
)

b2_int <- coef(m_interact)["diff_log_yield_hat:meal_typeTwo"]
b3_int <- coef(m_interact)["diff_log_yield_hat:meal_typeThree"]

test_row <- data.frame(
  comparison = "3-meal interaction minus 2-meal interaction",
  coef_2meal = b2_int,
  coef_3meal = b3_int,
  difference = test_3v2$estimate,
  std_error = test_3v2$se,
  p_value = test_3v2$p,
  interpretation = ifelse(
    test_3v2$p < 0.05,
    ifelse(test_3v2$estimate < 0,
           "3-meal significantly more negative than 2-meal — justifies 3-meal proxy",
           "3-meal significantly different from 2-meal (unexpected sign)"),
    "Difference not significant — choice between 2-meal and 3-meal proxy does not substantially affect estimate (reassuring)"
  ),
  stringsAsFactors = FALSE
)

write_csv(test_row, file.path(OUT, "tables/test_3meal_vs_2meal.csv"))
log_msg(sprintf("  2-meal interaction: %.1f | 3-meal interaction: %.1f",
                b2_int, b3_int))
log_msg(sprintf("  Difference (3 − 2): %.1f (SE %.1f, p=%s)",
                test_3v2$estimate, test_3v2$se, fmt_p(test_3v2$p)))
log_msg("  ", test_row$interpretation)

## ── Coefficient plot: four specifications ───────────────────────────────────
plot_specs <- spec_results %>%
  mutate(
    spec_id = factor(spec_id, levels = c("A", "B", "C", "D")),
    ci_lo = differential - 1.96 * differential_se,
    ci_hi = differential + 1.96 * differential_se,
    sig_simple = ifelse(differential_p < 0.05, "p < 0.05",
                        ifelse(differential_p < 0.10, "p < 0.10", "n.s."))
  )

p_specs <- ggplot(plot_specs, aes(x = differential, y = spec_label, color = spec_id)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = PAPER_DIFF, linetype = "dotted", color = "#0072B2") +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.15) +
  labs(
    title = "Attached/casual differential under alternative meal definitions",
    subtitle = "Dotted line = paper benchmark (−222.6 BDT); theory predicts attenuation A → D",
    x = "Differential (attached − casual pass-through, BDT)",
    y = NULL,
    color = "Spec"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

ggsave(file.path(OUT, "figures/coefplot_four_specifications.png"),
       p_specs, width = 9, height = 5, dpi = 150)

## ── Summary flags ─────────────────────────────────────────────────────────────
summary_flags <- data.frame(
  item = c(
    "All specs directionally consistent (diff < 0)",
    "Any reversed sign",
    "Spec A significant at 5%",
    "Attenuation pattern (|diff| A >= B,C,D)",
    "3-meal vs 2-meal interaction significant"
  ),
  value = c(
    all(spec_results$differential < 0),
    any(spec_results$differential > 0),
    spec_results$differential_p[spec_results$spec_id == "A"] < 0.05,
    all(abs(spec_results$differential[1]) >= abs(spec_results$differential[-1])),
    test_3v2$p < 0.05
  ),
  stringsAsFactors = FALSE
)
write_csv(summary_flags, file.path(OUT, "tables/summary_flags.csv"))

log_msg("\n=== SUMMARY FLAGS ===")
for (i in seq_len(nrow(summary_flags))) {
  log_msg("  ", summary_flags$item[i], ": ", summary_flags$value[i])
}
log_msg("\nDone. Outputs: ", OUT)
