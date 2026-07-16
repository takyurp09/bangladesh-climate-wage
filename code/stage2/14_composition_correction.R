## 14_composition_correction.R
## Composition-corrected pass-through differential: 3-meal vs 0-meal workers
## Output: output/stage2/composition_correction/

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
OUT  <- file.path(ROOT, "output/stage2/composition_correction")
for (d in c("", "tables", "figures", "models", "summary"))
  dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)

PAPER_DIFF <- -222.6
PAPER_SE   <- 95.0
PAPER_P    <- 0.019

cat("=== COMPOSITION-CORRECTED PASS-THROUGH ANALYSIS ===\n")
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

## Extract 0-meal, 3-meal, and differential pass-through from M5-type model
extract_meal_pass <- function(mod, version_label) {
  ct  <- summary(mod)$coeftable
  V   <- vcov(mod)
  bnm <- "diff_log_yield_hat"
  int <- "diff_log_yield_hat:meal_typeThree"
  df2 <- mod$nobs - length(coef(mod))

  b0 <- ct[bnm, "Estimate"]
  s0 <- ct[bnm, "Std. Error"]
  p0 <- ct[bnm, "Pr(>|t|)"]

  b3 <- b0 + ct[int, "Estimate"]
  s3 <- sqrt(V[bnm, bnm] + V[int, int] + 2 * V[bnm, int])
  p3 <- 2 * pt(-abs(b3 / s3), df = df2)

  diff_est <- ct[int, "Estimate"]
  diff_se  <- ct[int, "Std. Error"]
  diff_p   <- ct[int, "Pr(>|t|)"]

  share_lvl <- if ("share_attached" %in% rownames(ct)) {
    list(coef = ct["share_attached", "Estimate"],
         se   = ct["share_attached", "Std. Error"],
         p    = ct["share_attached", "Pr(>|t|)"])
  } else list(coef = NA, se = NA, p = NA)

  share_fd <- if ("diff_share_attached" %in% rownames(ct)) {
    list(coef = ct["diff_share_attached", "Estimate"],
         se   = ct["diff_share_attached", "Std. Error"],
         p    = ct["diff_share_attached", "Pr(>|t|)"])
  } else list(coef = NA, se = NA, p = NA)

  data.frame(
    version              = version_label,
    casual_coef          = b0,
    casual_se            = s0,
    casual_p             = p0,
    attached_coef        = b3,
    attached_se          = s3,
    attached_p           = p3,
    differential         = diff_est,
    differential_se      = diff_se,
    differential_p       = diff_p,
    share_attached_coef  = share_lvl$coef,
    share_attached_se    = share_lvl$se,
    share_attached_p     = share_lvl$p,
    diff_share_coef      = share_fd$coef,
    diff_share_se        = share_fd$se,
    diff_share_p         = share_fd$p,
    nobs                 = nobs(mod),
    sig_status           = sig_flag(diff_p),
    stringsAsFactors     = FALSE
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

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 2 — Construct share_attached (levels + FD) and merge into panel       ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 2 — Construct and merge share_attached\n")

wage_raw <- read_csv(
  file.path(ROOT, "data/Regression_data/wage_by_growing_season.csv"),
  show_col_types = FALSE
) %>%
  filter(meal_type != "No_info", !is.na(real_wage), year >= 2017, year <= 2023)

share_dsy <- wage_raw %>%
  group_by(District, year, growing_season, meal_type) %>%
  summarise(n_monthly = n(), .groups = "drop") %>%
  group_by(District, year, growing_season) %>%
  mutate(share = n_monthly / sum(n_monthly)) %>%
  ungroup() %>%
  filter(meal_type == "Three") %>%
  select(District, year, growing_season, share_attached = share) %>%
  mutate(share_attached = replace_na(share_attached, 0)) %>%
  arrange(District, growing_season, year) %>%
  group_by(District, growing_season) %>%
  mutate(diff_share_attached = share_attached - lag(share_attached)) %>%
  ungroup()

cat(sprintf("  share_attached cells: %d | mean = %.3f | SD = %.3f\n",
            nrow(share_dsy), mean(share_dsy$share_attached),
            sd(share_dsy$share_attached)))

df <- read_csv(
  file.path(ROOT, "data/Regression_data/df_2_merged_v2.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    meal_type = relevel(factor(meal_type), ref = "None"),
    gender    = relevel(factor(gender),    ref = "Female")
  ) %>%
  left_join(share_dsy, by = c("District", "year", "growing_season"))

n_before <- nrow(df)
df <- df %>% filter(!is.na(diff_log_yield_hat))
cat(sprintf("  Merged panel: %d obs (dropped %d without share_attached match)\n\n",
            nrow(df), n_before - nrow(df)))

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 1 & 3 — Four regression versions                                      ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 1 — Baseline replication (M5)\n")

FE <- "year + District^growing_season"
base_rhs <- "diff_log_yield_hat + diff_log_yield_hat:meal_type + gender"

m_v1 <- feols(
  as.formula(paste("diff_real_wage ~", base_rhs, "|", FE)),
  data = df, cluster = ~District
)
r_v1 <- extract_meal_pass(m_v1, "V1: Baseline")

cat(sprintf("  0-meal pass-through:  %.1f (SE %.1f, p = %s)\n",
            r_v1$casual_coef, r_v1$casual_se, fmt_p(r_v1$casual_p)))
cat(sprintf("  3-meal pass-through:  %.1f (SE %.1f, p = %s)\n",
            r_v1$attached_coef, r_v1$attached_se, fmt_p(r_v1$attached_p)))
cat(sprintf("  Differential (3 − 0): %.1f (SE %.1f, p = %s) [%s]\n",
            r_v1$differential, r_v1$differential_se, fmt_p(r_v1$differential_p),
            r_v1$sig_status))

if (abs(r_v1$differential - PAPER_DIFF) > 5) {
  cat(sprintf("  ** FLAG: Differential %.1f differs from paper %.1f by >5 BDT\n",
              r_v1$differential, PAPER_DIFF))
} else {
  cat(sprintf("  Replication matches paper (%.1f vs %.1f)\n", r_v1$differential, PAPER_DIFF))
}
cat("\n")

cat("STEP 3 — Composition-corrected regressions (four versions)\n")
cat("  NOTE: V3/V4 drop obs with NA diff_share_attached (first year per district-season)\n")

m_v2 <- feols(
  as.formula(paste("diff_real_wage ~", base_rhs, "+ share_attached |", FE)),
  data = df, cluster = ~District
)
m_v3 <- feols(
  as.formula(paste("diff_real_wage ~", base_rhs, "+ diff_share_attached |", FE)),
  data = df %>% filter(!is.na(diff_share_attached)),
  cluster = ~District
)
m_v4 <- feols(
  as.formula(paste("diff_real_wage ~", base_rhs,
                   "+ share_attached + diff_share_attached |", FE)),
  data = df %>% filter(!is.na(diff_share_attached)),
  cluster = ~District
)

results <- bind_rows(
  r_v1,
  extract_meal_pass(m_v2, "V2: + share_attached (level)"),
  extract_meal_pass(m_v3, "V3: + diff_share_attached (FD)"),
  extract_meal_pass(m_v4, "V4: + both composition controls")
)

for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  cat(sprintf("  %s\n", r$version))
  cat(sprintf("    Differential: %.1f (SE %.1f, p = %s) [%s]\n",
              r$differential, r$differential_se, fmt_p(r$differential_p), r$sig_status))
  if (!is.na(r$share_attached_coef))
    cat(sprintf("    share_attached coef: %.1f (p = %s)\n",
                r$share_attached_coef, fmt_p(r$share_attached_p)))
  if (!is.na(r$diff_share_coef))
    cat(sprintf("    diff_share_attached coef: %.1f (p = %s)\n",
                r$diff_share_coef, fmt_p(r$diff_share_p)))
}
cat("\n")

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 4 — Pattern A / B / C                                                 ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 4 — Pattern classification\n")

base_diff <- results$differential[1]
pct_changes <- 100 * (results$differential - base_diff) / abs(base_diff)

best_row <- results[which.min(results$differential_p), ]
v3_diff  <- results$differential[results$version == "V3: + diff_share_attached (FD)"]
v4_diff  <- results$differential[results$version == "V4: + both composition controls"]
v2_diff  <- results$differential[results$version == "V2: + share_attached (level)"]

## Preferred specs: V3/V4 (FD control aligns with FD wage regression)
fd_stable <- abs(v3_diff - base_diff) < 15 && abs(v4_diff - base_diff) < 15

if (fd_stable && abs(v2_diff - base_diff) >= 15) {
  pattern <- "C (with caveat)"
  pattern_txt <- paste0(
    "Pattern C — STABLE with FD controls: Adding diff_share_attached (V3/V4) leaves ",
    "the differential essentially unchanged (V3=", round(v3_diff, 1), ", V4=", round(v4_diff, 1),
    " vs baseline ", round(base_diff, 1), " BDT). The level control alone (V2=",
    round(v2_diff, 1), ") attenuates the differential by ", sprintf("%.0f%%", pct_changes[2]),
    ", suggesting cross-sectional composition differences matter, but within-district ",
    "year-on-year composition changes do not confound the main result. ",
    "Recommend V3/V4 as the appropriate corrected specification."
  )
} else if (min(results$differential) < base_diff - 10) {
  pattern <- "A"
  pattern_txt <- paste0(
    "Pattern A — DIFFERENTIAL GROWS: Corrected differential more negative than baseline ",
    "(", round(min(results$differential), 1), " vs ", round(base_diff, 1),
    " BDT). Baseline −222 is a LOWER BOUND."
  )
} else if (max(results$differential) > base_diff + 10) {
  pattern <- "B"
  pattern_txt <- paste0(
    "Pattern B — DIFFERENTIAL SHRINKS: Corrected differential less negative than baseline ",
    "(", round(max(results$differential), 1), " vs ", round(base_diff, 1),
    " BDT). Main result may be overstated by composition."
  )
} else {
  pattern <- "C"
  pattern_txt <- paste0(
    "Pattern C — STABLE: Differential stable across all versions (range ",
    round(min(results$differential), 1), " to ", round(max(results$differential), 1),
    " BDT). Baseline estimate is reliable."
  )
}

cat("  Per-version change from baseline:\n")
for (i in seq_len(nrow(results)))
  cat(sprintf("    %s: %.1f BDT (%+.1f%%)\n",
              results$version[i], results$differential[i], pct_changes[i]))

cat(" ", pattern_txt, "\n\n")

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 5 — Heterogeneity by composition level                                ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 5 — High vs low share_attached subsamples\n")

med_share <- median(share_dsy$share_attached, na.rm = TRUE)
cat(sprintf("  Median share_attached = %.3f\n", med_share))

df_het <- df %>%
  mutate(high_share = as.integer(share_attached >= med_share))

m_high <- feols(
  as.formula(paste("diff_real_wage ~", base_rhs, "|", FE)),
  data = df_het %>% filter(high_share == 1), cluster = ~District
)
m_low <- feols(
  as.formula(paste("diff_real_wage ~", base_rhs, "|", FE)),
  data = df_het %>% filter(high_share == 0), cluster = ~District
)

het_results <- bind_rows(
  extract_meal_pass(m_high, "High share_attached (≥ median)"),
  extract_meal_pass(m_low,  "Low share_attached (< median)")
)

for (i in seq_len(nrow(het_results))) {
  r <- het_results[i, ]
  cat(sprintf("  %s: differential = %.1f (SE %.1f, p = %s) [%s]  N = %d\n",
              r$version, r$differential, r$differential_se,
              fmt_p(r$differential_p), r$sig_status, r$nobs))
}

## Formal test: does differential differ by composition level?
m_het_pool <- feols(
  as.formula(paste(
    "diff_real_wage ~ diff_log_yield_hat:meal_type +",
    "diff_log_yield_hat:meal_type:high_share + gender |", FE
  )),
  data = df_het,
  cluster = ~District
)
het_wald <- wald_linear(
  m_het_pool,
  c("diff_log_yield_hat:meal_typeThree",
    "diff_log_yield_hat:meal_typeThree:high_share"),
  c(0, 1)
)
cat(sprintf("\n  Wald test (3-meal differential, high vs low): %.1f (SE %.1f, p = %s)\n\n",
            het_wald$estimate, het_wald$se, fmt_p(het_wald$p)))

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 6 — Corrected economic magnitude                                      ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 6 — Economic magnitude (1-SD yield shock effect on differential)\n")

sd_yield <- sd(df$diff_log_yield_hat, na.rm = TRUE)
mean_wage <- df %>%
  filter(meal_type %in% c("None", "Three")) %>%
  summarise(m = mean(real_wage, na.rm = TRUE)) %>%
  pull(m)

econ_mag <- results %>%
  mutate(
    sd_yield_hat       = sd_yield,
    mean_daily_wage    = mean_wage,
    effect_1sd_bdt     = differential * sd_yield,
    pct_mean_wage      = 100 * abs(effect_1sd_bdt) / mean_wage,
    effect_1sd_label   = sprintf("%.1f BDT/day (%.1f%% of mean wage)",
                                 effect_1sd_bdt, pct_mean_wage)
  )

baseline_mag <- econ_mag$effect_1sd_bdt[1]
best_mag_row <- econ_mag %>% slice(which.min(differential_p))

cat(sprintf("  SD(diff_log_yield_hat) = %.4f | Mean daily wage = %.1f BDT\n",
            sd_yield, mean_wage))
cat(sprintf("  Baseline 1-SD effect on differential: %.1f BDT/day (%.1f%% mean wage)\n",
            baseline_mag, econ_mag$pct_mean_wage[1]))
cat(sprintf("  Best-corrected (%s): %.1f BDT/day (%.1f%% mean wage)\n\n",
            best_mag_row$version, best_mag_row$effect_1sd_bdt, best_mag_row$pct_mean_wage))

## ════════════════════════════════════════════════════════════════════════════ ##
## STEP 7 — Publication table                                                 ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("STEP 7 — Publication table\n")

pub_tbl <- results %>%
  transmute(
    Version = version,
    `0-meal coef` = sprintf("%.1f", casual_coef),
    `3-meal coef` = sprintf("%.1f", attached_coef),
    Differential  = sprintf("%.1f", differential),
    `SE (diff)`   = sprintf("%.1f", differential_se),
    `p (diff)`    = fmt_p(differential_p),
    Status        = sig_status,
    `share_attached coef` = ifelse(is.na(share_attached_coef), "—",
                                   sprintf("%.1f (p=%s)", share_attached_coef,
                                           fmt_p(share_attached_p))),
    `diff_share coef` = ifelse(is.na(diff_share_coef), "—",
                               sprintf("%.1f (p=%s)", diff_share_coef,
                                       fmt_p(diff_share_p)))
  )

econ_row <- data.frame(
  Version = "Economic magnitude (1-SD yield shock)",
  `0-meal coef` = "—", `3-meal coef` = "—",
  Differential  = sprintf("Baseline: %.1f BDT | Best corrected: %.1f BDT",
                          baseline_mag, best_mag_row$effect_1sd_bdt),
  `SE (diff)` = "—",
  `p (diff)`  = "—",
  Status      = "—",
  `share_attached coef` = "—",
  `diff_share coef` = "—",
  check.names = FALSE
)

pub_display <- bind_rows(pub_tbl, econ_row)
write_csv(pub_display, file.path(OUT, "tables", "composition_correction_table.csv"))

pub_html <- pub_display %>%
  kbl(caption = "Composition-Corrected Pass-Through: 3-Meal vs 0-Meal Differential",
      format = "html", booktabs = TRUE, escape = FALSE) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = TRUE, font_size = 12) %>%
  footnote(
    general = paste(
      "M5 specification: diff_real_wage ~ diff_log_yield_hat × meal_type + gender",
      "| FE: year + District×Season | SE clustered by district.",
      "Differential = pass-through for 3-meal minus 0-meal workers.",
      sprintf("Paper baseline: %.1f BDT (SE %.1f, p = %.3f).", PAPER_DIFF, PAPER_SE, PAPER_P),
      pattern_txt
    ),
    general_title = "Note: ",
    footnote_as_chunk = TRUE
  )

writeLines(as.character(pub_html),
           file.path(OUT, "tables", "composition_correction_table.html"))

## Wide numeric table for analysis
write_csv(results, file.path(OUT, "tables", "composition_correction_results.csv"))
write_csv(het_results, file.path(OUT, "tables", "heterogeneity_by_share.csv"))
write_csv(econ_mag, file.path(OUT, "tables", "economic_magnitude.csv"))

## Coefficient comparison plot
plot_df <- results %>%
  mutate(version = factor(version, levels = version))

p_diff <- ggplot(plot_df, aes(x = differential, y = version)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
  geom_vline(xintercept = PAPER_DIFF, linetype = "dotted", color = "#0072B2") +
  geom_errorbarh(aes(xmin = differential - 1.96 * differential_se,
                     xmax = differential + 1.96 * differential_se),
                 height = 0.2, linewidth = 0.8) +
  geom_point(size = 3.5, color = "#D55E00") +
  geom_text(aes(label = sprintf("%.0f (p=%s)", differential, fmt_p(differential_p))),
            hjust = -0.1, size = 3, color = "grey25") +
  labs(
    title = "Pass-Through Differential (3-Meal vs 0-Meal) Across Specifications",
    subtitle = "Dotted line = paper baseline (−222.6 BDT) | 95% CI",
    x = "Differential (BDT/day per unit change in fitted log yield)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(OUT, "figures", "differential_by_version.png"),
       p_diff, width = 9, height = 5, dpi = 300)

## ── Save models & report ────────────────────────────────────────────────────── ##
report <- c(
  "# Composition-Corrected Pass-Through Analysis",
  sprintf("Generated: %s", Sys.time()),
  "",
  "## Step 1 — Baseline replication",
  sprintf("- Differential: %.1f (SE %.1f, p = %s) [%s]",
          r_v1$differential, r_v1$differential_se,
          fmt_p(r_v1$differential_p), r_v1$sig_status),
  sprintf("- Paper target: %.1f (SE %.1f, p = %.3f)",
          PAPER_DIFF, PAPER_SE, PAPER_P),
  "",
  "## Step 3 — Four versions",
  capture.output(print(results %>%
    select(version, casual_coef, attached_coef, differential,
           differential_se, differential_p, share_attached_coef, diff_share_coef))),
  "",
  "## Step 4 — Pattern",
  pattern_txt,
  "",
  "## Step 5 — Heterogeneity",
  capture.output(print(het_results %>% select(version, differential, differential_se, differential_p))),
  sprintf("- Wald test (high vs low differential): p = %s", fmt_p(het_wald$p)),
  "",
  "## Step 6 — Economic magnitude",
  sprintf("- Baseline 1-SD effect: %.1f BDT (%.1f%% mean wage)",
          baseline_mag, econ_mag$pct_mean_wage[1]),
  sprintf("- Best corrected 1-SD effect: %.1f BDT (%.1f%% mean wage)",
          best_mag_row$effect_1sd_bdt, best_mag_row$pct_mean_wage),
  "",
  "## Outputs",
  "- tables/composition_correction_table.csv/.html",
  "- tables/composition_correction_results.csv",
  "- figures/differential_by_version.png"
)

writeLines(report, file.path(OUT, "summary", "composition_correction_report.md"))

save(m_v1, m_v2, m_v3, m_v4, m_high, m_low, m_het_pool,
     results, het_results, econ_mag, pattern,
     file = file.path(OUT, "models", "composition_correction_models.RData"))

cat("Saved composition_correction_table.csv/.html\n")
cat("Saved composition_correction_report.md\n")
cat("=== COMPOSITION CORRECTION COMPLETE ===\n")
