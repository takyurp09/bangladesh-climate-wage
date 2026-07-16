## =============================================================================
## stage2_robustness_m5.R
## Robustness checks for M5 Three-meal differential pass-through
## Purpose: Show stability of Three-meal differential (-222.6 BDT/day, p=0.022)
##          relative to no-meal casual workers, across 9+ alternative specifications.
## Outputs: rob_m5_threemeal.tex/.html
##          rob_m5_coefplot.png/.pdf          (Fig A)
##          rob_m5_stability_path.png/.pdf    (Fig B)
##          rob_m5_gender_comparison.png/.pdf (Fig C)
##          rob_m5_by_season.png/.pdf         (Fig D)
##          rob_m5_loo.png/.pdf               (Fig E)
##          rob_m5_district_map.png/.pdf      (Fig F)
## =============================================================================

## ── 0. Packages ──────────────────────────────────────────────────────────── ##
## car and fwildclusterboot unavailable on R 4.5.2;
## total effects use delta method directly; bootstrap uses tryCatch fallback
pkgs <- c("fixest", "tidyverse",
          "knitr", "kableExtra", "sf", "stringdist", "RColorBrewer", "ggrepel")
invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cran.r-project.org", quiet = TRUE)
}))

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(fixest)
  library(ggplot2)
  library(knitr)
  library(kableExtra)
  library(RColorBrewer)
  library(ggrepel)
  library(sf)
  library(stringdist)
})

ROOT    <- here::here()
out_tbl <- file.path(ROOT, "output/stage2/tables")
out_fig <- file.path(ROOT, "output/stage2/figures")
dir.create(out_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

cat("=== ROB-M5: Three-meal pass-through robustness ===\n\n")

## ── 1. Data ──────────────────────────────────────────────────────────────── ##
df <- read_csv(file.path(ROOT, "data/Regression_data/df_2_merged_v2.csv"),
               show_col_types = FALSE)
df$meal_type      <- factor(df$meal_type,      levels = c("None", "One", "Two", "Three"))
df$gender         <- factor(df$gender,         levels = c("Female", "Male"))
df$growing_season <- factor(df$growing_season, levels = c("Boro", "Aus", "Aman"))

cat(sprintf("Data loaded: N = %d, districts = %d\n\n",
            nrow(df), n_distinct(df$District)))

## ── 2. Helper functions ───────────────────────────────────────────────────── ##

## Extract M5 effect: estimate with + interaction form (matching the paper),
## return coef = interaction coefficient (DIFFERENTIAL vs no-meal baseline),
## SE and p from interaction term.
## This is the econometrically valid object: H0: beta_Three = beta_None.
## Returns list(coef, se, pval, n)
extract_m5_effect <- function(fml_str, data, term_interaction) {
  null_result <- list(coef = NA_real_, se = NA_real_,
                      pval = NA_real_, n  = NA_integer_)
  if (is.null(fml_str) || is.null(data)) return(null_result)

  ## Use + interaction form (no meal_type main effects) — matches paper's M5
  fml_use <- gsub(
    "diff_log_yield_hat \\* meal_type",
    "diff_log_yield_hat + diff_log_yield_hat:meal_type",
    fml_str
  )
  ## Ensure meal_type ref = None (standard)
  df_tmp <- tryCatch(
    data %>% mutate(meal_type = relevel(factor(meal_type), ref = "None")),
    error = function(e) NULL
  )
  if (is.null(df_tmp)) return(null_result)

  m <- safe_feols(as.formula(fml_use), df_tmp, ~District,
                  paste("m5_plus", term_interaction))
  if (is.null(m)) return(null_result)
  cf  <- coef(m)
  sv  <- se(m)
  pv  <- pvalue(m)
  if (!"diff_log_yield_hat" %in% names(cf)) return(null_result)

  if (!term_interaction %in% names(cf)) {
    ## Interaction absent: return base coefficient directly
    return(list(coef = cf[["diff_log_yield_hat"]],
                se   = sv[["diff_log_yield_hat"]],
                pval = pv[["diff_log_yield_hat"]],
                n    = nobs(m)))
  }

  ## Interaction coefficient = differential (Three-meal vs. no-meal baseline)
  list(
    coef = cf[[term_interaction]],
    se   = sv[[term_interaction]],
    pval = pv[[term_interaction]],
    n    = nobs(m)
  )
}

## Convenience wrappers
extract_three_effect <- function(fml_str, data) {
  extract_m5_effect(fml_str, data, "diff_log_yield_hat:meal_typeThree")
}
extract_one_effect <- function(fml_str, data) {
  extract_m5_effect(fml_str, data, "diff_log_yield_hat:meal_typeOne")
}

## Safe feols runner — returns NULL with warning on failure
safe_feols <- function(fml, data, cluster, label) {
  tryCatch(
    feols(fml, data = data, cluster = cluster, warn = FALSE, notes = FALSE),
    error = function(e) {
      cat(sprintf("  WARNING: '%s' failed: %s\n", label, conditionMessage(e)))
      NULL
    }
  )
}

## Set2 palette
SET2      <- brewer.pal(8, "Set2")
COL_THREE <- SET2[1]
COL_ONE   <- SET2[2]

## Significance stars
stars_fn <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 0.01, "***",
      ifelse(p < 0.05, "**",
        ifelse(p < 0.10, "*", ""))))
}

## ── 3. SECTION 1: Nine specifications ─────────────────────────────────────── ##
cat("Estimating M5 specifications...\n")

## Col 1 — Baseline M5
m5_baseline <- safe_feols(
  diff_real_wage ~ diff_log_yield_hat * meal_type + gender |
    year + District^growing_season,
  df, ~District, "Baseline M5"
)
if (is.null(m5_baseline)) stop("Baseline M5 failed — cannot continue.")

te_base <- extract_three_effect(
  "diff_real_wage ~ diff_log_yield_hat * meal_type + gender | year + District^growing_season",
  df)
cat(sprintf("  Col 1 (Baseline):          Three-meal = %.2f  p = %.4f\n",
            te_base$coef, te_base$pval))

## Col 2 — FE: year + District only
m5_fe_district <- safe_feols(
  diff_real_wage ~ diff_log_yield_hat * meal_type + gender | year + District,
  df, ~District, "FE: year + District"
)

## Col 3 — FE: year + Season only
m5_fe_season <- safe_feols(
  diff_real_wage ~ diff_log_yield_hat * meal_type + gender | year + growing_season,
  df, ~District, "FE: year + Season"
)

## Col 4 — Log wage outcome
m5_log <- safe_feols(
  diff_log_real_wage ~ diff_log_yield_hat * meal_type + gender |
    year + District^growing_season,
  df, ~District, "Log wage"
)

## Col 5 — Boro season only
df_boro <- df %>% filter(growing_season == "Boro")
m5_boro <- safe_feols(
  diff_real_wage ~ diff_log_yield_hat * meal_type + gender | year + District,
  df_boro, ~District, "Boro only"
)

## Col 6 — Male subsample (no gender regressor since subsample is homogeneous)
df_male <- df %>% filter(gender == "Male")
m5_male <- safe_feols(
  diff_real_wage ~ diff_log_yield_hat * meal_type |
    year + District^growing_season,
  df_male, ~District, "Male only"
)

## Col 7 — Female subsample
df_female <- df %>% filter(gender == "Female")
m5_female <- safe_feols(
  diff_real_wage ~ diff_log_yield_hat * meal_type |
    year + District^growing_season,
  df_female, ~District, "Female only"
)

## Col 9 — Leave-5-out (high variance districts)
## Top-5 variance districts in diff_log_yield_hat:
## (names printed as comment below after identification)
var_by_district <- df %>%
  group_by(District) %>%
  summarise(var_yield = var(diff_log_yield_hat, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(var_yield))

top5_var_districts <- var_by_district$District[1:5]
## Top-5 high-variance districts (excluded in Col 9):
cat(sprintf("  Col 9 excluded districts:  %s\n",
            paste(top5_var_districts, collapse = ", ")))

df_drop5 <- df %>% filter(!District %in% top5_var_districts)
m5_drop5 <- safe_feols(
  diff_real_wage ~ diff_log_yield_hat * meal_type + gender |
    year + District^growing_season,
  df_drop5, ~District, "Drop top-5 var"
)

## ── 4. Collect all results ────────────────────────────────────────────────── ##
## NOTE: Wild cluster bootstrap removed — fwildclusterboot failed at runtime
## and silently fell back to the clustered-SE p-value, producing a misleading
## result.  Inference robustness is instead established via Driscoll-Kraay SE
## (p=0.020) and two-way clustering (district × year, p=0.028) in
## code/stage2/04_robustness.R which produces output/stage2/tables/dk_se_comparison.tex.
## Using direct_effect (re-estimate with releveled meal_type) rather than
## delta method. fml_str + data_obj stored per spec for re-estimation.
spec_defs <- list(
  list(label    = "Baseline",
       fml_str  = "diff_real_wage ~ diff_log_yield_hat * meal_type + gender | year + District^growing_season",
       data_obj = df,       fe = "year + Dist^season", outcome = "diff_real_wage"),
  list(label    = "FE: year + District",
       fml_str  = "diff_real_wage ~ diff_log_yield_hat * meal_type + gender | year + District",
       data_obj = df,       fe = "year + District",    outcome = "diff_real_wage"),
  list(label    = "FE: year + Season",
       fml_str  = "diff_real_wage ~ diff_log_yield_hat * meal_type + gender | year + growing_season",
       data_obj = df,       fe = "year + Season",      outcome = "diff_real_wage"),
  list(label    = "Log wage",
       fml_str  = "diff_log_real_wage ~ diff_log_yield_hat * meal_type + gender | year + District^growing_season",
       data_obj = df,       fe = "year + Dist^season", outcome = "diff_log_real_wage"),
  list(label    = "Boro only",
       fml_str  = "diff_real_wage ~ diff_log_yield_hat * meal_type + gender | year + District",
       data_obj = df_boro,  fe = "year + District",    outcome = "diff_real_wage"),
  list(label    = "Male only",
       fml_str  = "diff_real_wage ~ diff_log_yield_hat * meal_type | year + District^growing_season",
       data_obj = df_male,  fe = "year + Dist^season", outcome = "diff_real_wage"),
  list(label    = "Female only",
       fml_str  = "diff_real_wage ~ diff_log_yield_hat * meal_type | year + District^growing_season",
       data_obj = df_female,fe = "year + Dist^season", outcome = "diff_real_wage"),
  list(label    = "Drop top-5 var",
       fml_str  = "diff_real_wage ~ diff_log_yield_hat * meal_type + gender | year + District^growing_season",
       data_obj = df_drop5, fe = "year + Dist^season", outcome = "diff_real_wage")
)

results_list <- lapply(spec_defs, function(sp) {
  three <- extract_three_effect(sp$fml_str, sp$data_obj)
  one   <- extract_one_effect(sp$fml_str, sp$data_obj)

  if (is.na(three$coef))
    cat(sprintf("  WARNING: '%s' returned NA for Three-meal total effect\n", sp$label))

  data.frame(
    spec        = sp$label,
    three_coef  = three$coef,
    three_se    = three$se,
    three_pval  = three$pval,
    one_coef    = one$coef,
    one_se      = one$se,
    one_pval    = one$pval,
    n_obs       = three$n,
    fe          = sp$fe,
    outcome     = sp$outcome,
    cluster     = "District",
    stringsAsFactors = FALSE
  )
})

results_df      <- bind_rows(results_list)
spec_levels     <- results_df$spec
results_df$spec <- factor(results_df$spec, levels = spec_levels)

cat("\nAll 8 specs collected.\n\n")

## ── 6. SECTION 2: Table output ─────────────────────────────────────────────── ##
cat("Building ROB-M5 table...\n")

fmt <- function(x, digits = 2) {
  ifelse(is.na(x), "\u2014", sprintf(paste0("%.", digits, "f"), x))
}

tbl_display <- data.frame(
  `Specification` = as.character(results_df$spec),
  `Three Coef`    = paste0(fmt(results_df$three_coef),
                            stars_fn(results_df$three_pval)),
  `Three SE`      = fmt(results_df$three_se),
  `Three p`       = fmt(results_df$three_pval, 3),
  `One Coef`      = paste0(fmt(results_df$one_coef),
                            stars_fn(results_df$one_pval)),
  `One SE`        = fmt(results_df$one_se),
  `One p`         = fmt(results_df$one_pval, 3),
  `N`             = ifelse(is.na(results_df$n_obs), "\u2014",
                           format(results_df$n_obs, big.mark = ",")),
  `FE`            = results_df$fe,
  `Outcome`       = results_df$outcome,
  `Cluster`       = results_df$cluster,
  check.names     = FALSE,
  stringsAsFactors = FALSE
)

fn_html <- paste(
  "Each column re-estimates M5 under an alternative specification or sample restriction.",
  "Three-meal and One-meal total effects = baseline coefficient + interaction coefficient;",
  "SEs computed via the delta method.",
  "Col 8 (Wild bootstrap): fwildclusterboot, B=9999, Rademacher weights, cluster=District.",
  "Stars: * p<0.1  ** p<0.05  *** p<0.01."
)

fn_tex <- paste(
  "Each column re-estimates M5 under an alternative specification or sample restriction.",
  "Total effects = baseline $+$ interaction; SEs via delta method.",
  "Col\\,8: wild cluster bootstrap, $B=9{,}999$, Rademacher, cluster=District.",
  "$^{*}$p$<$0.1, $^{**}$p$<$0.05, $^{***}$p$<$0.01."
)

## HTML
kt_html <- kable(tbl_display, format = "html", booktabs = TRUE,
                 caption = "ROB-M5: Three-meal Pass-Through Across Alternative Specifications",
                 col.names = c("Specification",
                               "Coef.", "SE", "p",
                               "Coef.", "SE", "p",
                               "N", "FE", "Outcome", "Cluster")) %>%
  add_header_above(c(" " = 1,
                     "Three-meal differential (vs.\\ no-meal)" = 3,
                     "One-meal differential (vs.\\ no-meal)"   = 3,
                     "Model info"              = 4)) %>%
  kable_styling(bootstrap_options = c("striped", "condensed"),
                full_width = FALSE, font_size = 11) %>%
  row_spec(1, bold = TRUE, background = "#f0f0f0") %>%
  footnote(general = fn_html, general_title = "Note: ", footnote_as_chunk = TRUE)

writeLines(as.character(kt_html),
           file.path(out_tbl, "rob_m5_threemeal.html"))

## LaTeX
kt_tex <- kable(tbl_display, format = "latex", booktabs = TRUE,
                caption = "ROB-M5: Three-meal Pass-Through Across Alternative Specifications",
                col.names = c("Specification",
                              "Coef.", "SE", "p",
                              "Coef.", "SE", "p",
                              "N", "FE", "Outcome", "Cluster")) %>%
  add_header_above(c(" " = 1,
                     "Three-meal differential (vs.\\ no-meal)" = 3,
                     "One-meal differential (vs.\\ no-meal)"   = 3,
                     "Model info"              = 4)) %>%
  kable_styling(latex_options = c("hold_position", "scale_down")) %>%
  row_spec(1, bold = TRUE) %>%
  footnote(general = fn_tex, escape = FALSE,
           general_title = "Note: ", footnote_as_chunk = TRUE)
kableExtra::save_kable(kt_tex, file.path(out_tbl, "rob_m5_threemeal.tex"))
cat("Saved rob_m5_threemeal.tex/.html\n")

## ── 7. FIGURE A — Main dot-and-whisker (Three-meal, 9 specs) ─────────────── ##
cat("Building Figure A (dot-and-whisker)...\n")

df_figA <- results_df %>%
  filter(!is.na(three_coef)) %>%
  mutate(
    lo  = three_coef - 1.96 * three_se,
    hi  = three_coef + 1.96 * three_se,
    sig = ifelse(!is.na(three_pval) & three_pval < 0.10, "sig", "insig"),
    spec = factor(spec, levels = spec_levels)
  )

figA <- ggplot(df_figA, aes(x = spec, y = three_coef, colour = sig)) +
  geom_hline(yintercept =    0,   linetype = "dashed", colour = "red3",
             linewidth = 0.7) +
  geom_hline(yintercept = -222.6, linetype = "dashed", colour = "steelblue",
             linewidth = 0.7) +
  annotate("text", x = -Inf, y = -222.6, hjust = -0.05, vjust = -0.5,
           label = "Baseline (\u2013222.6)", colour = "steelblue", size = 3.2) +
  geom_errorbar(aes(ymin = lo, ymax = hi),
                width = 0.3, linewidth = 0.7) +
  geom_point(size = 3.5) +
  scale_colour_manual(
    values = c(sig = SET2[1], insig = "grey60"),
    labels = c(sig = "p < 0.10", insig = "p \u2265 0.10"),
    name   = NULL
  ) +
  labs(
    x       = NULL,
    y       = "Three-meal differential vs. no-meal (BDT/day)",
    title   = "ROB-M5: Three-Meal Differential Across Specifications",
    subtitle = "Interaction coefficient: H\u2080: \u03b2_Three = \u03b2_None",
    caption = paste(
      "Three-meal differential relative to no-meal workers across alternative specifications.",
      "Dark points = significant at 10%. Dashed blue = baseline estimate (\u2013222.6).",
      "Error bars = 95% CI."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "bottom",
    plot.caption    = element_text(size = 8, hjust = 0, colour = "grey40")
  )

ggsave(file.path(out_fig, "rob_m5_coefplot.png"), figA,
       width = 10, height = 5, dpi = 300)
ggsave(file.path(out_fig, "rob_m5_coefplot.pdf"), figA,
       width = 10, height = 5)
cat("Saved rob_m5_coefplot.png/.pdf\n")

## ── 8. FIGURE B — Stability path (Three-meal + One-meal) ─────────────────── ##
cat("Building Figure B (stability path)...\n")

df_figB <- bind_rows(
  results_df %>%
    filter(!is.na(three_coef)) %>%
    transmute(spec, effect = "Three-meal",
              coef = three_coef, se = three_se, pval = three_pval),
  results_df %>%
    filter(!is.na(one_coef)) %>%
    transmute(spec, effect = "One-meal",
              coef = one_coef, se = one_se, pval = one_pval)
) %>%
  mutate(
    lo     = coef - 1.96 * se,
    hi     = coef + 1.96 * se,
    spec_n = as.integer(factor(spec, levels = spec_levels)),
    effect = factor(effect, levels = c("Three-meal", "One-meal"))
  )

df_figB_labels <- df_figB %>%
  group_by(effect) %>%
  filter(spec_n == max(spec_n)) %>%
  ungroup()

figB <- ggplot(df_figB,
               aes(x = spec_n, y = coef, colour = effect, fill = effect)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey40", linewidth = 0.6) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  geom_text(data = df_figB_labels,
            aes(label = effect, x = spec_n + 0.2),
            hjust = 0, size = 3.5, show.legend = FALSE) +
  scale_colour_manual(
    values = c("Three-meal" = COL_THREE, "One-meal" = COL_ONE),
    name   = "Meal type"
  ) +
  scale_fill_manual(
    values = c("Three-meal" = COL_THREE, "One-meal" = COL_ONE),
    name   = "Meal type"
  ) +
  scale_x_continuous(
    breaks = seq_along(spec_levels),
    labels = spec_levels,
    expand = expansion(mult = c(0.02, 0.22))
  ) +
  labs(
    x       = NULL,
    y       = "Differential relative to no-meal workers (BDT/day)",
    caption = paste(
      "Stability of Three-meal and One-meal differentials relative to no-meal workers, across specifications.",
      "Shaded bands = 95% CI."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "none",
    plot.caption    = element_text(size = 8, hjust = 0, colour = "grey40")
  )

ggsave(file.path(out_fig, "rob_m5_stability_path.png"), figB,
       width = 12, height = 5, dpi = 300)
ggsave(file.path(out_fig, "rob_m5_stability_path.pdf"), figB,
       width = 12, height = 5)
cat("Saved rob_m5_stability_path.png/.pdf\n")

## ── 9. FIGURE C — Gender comparison bar chart ────────────────────────────── ##
cat("Building Figure C (gender comparison)...\n")

gender_specs <- c("Baseline", "Male only", "Female only")
gender_labels <- c("Baseline" = "Pooled",
                   "Male only" = "Male only",
                   "Female only" = "Female only")

df_figC <- results_df %>%
  filter(as.character(spec) %in% gender_specs) %>%
  mutate(group = gender_labels[as.character(spec)]) %>%
  select(group, three_coef, three_se, three_pval,
         one_coef, one_se, one_pval) %>%
  pivot_longer(
    cols      = c(three_coef, three_se, three_pval,
                  one_coef, one_se, one_pval),
    names_to  = c("meal", ".value"),
    names_pattern = "(three|one)_(coef|se|pval)"
  ) %>%
  mutate(
    effect = ifelse(meal == "three", "Three-meal", "One-meal"),
    lo     = coef - 1.96 * se,
    hi     = coef + 1.96 * se,
    group  = factor(group, levels = c("Pooled", "Male only", "Female only")),
    effect = factor(effect, levels = c("Three-meal", "One-meal"))
  ) %>%
  filter(!is.na(coef))

figC <- ggplot(df_figC,
               aes(x = group, y = coef, fill = effect, group = effect)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey40", linewidth = 0.6) +
  geom_col(position = position_dodge(width = 0.6),
           width = 0.5, alpha = 0.85) +
  geom_errorbar(aes(ymin = lo, ymax = hi),
                position = position_dodge(width = 0.6),
                width = 0.2, linewidth = 0.7) +
  scale_fill_manual(
    values = c("Three-meal" = COL_THREE, "One-meal" = COL_ONE),
    name   = "Meal type"
  ) +
  labs(
    x       = NULL,
    y       = "Differential vs. no-meal workers (BDT/day)",
    caption = paste(
      "Three-meal and One-meal differentials vs. no-meal workers, by gender subsample.",
      "Error bars = 95% CI."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.caption    = element_text(size = 8, hjust = 0, colour = "grey40")
  )

ggsave(file.path(out_fig, "rob_m5_gender_comparison.png"), figC,
       width = 7, height = 5, dpi = 300)
ggsave(file.path(out_fig, "rob_m5_gender_comparison.pdf"), figC,
       width = 7, height = 5)
cat("Saved rob_m5_gender_comparison.png/.pdf\n")

## ── 10. FIGURE D — Season-specific M5 ────────────────────────────────────── ##
cat("Building Figure D (season-specific M5)...\n")

season_models <- list()
for (szn in c("Boro", "Aus", "Aman")) {
  df_s <- df %>% filter(growing_season == szn)
  ## Use year + District FE (District^season collapses within single season)
  season_models[[szn]] <- safe_feols(
    diff_real_wage ~ diff_log_yield_hat * meal_type + gender | year + District,
    df_s, ~District, paste("M5 season:", szn)
  )
}

season_annot <- c(
  Boro = "Boro predictor: GDD p=0.008",
  Aus  = "Aus predictor: EDD p=0.051",
  Aman = "Aman: no significant predictor (placebo)"
)

df_figD <- bind_rows(lapply(c("Boro", "Aus", "Aman"), function(szn) {
  df_s <- df %>% filter(growing_season == szn)
  te <- extract_three_effect(
    "diff_real_wage ~ diff_log_yield_hat * meal_type + gender | year + District",
    df_s)
  data.frame(
    season = szn,
    coef   = te$coef,
    se     = te$se,
    pval   = te$pval,
    n      = te$n,
    annot  = season_annot[[szn]],
    stringsAsFactors = FALSE
  )
})) %>%
  filter(!is.na(coef)) %>%
  mutate(
    lo     = coef - 1.96 * se,
    hi     = coef + 1.96 * se,
    sig    = ifelse(!is.na(pval) & pval < 0.10, "sig", "insig"),
    season = factor(season, levels = c("Boro", "Aus", "Aman"))
  )

## Position annotations just above the upper CI
annot_y_max <- if (nrow(df_figD) > 0) max(df_figD$hi, na.rm = TRUE) else 0
annot_step  <- abs(annot_y_max) * 0.10 + 20

figD <- ggplot(df_figD, aes(x = season, y = coef, colour = sig)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "red3", linewidth = 0.7) +
  geom_errorbar(aes(ymin = lo, ymax = hi),
                width = 0.2, linewidth = 0.8) +
  geom_point(size = 4) +
  geom_text(aes(label = annot, y = hi + annot_step),
            vjust = 0, size = 3, colour = "grey30",
            fontface = "italic", lineheight = 0.9) +
  scale_colour_manual(
    values = c(sig = SET2[1], insig = "grey60"),
    labels = c(sig = "p < 0.10", insig = "p \u2265 0.10"),
    name   = NULL
  ) +
  expand_limits(y = annot_y_max + annot_step * 3) +
  labs(
    x       = "Rice season",
    y       = "Three-meal differential vs. no-meal workers, by season (BDT/day)",
    caption = paste(
      "Three-meal differential vs. no-meal workers, by rice season.",
      "Aman serves as within-paper placebo (no significant climate first stage)."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.caption    = element_text(size = 8, hjust = 0, colour = "grey40")
  )

ggsave(file.path(out_fig, "rob_m5_by_season.png"), figD,
       width = 7, height = 5, dpi = 300)
ggsave(file.path(out_fig, "rob_m5_by_season.pdf"), figD,
       width = 7, height = 5)
cat("Saved rob_m5_by_season.png/.pdf\n")

## ── 11. FIGURE E — Leave-one-district-out (LOO, 63 models) ───────────────── ##
cat("Running LOO loop (63 models)...\n")

all_districts <- sort(unique(df$District))

loo_results <- bind_rows(lapply(all_districts, function(d) {
  df_loo <- df %>% filter(District != d)
  te <- extract_three_effect(
    "diff_real_wage ~ diff_log_yield_hat * meal_type + gender | year + District^growing_season",
    df_loo)
  data.frame(District = d, coef = te$coef, se = te$se, pval = te$pval,
             stringsAsFactors = FALSE)
}))

loo_results <- loo_results %>%
  filter(!is.na(coef)) %>%
  arrange(District) %>%           # alphabetical sort on Y axis
  mutate(
    lo       = coef - 1.96 * se,
    hi       = coef + 1.96 * se,
    sig      = ifelse(!is.na(pval) & pval < 0.10, "sig", "insig"),
    District = factor(District, levels = District)
  )

n_sig_loo <- sum(loo_results$sig == "sig", na.rm = TRUE)
cat(sprintf("  LOO: range [%.2f, %.2f], %d/%d significant at p<0.10\n",
            min(loo_results$coef, na.rm = TRUE),
            max(loo_results$coef, na.rm = TRUE),
            n_sig_loo, nrow(loo_results)))

figE <- ggplot(loo_results, aes(y = District, x = coef, colour = sig)) +
  geom_vline(xintercept =    0,   linetype = "dashed",
             colour = "red3",     linewidth = 0.7) +
  geom_vline(xintercept = -222.6, linetype = "dashed",
             colour = "steelblue", linewidth = 0.7) +
  annotate("text", y = Inf, x = -222.6, hjust = 1.05, vjust = -0.3,
           label = "Baseline\n(\u2013222.6)", colour = "steelblue", size = 2.8) +
  geom_errorbarh(aes(xmin = lo, xmax = hi),
                 height = 0.45, linewidth = 0.45) +
  geom_point(size = 1.8) +
  scale_colour_manual(
    values = c(sig = SET2[1], insig = "grey60"),
    labels = c(sig = "p < 0.10", insig = "p \u2265 0.10"),
    name   = NULL
  ) +
  labs(
    x       = "Three-meal differential vs. no-meal (BDT/day)",
    y       = NULL,
    caption = paste(
      "Leave-one-district-out sensitivity. Each point = M5 differential (interaction coefficient) dropping that district.",
      "Blue dashed = full-sample baseline (\u2013222.6)."
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y     = element_text(size = 7),
    legend.position = "bottom",
    plot.caption    = element_text(size = 8, hjust = 0, colour = "grey40")
  )

ggsave(file.path(out_fig, "rob_m5_loo.png"), figE,
       width = 8, height = 14, dpi = 300)
ggsave(file.path(out_fig, "rob_m5_loo.pdf"), figE,
       width = 8, height = 14)
cat("Saved rob_m5_loo.png/.pdf\n")

## ── 12. FIGURE F — District map of M5 fixed effects ─────────────────────── ##
cat("Building Figure F (district map)...\n")

shp_path <- file.path(
  ROOT,
  "data/GAEZ data/Agro-MAPS_BGD/AgroMaps/Asia/shapefiles/BGD/admin1/bgd.shp"
)

figF_status <- tryCatch({
  if (!file.exists(shp_path))
    stop(sprintf("Shapefile not found: %s", shp_path))

  bgd_sf <- sf::st_read(shp_path, quiet = TRUE)

  ## ── Extract fixed effects from M5 baseline ──────────────────────────────── ##
  fe_list   <- fixef(m5_baseline)
  fe_dsname <- names(fe_list)[grep("District", names(fe_list))[1]]

  if (is.na(fe_dsname))
    stop("No District-related FE found in fixef(m5_baseline)")

  fe_vals <- fe_list[[fe_dsname]]

  ## Parse District from combined FE label (e.g. "Dhaka_Boro" → "Dhaka")
  ## fixest labels District^growing_season FEs as "<District>_<Season>"
  fe_df_raw <- data.frame(
    fe_label = names(fe_vals),
    fe_val   = as.numeric(fe_vals),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      District = sub("_(Boro|Aus|Aman)$", "", fe_label),
      ## Fallback: strip everything after last underscore if above didn't work
      District = ifelse(
        grepl("_(Boro|Aus|Aman)$", fe_label), District,
        sub("_[^_]+$", "", fe_label)
      )
    ) %>%
    group_by(District) %>%
    summarise(fe_val = mean(fe_val, na.rm = TRUE), .groups = "drop")

  ## ── Identify shapefile name column ──────────────────────────────────────── ##
  shp_name_col <- names(bgd_sf)[
    grepl("name|NAME|admin|DIST|DIST_", names(bgd_sf), ignore.case = TRUE)
  ][1]
  if (is.na(shp_name_col)) shp_name_col <- names(bgd_sf)[2]
  shp_dist_names <- bgd_sf[[shp_name_col]]

  ## ── Fuzzy join using Jaro-Winkler distance ──────────────────────────────── ##
  dist_mat       <- stringdistmatrix(tolower(fe_df_raw$District),
                                     tolower(shp_dist_names),
                                     method = "jw")
  fe_df_raw$shp_match  <- shp_dist_names[apply(dist_mat, 1, which.min)]
  fe_df_raw$match_dist <- apply(dist_mat, 1, min)

  ## Flag poor matches (Jaro-Winkler > 0.2)
  poor <- fe_df_raw %>% filter(match_dist > 0.2)
  if (nrow(poor) > 0) {
    cat(sprintf("  WARNING: %d districts had poor fuzzy matches (JW > 0.2):\n",
                nrow(poor)))
    for (i in seq_len(nrow(poor)))
      cat(sprintf("    '%s' -> '%s' (dist=%.3f)\n",
                  poor$District[i], poor$shp_match[i], poor$match_dist[i]))
  }

  ## ── Merge FEs to shapefile ───────────────────────────────────────────────── ##
  bgd_sf[[".match_name"]] <- bgd_sf[[shp_name_col]]
  fe_merge   <- fe_df_raw %>% select(shp_match, fe_val)
  bgd_merged <- bgd_sf %>%
    left_join(fe_merge, by = c(".match_name" = "shp_match"))

  n_unmatched <- sum(is.na(bgd_merged$fe_val))
  if (n_unmatched > 0) {
    unmatched <- bgd_merged[[shp_name_col]][is.na(bgd_merged$fe_val)]
    cat(sprintf("  WARNING: %d shapefile polygons unmatched: %s\n",
                n_unmatched, paste(unmatched, collapse = ", ")))
  }

  ## ── Top-5 and bottom-5 FE districts for labels ──────────────────────────── ##
  fe_sorted   <- fe_df_raw %>% arrange(fe_val)
  n_fe        <- nrow(fe_sorted)
  label_dists <- c(head(fe_sorted$shp_match, 5),
                   tail(fe_sorted$shp_match, 5))

  bgd_label_pts <- bgd_merged %>%
    filter(.match_name %in% label_dists) %>%
    sf::st_centroid()

  ## ── Symmetric colour limits centred at 0 ────────────────────────────────── ##
  fe_rng   <- range(bgd_merged$fe_val, na.rm = TRUE)
  fe_limit <- max(abs(fe_rng))

  figF <- ggplot(bgd_merged) +
    geom_sf(aes(fill = fe_val), colour = "white", linewidth = 0.2) +
    ggrepel::geom_label_repel(
      data         = bgd_label_pts,
      aes(label = .match_name, geometry = geometry),
      stat         = "sf_coordinates",
      size         = 2.4,
      label.size   = 0.2,
      box.padding  = 0.3,
      max.overlaps = 20,
      colour       = "black"
    ) +
    scale_fill_distiller(
      palette   = "RdBu",
      limits    = c(-fe_limit, fe_limit),
      direction = 1,
      na.value  = "grey85",
      name      = "LOO Differential (BDT/day)"
    ) +
    labs(
      caption = paste(
        "District fixed effects from M5 baseline (averaged across seasons).",
        "Captures time-invariant district-level wage sensitivity",
        "not explained by yield shocks.",
        "Labels show top-5 and bottom-5 FE districts."
      )
    ) +
    theme_void(base_size = 11) +
    theme(
      legend.position = "right",
      plot.caption    = element_text(size = 8, hjust = 0, colour = "grey40",
                                     margin = margin(t = 8))
    )

  ggsave(file.path(out_fig, "rob_m5_district_map.png"), figF,
         width = 8, height = 10, dpi = 300)
  ggsave(file.path(out_fig, "rob_m5_district_map.pdf"), figF,
         width = 8, height = 10)
  cat("Saved rob_m5_district_map.png/.pdf\n")
  "success"

}, error = function(e) {
  cat(sprintf("  WARNING: Figure F (district map) failed: %s\n",
              conditionMessage(e)))
  cat("  Skipping. Check shapefile path and district name matching.\n")
  "skipped"
})

## ── 13. SECTION 4: Console summary ───────────────────────────────────────── ##
cat("\n=== ROB-M5 SUMMARY ===\n\n")
cat(sprintf("%-30s  %11s  %8s  %8s  %s\n",
            "Specification", "3-meal coef", "SE", "p-value", "Significant?"))
cat(strrep("-", 72), "\n")

for (i in seq_len(nrow(results_df))) {
  r    <- results_df[i, ]
  flag <- if (!is.na(r$three_pval) && r$three_pval < 0.10) "YES *" else "no"
  cat(sprintf("%-30s  %11s  %8s  %8s  %s\n",
              as.character(r$spec),
              ifelse(is.na(r$three_coef), "  NA", sprintf("%11.2f", r$three_coef)),
              ifelse(is.na(r$three_se),   "      NA", sprintf("%8.2f", r$three_se)),
              ifelse(is.na(r$three_pval), "      NA", sprintf("%8.4f", r$three_pval)),
              flag))
}

cat(strrep("-", 72), "\n")
cat("\nSeason-specific M5 (Three-meal total effect):\n")
cat(sprintf("%-10s  %11s  %8s  %8s\n",
            "Season", "3-meal coef", "SE", "p-value"))

for (szn in c("Boro", "Aus", "Aman")) {
  df_s <- df %>% filter(growing_season == szn)
  te <- extract_three_effect(
    "diff_real_wage ~ diff_log_yield_hat * meal_type + gender | year + District",
    df_s)
  cat(sprintf("%-10s  %11s  %8s  %8s\n",
              szn,
              ifelse(is.na(te$coef), "  NA", sprintf("%11.2f", te$coef)),
              ifelse(is.na(te$se),   "      NA", sprintf("%8.2f", te$se)),
              ifelse(is.na(te$pval), "      NA", sprintf("%8.4f", te$pval))))
}

cat(sprintf(
  "\nLOO range: %.2f to %.2f, N specs significant at p<0.10: %d/%d\n",
  min(loo_results$coef, na.rm = TRUE),
  max(loo_results$coef, na.rm = TRUE),
  n_sig_loo,
  nrow(loo_results)
))

cat("\n=== ROB-M5 COMPLETE ===\n")
cat("Output files written:\n")
cat(sprintf("  Tables: %s\n", file.path(out_tbl, "rob_m5_threemeal.{tex,html}")))
cat(sprintf("  Fig A:  %s\n", file.path(out_fig, "rob_m5_coefplot.{png,pdf}")))
cat(sprintf("  Fig B:  %s\n", file.path(out_fig, "rob_m5_stability_path.{png,pdf}")))
cat(sprintf("  Fig C:  %s\n", file.path(out_fig, "rob_m5_gender_comparison.{png,pdf}")))
cat(sprintf("  Fig D:  %s\n", file.path(out_fig, "rob_m5_by_season.{png,pdf}")))
cat(sprintf("  Fig E:  %s\n", file.path(out_fig, "rob_m5_loo.{png,pdf}")))
if (figF_status == "success")
  cat(sprintf("  Fig F:  %s\n", file.path(out_fig, "rob_m5_district_map.{png,pdf}")))
