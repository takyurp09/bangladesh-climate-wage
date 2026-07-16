## stage2_by_season.R
## Second-stage regressions run separately by growing season
## Instrument strength: Boro (strong, p=0.008), Aus (weak-moderate, p=0.051), Aman (none → placebo)
## Output: output/stage2/tables/, output/stage2/figures/, output/stage2/summary/

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(fixest)
  library(modelsummary)
  library(ggplot2)
  library(patchwork)
})

ROOT    <- here::here()
out_tbl <- file.path(ROOT, "output/stage2/tables")
out_fig <- file.path(ROOT, "output/stage2/figures")
out_sum <- file.path(ROOT, "output/stage2/summary")
dir.create(out_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_sum, recursive = TRUE, showWarnings = FALSE)

cat("=== STAGE 2 BY SEASON ===\n")

## ── Load wage panel (already merged with levels yield_hat) ───────────────── ##
wage_raw <- read_csv(file.path(ROOT, "data/Regression_data/df_2_merged_v2.csv"),
                     show_col_types = FALSE) %>%
  filter(!is.na(log_real_wage), !is.na(log_yield_hat)) %>%
  mutate(
    gender    = relevel(factor(gender),    ref = "Female"),
    meal_type = relevel(factor(meal_type), ref = "None")
  )

cat(sprintf("Wage panel rows (levels): %d\n", nrow(wage_raw)))

## Stage 1 instrument strength labels
s1_strength <- c(
  Boro = "Boro (GDD p=0.008 **)",
  Aus  = "Aus (EDD p=0.051 *)",
  Aman = "Aman (no sig. instrument)"
)

SEASONS <- c("Boro", "Aus", "Aman")

## ── Build season-specific datasets ───────────────────────────────────────── ##
season_data <- lapply(SEASONS, function(s) {
  df_s <- wage_raw %>% filter(growing_season == s)

  cat(sprintf("%s: N=%d | districts=%d | years=%s\n",
              s, nrow(df_s), length(unique(df_s$District)),
              paste(sort(unique(df_s$year)), collapse = ",")))
  df_s
})
names(season_data) <- SEASONS

## ── Fit models ────────────────────────────────────────────────────────────── ##
fit_season <- function(df_s, season) {
  fe_spec <- "year + District"

  m3 <- feols(log_real_wage ~ log_yield_hat + gender + meal_type |
                year + District,
              data = df_s, cluster = ~District)

  m5 <- tryCatch(
    feols(log_real_wage ~ log_yield_hat +
            log_yield_hat:meal_type + gender |
            year + District,
          data = df_s, cluster = ~District),
    error = function(e) {
      cat(sprintf("  M5 failed for %s: %s\n", season, conditionMessage(e)))
      NULL
    }
  )

  list(m3 = m3, m5 = m5)
}

models <- lapply(SEASONS, function(s) fit_season(season_data[[s]], s))
names(models) <- SEASONS

## ── Print key results ─────────────────────────────────────────────────────── ##
cat("\n=== MAIN COEFFICIENTS (M3 per season) ===\n")
for (s in SEASONS) {
  m <- models[[s]]$m3
  cat(sprintf("%-5s coef=%9.4f SE=%8.4f p=%.4f N=%d\n",
              s,
              coef(m)["log_yield_hat"],
              se(m)["log_yield_hat"],
              pvalue(m)["log_yield_hat"],
              nobs(m)))
}

## ── Helper ────────────────────────────────────────────────────────────────── ##
coef_rename_map <- c(
  "log_yield_hat"                       = "log(YieldHat)",
  "genderMale"                          = "Male",
  "meal_typeOne"                        = "Meal: One",
  "meal_typeTwo"                        = "Meal: Two",
  "meal_typeThree"                      = "Meal: Three",
  "log_yield_hat:meal_typeOne"       = "YieldHat x Meal: One",
  "log_yield_hat:meal_typeTwo"       = "YieldHat x Meal: Two",
  "log_yield_hat:meal_typeThree"     = "YieldHat x Meal: Three"
)

save_table <- function(models_list, stem, title, add_rows = NULL) {
  args <- list(
    models_list,
    title       = title,
    stars       = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
    gof_omit    = "AIC|BIC|Log|Std|RMSE|F$",
    coef_rename = coef_rename_map,
    notes       = "* p<0.1 ** p<0.05 *** p<0.01. SE clustered by district. FE: year + District."
  )
  if (!is.null(add_rows)) args$add_rows <- add_rows
  do.call(modelsummary, c(args, list(output = file.path(out_tbl, paste0(stem, ".tex")))))
  do.call(modelsummary, c(args, list(output = file.path(out_tbl, paste0(stem, ".html")))))
  cat(sprintf("Saved %s.tex/.html\n", stem))
}

## ── Table: season_main (M3 per season) ───────────────────────────────────── ##
fe_rows_main <- data.frame(
  term       = c("Year FE", "District FE", "Stage 1 instrument"),
  M3_Boro    = c("Yes", "Yes", "GDD (p=0.008 **)"),
  M3_Aus     = c("Yes", "Yes", "EDD (p=0.051 *)"),
  M3_Aman    = c("Yes", "Yes", "None (placebo)"),
  stringsAsFactors = FALSE
)
save_table(
  models_list = list("M3 Boro" = models$Boro$m3,
                     "M3 Aus"  = models$Aus$m3,
                     "M3 Aman" = models$Aman$m3),
  stem        = "season_main",
  title       = "Season-Specific Yield Pass-Through (M3 per Season)",
  add_rows    = fe_rows_main
)

## ── Table: season_meal_interaction (M5 per season) ───────────────────────── ##
m5_list <- Filter(Negate(is.null), lapply(SEASONS, function(s) models[[s]]$m5))
names(m5_list) <- paste0("M5 ", SEASONS[!sapply(SEASONS, function(s) is.null(models[[s]]$m5))])

if (length(m5_list) > 0) {
  fe_rows_m5 <- data.frame(
    term = c("Year FE", "District FE"),
    stringsAsFactors = FALSE
  )
  for (nm in names(m5_list)) fe_rows_m5[[nm]] <- c("Yes", "Yes")
  save_table(
    models_list = m5_list,
    stem        = "season_meal_interaction",
    title       = "Season-Specific Yield × Meal-Type Interaction (M5 per Season)",
    add_rows    = fe_rows_m5
  )
}

## ── Figure: season_coefplot.png ───────────────────────────────────────────── ##
oi <- c(Boro = "#0072B2", Aus = "#E69F00", Aman = "#CC79A7")
season_labels <- c(
  Boro = "Boro\n(GDD p=0.008**)",
  Aus  = "Aus\n(EDD p=0.051*)",
  Aman = "Aman\n(no first-stage predictor)"
)

## Row 1: M3 main coef per season
m3_df <- bind_rows(lapply(SEASONS, function(s) {
  m <- models[[s]]$m3
  data.frame(
    season = s,
    coef   = coef(m)["log_yield_hat"],
    se     = se(m)["log_yield_hat"],
    p      = pvalue(m)["log_yield_hat"],
    N      = nobs(m),
    stringsAsFactors = FALSE
  )
})) %>%
  mutate(
    season = factor(season, levels = SEASONS),
    lo90 = coef - 1.645 * se, hi90 = coef + 1.645 * se,
    lo95 = coef - 1.96  * se, hi95 = coef + 1.96  * se,
    sig  = ifelse(p < 0.05, "p<0.05", ifelse(p < 0.1, "p<0.1", "n.s."))
  )

p_m3 <- ggplot(m3_df, aes(x = coef, y = season, color = season)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_linerange(aes(xmin = lo95, xmax = hi95), linewidth = 0.8) +
  geom_linerange(aes(xmin = lo90, xmax = hi90), linewidth = 1.8) +
  geom_point(size = 4) +
  geom_text(aes(label = sprintf("p=%.3f\nN=%d", p, N)),
            hjust = -0.15, vjust = 0.5, size = 3, color = "grey30") +
  scale_y_discrete(labels = season_labels) +
  scale_color_manual(values = oi, guide = "none") +
  labs(title = "Row 1: M3 — Main Pass-Through by Season",
       subtitle = "Thick = 90% CI, Thin = 95% CI. FE: year + District.",
       x = "Coefficient on ΔYieldHat", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

## Row 2: M5 interaction terms per season
int_df_list <- lapply(SEASONS, function(s) {
  m <- models[[s]]$m5
  if (is.null(m)) return(NULL)
  ct        <- summary(m)$coeftable
  base_row  <- "log_yield_hat"
  int_terms <- rownames(ct)[grepl("log_yield_hat:meal_type", rownames(ct))]
  if (length(int_terms) == 0) return(NULL)
  bind_rows(
    data.frame(season = s, meal = "None (baseline)",
               coef = ct[base_row, "Estimate"],
               se   = ct[base_row, "Std. Error"],
               p    = ct[base_row, "Pr(>|t|)"],
               stringsAsFactors = FALSE),
    data.frame(
      season = s,
      meal   = sub("log_yield_hat:meal_type", "", int_terms),
      coef   = ct[int_terms, "Estimate"],
      se     = ct[int_terms, "Std. Error"],
      p      = ct[int_terms, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
  )
})
int_df <- bind_rows(int_df_list) %>%
  mutate(
    season = factor(season, levels = SEASONS),
    meal   = factor(meal, levels = c("None (baseline)", "One", "Two", "Three")),
    lo90 = coef - 1.645 * se, hi90 = coef + 1.645 * se,
    lo95 = coef - 1.96  * se, hi95 = coef + 1.96  * se
  )

meal_colors <- c("None (baseline)" = "#999999", "One" = "#56B4E9",
                 "Two" = "#009E73", "Three" = "#D55E00")

p_m5 <- ggplot(int_df, aes(x = coef, y = meal, color = meal)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_linerange(aes(xmin = lo95, xmax = hi95), linewidth = 0.7) +
  geom_linerange(aes(xmin = lo90, xmax = hi90), linewidth = 1.5) +
  geom_point(size = 3) +
  facet_wrap(~ season, nrow = 1,
             labeller = labeller(season = season_labels)) +
  scale_color_manual(values = meal_colors, guide = "none") +
  labs(title = "Row 2: M5 — Three-meal Differential vs. No-meal Workers, by Season",
       subtitle = "Differential = interaction coefficient only. Thick = 90% CI, Thin = 95% CI.",
       x = "Three-meal differential vs. no-meal workers, by season (BDT/day)", y = "Meal type") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(size = 9))

p_combined <- p_m3 / p_m5 +
  plot_layout(heights = c(1, 1.5)) +
  plot_annotation(
    title    = "Season-Specific Yield Pass-Through to Agricultural Wages",
    subtitle = "Stage 2 | FE: year + District | Clustered SE by district",
    theme    = theme(plot.title    = element_text(size = 13, face = "bold"),
                     plot.subtitle = element_text(size = 10, color = "grey40"))
  )

ggsave(file.path(out_fig, "season_coefplot.png"), p_combined,
       width = 10, height = 8, dpi = 300)
cat("Saved season_coefplot.png\n")

## ── Summary report ────────────────────────────────────────────────────────── ##
# Collect M5 interaction results per season
m5_summary <- lapply(SEASONS, function(s) {
  m <- models[[s]]$m5
  if (is.null(m)) return(list(sig_meals = "M5 failed", strongest = "n/a"))
  ct        <- summary(m)$coeftable
  base_row  <- "log_yield_hat"
  int_terms <- rownames(ct)[grepl("log_yield_hat:meal_type", rownames(ct))]
  sig_meals <- if (length(int_terms) == 0) "none" else {
    sig <- int_terms[ct[int_terms, "Pr(>|t|)"] < 0.1]
    if (length(sig) == 0) "none at p<0.1" else
      paste(sub("log_yield_hat:meal_type", "", sig), collapse = ", ")
  }
  strongest <- if (length(int_terms) == 0) "n/a" else {
    differentials <- ct[int_terms, "Estimate"]
    sub("log_yield_hat:meal_type", "", int_terms[which.max(abs(differentials))])
  }
  list(sig_meals = sig_meals, strongest = strongest)
})
names(m5_summary) <- SEASONS

red_flags <- character(0)
for (s in SEASONS) {
  m <- models[[s]]$m3
  if (s == "Aman" && pvalue(m)["log_yield_hat"] < 0.1)
    red_flags <- c(red_flags, sprintf("Aman M3 significant (p=%.3f) — placebo concern", pvalue(m)["log_yield_hat"]))
  if (nobs(m) < 500)
    red_flags <- c(red_flags, sprintf("%s N=%d — small sample", s, nobs(m)))
}
if (length(red_flags) == 0) red_flags <- "None"

report <- c(
  "# Season-Specific Stage 2 Summary",
  sprintf("Generated: %s", Sys.time()),
  "",
  "## Main pass-through coefficient (M3) per season",
  "",
  "| Season | Coef | SE | p-value | N | Stage 1 instrument |",
  "|--------|------|----|---------|---|--------------------|",
  paste0(sapply(SEASONS, function(s) {
    m <- models[[s]]$m3
    sprintf("| %s | %.4f | %.4f | %.4f | %d | %s |",
            s,
            coef(m)["log_yield_hat"],
            se(m)["log_yield_hat"],
            pvalue(m)["log_yield_hat"],
            nobs(m),
            s1_strength[s])
  }), collapse = "\n"),
  "",
  "## Does Boro (strong instrument) show cleaner pass-through?",
  sprintf("- Boro: coef=%.4f, p=%.4f — %s",
          coef(models$Boro$m3)["log_yield_hat"],
          pvalue(models$Boro$m3)["log_yield_hat"],
          if (pvalue(models$Boro$m3)["log_yield_hat"] < 0.1) "YES, significant" else "No, still insignificant"),
  sprintf("- Aus:  coef=%.4f, p=%.4f",
          coef(models$Aus$m3)["log_yield_hat"],
          pvalue(models$Aus$m3)["log_yield_hat"]),
  sprintf("- Aman: coef=%.4f, p=%.4f — %s",
          coef(models$Aman$m3)["log_yield_hat"],
          pvalue(models$Aman$m3)["log_yield_hat"],
          if (pvalue(models$Aman$m3)["log_yield_hat"] >= 0.1) "consistent with placebo (n.s.)" else "SIGNIFICANT — placebo concern"),
  "",
  "## M5 meal-type interactions per season",
  "",
  "| Season | Significant meal types (p<0.1) | Strongest differential |",
  "|--------|-------------------------------|------------------------------|",
  paste0(sapply(SEASONS, function(s)
    sprintf("| %s | %s | %s |", s, m5_summary[[s]]$sig_meals, m5_summary[[s]]$strongest)
  ), collapse = "\n"),
  "",
  "## Red flags",
  paste("-", red_flags, collapse = "\n"),
  "",
  "## Outputs",
  "- output/stage2/tables/season_main.tex/.html",
  "- output/stage2/tables/season_meal_interaction.tex/.html",
  "- output/stage2/figures/season_coefplot.png"
)

writeLines(report, file.path(out_sum, "SEASON_SPECIFIC_SUMMARY.md"))
cat("Saved SEASON_SPECIFIC_SUMMARY.md\n")
cat("=== STAGE 2 BY SEASON COMPLETE ===\n")
