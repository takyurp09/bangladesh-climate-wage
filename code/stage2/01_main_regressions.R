## stage2_main.R
## Second-stage: yield pass-through to agricultural wages (LEVELS spec)
## Spec: log_real_wage ~ log_yield_hat | year + District^growing_season
## Output: output/stage2/tables/, output/stage2/models/

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(fixest)
  library(modelsummary)
  library(ggplot2)
})

ROOT     <- here::here()
out_tbl  <- file.path(ROOT, "output/stage2/tables")
out_mdl  <- file.path(ROOT, "output/stage2/models")
out_fig  <- file.path(ROOT, "output/stage2/figures")
out_sum  <- file.path(ROOT, "output/stage2/summary")
dir.create(out_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(out_mdl, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_sum, recursive = TRUE, showWarnings = FALSE)

YHAT <- "log_yield_hat"
cat("=== STAGE 2 MAIN REGRESSION (LEVELS) ===\n")

## ── Load data ─────────────────────────────────────────────────────────────── ##
df <- read_csv(file.path(ROOT, "data/Regression_data/df_2_merged_v2.csv"),
               show_col_types = FALSE) %>%
  mutate(
    gender    = relevel(factor(gender),    ref = "Female"),
    meal_type = relevel(factor(meal_type), ref = "None")
  )

cat(sprintf("N = %d | districts = %d | years = %s\n",
            nrow(df), length(unique(df$District)),
            paste(sort(unique(df$year)), collapse = ",")))

## ── Models M1–M4 ──────────────────────────────────────────────────────────── ##
m1 <- feols(as.formula(paste("log_real_wage ~", YHAT, "| year + District^growing_season")),
            data = df, cluster = ~District)

m2 <- feols(as.formula(paste("log_real_wage ~", YHAT, "+ gender | year + District^growing_season")),
            data = df, cluster = ~District)

m3 <- feols(as.formula(paste("log_real_wage ~", YHAT, "+ gender + meal_type | year + District^growing_season")),
            data = df, cluster = ~District)

m4 <- m3  # levels spec: M3 and M4 both log wages

## ── M3_extended: + controls (irrigation + land + pop) ─────────────────────── ##
m3_ext <- feols(as.formula(paste(
  "log_real_wage ~", YHAT, "+ gender + meal_type +",
  "share_Boro + share_Aus + pop_density + log_irrigation_total + cropping_intensity |",
  "year + District^growing_season"
)), data = df, cluster = ~District)

## ── M5: yield × meal_type interaction ─────────────────────────────────────── ##
m5 <- feols(as.formula(paste(
  "log_real_wage ~", YHAT, "+", YHAT, ":meal_type + gender |",
  "year + District^growing_season"
)), data = df, cluster = ~District)

cat("\n=== M3 RESULTS ===\n")
cat(sprintf("coef=%.4f SE=%.4f p=%.4f N=%d\n",
            coef(m3)[YHAT], se(m3)[YHAT],
            pvalue(m3)[YHAT], nobs(m3)))

cat("\n=== M3_extended RESULTS ===\n")
cat(sprintf("coef=%.4f SE=%.4f p=%.4f N=%d\n",
            coef(m3_ext)[YHAT], se(m3_ext)[YHAT],
            pvalue(m3_ext)[YHAT], nobs(m3_ext)))

cat("\n=== M5 MEAL INTERACTION RESULTS ===\n")
print(summary(m5)$coeftable)

## ── Wild cluster bootstrap ────────────────────────────────────────────────── ##
run_wild_boot <- function(model_full, model_restricted, df_b, B = 9999, label = "") {
  cat(sprintf("Running wild bootstrap on %s (B=%d)...\n", label, B))
  tryCatch({
    set.seed(42)
    clusters <- unique(df_b$District)
    G        <- length(clusters)
    resids_r <- residuals(model_restricted)
    t_actual <- coef(model_full)[YHAT] / se(model_full)[YHAT]
    t_boot <- numeric(B)
    for (b in seq_len(B)) {
      w_map <- setNames(sample(c(-1, 1), G, replace = TRUE), clusters)
      w_vec <- w_map[df_b$District]
      df_b$y_star <- fitted(model_restricted) + resids_r * w_vec
      m_b <- feols(as.formula(paste(
        "y_star ~", YHAT, "+ gender + meal_type | year + District^growing_season"
      )), data = df_b, cluster = ~District, warn = FALSE, notes = FALSE)
      t_boot[b] <- coef(m_b)[YHAT] / se(m_b)[YHAT]
    }
    p <- mean(abs(t_boot) >= abs(t_actual))
    cat(sprintf("  Wild bootstrap p (%s): %.4f\n", label, p))
    p
  }, error = function(e) {
    cat("  Wild bootstrap failed:", conditionMessage(e), "\n")
    NA_real_
  })
}

ctrl_vars <- c(YHAT, "gender", "meal_type",
               "share_Boro", "share_Aus", "pop_density",
               "log_irrigation_total", "cropping_intensity")

df_m3    <- df[complete.cases(df[, c("log_real_wage", YHAT, "gender", "meal_type")]), ]
df_m3ext <- df[complete.cases(df[, c("log_real_wage", ctrl_vars)]), ]

m3_r <- feols(as.formula(paste(
  "log_real_wage ~ gender + meal_type | year + District^growing_season"
)), data = df_m3, cluster = ~District)

m3ext_r <- feols(as.formula(paste(
  "log_real_wage ~ gender + meal_type + share_Boro + share_Aus +",
  "pop_density + log_irrigation_total + cropping_intensity |",
  "year + District^growing_season"
)), data = df_m3ext, cluster = ~District)

wild_p_m3    <- run_wild_boot(m3,     m3_r,    df_m3,    label = "M3")
wild_p_m3ext <- run_wild_boot(m3_ext, m3ext_r, df_m3ext, label = "M3_extended")

## ── Helper: save tex + html ───────────────────────────────────────────────── ##
save_table <- function(models_list, filename_stem, title, add_rows = NULL) {
  int_pat <- paste0(YHAT, ":meal_type")
  common_args <- list(
    models_list,
    title        = title,
    stars        = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
    gof_omit     = "AIC|BIC|Log|Std|RMSE|F$",
    coef_rename  = c(
      "log_yield_hat"                         = "log(Yield Hat)",
      "genderMale"                            = "Male",
      "meal_typeOne"                          = "Meal: One",
      "meal_typeTwo"                          = "Meal: Two",
      "meal_typeThree"                        = "Meal: Three",
      "share_Boro"                            = "Share Boro",
      "share_Aus"                             = "Share Aus",
      "pop_density"                           = "Pop. Density",
      "log_irrigation_total"                  = "log(Irrigation)",
      "cropping_intensity"                    = "Cropping Intensity",
      "log_yield_hat:meal_typeOne"            = "Yield x Meal: One",
      "log_yield_hat:meal_typeTwo"            = "Yield x Meal: Two",
      "log_yield_hat:meal_typeThree"          = "Yield x Meal: Three"
    ),
    notes        = "* p<0.1 ** p<0.05 *** p<0.01. Levels spec. SE clustered by district."
  )
  if (!is.null(add_rows)) common_args$add_rows <- add_rows
  do.call(modelsummary, c(common_args, list(
    output = file.path(out_tbl, paste0(filename_stem, ".tex"))
  )))
  do.call(modelsummary, c(common_args, list(
    output = file.path(out_tbl, paste0(filename_stem, ".html"))
  )))
  cat(sprintf("Saved %s.tex/.html\n", filename_stem))
}

## ── Table: main_table_v2 ──────────────────────────────────────────────────── ##
fe_rows_v2 <- data.frame(
  term         = c("Year FE", "District x Season FE", "Specification",
                   "Wild Bootstrap p (M3)", "Wild Bootstrap p (M3 extended)"),
  M1           = c("Yes", "Yes", "Levels", "", ""),
  M2           = c("Yes", "Yes", "Levels", "", ""),
  M3           = c("Yes", "Yes", "Levels",
                   ifelse(is.na(wild_p_m3),    "failed", sprintf("%.3f", wild_p_m3)), ""),
  M3_extended  = c("Yes", "Yes", "Levels", "",
                   ifelse(is.na(wild_p_m3ext), "failed", sprintf("%.3f", wild_p_m3ext))),
  M4           = c("Yes", "Yes", "Levels", "", ""),
  stringsAsFactors = FALSE
)

save_table(
  models_list   = list(M1 = m1, M2 = m2, M3 = m3, "M3 Extended" = m3_ext, M4 = m4),
  filename_stem = "main_table_v2",
  title         = "Yield Pass-Through to Agricultural Wages (Levels, Second Stage)",
  add_rows      = fe_rows_v2
)

## ── Table: meal_interaction_table (M5) ───────────────────────────────────── ##
fe_rows_m5 <- data.frame(
  term = c("Year FE", "District x Season FE", "Specification"),
  M5   = c("Yes", "Yes", "Levels"),
  stringsAsFactors = FALSE
)

save_table(
  models_list   = list(M5 = m5),
  filename_stem = "meal_interaction_table",
  title         = "Yield Pass-Through by Meal Type (Levels, M5)",
  add_rows      = fe_rows_m5
)

## ── Figures ───────────────────────────────────────────────────────────────── ##
oi_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")

extract_coef <- function(mod, name) {
  ct <- summary(mod)$coeftable
  if (!YHAT %in% rownames(ct)) return(NULL)
  data.frame(model = name, coef = ct[YHAT, "Estimate"],
             se = ct[YHAT, "Std. Error"], stringsAsFactors = FALSE)
}

coef_df <- bind_rows(
  extract_coef(m1,     "M1"),
  extract_coef(m2,     "M2"),
  extract_coef(m3,     "M3"),
  extract_coef(m3_ext, "M3 Extended"),
  extract_coef(m5,     "M5 (yield×meal)")
) %>% mutate(
  model = factor(model, levels = rev(c("M1","M2","M3","M3 Extended","M5 (yield×meal)"))),
  lo90  = coef - 1.645 * se, hi90  = coef + 1.645 * se,
  lo95  = coef - 1.96  * se, hi95  = coef + 1.96  * se
)

p_coef <- ggplot(coef_df, aes(x = coef, y = model, color = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_linerange(aes(xmin = lo95, xmax = hi95), linewidth = 0.8) +
  geom_linerange(aes(xmin = lo90, xmax = hi90), linewidth = 1.6) +
  geom_point(size = 3) +
  scale_color_manual(values = oi_cols, guide = "none") +
  labs(
    title    = "Yield Pass-Through Coefficient (log Yield Hat)",
    subtitle = "Levels spec. Thick = 90% CI, thin = 95% CI. Clustered SE by district.",
    x        = "Coefficient estimate",
    y        = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_fig, "coefplot_main.png"), p_coef, width = 8, height = 5, dpi = 300)
cat("Saved coefplot_main.png\n")

m5_ct <- summary(m5)$coeftable
int_terms <- rownames(m5_ct)[grepl(paste0(YHAT, ":meal_type"), rownames(m5_ct))]

meal_df <- bind_rows(
  data.frame(term = "None (baseline)", coef = m5_ct[YHAT, "Estimate"],
             se = m5_ct[YHAT, "Std. Error"], stringsAsFactors = FALSE),
  data.frame(
    term = sub(paste0(YHAT, ":meal_type"), "", int_terms),
    coef = m5_ct[int_terms, "Estimate"] + m5_ct[YHAT, "Estimate"],
    se   = sqrt(m5_ct[int_terms, "Std. Error"]^2 + m5_ct[YHAT, "Std. Error"]^2),
    stringsAsFactors = FALSE
  )
) %>% mutate(
  term = factor(term, levels = c("None (baseline)", "One", "Two", "Three")),
  lo90 = coef - 1.645 * se, hi90 = coef + 1.645 * se,
  lo95 = coef - 1.96  * se, hi95 = coef + 1.96  * se
)

p_meal <- ggplot(meal_df, aes(x = coef, y = term, color = term)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_linerange(aes(xmin = lo95, xmax = hi95), linewidth = 0.8) +
  geom_linerange(aes(xmin = lo90, xmax = hi90), linewidth = 1.6) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("None (baseline)" = "#E69F00", "One" = "#56B4E9",
               "Two" = "#009E73", "Three" = "#D55E00"),
    guide = "none"
  ) +
  labs(
    title    = "Pass-Through by Meal Type (M5, Levels)",
    subtitle = "Total effect = baseline + interaction.",
    x        = "Total coefficient on log(Yield Hat)",
    y        = "Meal type"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_fig, "meal_interaction_coefplot.png"), p_meal, width = 7, height = 4, dpi = 300)
cat("Saved meal_interaction_coefplot.png\n")

## ── Save models + summary ─────────────────────────────────────────────────── ##
save(m1, m2, m3, m3_ext, m4, m5, wild_p_m3, wild_p_m3ext,
     file = file.path(out_mdl, "stage2_main_models_v2.RData"))

m5_ct <- summary(m5)$coeftable
int_rows <- rownames(m5_ct)[grepl(paste0(YHAT, ":meal_type"), rownames(m5_ct))]

summary_lines <- c(
  "# Stage 2 Final Summary (Levels Spec)",
  sprintf("Generated: %s", Sys.time()),
  "",
  "## Specification",
  "- Stage 1: levels (log_yield ~ GDD/EDD | year + district), season-specific",
  "- Stage 2: log_real_wage ~ log_yield_hat | year + District^growing_season",
  sprintf("- Data: df_2_merged_levels.csv (N=%d)", nrow(df)),
  sprintf("- M3_extended obs: %d", nobs(m3_ext)),
  "",
  "## M3 vs M3_extended",
  sprintf("| Model | Coef | SE | p-value | Wild-boot p |"),
  sprintf("|-------|------|----|---------|--------------|"),
  sprintf("| M3 | %.4f | %.4f | %.4f | %s |",
          coef(m3)[YHAT], se(m3)[YHAT], pvalue(m3)[YHAT],
          ifelse(is.na(wild_p_m3), "failed", sprintf("%.3f", wild_p_m3))),
  sprintf("| M3_extended | %.4f | %.4f | %.4f | %s |",
          coef(m3_ext)[YHAT], se(m3_ext)[YHAT], pvalue(m3_ext)[YHAT],
          ifelse(is.na(wild_p_m3ext), "failed", sprintf("%.3f", wild_p_m3ext))),
  "",
  "## M5: Meal-type interaction (Three-meal differential)",
  if (length(int_rows) > 0) {
    r3 <- grep("Three", int_rows, value = TRUE)[1]
    sprintf("- Three-meal interaction: coef=%.4f, p=%.4f",
            m5_ct[r3, "Estimate"], m5_ct[r3, "Pr(>|t|)"])
  } else "- No interaction terms found.",
  sprintf("- Baseline (None): coef=%.4f, p=%.4f",
          m5_ct[YHAT, "Estimate"], m5_ct[YHAT, "Pr(>|t|)"]),
  "",
  "## Outputs",
  "- output/stage2/tables/main_table_v2.tex/.html",
  "- output/stage2/tables/meal_interaction_table.tex/.html",
  "- output/stage2/figures/coefplot_main.png",
  "- output/stage2/figures/meal_interaction_coefplot.png"
)

writeLines(summary_lines, file.path(out_sum, "STAGE2_FINAL_SUMMARY_v2.md"))
cat("Saved STAGE2_FINAL_SUMMARY_v2.md\n")
cat("=== STAGE 2 MAIN (LEVELS) COMPLETE ===\n")
