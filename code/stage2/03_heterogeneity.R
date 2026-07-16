## stage2_extend.R
## Heterogeneity analysis: structural vars × yield pass-through
## Output: output/stage2/tables/heterogeneity_*.tex/.html

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(fixest)
  library(modelsummary)
})

ROOT    <- here::here()
out_tbl <- file.path(ROOT, "output/stage2/tables")
dir.create(out_tbl, recursive = TRUE, showWarnings = FALSE)

cat("=== STAGE 2 HETEROGENEITY ===\n")

## ── Load data ─────────────────────────────────────────────────────────────── ##
df <- read_csv(file.path(ROOT, "data/Regression_data/df_2_merged_struct.csv"),
               show_col_types = FALSE) %>%
  mutate(
    gender    = relevel(factor(gender),    ref = "Female"),
    meal_type = relevel(factor(meal_type), ref = "None")
  )

cat(sprintf("N = %d | struct N (non-NA irrigation) = %d\n",
            nrow(df), sum(!is.na(df$z_irrigation_2019))))

# Structural variables (2019 levels, z-scored)
struct_vars <- c(
  irrigation = "z_irrigation_2019",
  intensity  = "z_intensity_2019",
  holdings   = "z_holdings_2019",
  gca        = "z_gca_2019"
)
struct_labels <- c(
  z_irrigation_2019 = "Irrigation (z)",
  z_intensity_2019  = "Intensity (z)",
  z_holdings_2019   = "Holdings (z)",
  z_gca_2019        = "GCA (z)"
)

## ── Helper ────────────────────────────────────────────────────────────────── ##
save_table_pair <- function(models, stem, title, coef_rename = NULL, ...) {
  args <- list(
    models,
    title    = title,
    stars    = c("*" = 0.1, "**" = 0.05, "***" = 0.01),
    gof_omit = "AIC|BIC|Log|Std|RMSE|F$",
    notes    = "* p<0.1 ** p<0.05 *** p<0.01. SE clustered by district.",
    ...
  )
  if (!is.null(coef_rename)) args$coef_rename <- coef_rename
  do.call(modelsummary, c(args, list(output = file.path(out_tbl, paste0(stem, ".tex")))))
  do.call(modelsummary, c(args, list(output = file.path(out_tbl, paste0(stem, ".html")))))
  cat(sprintf("Saved %s.tex/.html\n", stem))
}

## ── 1. Interaction models ─────────────────────────────────────────────────── ##
int_models <- list()
for (v in names(struct_vars)) {
  zvar <- struct_vars[[v]]
  fml  <- as.formula(sprintf(
    "log_real_wage ~ log_yield_hat * %s + gender + meal_type | year + District^growing_season",
    zvar))
  int_models[[v]] <- feols(fml, data = df, cluster = ~District)
  cat(sprintf("Interaction model (%s): N=%d\n",
              v, nobs(int_models[[v]])))
}

cr <- c(
  "log_yield_hat"                = "Log Yield Hat",
  "log_yield_hat:z_irrigation_2019" = "Yield x Irrigation",
  "log_yield_hat:z_intensity_2019"  = "Yield x Intensity",
  "log_yield_hat:z_holdings_2019"   = "Yield x Holdings",
  "log_yield_hat:z_gca_2019"        = "Yield x GCA",
  "z_irrigation_2019"   = "Irrigation (z)",
  "z_intensity_2019"    = "Intensity (z)",
  "z_holdings_2019"     = "Holdings (z)",
  "z_gca_2019"          = "GCA (z)",
  "genderMale"          = "Male",
  "meal_typeOne"        = "Meal: One",
  "meal_typeTwo"        = "Meal: Two",
  "meal_typeThree"      = "Meal: Three"
)

save_table_pair(
  int_models,
  stem         = "heterogeneity_interaction",
  title        = "Heterogeneity: Yield Pass-Through x Structural Variables",
  coef_rename  = cr
)

## ── 2. Gender triple interaction ──────────────────────────────────────────── ##
triple_models <- list()
for (v in names(struct_vars)) {
  zvar <- struct_vars[[v]]
  fml  <- as.formula(sprintf(
    "log_real_wage ~ log_yield_hat * %s * gender + meal_type | year + District^growing_season",
    zvar))
  triple_models[[v]] <- feols(fml, data = df, cluster = ~District)
  cat(sprintf("Triple interaction model (%s): N=%d\n",
              v, nobs(triple_models[[v]])))
}

cr_triple <- c(cr,
  "log_yield_hat:z_irrigation_2019:genderMale" = "Yield x Irrig x Male",
  "log_yield_hat:z_intensity_2019:genderMale"  = "Yield x Intens x Male",
  "log_yield_hat:z_holdings_2019:genderMale"   = "Yield x Holdings x Male",
  "log_yield_hat:z_gca_2019:genderMale"        = "Yield x GCA x Male"
)

save_table_pair(
  triple_models,
  stem        = "heterogeneity_gender_triple",
  title       = "Heterogeneity: Yield Pass-Through x Structural x Gender",
  coef_rename = cr_triple
)

## ── 3. Quartile pass-through ──────────────────────────────────────────────── ##
quartile_results <- list()

for (v in names(struct_vars)) {
  zvar <- struct_vars[[v]]
  # quartile variable for each district (constant across rows per district)
  df_q <- df %>%
    group_by(District) %>%
    mutate(q_var = first(!!sym(zvar))) %>%
    ungroup() %>%
    filter(!is.na(q_var)) %>%
    mutate(quartile = ntile(q_var, 4))

  q_models <- list()
  for (q in 1:4) {
    df_sub <- df_q %>% filter(quartile == q)
    if (nrow(df_sub) < 30) {
      q_models[[paste0("Q", q)]] <- NULL
      next
    }
    q_models[[paste0("Q", q)]] <- tryCatch(
      feols(log_real_wage ~ log_yield_hat + gender + meal_type |
              year + District^growing_season,
            data = df_sub, cluster = ~District),
      error = function(e) NULL
    )
  }
  q_models <- Filter(Negate(is.null), q_models)
  quartile_results[[v]] <- q_models
  cat(sprintf("Quartile models (%s): %d quartiles fit\n", v, length(q_models)))
}

# Flatten into one named list for table
all_q <- unlist(quartile_results, recursive = FALSE)
names(all_q) <- gsub("\\.", " ", names(all_q))

if (length(all_q) > 0) {
  save_table_pair(
    all_q,
    stem        = "heterogeneity_quartile",
    title       = "Quartile Pass-Through by Structural Variable",
    coef_rename = c("log_yield_hat" = "Log Yield Hat",
                    "genderMale" = "Male",
                    "meal_typeOne" = "Meal: One",
                    "meal_typeTwo" = "Meal: Two",
                    "meal_typeThree" = "Meal: Three")
  )
}

## ── Save for figures ──────────────────────────────────────────────────────── ##
save(int_models, triple_models, quartile_results,
     file = file.path(ROOT, "output/stage2/models/stage2_extend_models.RData"))

cat("=== HETEROGENEITY COMPLETE ===\n")
