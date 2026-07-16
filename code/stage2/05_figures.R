## stage2_figures.R
## Publication figures for stage 2
## FIG1: Coefplot main models, FIG2: Quartile pass-through
## FIG3: Randomization distribution, FIG4: Gender × meal_type

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(fixest)
  library(ggplot2)
  library(patchwork)
})

ROOT    <- here::here()
out_fig <- file.path(ROOT, "output/stage2/figures")
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

# Okabe-Ito palette
OI <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
        "#0072B2", "#D55E00", "#CC79A7", "#000000")

cat("=== STAGE 2 FIGURES ===\n")

## ── Load models (LEVELS spec) ──────────────────────────────────────────────── ##
load(file.path(ROOT, "output/stage2/models/stage2_main_models_v2.RData"))

extend_path <- file.path(ROOT, "output/stage2/models/stage2_extend_models.RData")
use_cached <- FALSE
if (file.exists(extend_path)) {
  tmp <- new.env()
  load(extend_path, envir = tmp)
  if (length(tmp$quartile_results) > 0) {
    qm <- tmp$quartile_results[[1]][[1]]
    use_cached <- !is.null(qm) && "log_yield_hat" %in% names(coef(qm))
  }
}
if (use_cached) {
  load(extend_path)
} else {
  cat("Note: refitting quartile models (levels) for FIG2...\n")
  struct_vars <- c(
    irrigation = "z_irrigation_2019",
    intensity  = "z_intensity_2019",
    holdings   = "z_holdings_2019",
    gca        = "z_gca_2019"
  )
  df_qsrc <- read_csv(file.path(ROOT, "data/Regression_data/df_2_merged_struct.csv"),
                      show_col_types = FALSE) %>%
    mutate(
      gender    = relevel(factor(gender), ref = "Female"),
      meal_type = relevel(factor(meal_type), ref = "None")
    )
  quartile_results <- list()
  for (v in names(struct_vars)) {
    zvar <- struct_vars[[v]]
    df_q <- df_qsrc %>%
      group_by(District) %>%
      mutate(q_var = first(.data[[zvar]])) %>%
      ungroup() %>%
      filter(!is.na(q_var)) %>%
      mutate(quartile = ntile(q_var, 4))
    q_models <- list()
    for (q in 1:4) {
      df_sub <- df_q %>% filter(quartile == q)
      if (nrow(df_sub) < 30) next
      q_models[[paste0("Q", q)]] <- tryCatch(
        feols(log_real_wage ~ log_yield_hat + gender + meal_type |
                year + District^growing_season,
              data = df_sub, cluster = ~District),
        error = function(e) NULL
      )
    }
    quartile_results[[v]] <- Filter(Negate(is.null), q_models)
  }
}

## ── Helper: extract coef + CI from fixest model ───────────────────────────── ##
extract_coef <- function(model, param, label) {
  if (is.null(model)) return(NULL)
  cf  <- coef(model)[param]
  se_ <- se(model)[param]
  if (is.na(cf) || is.na(se_)) return(NULL)
  data.frame(
    label = label,
    est   = cf,
    lo    = cf - 1.96 * se_,
    hi    = cf + 1.96 * se_,
    p     = pvalue(model)[param],
    stringsAsFactors = FALSE
  )
}

## ── FIG1: Coefplot across M1–M4 + key robustness specs ───────────────────── ##
cat("Building FIG1...\n")

# Load robustness models inline for coefplot (LEVELS)
df_main <- read_csv(file.path(ROOT, "data/Regression_data/df_2_merged_v2.csv"),
                    show_col_types = FALSE) %>%
  mutate(gender = relevel(factor(gender), ref = "Female"),
         meal_type = relevel(factor(meal_type), ref = "None"))

rob5_m <- tryCatch(
  feols(log_real_wage ~ log_yield_hat + gender + meal_type |
          year + District^growing_season + District[year],
        data = df_main, cluster = ~District),
  error = function(e) NULL
)
rob6_m <- {
  q01 <- quantile(df_main$log_real_wage, 0.01, na.rm = TRUE)
  q99 <- quantile(df_main$log_real_wage, 0.99, na.rm = TRUE)
  df_w <- df_main %>% mutate(log_real_wage = pmax(pmin(log_real_wage, q99), q01))
  feols(log_real_wage ~ log_yield_hat + gender + meal_type |
          year + District^growing_season, data = df_w, cluster = ~District)
}

coef_list <- list(
  extract_coef(m1,    "log_yield_hat", "M1: Yield only"),
  extract_coef(m2,    "log_yield_hat", "M2: + Gender"),
  extract_coef(m3,    "log_yield_hat", "M3: Main"),
  extract_coef(m4,    "log_yield_hat", "M4: Log wage"),
  extract_coef(rob5_m,"log_yield_hat", "ROB5: + Trend"),
  extract_coef(rob6_m,"log_yield_hat", "ROB6: Winsorize")
)
df_coef <- dplyr::bind_rows(Filter(Negate(is.null), coef_list))
df_coef$label <- factor(df_coef$label, levels = rev(df_coef$label))
df_coef$sig   <- ifelse(df_coef$p < 0.1, "sig", "insig")

fig1 <- ggplot(df_coef, aes(x = est, y = label, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.7) +
  geom_point(size = 3) +
  scale_color_manual(values = c(sig = OI[6], insig = OI[2]),
                     labels = c(sig = "p<0.1", insig = "p>=0.1"),
                     name   = NULL) +
  labs(x = "Coefficient on log(Yield Hat)",
       y = NULL,
       title = "Wage Pass-Through: Main and Robustness Specs") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(out_fig, "coefplot_main.png"), fig1,
       width = 8, height = 5, dpi = 300)
cat("Saved coefplot_main.png\n")

## ── FIG2: Quartile pass-through (2x2 grid) ────────────────────────────────── ##
cat("Building FIG2...\n")
struct_labels <- c(
  irrigation = "Irrigation",
  intensity  = "Crop Intensity",
  holdings   = "Holdings",
  gca        = "Gross Cropped Area"
)

panels_q <- list()
for (v in names(quartile_results)) {
  q_models <- quartile_results[[v]]
  rows <- lapply(names(q_models), function(qn) {
    m <- q_models[[qn]]
    if (is.null(m)) return(NULL)
    extract_coef(m, "log_yield_hat", qn)
  })
  rows <- dplyr::bind_rows(Filter(Negate(is.null), rows))
  if (nrow(rows) == 0) next
  rows$quartile <- as.integer(gsub("Q", "", rows$label))
  rows$sig      <- ifelse(rows$p < 0.1, "sig", "insig")
  panels_q[[v]] <- ggplot(rows, aes(x = quartile, y = est, color = sig)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) +
    geom_point(size = 3) +
    scale_color_manual(values = c(sig = OI[6], insig = OI[2]), name = NULL) +
    scale_x_continuous(breaks = 1:4) +
    labs(x = "Quartile", y = "Coef (log Yield Hat)",
         title = struct_labels[v]) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none")
}

if (length(panels_q) >= 1) {
  n_panels <- length(panels_q)
  ncols    <- min(2, n_panels)
  fig2 <- wrap_plots(panels_q, ncol = ncols) +
    plot_annotation(title = "Quartile Pass-Through by Structural Variable")
  ggsave(file.path(out_fig, "quartile_passthrough.png"), fig2,
         width = 10, height = 8, dpi = 300)
  cat("Saved quartile_passthrough.png\n")
} else {
  cat("FIG2 skipped: no quartile results\n")
}

## ── FIG3: Randomization distribution ─────────────────────────────────────── ##
cat("Building FIG3...\n")
rand_file <- file.path(ROOT, "output/stage2/models/rob4_randomization.RData")
if (file.exists(rand_file)) {
  load(rand_file)  # perm_coefs, actual_coef, rand_p
  df_rand <- data.frame(coef = perm_coefs[!is.na(perm_coefs)])
  fig3 <- ggplot(df_rand, aes(x = coef)) +
    geom_histogram(fill = OI[2], color = "white", bins = 40) +
    geom_vline(xintercept = actual_coef, color = OI[6], linewidth = 1.2, linetype = "solid") +
    annotate("text", x = actual_coef * 1.05, y = Inf, vjust = 1.5,
             label = sprintf("Actual coef\n%.2f\np=%.3f", actual_coef, rand_p),
             color = OI[6], size = 3.5, hjust = 0) +
    labs(x = "Permuted Coefficient",
         y = "Count",
         title = "Randomization Test: Distribution of Permuted Coefficients") +
    theme_minimal(base_size = 12)
  ggsave(file.path(out_fig, "randomization_test.png"), fig3,
         width = 8, height = 5, dpi = 300)
  cat("Saved randomization_test.png\n")
} else {
  cat("FIG3 skipped: randomization RData not found\n")
}

## ── FIG4: Gender × meal_type coefplot ────────────────────────────────────── ##
cat("Building FIG4...\n")
# Extract gender and meal_type coefficients from M3
coef_m3 <- coef(m3)
se_m3   <- se(m3)
p_m3    <- pvalue(m3)

gm_terms <- names(coef_m3)[grepl("gender|meal_type", names(coef_m3))]
df_gm <- data.frame(
  term  = gm_terms,
  est   = coef_m3[gm_terms],
  lo    = coef_m3[gm_terms] - 1.96 * se_m3[gm_terms],
  hi    = coef_m3[gm_terms] + 1.96 * se_m3[gm_terms],
  p     = p_m3[gm_terms],
  label = c("Male" = "Male", "meal_typeOne" = "Meal: One",
            "meal_typeThree" = "Meal: Three", "meal_typeTwo" = "Meal: Two")[gm_terms],
  stringsAsFactors = FALSE
)
df_gm$label[is.na(df_gm$label)] <- df_gm$term[is.na(df_gm$label)]
df_gm$sig   <- ifelse(df_gm$p < 0.1, "sig", "insig")
df_gm$label <- factor(df_gm$label, levels = rev(df_gm$label))

fig4 <- ggplot(df_gm, aes(x = est, y = label, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.3, linewidth = 0.7) +
  geom_point(size = 3) +
  scale_color_manual(values = c(sig = OI[6], insig = OI[2]),
                     labels = c(sig = "p<0.1", insig = "p>=0.1"), name = NULL) +
  labs(x = "Coefficient (log wages vs reference)",
       y = NULL,
       title = "Gender and Meal-Type Wage Differentials (M3 Main Spec)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(out_fig, "gender_meal_coefplot.png"), fig4,
       width = 8, height = 4, dpi = 300)
cat("Saved gender_meal_coefplot.png\n")

# ============================================================
# FIGURE 6: Raw wage trends by meal type
# ============================================================
cat("Building FIG6...\n")

suppressPackageStartupMessages({
  library(RColorBrewer)
  library(ggrepel)
})

## Load level wages from v2 panel (real_wage, not diff)
df_fig6 <- tryCatch(
  read_csv(file.path(ROOT, "data/Regression_data/df_2_merged_v2.csv"),
           show_col_types = FALSE),
  error = function(e) {
    cat("FIG6: df_2_merged_v2.csv not found, trying df_2_merged.csv\n")
    read_csv(file.path(ROOT, "data/Regression_data/df_2_merged.csv"),
             show_col_types = FALSE)
  }
) %>%
  mutate(
    gender    = relevel(factor(gender),    ref = "Female"),
    meal_type = relevel(factor(meal_type), ref = "None")
  )

## Compute group means and SEs
df_fig6_agg <- df_fig6 %>%
  group_by(year, gender, meal_type) %>%
  summarise(
    mean_wage = mean(real_wage, na.rm = TRUE),
    se_wage   = sd(real_wage,   na.rm = TRUE) / sqrt(sum(!is.na(real_wage))),
    .groups   = "drop"
  )

## Sanity check: flag if grand mean is outside expected range
grand_mean <- mean(df_fig6_agg$mean_wage, na.rm = TRUE)
if (grand_mean < 300 || grand_mean > 500) {
  # WARNING: grand mean real wage (grand_mean BDT/day) is outside the expected
  # 300–500 BDT/day range. Check CPI deflation or data source before publishing.
  cat(sprintf("FIG6 WARNING: grand mean wage = %.1f BDT/day (expected 300-500).\n",
              grand_mean))
}

## Rightmost-year labels for direct line labelling
df_fig6_labels <- df_fig6_agg %>%
  filter(year == max(year))

## Set2 palette (4 meal types)
set2_cols <- brewer.pal(4, "Set2")
names(set2_cols) <- c("None", "One", "Two", "Three")

fig6 <- ggplot(df_fig6_agg,
               aes(x = year, y = mean_wage,
                   colour = meal_type, fill = meal_type, group = meal_type)) +
  geom_ribbon(aes(ymin = mean_wage - se_wage,
                  ymax = mean_wage + se_wage),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_text(
    data    = df_fig6_labels,
    aes(label = meal_type),
    hjust   = -0.15,
    size    = 3.5,
    fontface = "plain",
    show.legend = FALSE
  ) +
  facet_wrap(~gender, ncol = 2) +
  scale_colour_manual(values = set2_cols, name = "Meal type") +
  scale_fill_manual(  values = set2_cols, name = "Meal type") +
  scale_x_continuous(
    breaks = 2017:2025,
    expand = expansion(mult = c(0.02, 0.18))   # right margin for labels
  ) +
  labs(
    x       = "Year",
    y       = "Mean real wage (BDT/day)",
    caption = paste(
      "Note: Lines show mean real wages (BDT/day) deflated by CPI.",
      "Shaded bands = \u00b11 SE. Source: BBS agricultural wage survey 2017\u20132025."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text        = element_text(size = 11),
    axis.title       = element_text(size = 11),
    strip.text       = element_text(size = 12, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90"),
    panel.grid.minor.y = element_blank(),
    legend.position  = "none",
    plot.caption     = element_text(size = 9, hjust = 0, colour = "grey40")
  )

ggsave(file.path(out_fig, "fig6_wage_trends_by_meal.png"), fig6,
       width = 10, height = 5, dpi = 300)
ggsave(file.path(out_fig, "fig6_wage_trends_by_meal.pdf"), fig6,
       width = 10, height = 5)
cat("Saved fig6_wage_trends_by_meal.png/.pdf\n")

cat("=== FIGURES COMPLETE ===\n")
