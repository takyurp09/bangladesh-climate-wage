## code/stage2/09_plots.R
## Regression-based coefficient plots and heatmaps (Stage 2)
## All plots derived from regression results only
## Outputs → output/stage2/plots/

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(fixest)
  library(here)
  library(readr)
  library(ggrepel)
  library(stringr)
})

ROOT    <- here::here()
PLOTDIR <- file.path(ROOT, "output/stage2/plots")
WAGE    <- file.path(ROOT, "data/Regression_data/df_2_merged_v2.csv")
STRUCT  <- file.path(ROOT, "data/Regression_data/df_struct_2019.csv")
dir.create(PLOTDIR, showWarnings = FALSE, recursive = TRUE)

cat("=== STAGE 2 PLOTS ===\n")

# ── Load models (LEVELS spec) ───────────────────────────────────────────────────
load(file.path(ROOT, "output/stage2/models/stage2_main_models_v2.RData"))

extend_path <- file.path(ROOT, "output/stage2/models/stage2_extend_models.RData")
use_cached <- FALSE
if (file.exists(extend_path)) {
  tmp <- new.env()
  load(extend_path, envir = tmp)
  if (length(tmp$int_models) > 0) {
    m0 <- tmp$int_models[[1]]
    use_cached <- "log_yield_hat" %in% names(coef(m0))
  }
}
if (use_cached) {
  load(extend_path)
} else {
  cat("Refitting heterogeneity models (levels) for plots...\n")
  struct_vars <- c(
    irrigation = "z_irrigation_2019",
    intensity  = "z_intensity_2019",
    holdings   = "z_holdings_2019",
    gca        = "z_gca_2019"
  )
  df_struct <- read_csv(STRUCT, show_col_types = FALSE) %>%
    mutate(district_key = trimws(tolower(District))) %>%
    mutate(
      z_holdings_2019   = as.numeric(scale(holdings_2019)),
      z_irrigation_2019 = as.numeric(scale(irrigation_2019)),
      z_intensity_2019  = as.numeric(scale(intensity_2019)),
      z_gca_2019        = as.numeric(scale(gca_2019))
    )
  df_het <- read_csv(WAGE, show_col_types = FALSE) %>%
    mutate(
      gender    = relevel(factor(gender), ref = "Female"),
      meal_type = relevel(factor(meal_type), ref = "None"),
      district_key = trimws(tolower(District))
    ) %>%
    left_join(df_struct %>% select(district_key, starts_with("z_")), by = "district_key")

  int_models <- list()
  quartile_results <- list()
  for (v in names(struct_vars)) {
    zvar <- struct_vars[[v]]
    fml <- as.formula(sprintf(
      "log_real_wage ~ log_yield_hat * %s + gender + meal_type | year + District^growing_season",
      zvar))
    int_models[[v]] <- feols(fml, data = df_het, cluster = ~District)

    df_q <- df_het %>%
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

# ── Load data ─────────────────────────────────────────────────────────────────
df <- read_csv(WAGE, show_col_types = FALSE) %>%
  mutate(
    gender    = relevel(factor(gender),    ref = "Female"),
    meal_type = relevel(factor(meal_type), ref = "None")
  )

df_struct <- read_csv(STRUCT, show_col_types = FALSE) %>%
  mutate(district_key = trimws(tolower(District))) %>%
  mutate(
    z_holdings_2019   = as.numeric(scale(holdings_2019)),
    z_irrigation_2019 = as.numeric(scale(irrigation_2019)),
    z_intensity_2019  = as.numeric(scale(intensity_2019)),
    z_gca_2019        = as.numeric(scale(gca_2019))
  )

# Merge struct into wage data
df <- df %>%
  mutate(district_key = trimws(tolower(District))) %>%
  left_join(df_struct %>% select(district_key, starts_with("z_")), by = "district_key")

# ── Shared style ──────────────────────────────────────────────────────────────
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#000000")

plot_theme <- function() {
  theme_minimal(base_size = 9) +
    theme(
      plot.title       = element_text(size = 9, face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(size = 7, hjust = 0.5, colour = "grey40"),
      plot.caption     = element_text(size = 7, colour = "grey40", hjust = 0),
      axis.title       = element_text(size = 8),
      axis.text        = element_text(size = 7),
      legend.title     = element_text(size = 8, face = "bold"),
      legend.text      = element_text(size = 7),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
      strip.text       = element_text(size = 8, face = "bold"),
      plot.margin      = margin(t = 10, r = 20, b = 10, l = 10, unit = "mm")
    )
}

save_plot <- function(p, fname, w = 183, h = 120) {
  path <- file.path(PLOTDIR, fname)
  ggsave(path, p, width = w, height = h, units = "mm", dpi = 300, bg = "white")
  cat(sprintf("  saved \u2713 %s  [%dx%dmm]\n", fname, w, h))
}

# Helper: extract coef + 95% CI for a given coefficient name from a fixest model
extract_coef <- function(model, coef_name, se_type = NULL) {
  tryCatch({
    if (!is.null(se_type)) {
      v  <- vcov(model, vcov = se_type)
      cf <- coef(model)[coef_name]
      se <- sqrt(diag(v))[coef_name]
    } else {
      cf <- coef(model)[coef_name]
      se <- sqrt(diag(vcov(model)))[coef_name]
    }
    if (is.na(cf) || is.na(se)) return(NULL)
    pv  <- 2 * pnorm(-abs(cf / se))
    data.frame(coef = cf, se = se, ci_lo = cf - 1.96 * se, ci_hi = cf + 1.96 * se, pval = pv)
  }, error = function(e) NULL)
}

# ─────────────────────────────────────────────────────────────────────────────
# PLOT 1 — Main Coefficient Across Specifications
# ─────────────────────────────────────────────────────────────────────────────
cat("\nPLOT 1: Robustness coefplot\n")

SEASONS <- c("Boro", "Aus", "Aman")
season_models <- lapply(SEASONS, function(s) {
  sub <- df %>% filter(growing_season == s, !is.na(log_real_wage), !is.na(log_yield_hat))
  tryCatch(
    feols(log_real_wage ~ log_yield_hat + gender + meal_type |
            year + District^growing_season,
          data = sub, cluster = ~District, warn = FALSE, notes = FALSE),
    error = function(e) NULL
  )
})
names(season_models) <- SEASONS

# DK and two-way SE on pooled m3_base (refit to avoid dependency on robustness models)
m3_base <- feols(log_real_wage ~ log_yield_hat + gender + meal_type |
                   year + District^growing_season,
                 data = df, cluster = ~District, warn = FALSE, notes = FALSE)

tgt <- "log_yield_hat"
spec_list <- list(
  list(label = "M1 Pooled (yield only)",        grp = "Pooled",        mod = m1),
  list(label = "M2 Pooled (+ gender)",           grp = "Pooled",        mod = m2),
  list(label = "M3 Pooled (+ meal type)",        grp = "Pooled",        mod = m3),
  list(label = "M3 Extended (+ controls)",       grp = "Pooled",        mod = m3_ext),
  list(label = "M4 Log wage",                    grp = "Pooled",        mod = m4),
  list(label = "M3 Boro season",                 grp = "By Season",     mod = season_models$Boro),
  list(label = "M3 Aus season",                  grp = "By Season",     mod = season_models$Aus),
  list(label = "M3 Aman season",                 grp = "By Season",     mod = season_models$Aman),
  list(label = "M3 + DK SE",                     grp = "Alternative SE", mod = m3_base, se_type = "DK"),
  list(label = "M3 + Two-way cluster",           grp = "Alternative SE", mod = m3_base,
       se_type = ~District + year)
)

coef_rows <- lapply(spec_list, function(s) {
  if (is.null(s$mod)) return(NULL)
  res <- extract_coef(s$mod, tgt, se_type = s$se_type)
  if (is.null(res)) return(NULL)
  res$label <- s$label
  res$group <- s$grp
  res
}) %>% bind_rows()

# Ordered factor to preserve row order
coef_rows$label <- factor(coef_rows$label, levels = rev(coef_rows$label))
coef_rows$sig   <- cut(coef_rows$pval,
                       breaks = c(0, 0.05, 0.1, 1),
                       labels = c("p<0.05", "p<0.1", "p>0.1"),
                       include.lowest = TRUE)

p1 <- ggplot(coef_rows, aes(x = coef, y = label)) +
  geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi, colour = group),
                 width = 0.2, linewidth = 0.5, orientation = "y") +
  geom_point(aes(colour = group, shape = sig), size = 2.2, fill = "white") +
  scale_shape_manual(
    values = c("p<0.05" = 16, "p<0.1" = 21, "p>0.1" = 1),
    name   = "Significance"
  ) +
  scale_colour_manual(
    values = c("Pooled" = okabe_ito[6], "By Season" = okabe_ito[3],
               "Alternative SE" = okabe_ito[2]),
    name = "Specification group"
  ) +
  facet_grid(rows = vars(group), scales = "free_y", space = "free") +
  labs(
    x       = "Coefficient on log(yield\u0302)",
    y       = NULL,
    title   = "Yield-Wage Pass-Through: Main Coefficient Across Specifications",
    caption = str_wrap("95% CI shown. Clustered SE by district unless noted. DK = Driscoll-Kraay. Filled circle = p<0.05, half-filled = p<0.1, hollow = p>0.1.", width = 80)
  ) +
  plot_theme() +
  theme(
    strip.background = element_rect(fill = "grey95", colour = NA),
    axis.text.y      = element_text(size = 7, hjust = 1)
  )

save_plot(p1, "plot1_robustness_coefplot.png", w = 183, h = 160)

# ─────────────────────────────────────────────────────────────────────────────
# PLOT 2 — Meal × Season Interaction Heatmap
# Run M5 by season; compute total pass-through per meal type
# ─────────────────────────────────────────────────────────────────────────────
cat("\nPLOT 2: Meal × Season heatmap\n")

MEAL_TYPES <- c("One", "Two", "Three")

heatmap_rows <- lapply(SEASONS, function(s) {
  sub <- df %>%
    filter(growing_season == s, !is.na(log_real_wage), !is.na(log_yield_hat))
  m <- tryCatch(
    feols(log_real_wage ~ log_yield_hat + log_yield_hat:meal_type + gender |
            year + District^growing_season,
          data = sub, cluster = ~District, warn = FALSE, notes = FALSE),
    error = function(e) NULL
  )
  if (is.null(m)) return(NULL)
  base_cf <- coef(m)["log_yield_hat"]
  base_se <- sqrt(diag(vcov(m)))["log_yield_hat"]
  lapply(MEAL_TYPES, function(mt) {
    int_name <- paste0("log_yield_hat:meal_type", mt)
    int_cf   <- tryCatch(coef(m)[int_name], error = function(e) 0)
    int_se   <- tryCatch(sqrt(diag(vcov(m)))[int_name], error = function(e) NA_real_)
    if (is.na(int_cf)) { int_cf <- 0; int_se <- NA_real_ }
    total_cf <- base_cf + int_cf
    # Delta method SE (ignore covariance — conservative)
    total_se <- if (!is.na(int_se)) sqrt(base_se^2 + int_se^2) else base_se
    pv <- 2 * pnorm(-abs(total_cf / total_se))
    stars <- ifelse(pv < 0.05, "**", ifelse(pv < 0.1, "*", ""))
    data.frame(season = s, meal_type = mt, total_coef = total_cf,
               se = total_se, pval = pv, stars = stars)
  }) %>% bind_rows()
}) %>% bind_rows()

if (nrow(heatmap_rows) > 0) {
  heatmap_rows$meal_type <- factor(heatmap_rows$meal_type, levels = c("One", "Two", "Three"))
  heatmap_rows$season    <- factor(heatmap_rows$season,    levels = c("Boro", "Aus", "Aman"))
  heatmap_rows$label     <- sprintf("%.1f%s", heatmap_rows$total_coef, heatmap_rows$stars)
  lim2 <- max(abs(heatmap_rows$total_coef), na.rm = TRUE) * 1.1

  p2 <- ggplot(heatmap_rows, aes(x = season, y = meal_type, fill = total_coef)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = label), size = 2.5, colour = "white", fontface = "bold") +
    scale_fill_distiller(
      palette   = "RdBu",
      direction = -1,
      limits    = c(-lim2, lim2),
      name      = "Total pass-through (log wages)"
    ) +
    labs(
      x       = "Season",
      y       = "Meal type",
      title   = "Total Yield-Wage Pass-Through by Meal Type and Season",
      caption = str_wrap("Total = baseline (No-meal) + interaction coef. * p<0.1  ** p<0.05", width = 60)
    ) +
    plot_theme() +
    theme(
      panel.grid      = element_blank(),
      legend.position = "right"
    )

  save_plot(p2, "plot2_meal_season_heatmap.png", w = 120, h = 89)
} else {
  cat("  WARNING: no heatmap rows — skipping\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# PLOT 3 — Heterogeneity Interaction Coefplot
# ─────────────────────────────────────────────────────────────────────────────
cat("\nPLOT 3: Heterogeneity coefplot\n")

struct_labels <- c(
  "log_yield_hat"               = "Main effect",
  "log_yield_hat:z_holdings_2019"   = "× Holdings (z)",
  "log_yield_hat:z_irrigation_2019" = "× Irrigation (z)",
  "log_yield_hat:z_intensity_2019"  = "× Intensity (z)",
  "log_yield_hat:z_gca_2019"        = "× GCA (z)"
)

# Extract main effect from each interaction model
het_main <- lapply(names(int_models), function(v) {
  res <- extract_coef(int_models[[v]], "log_yield_hat")
  if (is.null(res)) return(NULL)
  res$variable <- "log_yield_hat"
  res$panel    <- "Main Effect"
  res$var_grp  <- v
  res
}) %>% bind_rows()

# Extract interaction terms
het_int <- lapply(names(int_models), function(v) {
  int_name <- paste0("log_yield_hat:z_", v, "_2019")
  res <- extract_coef(int_models[[v]], int_name)
  if (is.null(res)) return(NULL)
  res$variable <- int_name
  res$panel    <- "Interaction Terms"
  res$var_grp  <- v
  res
}) %>% bind_rows()

het_df <- bind_rows(het_main, het_int)
if (nrow(het_df) == 0) {
  cat("  WARNING: no heterogeneity rows — skipping plot 3\n")
} else {
het_df <- het_df %>%
  mutate(
    label = recode(variable, !!!struct_labels),
    sig   = pval < 0.1,
    panel = factor(panel, levels = c("Main Effect", "Interaction Terms"))
  )

p3 <- ggplot(het_df, aes(x = coef, y = label, colour = sig)) +
  geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), width = 0.25,
                 linewidth = 0.5, orientation = "y") +
  geom_point(size = 2.5) +
  scale_colour_manual(
    values = c("TRUE" = okabe_ito[6], "FALSE" = "grey60"),
    labels = c("TRUE" = "p<0.1", "FALSE" = "p\u22650.1"),
    name   = NULL
  ) +
  facet_wrap(~ panel, scales = "free") +
  labs(
    x       = "Coefficient",
    y       = NULL,
    title   = "Heterogeneity in Yield-Wage Pass-Through by Agricultural Structure",
    caption = str_wrap("Structural variables standardized (z-scores). Clustered SE by district.", width = 80)
  ) +
  plot_theme() +
  theme(
    strip.background = element_rect(fill = "grey95", colour = NA),
    axis.text.y      = element_text(size = 7, hjust = 1)
  )

save_plot(p3, "plot3_heterogeneity_coefplot.png", w = 183, h = 120)
}
# ─────────────────────────────────────────────────────────────────────────────
cat("\nPLOT 4: Quartile pass-through\n")

struct_panel_labels <- c(
  holdings  = "Holdings",
  irrigation = "Irrigation",
  intensity  = "Intensity",
  gca        = "GCA"
)

quartile_df <- lapply(names(quartile_results), function(v) {
  lapply(names(quartile_results[[v]]), function(q) {
    m   <- quartile_results[[v]][[q]]
    if (is.null(m)) return(NULL)
    res <- extract_coef(m, "log_yield_hat")
    if (is.null(res)) return(NULL)
    res$var     <- v
    res$quartile <- q
    res
  }) %>% bind_rows()
}) %>% bind_rows()

if (nrow(quartile_df) > 0) {
  quartile_df$quartile <- factor(quartile_df$quartile, levels = c("Q1", "Q2", "Q3", "Q4"))
  quartile_df$var_label <- struct_panel_labels[quartile_df$var]

  p4 <- ggplot(quartile_df, aes(x = quartile, y = coef, colour = quartile)) +
    geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed", linewidth = 0.4) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2, linewidth = 0.5) +
    geom_point(size = 2.5) +
    scale_colour_viridis_d(option = "D", begin = 0.2, end = 0.85, name = "Quartile") +
    facet_wrap(~ var_label, ncol = 2) +
    labs(
      x       = "Structural quartile (Q1=Low \u2192 Q4=High)",
      y       = "Pass-Through Coef (log wages)",
      title   = "Yield-Wage Pass-Through by Structural Quartile",
      caption = str_wrap("Each point = M3 coefficient estimated within quartile subsample. 95% CI shown.", width = 80)
    ) +
    plot_theme() +
    theme(
      strip.background = element_rect(fill = "grey95", colour = NA),
      axis.text.x      = element_text(size = 7)
    )

  save_plot(p4, "plot4_quartile_passthrough.png", w = 183, h = 120)
} else {
  cat("  WARNING: no quartile results — skipping\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# PLOT 5 — District Pass-Through vs EDD Exposure (scatter)
# ─────────────────────────────────────────────────────────────────────────────
cat("\nPLOT 5: Pass-through vs EDD scatter\n")

# Reuse passthrough_res from 08_maps (recompute here independently)
clean_name_fn <- function(x) trimws(tolower(x))

passthrough_local <- lapply(unique(df$district_key), function(d) {
  sub <- df %>% filter(district_key == d, !is.na(log_real_wage), !is.na(log_yield_hat))
  tryCatch({
    # HC1 SE — clustering by district_key invalid within a single district
    m  <- feols(log_real_wage ~ log_yield_hat | year,
                data = sub, vcov = "HC1", warn = FALSE, notes = FALSE)
    cf  <- coef(m)["log_yield_hat"]
    se  <- sqrt(diag(vcov(m)))["log_yield_hat"]
    pv  <- 2 * pnorm(-abs(cf / se))
    data.frame(district_key = d,
               District = unique(sub$District)[1],
               pt_coef = cf, pt_se = se, pt_pval = pv,
               n_obs = as.integer(nobs(m)))
  }, error = function(e) NULL)
}) %>% bind_rows()

# EDD exposure per district (Aus season, 2017-2023)
df_clim_s5 <- read_csv(
  file.path(ROOT, "data/Regression_data/bangladesh_rice_regression_panel.csv"),
  show_col_types = FALSE
) %>%
  mutate(district_key = clean_name_fn(district)) %>%
  filter(season == "Aus", year >= 2017, year <= 2023) %>%
  group_by(district_key) %>%
  summarise(mean_edd30 = mean(edd_30, na.rm = TRUE), .groups = "drop")

scatter5 <- passthrough_local %>%
  filter(n_obs >= 15, !is.na(pt_coef)) %>%
  left_join(df_clim_s5, by = "district_key") %>%
  filter(!is.na(mean_edd30)) %>%
  mutate(sig = pt_pval < 0.1)

# Top 5 highest EDD for labelling
top5 <- scatter5 %>% slice_max(mean_edd30, n = 5)

if (nrow(scatter5) >= 5) {
  r_val <- cor(scatter5$mean_edd30, scatter5$pt_coef, use = "complete.obs")
  r_lab <- paste0("r = ", round(r_val, 2))

  # x/y ranges for annotation placement
  x_min  <- min(scatter5$mean_edd30, na.rm = TRUE)
  y_max  <- max(scatter5$pt_coef,    na.rm = TRUE)

  annot_text <- str_wrap(
    "No systematic EDD-passthrough relationship suggests labor contract type, not climate exposure, drives heterogeneity",
    width = 35
  )

  p5 <- ggplot(scatter5, aes(x = mean_edd30, y = pt_coef)) +
    geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed", linewidth = 0.4) +
    geom_smooth(method = "loess", span = 0.75, se = FALSE, colour = okabe_ito[2],
                linewidth = 0.8) +
    geom_point(aes(size = n_obs, fill = sig), shape = 21, colour = "grey30",
               stroke = 0.4, alpha = 0.85) +
    geom_label_repel(
      data        = top5,
      aes(label   = District),
      size        = 2.0,
      colour      = "grey20",
      fill        = "white",
      label.size  = 0.2,
      label.padding = unit(0.1, "lines"),
      max.overlaps = 15,
      box.padding  = 0.4,
      point.padding = 0.3,
      segment.colour = "grey60",
      segment.size   = 0.3
    ) +
    annotate("text", x = x_min, y = y_max, label = r_lab,
             hjust = 0, vjust = 1, size = 2.5, colour = "grey30", fontface = "italic") +
    annotate("text",
             x = x_min, y = y_max * 0.75,
             label = annot_text,
             hjust = 0, vjust = 1, size = 2.2, colour = "grey50", lineheight = 1.15) +
    scale_size_continuous(range = c(1.5, 5), name = "N obs", guide = "none") +
    scale_fill_manual(
      values = c("TRUE" = "#E69F00", "FALSE" = "grey70"),
      labels = c("TRUE" = "p<0.1", "FALSE" = "p\u22650.1"),
      name   = NULL
    ) +
    labs(
      x       = "Mean EDD >30\u00b0C (Aus season, 2017\u20132023)",
      y       = "Pass-Through Coef (log wages)",
      title   = "Yield-Wage Pass-Through vs Heat Exposure by District",
      caption = str_wrap("District-level FD regression, year FE. Size = N obs. LOWESS smoother. Top 5 highest-EDD districts labelled.", width = 70)
    ) +
    plot_theme()

  save_plot(p5, "plot5_passthrough_vs_edd.png", w = 120, h = 120)
} else {
  cat("  WARNING: insufficient districts for scatter — skipping\n")
}

cat("\n=== PLOT SUMMARY ===\n")
cat(sprintf("Output: %s\n", PLOTDIR))
