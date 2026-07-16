# ============================================================
# fig_stage1_visualizations.R
# Four Stage 1 publication-quality figures for Appendix A
#
# Figure 1: GDD/EDD Construction Diagram (synthetic curve)
# Figure 2: Nonparametric Temperature Bin Coefficients (2-panel)
# Figure 3: Threshold Sensitivity Coefficient Plot (2-panel stacked)
# Figure 4: Main First-Stage Coefficient Plot (4-panel)
#
# Output: output/stage1/figures/*.pdf
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(fixest)
  library(scales)
  library(arrow)
  library(here)
})

ROOT    <- here::here()
FIGDIR  <- file.path(ROOT, "output/stage1/figures")
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

# ── Global theme (matches existing paper figures) ─────────────────────────────
paper_theme <- function(base_size = 11) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      plot.title       = element_text(face = "bold", size = base_size, hjust = 0),
      plot.subtitle    = element_text(size = base_size - 1, hjust = 0, color = "grey40"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      axis.text        = element_text(size = base_size - 1, color = "grey30"),
      axis.title       = element_text(size = base_size - 0.5),
      strip.text       = element_text(face = "bold", size = base_size),
      legend.position  = "bottom",
      legend.text      = element_text(size = base_size - 1),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

theme_set(paper_theme())

# Okabe-Ito palette (colorblind-safe) — consistent with existing figures
OI <- c(Boro = "#E69F00", Aus = "#56B4E9", Aman = "#009E73",
        Pooled = "#999999")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1: GDD/EDD Construction Diagram
# ══════════════════════════════════════════════════════════════════════════════
cat("=== Figure 1: GDD/EDD Construction Diagram ===\n")

# Synthetic daily temperature via sinusoidal interpolation
T_mean <- 27; T_max <- 35; T_min <- 19
h_grid <- seq(0, 24, by = 0.1)
T_h    <- T_mean + (T_max - T_min) / 2 * cos(pi * (h_grid - 14) / 12)

df_curve <- data.frame(hour = h_grid, temp = T_h)

# Shading data: GDD zone [10, 30] and EDD zone [>30]
df_gdd  <- df_curve %>% mutate(
  y_lo = pmax(pmin(temp, 30), 10),
  y_hi = pmin(temp, 30),
  fill_gdd = ifelse(temp >= 10, pmin(temp, 30), NA)
)

# Build ribbon data
df_ribbon <- df_curve %>%
  mutate(
    gdd_lo = 10,
    gdd_hi = ifelse(temp >= 10, pmin(temp, 30), 10),
    edd_lo = 30,
    edd_hi = ifelse(temp > 30, temp, 30)
  )

fig1 <- ggplot(df_curve, aes(x = hour, y = temp)) +
  # EDD shading (above 30°C) — only where temp > 30
  geom_ribbon(data = df_ribbon %>% filter(temp > 30),
              aes(x = hour, ymin = edd_lo, ymax = edd_hi),
              fill = "#B71C1C", alpha = 0.30, inherit.aes = FALSE) +
  # GDD shading (10–30°C) — only where temp >= 10
  geom_ribbon(data = df_ribbon %>% filter(temp >= 10),
              aes(x = hour, ymin = gdd_lo, ymax = gdd_hi),
              fill = "#2E7D32", alpha = 0.30, inherit.aes = FALSE) +
  # Temperature curve
  geom_line(color = "#1565C0", linewidth = 1.1) +
  # Reference lines
  geom_hline(yintercept = 10, linetype = "dashed", color = "#2E7D32",
             linewidth = 0.7) +
  geom_hline(yintercept = 30, linetype = "solid", color = "#B71C1C",
             linewidth = 0.7) +
  # Zone labels
  annotate("text", x = 2, y = 21, label = "GDD accumulation [10, 30°C]",
           color = "#2E7D32", size = 3.2, fontface = "bold", hjust = 0) +
  annotate("text", x = 2, y = 33.5, label = "EDD accumulation (>30°C)",
           color = "#B71C1C", size = 3.2, fontface = "bold", hjust = 0) +
  # Reference line labels
  annotate("text", x = 23.5, y = 11.2, label = "GDD lower bound",
           color = "#2E7D32", size = 2.8, hjust = 1) +
  annotate("text", x = 23.5, y = 31.2, label = "GDD upper / EDD threshold",
           color = "#B71C1C", size = 2.8, hjust = 1) +
  scale_x_continuous(breaks = seq(0, 24, by = 4), limits = c(0, 24)) +
  scale_y_continuous(breaks = seq(10, 40, by = 5)) +
  labs(
    title = "Daily Temperature Profile and Degree-Day Construction",
    x     = "Hour of day",
    y     = "Temperature (°C)"
  ) +
  paper_theme(base_size = 11)

ggsave(
  file.path(FIGDIR, "fig_gdd_edd_definition.pdf"),
  fig1, width = 6, height = 4, device = "pdf", dpi = 300, bg = "white"
)
cat("✓ fig_gdd_edd_definition.pdf saved\n")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2: Nonparametric Temperature Bin Coefficients
# ══════════════════════════════════════════════════════════════════════════════
cat("\n=== Figure 2: Nonparametric Temperature Bin Coefficients ===\n")

bins_file <- file.path(ROOT, "data/Regression_data/bangladesh_rice_temperature_bins.parquet")
if (!file.exists(bins_file)) stop("Temperature bins parquet not found: ", bins_file)

df_bins <- read_parquet(bins_file) %>%
  mutate(
    district  = as.character(district),
    season    = as.character(season),
    year      = as.numeric(year),
    log_yield = log(yield_per_ha),
    temp_0_10  = temp_bin_0_1  + temp_bin_1_2  + temp_bin_2_3  + temp_bin_3_4  +
                 temp_bin_4_5  + temp_bin_5_6  + temp_bin_6_7  + temp_bin_7_8  +
                 temp_bin_8_9  + temp_bin_9_10,
    temp_10_15 = temp_bin_10_11 + temp_bin_11_12 + temp_bin_12_13 +
                 temp_bin_13_14 + temp_bin_14_15,
    temp_15_20 = temp_bin_15_16 + temp_bin_16_17 + temp_bin_17_18 +
                 temp_bin_18_19 + temp_bin_19_20,
    temp_20_25 = temp_bin_20_21 + temp_bin_21_22 + temp_bin_22_23 +
                 temp_bin_23_24 + temp_bin_24_25,
    temp_25_28 = temp_bin_25_26 + temp_bin_26_27 + temp_bin_27_28,
    temp_28_30 = temp_bin_28_29 + temp_bin_29_30,
    temp_30_32 = temp_bin_30_31 + temp_bin_31_32,
    temp_32_35 = temp_bin_32_33 + temp_bin_33_34 + temp_bin_34_35,
    temp_35_40 = temp_bin_35_36 + temp_bin_36_37 + temp_bin_37_38 +
                 temp_bin_38_39 + temp_bin_39_40
  ) %>%
  arrange(district, season, year) %>%
  group_by(district, season) %>%
  mutate(
    diff_log_yield  = log_yield - lag(log_yield),
    diff_temp_0_10  = temp_0_10  - lag(temp_0_10),
    diff_temp_10_15 = temp_10_15 - lag(temp_10_15),
    diff_temp_15_20 = temp_15_20 - lag(temp_15_20),
    diff_temp_20_25 = temp_20_25 - lag(temp_20_25),
    diff_temp_28_30 = temp_28_30 - lag(temp_28_30),
    diff_temp_30_32 = temp_30_32 - lag(temp_30_32),
    diff_temp_32_35 = temp_32_35 - lag(temp_32_35),
    diff_temp_35_40 = temp_35_40 - lag(temp_35_40)
  ) %>%
  ungroup() %>%
  filter(year >= 2013, !is.na(diff_log_yield))

df_boro <- df_bins %>% filter(season == "Boro")
df_aus  <- df_bins %>% filter(season == "Aus")
cat("After FD — Boro N:", nrow(df_boro), "| Aus N:", nrow(df_aus), "\n")

bin_formula <- diff_log_yield ~
  diff_temp_0_10 + diff_temp_10_15 + diff_temp_15_20 + diff_temp_20_25 +
  diff_temp_28_30 + diff_temp_30_32 + diff_temp_32_35 + diff_temp_35_40 | year

m_boro_bins <- feols(bin_formula, data = df_boro, cluster = ~district, warn = FALSE, notes = FALSE)
m_aus_bins  <- feols(bin_formula, data = df_aus,  cluster = ~district, warn = FALSE, notes = FALSE)
cat("Bin regressions complete\n")

bin_vars <- c(
  "diff_temp_0_10","diff_temp_10_15","diff_temp_15_20","diff_temp_20_25",
  "diff_temp_28_30","diff_temp_30_32","diff_temp_32_35","diff_temp_35_40"
)
bin_labels   <- c("[0,10)","[10,15)","[15,20)","[20,25)","[28,30)","[30,32)","[32,35)","[35,40]")
bin_midpoints <- c(5, 12.5, 17.5, 22.5, 29, 31, 33.5, 37.5)

extract_bins <- function(model, season_name) {
  ct <- coeftable(model)
  data.frame(
    season      = season_name,
    bin_label   = bin_labels,
    bin_mid     = bin_midpoints,
    coef        = ct[bin_vars, "Estimate"],
    se          = ct[bin_vars, "Std. Error"],
    pval        = ct[bin_vars, "Pr(>|t|)"],
    row.names   = NULL
  ) %>%
    mutate(
      ci_lo       = coef - 1.96 * se,
      ci_hi       = coef + 1.96 * se,
      significant = pval < 0.10,
      bin_label   = factor(bin_label, levels = bin_labels)
    )
}

df_bin_coefs <- bind_rows(
  extract_bins(m_boro_bins, "Panel A: Boro Season"),
  extract_bins(m_aus_bins,  "Panel B: Aus Season")
) %>%
  mutate(season = factor(season, levels = c("Panel A: Boro Season", "Panel B: Aus Season")))

# Y range: symmetric around data range
y_lim <- c(
  min(df_bin_coefs$ci_lo, na.rm = TRUE) * 1.15,
  max(df_bin_coefs$ci_hi, na.rm = TRUE) * 1.15
)

fig2 <- ggplot(df_bin_coefs, aes(x = bin_mid, y = coef)) +
  facet_wrap(~season, nrow = 1) +
  # Zone shading
  annotate("rect", xmin = 10, xmax = 30, ymin = -Inf, ymax = Inf,
           fill = "#2E7D32", alpha = 0.08) +
  annotate("rect", xmin = 30, xmax = 42, ymin = -Inf, ymax = Inf,
           fill = "#B71C1C", alpha = 0.08) +
  # Zone boundary lines
  geom_vline(xintercept = c(10, 30), linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  # Reference line
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
             linewidth = 0.5) +
  # Error bars
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.8,
                linewidth = 0.7, color = "grey50") +
  # Points — color by season significance
  geom_point(data = . %>% filter(!significant), size = 3, shape = 21,
             fill = "white", color = "gray50", stroke = 1) +
  geom_point(data = . %>% filter(significant, season == "Panel A: Boro Season"),
             size = 3, shape = 16, color = "#1B5E20") +
  geom_point(data = . %>% filter(significant, season == "Panel B: Aus Season"),
             size = 3, shape = 16, color = "#B71C1C") +
  # Zone text annotations using dedicated label layer
  geom_text(data = data.frame(
    season   = factor(c("Panel A: Boro Season","Panel A: Boro Season",
                        "Panel B: Aus Season","Panel B: Aus Season"),
                      levels = c("Panel A: Boro Season","Panel B: Aus Season")),
    bin_mid  = c(20, 35, 20, 35),
    coef     = y_lim[2] * 0.9,
    label    = c("GDD[10,30]","EDD>30","GDD[10,30]","EDD>30"),
    col_lab  = c("#2E7D32","#B71C1C","#2E7D32","#B71C1C")
  ),
  aes(x = bin_mid, y = coef, label = label, color = I(col_lab)),
  size = 2.9, fontface = "bold", inherit.aes = FALSE) +
  scale_x_continuous(
    name   = "Temperature bin midpoint (°C)",
    breaks = bin_midpoints,
    labels = bin_labels
  ) +
  scale_y_continuous(
    name   = expression("Coefficient on "*Delta*"temp bin (effect on "*Delta*"log yield)"),
    limits = y_lim
  ) +
  labs(title = "Nonparametric Temperature-Yield Response by Season") +
  paper_theme(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8))

ggsave(
  file.path(FIGDIR, "fig_temp_bins_byseason.pdf"),
  fig2, width = 9, height = 5, device = "pdf", dpi = 300, bg = "white"
)
cat("✓ fig_temp_bins_byseason.pdf saved\n")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3: Threshold Sensitivity Coefficient Plot
# ══════════════════════════════════════════════════════════════════════════════
cat("\n=== Figure 3: Threshold Sensitivity ===\n")

# Load data and run threshold regressions (same as rob1_thresholds_byseason.R)
df_raw <- read_csv(
  file.path(ROOT, "data/Regression_data/bangladesh_rice_regression_panel.csv"),
  show_col_types = FALSE
) %>%
  arrange(district, season, year) %>%
  group_by(district, season) %>%
  mutate(
    log_yield      = log(yield_per_ha),
    diff_log_yield = log_yield - lag(log_yield),
    diff_gdd_8_30  = gdd_8_30  - lag(gdd_8_30),
    diff_gdd_10_30 = gdd_10_30 - lag(gdd_10_30),
    diff_gdd_12_30 = gdd_12_30 - lag(gdd_12_30),
    diff_gdd_15_30 = gdd_15_30 - lag(gdd_15_30),
    diff_gdd_10_32 = gdd_10_32 - lag(gdd_10_32),
    diff_gdd_10_35 = gdd_10_35 - lag(gdd_10_35),
    diff_edd_28    = edd_28    - lag(edd_28),
    diff_edd_30    = edd_30    - lag(edd_30),
    diff_edd_32    = edd_32    - lag(edd_32),
    diff_edd_35    = edd_35    - lag(edd_35)
  ) %>%
  ungroup() %>%
  filter(year >= 2013, !is.na(diff_log_yield))

df_boro3 <- df_raw %>% filter(season == "Boro", !is.na(diff_gdd_10_30))
df_aus3  <- df_raw %>% filter(season == "Aus",  !is.na(diff_edd_30))

mb1 <- feols(diff_log_yield ~ diff_gdd_8_30  | year, data = df_boro3, cluster = ~district, warn = FALSE, notes = FALSE)
mb2 <- feols(diff_log_yield ~ diff_gdd_10_30 | year, data = df_boro3, cluster = ~district, warn = FALSE, notes = FALSE)
mb3 <- feols(diff_log_yield ~ diff_gdd_12_30 | year, data = df_boro3, cluster = ~district, warn = FALSE, notes = FALSE)
mb4 <- feols(diff_log_yield ~ diff_gdd_15_30 | year, data = df_boro3, cluster = ~district, warn = FALSE, notes = FALSE)
mb5 <- feols(diff_log_yield ~ diff_gdd_10_32 | year, data = df_boro3, cluster = ~district, warn = FALSE, notes = FALSE)
mb6 <- feols(diff_log_yield ~ diff_gdd_10_35 | year, data = df_boro3, cluster = ~district, warn = FALSE, notes = FALSE)

ma1 <- feols(diff_log_yield ~ diff_edd_28 | year, data = df_aus3, cluster = ~district, warn = FALSE, notes = FALSE)
ma2 <- feols(diff_log_yield ~ diff_edd_30 | year, data = df_aus3, cluster = ~district, warn = FALSE, notes = FALSE)
ma3 <- feols(diff_log_yield ~ diff_edd_32 | year, data = df_aus3, cluster = ~district, warn = FALSE, notes = FALSE)
ma4 <- feols(diff_log_yield ~ diff_edd_35 | year, data = df_aus3, cluster = ~district, warn = FALSE, notes = FALSE)

extract_thresh <- function(m, spec_label, is_main) {
  ct   <- coeftable(m)
  coef <- ct[1, "Estimate"]
  se   <- ct[1, "Std. Error"]
  pval <- ct[1, "Pr(>|t|)"]
  data.frame(
    spec    = spec_label,
    coef    = coef,
    se      = se,
    pval    = pval,
    ci_lo   = coef - 1.96 * se,
    ci_hi   = coef + 1.96 * se,
    is_main = is_main,
    sig     = pval < 0.10
  )
}

boro_labels <- c("GDD[8,30]","GDD[10,30]*","GDD[12,30]","GDD[15,30]","GDD[10,32]","GDD[10,35]")
aus_labels  <- c("EDD>28","EDD>30*","EDD>32","EDD>35")

df_boro_thresh <- bind_rows(
  extract_thresh(mb1, boro_labels[1], FALSE),
  extract_thresh(mb2, boro_labels[2], TRUE),
  extract_thresh(mb3, boro_labels[3], FALSE),
  extract_thresh(mb4, boro_labels[4], FALSE),
  extract_thresh(mb5, boro_labels[5], FALSE),
  extract_thresh(mb6, boro_labels[6], FALSE)
) %>% mutate(spec = factor(spec, levels = boro_labels))

df_aus_thresh <- bind_rows(
  extract_thresh(ma1, aus_labels[1], FALSE),
  extract_thresh(ma2, aus_labels[2], TRUE),
  extract_thresh(ma3, aus_labels[3], FALSE),
  extract_thresh(ma4, aus_labels[4], FALSE)
) %>% mutate(spec = factor(spec, levels = aus_labels))

cat("Threshold regressions complete\n")
cat("Boro results:\n")
print(df_boro_thresh[, c("spec","coef","se","pval")])
cat("Aus results:\n")
print(df_aus_thresh[, c("spec","coef","se","pval")])

# Panel A: Boro
# Y ceiling for annotation placement
y_top_boro <- max(df_boro_thresh$ci_hi) * 1.35
p3a <- ggplot(df_boro_thresh, aes(x = spec, y = coef)) +
  # Separator between lower-bound and upper-bound groups
  geom_vline(xintercept = 4.5, linetype = "dashed", color = "grey60",
             linewidth = 0.4) +
  # Group labels via discrete x geom_text
  geom_text(
    data = data.frame(spec = factor("GDD[12,30]", levels = boro_labels),
                      y = y_top_boro, label = "Lower bound varies"),
    aes(x = spec, y = y, label = label), inherit.aes = FALSE,
    size = 2.8, color = "grey40", nudge_x = -0.5
  ) +
  geom_text(
    data = data.frame(spec = factor("GDD[10,32]", levels = boro_labels),
                      y = y_top_boro, label = "Upper bound varies"),
    aes(x = spec, y = y, label = label), inherit.aes = FALSE,
    size = 2.8, color = "grey40"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
             linewidth = 0.5) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2,
                linewidth = 0.7, color = "grey50") +
  # Insignificant non-main
  geom_point(data = . %>% filter(!is_main, !sig), size = 3, shape = 21,
             fill = "white", color = "gray60", stroke = 1) +
  # Significant non-main
  geom_point(data = . %>% filter(!is_main, sig), size = 3, shape = 16,
             color = "#2E7D32") +
  # Main spec (larger, dark)
  geom_point(data = . %>% filter(is_main), size = 5, shape = 16,
             color = "#1B5E20") +
  geom_label(data = . %>% filter(is_main),
             aes(label = "Main spec"),
             nudge_y = max(df_boro_thresh$ci_hi) * 0.25,
             size = 2.8, color = "#1B5E20", fill = "white",
             linewidth = 0.2, fontface = "bold") +
  scale_y_continuous(labels = label_scientific(digits = 2),
                     expand = expansion(mult = c(0.05, 0.30))) +
  labs(
    title = "Panel A: Boro GDD Threshold Sensitivity",
    x     = NULL,
    y     = "Coefficient estimate"
  ) +
  paper_theme(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# Panel B: Aus
p3b <- ggplot(df_aus_thresh, aes(x = spec, y = coef)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
             linewidth = 0.5) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2,
                linewidth = 0.7, color = "grey50") +
  geom_point(data = . %>% filter(!is_main, !sig), size = 3, shape = 21,
             fill = "white", color = "gray60", stroke = 1) +
  geom_point(data = . %>% filter(!is_main, sig), size = 3, shape = 16,
             color = "#B71C1C") +
  geom_point(data = . %>% filter(is_main), size = 5, shape = 16,
             color = "#7B0000") +
  geom_label(data = . %>% filter(is_main),
             aes(label = "Main spec"),
             nudge_y = min(df_aus_thresh$ci_lo) * 0.25,
             size = 2.8, color = "#7B0000", fill = "white",
             linewidth = 0.2, fontface = "bold") +
  scale_y_continuous(labels = label_scientific(digits = 2)) +
  labs(
    title = "Panel B: Aus EDD Threshold Sensitivity",
    x     = "Threshold specification",
    y     = "Coefficient estimate"
  ) +
  paper_theme(base_size = 10) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

fig3 <- p3a / p3b +
  plot_annotation(
    title = "Stage 1 Threshold Sensitivity: Boro GDD and Aus EDD",
    theme = theme(plot.title = element_text(face = "bold", size = 11))
  )

ggsave(
  file.path(FIGDIR, "fig_threshold_sensitivity.pdf"),
  fig3, width = 7, height = 8, device = "pdf", dpi = 300, bg = "white"
)
cat("✓ fig_threshold_sensitivity.pdf saved\n")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4: Main First-Stage Coefficient Plot (4-panel)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n=== Figure 4: Main First-Stage Coefficient Plot ===\n")

# Load models (run from 01_main_regressions.R) or re-estimate
models_file <- file.path(ROOT, "output/stage1/models/stage1_main_models.RData")
if (file.exists(models_file)) {
  load(models_file)
  cat("Loaded existing stage1 models from RData\n")
  # We have: m_main, m_boro, m_aus, m_aman
} else {
  cat("RData not found — re-estimating from CSV\n")
  df_est4 <- read_csv(
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
    filter(!is.na(diff_log_yield), !is.na(diff_gdd_10_30), !is.na(diff_edd_30))
  m_main <- feols(diff_log_yield ~ diff_gdd_10_30 + diff_edd_30 | year + district^season,
                  data = df_est4, cluster = ~district, warn = FALSE, notes = FALSE)
  m_boro <- feols(diff_log_yield ~ diff_gdd_10_30 | year,
                  data = df_est4 %>% filter(season == "Boro"), cluster = ~district, warn = FALSE, notes = FALSE)
  m_aus  <- feols(diff_log_yield ~ diff_edd_30 | year,
                  data = df_est4 %>% filter(season == "Aus"), cluster = ~district, warn = FALSE, notes = FALSE)
  m_aman <- feols(diff_log_yield ~ diff_gdd_10_30 + diff_edd_30 | year,
                  data = df_est4 %>% filter(season == "Aman"), cluster = ~district, warn = FALSE, notes = FALSE)
}

# Build coefficient data frame
get_coef_row <- function(model, predictor, season_label, pred_type) {
  ct <- coeftable(model)
  if (!predictor %in% rownames(ct)) return(NULL)
  coef <- ct[predictor, "Estimate"]
  se   <- ct[predictor, "Std. Error"]
  pval <- ct[predictor, "Pr(>|t|)"]
  data.frame(
    season    = season_label,
    predictor = pred_type,
    coef      = coef,
    se        = se,
    pval      = pval,
    ci_lo     = coef - 1.96 * se,
    ci_hi     = coef + 1.96 * se,
    sig       = pval < 0.10
  )
}

df_fig4 <- bind_rows(
  get_coef_row(m_main, "diff_gdd_10_30", "Pooled\n(attenuated)", "GDD [10,30]"),
  get_coef_row(m_main, "diff_edd_30",    "Pooled\n(attenuated)", "EDD >30"),
  get_coef_row(m_boro, "diff_gdd_10_30", "Boro\n(GDD, p=0.008)**", "GDD [10,30]"),
  get_coef_row(m_aus,  "diff_edd_30",    "Aus\n(EDD, p=0.051)*",  "EDD >30"),
  get_coef_row(m_aman, "diff_gdd_10_30", "Aman\n(placebo, p>0.7)", "GDD [10,30]"),
  get_coef_row(m_aman, "diff_edd_30",    "Aman\n(placebo, p>0.7)", "EDD >30")
) %>%
  mutate(
    season    = factor(season, levels = c(
      "Pooled\n(attenuated)",
      "Boro\n(GDD, p=0.008)**",
      "Aus\n(EDD, p=0.051)*",
      "Aman\n(placebo, p>0.7)"
    )),
    predictor = factor(predictor, levels = c("GDD [10,30]", "EDD >30")),
    point_shape = ifelse(sig, 16, 21)
  )

cat("Figure 4 data:\n")
print(df_fig4[, c("season","predictor","coef","pval","sig")])

fig4 <- ggplot(df_fig4, aes(x = predictor, y = coef, color = predictor)) +
  facet_wrap(~season, nrow = 1, scales = "free_y") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
             linewidth = 0.5) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.18,
                linewidth = 0.9) +
  # Insignificant: hollow
  geom_point(data = . %>% filter(!sig), size = 3.5, shape = 21,
             aes(fill = predictor), stroke = 1.2) +
  # Significant: filled
  geom_point(data = . %>% filter(sig), size = 3.5, shape = 16) +
  scale_color_manual(
    values  = c("GDD [10,30]" = "#2E7D32", "EDD >30" = "#B71C1C"),
    guide   = "none"
  ) +
  scale_fill_manual(
    values  = c("GDD [10,30]" = "#2E7D32", "EDD >30" = "#B71C1C"),
    guide   = "none"
  ) +
  scale_y_continuous(labels = label_scientific(digits = 2)) +
  labs(
    title = "First-Stage Climate Predictors by Season",
    x     = "Climate predictor",
    y     = expression("Coefficient on "*Delta*"climate predictor")
  ) +
  paper_theme(base_size = 11) +
  theme(
    strip.text   = element_text(face = "bold", size = 8.5),
    axis.text.x  = element_text(size = 8.5)
  )

ggsave(
  file.path(FIGDIR, "fig_stage1_coefplot.pdf"),
  fig4, width = 10, height = 5, device = "pdf", dpi = 300, bg = "white"
)
cat("✓ fig_stage1_coefplot.pdf saved\n")

cat("\n=== All four figures saved to", FIGDIR, "===\n")
