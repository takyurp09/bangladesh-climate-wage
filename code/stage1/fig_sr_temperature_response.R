# ============================================================
# fig_sr_temperature_response.R
# Schlenker-Roberts style temperature-yield response figure
# for Boro and Aus rice seasons, Bangladesh 2013-2023
# Creates: output/stage1/figures/fig_sr_temperature_response.pdf
# Does NOT modify any existing script or data file
# ============================================================

# 1. Packages
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(fixest)
  library(arrow)
  library(here)
})

ROOT   <- here::here()
FIGDIR <- file.path(ROOT, "output/stage1/figures")
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

# ── Global theme (matches existing paper figures) ──────────────────────────
paper_theme <- function(base_size = 11) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      plot.title       = element_text(face = "bold", size = base_size, hjust = 0),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      axis.text        = element_text(size = base_size - 1.5, color = "grey30"),
      axis.title       = element_text(size = base_size - 0.5),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}
theme_set(paper_theme())

# ── 2. Load data (read-only) ───────────────────────────────────────────────
bins_file <- file.path(ROOT, "data/Regression_data/bangladesh_rice_temperature_bins.parquet")
if (!file.exists(bins_file)) stop("Not found: ", bins_file)

df_raw <- read_parquet(bins_file) %>%
  mutate(
    district  = as.character(district),
    season    = as.character(season),
    year      = as.numeric(year),
    log_yield = log(yield_per_ha)
  )
cat("Temperature bins loaded: N =", nrow(df_raw), "| 1-deg bins: temp_bin_0_1 to temp_bin_39_40\n")

# ── 3. Construct 3-degree bins (temporary in-memory objects only) ──────────
# 3-degree bins: [0,3), [3,6), ..., [39,42)
# temp_bin_X_Y is the 1-degree bin from X to X+1 (i.e., hours at [X, X+1) °C)
# Each 3-deg bin = sum of 3 consecutive 1-deg bins

make_3deg_bins <- function(df) {
  df %>%
    mutate(
      b_0_3   = temp_bin_0_1  + temp_bin_1_2  + temp_bin_2_3,
      b_3_6   = temp_bin_3_4  + temp_bin_4_5  + temp_bin_5_6,
      b_6_9   = temp_bin_6_7  + temp_bin_7_8  + temp_bin_8_9,
      b_9_12  = temp_bin_9_10 + temp_bin_10_11 + temp_bin_11_12,
      b_12_15 = temp_bin_12_13 + temp_bin_13_14 + temp_bin_14_15,
      b_15_18 = temp_bin_15_16 + temp_bin_16_17 + temp_bin_17_18,
      b_18_21 = temp_bin_18_19 + temp_bin_19_20 + temp_bin_20_21,
      b_21_24 = temp_bin_21_22 + temp_bin_22_23 + temp_bin_23_24,
      # b_24_27 is REFERENCE bin — still compute for exposure/centering
      b_24_27 = temp_bin_24_25 + temp_bin_25_26 + temp_bin_26_27,
      b_27_30 = temp_bin_27_28 + temp_bin_28_29 + temp_bin_29_30,
      b_30_33 = temp_bin_30_31 + temp_bin_31_32 + temp_bin_32_33,
      b_33_36 = temp_bin_33_34 + temp_bin_34_35 + temp_bin_35_36,
      b_36_39 = temp_bin_36_37 + temp_bin_37_38 + temp_bin_38_39,
      b_39_42 = temp_bin_39_40  # only one 1-deg bin at top end
    )
}

df_bins <- make_3deg_bins(df_raw)

# Bin metadata: all 14 bins (including reference)
bin_names <- c("b_0_3","b_3_6","b_6_9","b_9_12","b_12_15","b_15_18",
               "b_18_21","b_21_24","b_24_27","b_27_30","b_30_33",
               "b_33_36","b_36_39","b_39_42")
bin_mid   <- c(1.5, 4.5, 7.5, 10.5, 13.5, 16.5, 19.5, 22.5,
               25.5, 28.5, 31.5, 34.5, 37.5, 40.5)
ref_name  <- "b_24_27"   # reference bin (index 9, midpoint 25.5)
non_ref   <- setdiff(bin_names, ref_name)

# ── 4. Compute exposure (levels) BEFORE first-differencing ─────────────────
# Mean accumulated degree-days per bin per season across all district-year obs
exposure_boro <- df_bins %>%
  filter(season == "Boro") %>%
  summarise(across(all_of(bin_names), mean, na.rm = TRUE)) %>%
  tidyr::pivot_longer(everything(), names_to = "bin", values_to = "exposure") %>%
  mutate(midpoint = bin_mid[match(bin, bin_names)])

exposure_aus <- df_bins %>%
  filter(season == "Aus") %>%
  summarise(across(all_of(bin_names), mean, na.rm = TRUE)) %>%
  tidyr::pivot_longer(everything(), names_to = "bin", values_to = "exposure") %>%
  mutate(midpoint = bin_mid[match(bin, bin_names)])

cat("Exposure computed for Boro and Aus\n")

# ── 5. First-difference within district×season ─────────────────────────────
df_fd <- df_bins %>%
  arrange(district, season, year) %>%
  group_by(district, season) %>%
  mutate(
    diff_log_yield = log_yield - lag(log_yield),
    across(all_of(bin_names), ~ . - lag(.), .names = "diff_{.col}")
  ) %>%
  ungroup() %>%
  filter(year >= 2013, !is.na(diff_log_yield))

df_boro <- df_fd %>% filter(season == "Boro")
df_aus  <- df_fd %>% filter(season == "Aus")
cat("After FD — Boro N:", nrow(df_boro), "| Aus N:", nrow(df_aus), "\n")

# ── 6. Run fine-bin regressions (reference bin = b_24_27 omitted) ──────────
non_ref_fd <- paste0("diff_", non_ref)

reg_formula <- as.formula(
  paste("diff_log_yield ~",
        paste(non_ref_fd, collapse = " + "),
        "| year")
)

m_boro <- feols(reg_formula, data = df_boro, cluster = ~district, warn = FALSE, notes = FALSE)
m_aus  <- feols(reg_formula, data = df_aus,  cluster = ~district, warn = FALSE, notes = FALSE)
cat("Fine-bin regressions complete\n")

# ── 7. Extract coefficients and add reference bin (coef=0, se=0) ───────────
extract_coefs <- function(model, season_label, exposure_df) {
  ct <- coeftable(model)
  df_coef <- data.frame(
    bin      = sub("^diff_", "", rownames(ct)),
    coef     = ct[, "Estimate"],
    se       = ct[, "Std. Error"],
    pval     = ct[, "Pr(>|t|)"],
    row.names = NULL
  ) %>%
    # Add reference bin row
    bind_rows(data.frame(bin = ref_name, coef = 0, se = 0, pval = 1)) %>%
    mutate(
      midpoint  = bin_mid[match(bin, bin_names)],
      ci_lo     = coef - 1.96 * se,
      ci_hi     = coef + 1.96 * se,
      sig       = pval < 0.10,
      season    = season_label
    ) %>%
    arrange(midpoint)

  # Join exposure weights
  df_coef <- df_coef %>%
    left_join(exposure_df %>% select(bin, exposure), by = "bin")

  # ── Exposure-weighted centering ───────────────────────────────────────────
  # Normalize weights to sum to 1
  total_exp <- sum(df_coef$exposure, na.rm = TRUE)
  df_coef <- df_coef %>%
    mutate(
      weight       = exposure / total_exp,
      center_val   = sum(coef * weight, na.rm = TRUE)
    ) %>%
    mutate(
      coef_c  = coef  - center_val,
      ci_lo_c = ci_lo - center_val,
      ci_hi_c = ci_hi - center_val
    )

  cat(season_label, "— center_val:", round(df_coef$center_val[1], 6), "\n")
  df_coef
}

df_boro_coef <- extract_coefs(m_boro, "Boro", exposure_boro)
df_aus_coef  <- extract_coefs(m_aus,  "Aus",  exposure_aus)

# Print key results
cat("\nBoro bin coefficients (centered):\n")
print(df_boro_coef[, c("bin","midpoint","coef_c","se","pval","sig","exposure")])
cat("\nAus bin coefficients (centered):\n")
print(df_aus_coef[, c("bin","midpoint","coef_c","se","pval","sig","exposure")])

# ── 8. Determine shared y-scale for both top panels ─────────────────────────
y_lo <- min(c(df_boro_coef$ci_lo_c, df_aus_coef$ci_lo_c), na.rm = TRUE)
y_hi <- max(c(df_boro_coef$ci_hi_c, df_aus_coef$ci_hi_c), na.rm = TRUE)
y_pad <- (y_hi - y_lo) * 0.12
y_range <- c(y_lo - y_pad, y_hi + y_pad)

# ── Helper: build top response-curve panel ──────────────────────────────────
make_top_panel <- function(df_coef, title_str, point_color) {
  # Filter sparse bins (exposure < 0.5 degree-days) to avoid tail instability
  df_step <- df_coef %>% arrange(midpoint) %>% filter(exposure >= 0.5)

  ggplot(df_step, aes(x = midpoint, y = coef_c)) +
    # Layer 1: background zone shading
    annotate("rect", xmin = 9, xmax = 30, ymin = -Inf, ymax = Inf,
             fill = "#E8F5E9", alpha = 1) +
    annotate("rect", xmin = 30, xmax = 43, ymin = -Inf, ymax = Inf,
             fill = "#FFEBEE", alpha = 1) +
    # Layer 2: CI ribbon — mapped to fill for legend
    geom_ribbon(aes(ymin = ci_lo_c, ymax = ci_hi_c,
                    fill = "95% Confidence Interval"),
                alpha = 0.35) +
    # Layer 3: Step function — mapped to color for legend
    geom_step(aes(color = "Bin Coefficients (Step Function)"),
              linewidth = 1.0) +
    # Layer 4: Loess smooth — mapped to color for legend
    geom_smooth(aes(weight = exposure,
                    color = "Loess Smooth (span = 0.75)"),
                method = "loess", span = 0.75,
                linewidth = 0.9, se = FALSE) +
    # Layer 5: Reference lines
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "gray40", linewidth = 0.5) +
    geom_vline(xintercept = 10, linetype = "dotted",
               color = "#2E7D32", linewidth = 0.6) +
    geom_vline(xintercept = 30, linetype = "dotted",
               color = "#B71C1C", linewidth = 0.6) +
    # Threshold text labels (y = Inf pins them inside the top of the panel)
    annotate("text", x = 10.3, y = Inf, vjust = 1.5,
             label = "10\u00b0C", size = 2.6, color = "#2E7D32",
             hjust = 0, fontface = "bold") +
    annotate("text", x = 30.3, y = Inf, vjust = 1.5,
             label = "30\u00b0C", size = 2.6, color = "#B71C1C",
             hjust = 0, fontface = "bold") +
    # Layer 6: Points (significant filled, insignificant hollow)
    geom_point(data = . %>% filter(!sig, bin != ref_name),
               size = 2.0, shape = 21, fill = "white",
               color = "gray50", stroke = 0.9) +
    geom_point(data = . %>% filter(sig),
               size = 2.0, shape = 16, color = point_color) +
    # Reference bin point (always open gray)
    geom_point(data = . %>% filter(bin == ref_name),
               size = 2.0, shape = 21, fill = "white",
               color = "gray40", stroke = 0.9) +
    # Legend scales
    scale_color_manual(
      name   = NULL,
      values = c(
        "Bin Coefficients (Step Function)" = "steelblue",
        "Loess Smooth (span = 0.75)"       = "black"
      )
    ) +
    scale_fill_manual(
      name   = NULL,
      values = c("95% Confidence Interval" = "gray70"),
      guide  = guide_legend(
        override.aes = list(alpha = 0.35, linetype = "blank", shape = NA)
      )
    ) +
    scale_x_continuous(limits = c(0, 43),
                       breaks  = seq(0, 42, by = 6),
                       labels  = NULL,   # no x labels on top panel
                       expand  = expansion(mult = c(0, 0))) +
    scale_y_continuous(breaks = seq(-0.02, 0.02, by = 0.005)) +
    coord_cartesian(ylim = c(-0.02, 0.02)) +
    labs(
      title = title_str,
      x     = NULL,
      y     = "Log Yield Effect"
    ) +
    paper_theme(base_size = 10) +
    theme(
      axis.text.x       = element_blank(),
      axis.ticks.x      = element_blank(),
      plot.title        = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.position   = c(0.85, 0.85),
      legend.text       = element_text(size = 7),
      legend.key.size   = unit(0.45, "cm"),
      legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
    )
}

# ── Helper: build bottom histogram panel ───────────────────────────────────
make_bot_panel <- function(exposure_df) {
  exposure_df <- exposure_df %>% filter(exposure >= 0.5)
  ggplot(exposure_df, aes(x = midpoint, y = exposure)) +
    geom_col(fill = "#4CAF50", color = "white", width = 2.5, alpha = 0.85) +
    geom_vline(xintercept = 10, linetype = "dotted",
               color = "#2E7D32", linewidth = 0.6) +
    geom_vline(xintercept = 30, linetype = "dotted",
               color = "#B71C1C", linewidth = 0.6) +
    scale_x_continuous(
      limits = c(0, 43),
      breaks  = seq(0, 42, by = 6),
      labels  = c("0","6","12","18","24","30","36","42"),
      expand  = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      name     = "",
      sec.axis = sec_axis(~ ., name = "Exposure (Degree-Days)")
    ) +
    labs(x = "Temperature (Celsius)") +
    paper_theme(base_size = 9) +
    theme(
      axis.title.y.right = element_text(color = "#2E7D32", size = 8,
                                        angle = 90, vjust = 0.5),
      axis.text.y.left   = element_text(size = 7),
      axis.text.y.right  = element_text(size = 7, color = "#2E7D32"),
      axis.text.x        = element_text(size = 8),
      axis.title.x       = element_text(size = 9)
    )
}

# ── 9–12. Build all four panels ─────────────────────────────────────────────
boro_top <- make_top_panel(df_boro_coef,
                           title_str    = "Boro Season (Jan\u2013May)",
                           point_color  = "#1565C0")
boro_bot <- make_bot_panel(exposure_boro)

aus_top  <- make_top_panel(df_aus_coef,
                           title_str    = "Aus Season (Apr\u2013Aug)",
                           point_color  = "#B71C1C")
aus_bot  <- make_bot_panel(exposure_aus)

# ── 13. Assemble with patchwork ─────────────────────────────────────────────
boro_combined <- boro_top / boro_bot +
  plot_layout(heights = c(3, 1))

aus_combined  <- aus_top / aus_bot +
  plot_layout(heights = c(3, 1))

final_figure <- boro_combined | aus_combined

# ── 14. Save ────────────────────────────────────────────────────────────────
out_pdf <- file.path(FIGDIR, "fig_sr_temperature_response.pdf")
out_png <- file.path(FIGDIR, "fig_sr_temperature_response.png")

ggsave(out_pdf, final_figure, width = 10, height = 7,
       device = "pdf", dpi = 300, bg = "white")
cat("Saved:", out_pdf, "\n")

ggsave(out_png, final_figure, width = 10, height = 7,
       dpi = 150, bg = "white")
cat("Saved:", out_png, "\n")

cat("\n=== Done ===\n")
cat("3-degree bins constructed from 1-degree bins (temp_bin_0_1 ... temp_bin_39_40)\n")
cat("Reference bin: b_24_27 (midpoint=25.5 deg C)\n")
cat(sprintf("Boro regression: N=%d district-season-year cells\n", nrow(df_boro)))
cat(sprintf("Aus  regression: N=%d district-season-year cells\n", nrow(df_aus)))
