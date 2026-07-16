# ============================================================
# ROBUSTNESS 4: TEMPERATURE BINS (Schlenker-Roberts Style)
# ============================================================

library(dplyr)
library(readr)
library(lfe)
library(arrow)
library(stargazer)
library(ggplot2)
library(tidyr)

library(here)

ROOT <- here::here()

cat("\n=== ROBUSTNESS 4: TEMPERATURE BINS ===\n\n")

# Check if temperature bins data exists
bins_file <- file.path(ROOT, "data/Regression_data/bangladesh_rice_temperature_bins.parquet")

if (!file.exists(bins_file)) {
  cat("ERROR: Temperature bins file not found!\n")
  cat("Expected:", bins_file, "\n")
  cat("Please run Script 7 first to generate temperature bin data.\n")
  quit(status = 1)
}

# Load temperature bins data
df_bins <- read_parquet(bins_file)

cat("Temperature bins data loaded\n")
cat("Sample size:", nrow(df_bins), "\n")
cat("Variables:", ncol(df_bins), "\n\n")

# Prepare data
df_bins <- df_bins %>%
  mutate(
    district = as.factor(district),
    season = as.factor(season),
    year = as.numeric(year),
    log_yield = log(yield_per_ha)
  ) %>%
  arrange(district, season, year)

# Create first differences for temperature bins
temp_bin_cols <- grep("^temp_bin_", names(df_bins), value = TRUE)
cat("Found", length(temp_bin_cols), "temperature bin variables\n\n")

# Create coarser bins (5°C intervals for easier interpretation)
df_bins <- df_bins %>%
  mutate(
    temp_0_10 = temp_bin_0_1 + temp_bin_1_2 + temp_bin_2_3 + temp_bin_3_4 + 
                temp_bin_4_5 + temp_bin_5_6 + temp_bin_6_7 + temp_bin_7_8 + 
                temp_bin_8_9 + temp_bin_9_10,
    temp_10_15 = temp_bin_10_11 + temp_bin_11_12 + temp_bin_12_13 + temp_bin_13_14 + temp_bin_14_15,
    temp_15_20 = temp_bin_15_16 + temp_bin_16_17 + temp_bin_17_18 + temp_bin_18_19 + temp_bin_19_20,
    temp_20_25 = temp_bin_20_21 + temp_bin_21_22 + temp_bin_22_23 + temp_bin_23_24 + temp_bin_24_25,
    temp_25_28 = temp_bin_25_26 + temp_bin_26_27 + temp_bin_27_28,  # Optimal range
    temp_28_30 = temp_bin_28_29 + temp_bin_29_30,
    temp_30_32 = temp_bin_30_31 + temp_bin_31_32,
    temp_32_35 = temp_bin_32_33 + temp_bin_33_34 + temp_bin_34_35,
    temp_35_40 = temp_bin_35_36 + temp_bin_36_37 + temp_bin_37_38 + temp_bin_38_39 + temp_bin_39_40
  )

# Create first differences
df_bins <- df_bins %>%
  group_by(district, season) %>%
  mutate(
    diff_log_yield = log_yield - lag(log_yield),
    diff_temp_0_10 = temp_0_10 - lag(temp_0_10),
    diff_temp_10_15 = temp_10_15 - lag(temp_10_15),
    diff_temp_15_20 = temp_15_20 - lag(temp_15_20),
    diff_temp_20_25 = temp_20_25 - lag(temp_20_25),
    diff_temp_25_28 = temp_25_28 - lag(temp_25_28),  # Omit as reference
    diff_temp_28_30 = temp_28_30 - lag(temp_28_30),
    diff_temp_30_32 = temp_30_32 - lag(temp_30_32),
    diff_temp_32_35 = temp_32_35 - lag(temp_32_35),
    diff_temp_35_40 = temp_35_40 - lag(temp_35_40)
  ) %>%
  ungroup() %>%
  filter(!is.na(diff_log_yield)) %>%
  filter(year >= 2013 & year <= 2023) %>%
  mutate(dist_season = interaction(district, season, drop = TRUE))

cat("After differencing, sample size:", nrow(df_bins), "\n\n")

# ============================================================
# REGRESSION: Temperature bins (omit 25-28°C as reference)
# ============================================================

cat("Running temperature bins regression...\n")

m_bins <- felm(diff_log_yield ~ diff_temp_0_10 + diff_temp_10_15 + diff_temp_15_20 + 
                diff_temp_20_25 + diff_temp_28_30 + diff_temp_30_32 + 
                diff_temp_32_35 + diff_temp_35_40 | 
                year + dist_season | 0 | district, 
               data = df_bins)

# Compare with baseline GDD/EDD model
df_bins_full <- read_csv(file.path(ROOT, "data/Regression_data/bangladesh_rice_regression_panel.csv")) %>%
  mutate(
    district = as.factor(district),
    season = as.factor(season),
    year = as.numeric(year),
    log_yield = log(yield_per_ha)
  ) %>%
  arrange(district, season, year) %>%
  group_by(district, season) %>%
  mutate(
    diff_log_yield = log_yield - lag(log_yield),
    diff_gdd_10_30 = gdd_10_30 - lag(gdd_10_30),
    diff_edd_30 = edd_30 - lag(edd_30),
    diff_pr1 = pr1 - lag(pr1),
    diff_pr2 = pr2 - lag(pr2)
  ) %>%
  ungroup() %>%
  filter(!is.na(diff_log_yield)) %>%
  filter(year >= 2013 & year <= 2023) %>%
  mutate(dist_season = interaction(district, season, drop = TRUE))

m_baseline <- felm(diff_log_yield ~ diff_gdd_10_30 + diff_edd_30 + diff_pr1 + diff_pr2 | 
                    year + dist_season | 0 | district, 
                   data = df_bins_full)

# TABLE
stargazer(m_baseline, m_bins,
          type = "latex",
          out = file.path(ROOT, "output/stage1/tables/rob4_temperature_bins.tex"),
          title = "Temperature Bins vs GDD/EDD Specification",
          column.labels = c("Baseline (GDD/EDD)", "Temperature Bins"),
          omit.stat = c("f", "ser"),
          star.cutoffs = c(0.1, 0.05, 0.01))
stargazer(m_baseline, m_bins,
          type = "html",
          out = file.path(ROOT, "output/stage1/tables/rob4_temperature_bins.html"),
          title = "Temperature Bins vs GDD/EDD Specification",
          column.labels = c("Baseline (GDD/EDD)", "Temperature Bins"),
          omit.stat = c("f", "ser"),
          star.cutoffs = c(0.1, 0.05, 0.01))

cat("✓ Table saved\n")

# ============================================================
# PLOT: Temperature bin coefficients
# ============================================================

cat("Creating temperature bins plot...\n")

# Extract coefficients
coef_summary <- summary(m_bins)$coefficients
bin_coefs <- data.frame(
  temp_bin = c("0-10", "10-15", "15-20", "20-25", "25-28", "28-30", "30-32", "32-35", "35-40"),
  temp_center = c(5, 12.5, 17.5, 22.5, 26.5, 29, 31, 33.5, 37.5),
  coef = c(
    coef_summary["diff_temp_0_10", "Estimate"],
    coef_summary["diff_temp_10_15", "Estimate"],
    coef_summary["diff_temp_15_20", "Estimate"],
    coef_summary["diff_temp_20_25", "Estimate"],
    0,  # Reference bin (25-28°C)
    coef_summary["diff_temp_28_30", "Estimate"],
    coef_summary["diff_temp_30_32", "Estimate"],
    coef_summary["diff_temp_32_35", "Estimate"],
    coef_summary["diff_temp_35_40", "Estimate"]
  ),
  se = c(
    coef_summary["diff_temp_0_10", "Cluster s.e."],
    coef_summary["diff_temp_10_15", "Cluster s.e."],
    coef_summary["diff_temp_15_20", "Cluster s.e."],
    coef_summary["diff_temp_20_25", "Cluster s.e."],
    0,
    coef_summary["diff_temp_28_30", "Cluster s.e."],
    coef_summary["diff_temp_30_32", "Cluster s.e."],
    coef_summary["diff_temp_32_35", "Cluster s.e."],
    coef_summary["diff_temp_35_40", "Cluster s.e."]
  )
) %>%
  mutate(
    lower = coef - 1.96 * se,
    upper = coef + 1.96 * se
  )

# Plot
p_bins <- ggplot(bin_coefs, aes(x = temp_center, y = coef)) +
  geom_line(size = 1.2, color = "blue") +
  geom_point(size = 3, color = "blue") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 26.5, linetype = "dotted", color = "green", size = 0.8) +
  annotate("text", x = 26.5, y = max(bin_coefs$upper), 
           label = "Reference\n(25-28°C)", vjust = -0.5, size = 3) +
  labs(
    title = "Temperature-Yield Relationship (Temperature Bins)",
    subtitle = "Coefficient of marginal effect of 1 day exposure to each temperature range",
    x = "Temperature (°C)",
    y = "Effect on Δlog(Yield)",
    caption = "95% confidence intervals shown. Reference category: 25-28°C (omitted)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12)
  )

ggsave(file.path(ROOT, "output/stage1/plots/rob4_temperature_bins.png"),
       p_bins, width = 10, height = 6, dpi = 300)

cat("✓ Plot saved\n\n")

# Print summary
cat("=== TEMPERATURE BIN COEFFICIENTS ===\n")
print(bin_coefs %>% select(temp_bin, coef, se))

cat("\n✓✓✓ TEMPERATURE BINS ROBUSTNESS COMPLETE ✓✓✓\n")
cat("\nAll robustness checks (1-4) are now complete!\n")
cat("Check output/stage1/ for all outputs.\n\n")

