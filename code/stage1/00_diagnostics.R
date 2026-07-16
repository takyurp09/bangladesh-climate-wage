suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(arrow)
  library(here); library(tidyr)
})

ROOT <- here::here()
DIAG <- file.path(ROOT, "output/stage1/diagnostics")
DATA <- file.path(ROOT, "data/Regression_data/bangladesh_rice_regression_panel.csv")
BINS <- file.path(ROOT, "data/Regression_data/bangladesh_rice_temperature_bins.parquet")

# ── Load & compute first differences ────────────────────────────────────────
df <- read_csv(DATA, show_col_types = FALSE) %>%
  arrange(district, season, year) %>%
  group_by(district, season) %>%
  mutate(
    log_yield      = log(yield_per_ha),
    diff_gdd_10_30 = gdd_10_30 - lag(gdd_10_30)
  ) %>%
  ungroup()

cat("\n========== 0. BASIC SHAPE ===========\n")
cat("Rows:", nrow(df), "| Years:", min(df$year), "-", max(df$year),
    "| Districts:", n_distinct(df$district), "| Seasons:", n_distinct(df$season), "\n")
cat("Seasons:", paste(sort(unique(df$season)), collapse=", "), "\n")

# ── CHECK 1: GDD/EDD collinearity ───────────────────────────────────────────
cat("\n========== 1. GDD / EDD COLLINEARITY ===========\n")

corr_df <- df %>%
  group_by(district, season) %>%
  summarise(
    r_levels = cor(gdd_10_30, edd_30, use = "complete.obs"),
    r_fd     = cor(diff_gdd_10_30, diff_edd_30, use = "complete.obs"),
    .groups  = "drop"
  )

cat("Levels GDD~EDD   : median r =", round(median(corr_df$r_levels, na.rm=TRUE), 3),
    " | >0.7:", sum(corr_df$r_levels > 0.7, na.rm=TRUE), "of", nrow(corr_df), "units\n")
cat("FD    ΔGDD~ΔEDD  : median r =", round(median(corr_df$r_fd, na.rm=TRUE), 3),
    " | >0.7:", sum(corr_df$r_fd > 0.7, na.rm=TRUE), "of", nrow(corr_df), "units\n")

cat("\nBy season (FD correlations):\n")
corr_df %>% group_by(season) %>%
  summarise(median_r = round(median(r_fd, na.rm=TRUE), 3),
            pct_gt07 = round(mean(abs(r_fd) > 0.7, na.rm=TRUE), 3)) %>%
  print()

# --- Plot 1a: ΔGDD vs Δlog(yield) ---
fd_df <- df %>% filter(!is.na(diff_gdd_10_30), !is.na(diff_log_yield))
p1a <- ggplot(fd_df, aes(x = diff_gdd_10_30, y = diff_log_yield, color = season)) +
  geom_point(alpha = 0.25, size = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
  facet_wrap(~season, scales = "free_x") +
  labs(title = "DIAG 1a: ΔGDD_10_30 vs Δlog(yield), by season",
       x = "ΔGDD (10–30°C)", y = "Δlog(yield)") +
  theme_minimal() + theme(legend.position = "none")
ggsave(file.path(DIAG, "diag1a_gdd_yield_scatter.png"), p1a, width = 10, height = 4, dpi = 150)
cat("✓ diag1a saved\n")

# --- Plot 1b: ΔEDD vs Δlog(yield) ---
p1b <- ggplot(fd_df, aes(x = diff_edd_30, y = diff_log_yield, color = season)) +
  geom_point(alpha = 0.25, size = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
  facet_wrap(~season, scales = "free_x") +
  labs(title = "DIAG 1b: ΔEDD_30 vs Δlog(yield), by season",
       x = "ΔEDD (>30°C)", y = "Δlog(yield)") +
  theme_minimal() + theme(legend.position = "none")
ggsave(file.path(DIAG, "diag1b_edd_yield_scatter.png"), p1b, width = 10, height = 4, dpi = 150)
cat("✓ diag1b saved\n")

# --- Plot 1c: GDD vs EDD (collinearity check) ---
p1c <- ggplot(df, aes(x = gdd_10_30, y = edd_30, color = season)) +
  geom_point(alpha = 0.15, size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  facet_wrap(~season, scales = "free") +
  labs(title = "DIAG 1c: GDD_10_30 vs EDD_30 (levels) — collinearity",
       x = "GDD (10–30°C)", y = "EDD (>30°C)") +
  theme_minimal() + theme(legend.position = "none")
ggsave(file.path(DIAG, "diag1c_gdd_edd_collinearity.png"), p1c, width = 10, height = 4, dpi = 150)
cat("✓ diag1c saved\n")

# ── CHECK 2: Temperature bin exposure ───────────────────────────────────────
cat("\n========== 2. TEMPERATURE BIN EXPOSURE ===========\n")
if (file.exists(BINS)) {
  bins <- read_parquet(BINS)
  cat("Rows:", nrow(bins), "| Cols:", ncol(bins), "\n")
  bin_cols <- grep("^temp_bin_", names(bins), value = TRUE)
  cat("Bin columns:", length(bin_cols), "\n")

  bin_means <- bins %>%
    group_by(season) %>%
    summarise(across(all_of(bin_cols), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_longer(-season, names_to = "bin", values_to = "mean_days") %>%
    mutate(bin_lo = as.numeric(gsub("temp_bin_(\\d+)_.*", "\\1", bin))) %>%
    arrange(season, bin_lo)

  cat("\nTop 5 bins by mean exposure (all seasons pooled):\n")
  bin_means %>% group_by(season) %>% slice_max(mean_days, n = 5) %>%
    mutate(mean_days = round(mean_days, 2)) %>% print(n = 20)

  cat("\n28–36°C zone:\n")
  bin_means %>% filter(bin_lo >= 28, bin_lo <= 36) %>%
    mutate(mean_days = round(mean_days, 3)) %>% print(n = 30)

  p2 <- ggplot(bin_means, aes(x = bin_lo, y = mean_days, color = season, group = season)) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.5) +
    geom_vline(xintercept = c(10, 30), linetype = "dashed", color = "grey50", linewidth = 0.6) +
    annotate("text", x = 10.3, y = max(bin_means$mean_days) * 0.97, label = "GDD base (10°C)",
             hjust = 0, size = 2.8, color = "grey40") +
    annotate("text", x = 30.3, y = max(bin_means$mean_days) * 0.97, label = "EDD threshold (30°C)",
             hjust = 0, size = 2.8, color = "grey40") +
    labs(title = "DIAG 2: Mean temperature bin exposure by season",
         x = "Temperature bin lower bound (°C)", y = "Mean days", color = "Season") +
    theme_minimal()
  ggsave(file.path(DIAG, "diag2_bin_exposure_profile.png"), p2, width = 9, height = 5, dpi = 150)
  cat("✓ diag2 saved\n")
} else {
  cat("WARNING: parquet file not found at", BINS, "\n")
}

# ── CHECK 3: Time-series break ───────────────────────────────────────────────
cat("\n========== 3. TIME-SERIES BREAK ===========\n")

yr_all <- df %>%
  group_by(year) %>%
  summarise(n = n(), mean_edd = mean(edd_30, na.rm=TRUE),
            mean_gdd = mean(gdd_10_30, na.rm=TRUE),
            mean_yield = mean(log_yield, na.rm=TRUE), .groups="drop")
cat("\nYear stats:\n"); print(yr_all, n = 20)

# Pre/post split
cat("\nPre-2019 vs 2019+:\n")
df %>%
  mutate(period = ifelse(year < 2019, "pre-2019", "2019+")) %>%
  group_by(season, period) %>%
  summarise(mean_edd   = round(mean(edd_30, na.rm=TRUE), 2),
            mean_yield = round(mean(log_yield, na.rm=TRUE), 4),
            n          = n(), .groups = "drop") %>%
  arrange(season, period) %>% print()

p3a <- ggplot(yr_all, aes(x = year, y = mean_edd)) +
  geom_line(color = "#E63946", linewidth = 1) + geom_point(color = "#E63946", size = 2.5) +
  geom_vline(xintercept = 2018.5, linetype = "dashed", color = "grey40") +
  annotate("text", x = 2018.7, y = max(yr_all$mean_edd),
           label = "2018|2019", hjust = 0, size = 3, color = "grey40") +
  labs(title = "DIAG 3a: Mean EDD_30 by year (all seasons)", x = "Year", y = "Mean EDD_30") +
  theme_minimal()
ggsave(file.path(DIAG, "diag3a_edd_by_year.png"), p3a, width = 8, height = 4, dpi = 150)

p3b <- ggplot(yr_all, aes(x = year, y = mean_yield)) +
  geom_line(color = "#457B9D", linewidth = 1) + geom_point(color = "#457B9D", size = 2.5) +
  geom_vline(xintercept = 2018.5, linetype = "dashed", color = "grey40") +
  labs(title = "DIAG 3b: Mean log(yield) by year (all seasons)", x = "Year", y = "Mean log(yield)") +
  theme_minimal()
ggsave(file.path(DIAG, "diag3b_yield_by_year.png"), p3b, width = 8, height = 4, dpi = 150)

p3c <- ggplot(yr_all, aes(x = year, y = n)) +
  geom_col(fill = "#2b8a3e") +
  labs(title = "DIAG 3c: Obs count per year", x = "Year", y = "N") +
  theme_minimal()
ggsave(file.path(DIAG, "diag3c_n_by_year.png"), p3c, width = 8, height = 4, dpi = 150)

yr_seas <- df %>%
  group_by(year, season) %>%
  summarise(mean_edd   = mean(edd_30, na.rm=TRUE),
            mean_yield = mean(log_yield, na.rm=TRUE), .groups="drop")

p3d <- ggplot(yr_seas, aes(x = year, y = mean_edd, color = season, group = season)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  geom_vline(xintercept = 2018.5, linetype = "dashed", color = "grey50") +
  facet_wrap(~season, scales = "free_y") +
  labs(title = "DIAG 3d: Mean EDD_30 by year and season", x = "Year", y = "Mean EDD_30") +
  theme_minimal() + theme(legend.position = "none")
ggsave(file.path(DIAG, "diag3d_edd_by_year_season.png"), p3d, width = 10, height = 4, dpi = 150)

p3e <- ggplot(yr_seas, aes(x = year, y = mean_yield, color = season, group = season)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  geom_vline(xintercept = 2018.5, linetype = "dashed", color = "grey50") +
  facet_wrap(~season, scales = "free_y") +
  labs(title = "DIAG 3e: Mean log(yield) by year and season", x = "Year", y = "Mean log(yield)") +
  theme_minimal() + theme(legend.position = "none")
ggsave(file.path(DIAG, "diag3e_yield_by_year_season.png"), p3e, width = 10, height = 4, dpi = 150)
cat("✓ Plots 3a–3e saved\n")

cat("\n========== DONE ===========\n")
