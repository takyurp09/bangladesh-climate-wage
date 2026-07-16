# ============================================================
# ROBUSTNESS CHECKS: TABLES + PLOTS
# ============================================================

library(dplyr)
library(readr)
library(lfe)
library(stargazer)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(here)

ROOT <- here::here()

# Load data
df <- read_csv(file.path(ROOT, "data/Regression_data/bangladesh_rice_regression_panel.csv"))

# Prepare data
df <- df %>%
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
    diff_gdd_10_32 = gdd_10_32 - lag(gdd_10_32),
    diff_gdd_10_35 = gdd_10_35 - lag(gdd_10_35),
    diff_gdd_8_30 = gdd_8_30 - lag(gdd_8_30),
    diff_gdd_12_30 = gdd_12_30 - lag(gdd_12_30),
    diff_gdd_15_30 = gdd_15_30 - lag(gdd_15_30),
    diff_edd_28 = edd_28 - lag(edd_28),
    diff_edd_30 = edd_30 - lag(edd_30),
    diff_edd_32 = edd_32 - lag(edd_32),
    diff_edd_35 = edd_35 - lag(edd_35),
    diff_pr1 = pr1 - lag(pr1),
    diff_pr2 = pr2 - lag(pr2)
  ) %>%
  ungroup() %>%
  filter(!is.na(diff_log_yield)) %>%
  filter(year >= 2013 & year <= 2023) %>%
  mutate(
    dist_season = interaction(district, season, drop = TRUE),
    diff_gdd_10_30_sq = diff_gdd_10_30^2,
    diff_edd_30_sq = diff_edd_30^2,
    gdd_pr1 = diff_gdd_10_30 * diff_pr1,
    edd_pr1 = diff_edd_30 * diff_pr1,
    gdd_edd = diff_gdd_10_30 * diff_edd_30
  )

cat("Sample size:", nrow(df), "\n")

# ============================================================
# ROBUSTNESS 1: ALTERNATIVE THRESHOLDS
# ============================================================

cat("\n=== ROBUSTNESS 1: ALTERNATIVE THRESHOLDS ===\n")

# Run models
m1_baseline <- felm(diff_log_yield ~ diff_gdd_10_30 + diff_edd_30 + diff_pr1 + diff_pr2 | 
                     year + dist_season | 0 | district, data = df)
m1_cap32 <- felm(diff_log_yield ~ diff_gdd_10_32 + diff_edd_32 + diff_pr1 + diff_pr2 | 
                  year + dist_season | 0 | district, data = df)
m1_cap35 <- felm(diff_log_yield ~ diff_gdd_10_35 + diff_edd_35 + diff_pr1 + diff_pr2 | 
                  year + dist_season | 0 | district, data = df)
m1_base8 <- felm(diff_log_yield ~ diff_gdd_8_30 + diff_edd_30 + diff_pr1 + diff_pr2 | 
                  year + dist_season | 0 | district, data = df)
m1_base12 <- felm(diff_log_yield ~ diff_gdd_12_30 + diff_edd_30 + diff_pr1 + diff_pr2 | 
                   year + dist_season | 0 | district, data = df)
m1_base15 <- felm(diff_log_yield ~ diff_gdd_15_30 + diff_edd_28 + diff_pr1 + diff_pr2 | 
                   year + dist_season | 0 | district, data = df)

# TABLE
stargazer(m1_baseline, m1_cap32, m1_cap35, m1_base8, m1_base12, m1_base15,
          type = "latex",
          out = file.path(ROOT, "output/stage1/tables/rob1_thresholds.tex"),
          title = "Alternative Temperature Thresholds",
          column.labels = c("10/30", "10/32", "10/35", "8/30", "12/30", "15/30"),
          omit.stat = c("f", "ser"),
          star.cutoffs = c(0.1, 0.05, 0.01))
stargazer(m1_baseline, m1_cap32, m1_cap35, m1_base8, m1_base12, m1_base15,
          type = "html",
          out = file.path(ROOT, "output/stage1/tables/rob1_thresholds.html"),
          title = "Alternative Temperature Thresholds",
          column.labels = c("10/30", "10/32", "10/35", "8/30", "12/30", "15/30"),
          omit.stat = c("f", "ser"),
          star.cutoffs = c(0.1, 0.05, 0.01))

# PLOT: Coefficient comparison
extract_coef <- function(model, var_pattern) {
  s <- summary(model)
  cf <- s$coefficients
  idx <- grep(var_pattern, rownames(cf))[1]
  if(length(idx) == 0) return(c(NA, NA, NA))
  c(cf[idx, "Estimate"], 
    cf[idx, "Estimate"] - 1.96*cf[idx, "Cluster s.e."],
    cf[idx, "Estimate"] + 1.96*cf[idx, "Cluster s.e."])
}

coef_data <- data.frame(
  model = c("Base10/Cap30", "Base10/Cap32", "Base10/Cap35", 
            "Base8/Cap30", "Base12/Cap30", "Base15/Cap30"),
  gdd_coef = c(extract_coef(m1_baseline, "gdd")[1],
               extract_coef(m1_cap32, "gdd")[1],
               extract_coef(m1_cap35, "gdd")[1],
               extract_coef(m1_base8, "gdd")[1],
               extract_coef(m1_base12, "gdd")[1],
               extract_coef(m1_base15, "gdd")[1]),
  gdd_lower = c(extract_coef(m1_baseline, "gdd")[2],
                extract_coef(m1_cap32, "gdd")[2],
                extract_coef(m1_cap35, "gdd")[2],
                extract_coef(m1_base8, "gdd")[2],
                extract_coef(m1_base12, "gdd")[2],
                extract_coef(m1_base15, "gdd")[2]),
  gdd_upper = c(extract_coef(m1_baseline, "gdd")[3],
                extract_coef(m1_cap32, "gdd")[3],
                extract_coef(m1_cap35, "gdd")[3],
                extract_coef(m1_base8, "gdd")[3],
                extract_coef(m1_base12, "gdd")[3],
                extract_coef(m1_base15, "gdd")[3]),
  edd_coef = c(extract_coef(m1_baseline, "edd")[1],
               extract_coef(m1_cap32, "edd")[1],
               extract_coef(m1_cap35, "edd")[1],
               extract_coef(m1_base8, "edd")[1],
               extract_coef(m1_base12, "edd")[1],
               extract_coef(m1_base15, "edd")[1]),
  edd_lower = c(extract_coef(m1_baseline, "edd")[2],
                extract_coef(m1_cap32, "edd")[2],
                extract_coef(m1_cap35, "edd")[2],
                extract_coef(m1_base8, "edd")[2],
                extract_coef(m1_base12, "edd")[2],
                extract_coef(m1_base15, "edd")[2]),
  edd_upper = c(extract_coef(m1_baseline, "edd")[3],
                extract_coef(m1_cap32, "edd")[3],
                extract_coef(m1_cap35, "edd")[3],
                extract_coef(m1_base8, "edd")[3],
                extract_coef(m1_base12, "edd")[3],
                extract_coef(m1_base15, "edd")[3])
)

p1 <- ggplot(coef_data, aes(x = model, y = gdd_coef)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = gdd_lower, ymax = gdd_upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(title = "GDD Coefficients Across Thresholds",
       x = "Specification", y = "Coefficient (95% CI)") +
  theme_minimal()

p2 <- ggplot(coef_data, aes(x = model, y = edd_coef)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = edd_lower, ymax = edd_upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(title = "EDD Coefficients Across Thresholds",
       x = "Specification", y = "Coefficient (95% CI)") +
  theme_minimal()

png(file.path(ROOT, "output/stage1/plots/rob1_thresholds.png"), width = 10, height = 6, units = "in", res = 300)
grid.arrange(p1, p2, ncol = 2)
dev.off()

cat("✓ Robustness 1 saved\n")

# ============================================================

# ============================================================
# ROBUSTNESS 3: SAMPLE SPLITS
# ============================================================

cat("\n=== ROBUSTNESS 3: SAMPLE SPLITS ===\n")

# By season
m3_boro <- felm(diff_log_yield ~ diff_gdd_10_30 + diff_edd_30 + diff_pr1 + diff_pr2 | 
                 year + district | 0 | district, data = df %>% filter(season == "Boro"))
m3_aus <- felm(diff_log_yield ~ diff_gdd_10_30 + diff_edd_30 + diff_pr1 + diff_pr2 | 
                year + district | 0 | district, data = df %>% filter(season == "Aus"))
m3_aman <- felm(diff_log_yield ~ diff_gdd_10_30 + diff_edd_30 + diff_pr1 + diff_pr2 | 
                 year + district | 0 | district, data = df %>% filter(season == "Aman"))

stargazer(m3_boro, m3_aus, m3_aman,
          type = "latex",
          out = file.path(ROOT, "output/stage1/tables/rob3_by_season.tex"),
          title = "Heterogeneous Effects by Season",
          column.labels = c("Boro", "Aus", "Aman"),
          omit.stat = c("f", "ser"),
          star.cutoffs = c(0.1, 0.05, 0.01))
stargazer(m3_boro, m3_aus, m3_aman,
          type = "html",
          out = file.path(ROOT, "output/stage1/tables/rob3_by_season.html"),
          title = "Heterogeneous Effects by Season",
          column.labels = c("Boro", "Aus", "Aman"),
          omit.stat = c("f", "ser"),
          star.cutoffs = c(0.1, 0.05, 0.01))

# PLOT: Coefficients by season
season_coefs <- data.frame(
  season = c("Boro", "Aus", "Aman"),
  gdd = c(extract_coef(m3_boro, "gdd")[1],
          extract_coef(m3_aus, "gdd")[1],
          extract_coef(m3_aman, "gdd")[1]),
  gdd_lower = c(extract_coef(m3_boro, "gdd")[2],
                extract_coef(m3_aus, "gdd")[2],
                extract_coef(m3_aman, "gdd")[2]),
  gdd_upper = c(extract_coef(m3_boro, "gdd")[3],
                extract_coef(m3_aus, "gdd")[3],
                extract_coef(m3_aman, "gdd")[3]),
  edd = c(extract_coef(m3_boro, "edd")[1],
          extract_coef(m3_aus, "edd")[1],
          extract_coef(m3_aman, "edd")[1]),
  edd_lower = c(extract_coef(m3_boro, "edd")[2],
                extract_coef(m3_aus, "edd")[2],
                extract_coef(m3_aman, "edd")[2]),
  edd_upper = c(extract_coef(m3_boro, "edd")[3],
                extract_coef(m3_aus, "edd")[3],
                extract_coef(m3_aman, "edd")[3])
)

p4 <- ggplot(season_coefs, aes(x = season, y = edd, color = season)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = edd_lower, ymax = edd_upper), width = 0.2, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Boro" = "#2E86AB", "Aus" = "#A23B72", "Aman" = "#F18F01")) +
  labs(title = "Heat Stress Effects by Season",
       subtitle = "EDD coefficients with 95% CI",
       x = "Season", y = "ΔEDD Coefficient") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(ROOT, "output/stage1/plots/rob3_by_season.png"), p4, width = 8, height = 6, dpi = 300)

# By time period
m3_early <- felm(diff_log_yield ~ diff_gdd_10_30 + diff_edd_30 + diff_pr1 + diff_pr2 | 
                  year + dist_season | 0 | district, 
                  data = df %>% filter(year <= 2017))
m3_late <- felm(diff_log_yield ~ diff_gdd_10_30 + diff_edd_30 + diff_pr1 + diff_pr2 | 
                 year + dist_season | 0 | district, 
                 data = df %>% filter(year >= 2018))

stargazer(m3_early, m3_late,
          type = "latex",
          out = file.path(ROOT, "output/stage1/tables/rob3_by_time.tex"),
          title = "Stability Over Time",
          column.labels = c("2013-2017", "2018-2023"),
          omit.stat = c("f", "ser"),
          star.cutoffs = c(0.1, 0.05, 0.01))
stargazer(m3_early, m3_late,
          type = "html",
          out = file.path(ROOT, "output/stage1/tables/rob3_by_time.html"),
          title = "Stability Over Time",
          column.labels = c("2013-2017", "2018-2023"),
          omit.stat = c("f", "ser"),
          star.cutoffs = c(0.1, 0.05, 0.01))

cat("✓ Robustness 3 saved\n")

cat("\n✓✓✓ ALL ROBUSTNESS CHECKS COMPLETE ✓✓✓\n")
cat("\nOutputs saved to:\n")
cat("  - output/stage1/tables/\n")
cat("  - output/stage1/plots/\n\n")

