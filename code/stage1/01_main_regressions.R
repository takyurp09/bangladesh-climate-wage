# stage1_main.R
# First-stage: climate effects on rice yield (Bangladesh, 2013-2025)
# MAIN SPEC: levels (log_yield ~ GDD/EDD | year + district), season-specific
# APPENDIX: first-differences (legacy comparison)

suppressPackageStartupMessages({
  library(fixest)
  library(modelsummary)
  library(dplyr)
  library(readr)
  library(here)
})
if (!requireNamespace("kableExtra", quietly = TRUE))
  warning("kableExtra not installed — tables will still export via modelsummary.")

ROOT <- here::here()
for (d in c("output/stage1/tables", "output/stage1/models", "output/stage1/fitted"))
  dir.create(file.path(ROOT, d), recursive = TRUE, showWarnings = FALSE)

# ── Data ──────────────────────────────────────────────────────────────────────
df_raw <- read_csv(
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
  mutate(
    precip         = pr1,
    precip_sq      = pr1^2,
    diff_precip    = diff_pr1,
    diff_precip_sq = diff_pr1^2
  )

df_est <- df_raw %>%
  filter(!is.na(log_yield), !is.na(gdd_10_30), !is.na(edd_30))

df_fd <- df_raw %>%
  filter(!is.na(diff_log_yield), !is.na(diff_gdd_10_30), !is.na(diff_edd_30))

cat("Levels estimation N:", nrow(df_est), "\n")
cat("FD estimation N:", nrow(df_fd), "\n")

# ── MAIN: Levels, season-specific ─────────────────────────────────────────────
m_lv_boro <- feols(
  log_yield ~ gdd_10_30 | year + district,
  data = df_est %>% filter(season == "Boro"), cluster = ~district
)
m_lv_aus <- feols(
  log_yield ~ edd_30 | year + district,
  data = df_est %>% filter(season == "Aus"), cluster = ~district
)
m_lv_aman <- feols(
  log_yield ~ gdd_10_30 + edd_30 | year + district,
  data = df_est %>% filter(season == "Aman"), cluster = ~district
)

# Pooled levels (appendix comparison)
m_lv_pooled <- feols(
  log_yield ~ gdd_10_30 + edd_30 + precip + precip_sq |
    year + district^season,
  data = df_est, cluster = ~district
)

# ── APPENDIX: FD season-specific ──────────────────────────────────────────────
m_fd_boro <- feols(
  diff_log_yield ~ diff_gdd_10_30 | year,
  data = df_fd %>% filter(season == "Boro"), cluster = ~district
)
m_fd_aus <- feols(
  diff_log_yield ~ diff_edd_30 | year,
  data = df_fd %>% filter(season == "Aus"), cluster = ~district
)
m_fd_aman <- feols(
  diff_log_yield ~ diff_gdd_10_30 + diff_edd_30 | year,
  data = df_fd %>% filter(season == "Aman"), cluster = ~district
)

m_fd_pooled <- feols(
  diff_log_yield ~ diff_gdd_10_30 + diff_edd_30 + diff_precip + diff_precip_sq |
    year + district^season,
  data = df_fd, cluster = ~district
)

# ── Table helpers ─────────────────────────────────────────────────────────────
stars <- c("*" = 0.1, "**" = 0.05, "***" = 0.01)

cm_lv <- c(
  "gdd_10_30" = "GDD (10--30$^\\circ$C)",
  "edd_30"    = "EDD ($>$30$^\\circ$C)",
  "precip"    = "Precipitation",
  "precip_sq" = "Precipitation$^2$"
)
cm_lv_html <- c(
  "gdd_10_30" = "GDD (10–30°C)",
  "edd_30"    = "EDD (>30°C)",
  "precip"    = "Precipitation",
  "precip_sq" = "Precipitation²"
)
cm_fd <- c(
  "diff_gdd_10_30" = "$\\Delta$GDD (10--30$^\\circ$C)",
  "diff_edd_30"    = "$\\Delta$EDD ($>$30$^\\circ$C)",
  "diff_precip"    = "$\\Delta$Precipitation",
  "diff_precip_sq" = "$\\Delta$Precipitation$^2$"
)
cm_fd_html <- c(
  "diff_gdd_10_30" = "ΔGDD (10–30°C)",
  "diff_edd_30"    = "ΔEDD (>30°C)",
  "diff_precip"    = "ΔPrecipitation",
  "diff_precip_sq" = "ΔPrecipitation²"
)

gof_sel <- data.frame(
  raw = c("nobs", "r.squared", "within.r.squared"),
  clean = c("Observations", "$R^2$", "Within $R^2$"),
  fmt = c(0, 3, 3)
)
gof_sel_html <- data.frame(
  raw = c("nobs", "r.squared", "within.r.squared"),
  clean = c("Observations", "R²", "Within R²"),
  fmt = c(0, 3, 3)
)

save_table <- function(models, cm, cm_h, gof, gof_h, title, title_h, note, base_path) {
  modelsummary(models, coef_map = cm, stars = stars, gof_map = gof,
               title = title, notes = note,
               output = paste0(base_path, ".tex"))
  modelsummary(models, coef_map = cm_h, stars = stars, gof_map = gof_h,
               title = title_h, notes = gsub("\\$|\\\\|\\{|\\}", "", note),
               output = paste0(base_path, ".html"))
}

note_lv <- "Levels spec. Clustered SE by district. FE: year + district (season-specific)."
note_fd <- "FD spec (appendix). Clustered SE by district. Year FE."

save_table(
  list("Pooled" = m_lv_pooled, "Boro (GDD)" = m_lv_boro,
       "Aus (EDD)" = m_lv_aus, "Aman (GDD+EDD)" = m_lv_aman),
  cm_lv, cm_lv_html, gof_sel, gof_sel_html,
  "First Stage: Climate Effects on Rice Yield (Levels, 2013--2025)",
  "First Stage: Climate Effects on Rice Yield (Levels, 2013–2025)",
  note_lv,
  file.path(ROOT, "output/stage1/tables/main_table")
)
cat("✓ main_table.tex + .html (levels)\n")

save_table(
  list("Pooled" = m_fd_pooled, "Boro (GDD)" = m_fd_boro,
       "Aus (EDD)" = m_fd_aus, "Aman (GDD+EDD)" = m_fd_aman),
  cm_fd, cm_fd_html, gof_sel, gof_sel_html,
  "Appendix: Climate--Yield Relationship (First Differences, 2013--2025)",
  "Appendix: Climate-Yield Relationship (First Differences, 2013–2025)",
  note_fd,
  file.path(ROOT, "output/stage1/tables/appendix_fd")
)
cat("✓ appendix_fd.tex + .html\n")

# ── Save models ───────────────────────────────────────────────────────────────
save(m_lv_boro, m_lv_aus, m_lv_aman, m_lv_pooled,
     m_fd_boro, m_fd_aus, m_fd_aman, m_fd_pooled,
     file = file.path(ROOT, "output/stage1/models/stage1_main_models.RData"))
cat("✓ models RData\n")

# ── Fitted values (levels, 2017+ for Stage 2) ───────────────────────────────
df_est$yield_hat <- NA_real_
season_models <- list(Boro = m_lv_boro, Aus = m_lv_aus, Aman = m_lv_aman)
for (season_name in names(season_models)) {
  model <- season_models[[season_name]]
  sub_idx <- which(df_est$season == season_name)
  global_idx <- sub_idx[obs(model)]
  df_est$yield_hat[global_idx] <- fitted(model)
}
df_est$residual <- df_est$log_yield - df_est$yield_hat

df_fitted <- df_est %>%
  filter(year >= 2017, !is.na(yield_hat)) %>%
  mutate(district = gsub("Cox.s bazar", "Cox's bazar", district, ignore.case = TRUE)) %>%
  distinct(district, season, year, .keep_all = TRUE) %>%
  select(district, season, year, log_yield, yield_hat, residual)

max_year <- max(df_fitted$year, na.rm = TRUE)
out_name <- sprintf("yield_hat_levels_2017_%d.csv", max_year)
write.csv(df_fitted,
          file.path(ROOT, "output/stage1/fitted", out_name),
          row.names = FALSE)
# Canonical filenames for Stage 2
write.csv(df_fitted,
          file.path(ROOT, "output/stage1/fitted/yield_hat_2017_2023.csv"),
          row.names = FALSE)
write.csv(df_fitted,
          file.path(ROOT, "output/stage1/fitted/yield_hat_2017_2025.csv"),
          row.names = FALSE)
cat("✓", out_name, "(N =", nrow(df_fitted), ")\n")

cat("\n=== STAGE1_MAIN COMPLETE (LEVELS) ===\n")
cat("Fitted N:", nrow(df_fitted), "| Years:", min(df_fitted$year), "-", max_year, "\n")
