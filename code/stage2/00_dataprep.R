## stage2_00_dataprep.R
## Build merged panel for stage 2 (LEVELS spec): wage × yield_hat × controls
## Main spec: log_real_wage ~ yield_hat | year + District^growing_season
## Output: data/Regression_data/df_2_merged_levels.csv
##         data/Regression_data/df_2_merged_struct.csv

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
})

ROOT <- here::here()
cat("=== STAGE 2 DATA PREP (LEVELS) ===\n")
cat("Root:", ROOT, "\n")

## ── 1a. Yield hat (LEVELS — fitted log yield) ─────────────────────────────── ##
yhat_files <- list.files(file.path(ROOT, "output/stage1/fitted"),
                          pattern = "yield_hat_(levels_)?2017_[0-9]+\\.csv", full.names = TRUE)
if (length(yhat_files) == 0) {
  yhat_files <- list.files(file.path(ROOT, "output/stage1/fitted"),
                           pattern = "yield_hat_2017_[0-9]+\\.csv", full.names = TRUE)
}
if (length(yhat_files) > 0) {
  yhat_file <- sort(yhat_files, decreasing = TRUE)[1]
} else {
  yhat_file <- file.path(ROOT, "output/stage1/fitted/yield_hat_2017_2023.csv")
}
cat("Yield hat file:", basename(yhat_file), "\n")
yield_raw <- read_csv(yhat_file, show_col_types = FALSE)

df_yield <- yield_raw %>%
  rename(District = district, growing_season = season) %>%
  mutate(District = gsub("Cox.s bazar", "Cox's bazar", District, ignore.case = TRUE)) %>%
  filter(year >= 2017, !is.na(yield_hat)) %>%
  distinct(District, growing_season, year, .keep_all = TRUE) %>%
  transmute(
    District, growing_season, year,
    log_yield_hat = yield_hat,
    # legacy aliases for scripts not yet migrated
    diff_log_yield_hat = yield_hat
  )

cat(sprintf("Yield hat rows (levels): %d  |  districts: %d  |  years: %s\n",
            nrow(df_yield),
            length(unique(df_yield$District)),
            paste(sort(unique(df_yield$year)), collapse = ",")))

## ── 1b. Wage ──────────────────────────────────────────────────────────────── ##
## Use extended wage panel (2024-2025 added); falls back to 2023 if not present
wage_file <- if (file.exists(file.path(ROOT, "data/Regression_data/wage_by_growing_season_2025.csv"))) {
  "data/Regression_data/wage_by_growing_season_2025.csv"
} else {
  "data/Regression_data/wage_by_growing_season.csv"
}
cat("Wage file:", wage_file, "\n")
wage_raw <- read_csv(file.path(ROOT, wage_file), show_col_types = FALSE)

df_wage <- wage_raw %>%
  mutate(District = gsub("Cox.s bazar", "Cox's bazar", District, ignore.case = TRUE)) %>%
  filter(meal_type != "No_info") %>%
  filter(!is.na(real_wage)) %>%
  group_by(District, year, growing_season, gender, meal_type) %>%
  summarise(real_wage = mean(real_wage, na.rm = TRUE), .groups = "drop") %>%
  arrange(District, growing_season, gender, meal_type, year) %>%
  group_by(District, growing_season, gender, meal_type) %>%
  mutate(
    log_real_wage       = log(real_wage),
    diff_real_wage      = real_wage      - lag(real_wage),
    diff_log_real_wage  = log_real_wage  - lag(log_real_wage)
  ) %>%
  ungroup() %>%
  filter(year >= 2017, !is.na(log_real_wage))

# Relevel factors
df_wage <- df_wage %>%
  mutate(
    gender    = relevel(factor(gender),    ref = "Female"),
    meal_type = relevel(factor(meal_type), ref = "None")
  )

cat(sprintf("Wage rows (levels): %d\n", nrow(df_wage)))
cat("Wage N by gender × meal_type:\n")
print(df_wage %>% count(gender, meal_type))

## ── 1c. Merge wage × yield ────────────────────────────────────────────────── ##
df_merged <- df_wage %>%
  left_join(df_yield, by = c("District", "growing_season", "year"))

n_before <- nrow(df_merged)
df_merged <- df_merged %>% na.omit()
n_after   <- nrow(df_merged)

cat(sprintf("Merged rows before na.omit: %d  |  after: %d  |  dropped: %d\n",
            n_before, n_after, n_before - n_after))
cat(sprintf("Merged: %d obs | %d districts | %d years | seasons: %s\n",
            nrow(df_merged),
            length(unique(df_merged$District)),
            length(unique(df_merged$year)),
            paste(sort(unique(df_merged$growing_season)), collapse = ",")))

## District name diagnostics
wage_districts  <- toupper(unique(df_wage$District))
yield_districts <- toupper(unique(df_yield$District))
unmatched <- setdiff(wage_districts, yield_districts)
cat(sprintf("District name mismatches (wage not in yield): %d\n", length(unmatched)))
if (length(unmatched) > 0) cat("  Unmatched:", paste(unmatched, collapse = ", "), "\n")
if (length(unmatched) > 10) stop("STOPPING: >10 district mismatches. Manual review required.")

## ── 1d. Structural data ───────────────────────────────────────────────────── ##
df_struct_2019 <- read_csv(file.path(ROOT, "data/Regression_data/df_struct_2019.csv"),
                           show_col_types = FALSE)

# Log-change 1996→2019 from panel
panel_raw <- read_csv(file.path(ROOT, "data/panel_1996_2008_2019.csv"),
                      show_col_types = FALSE)

df_struct_change <- panel_raw %>%
  filter(Year %in% c(1996, 2019)) %>%
  group_by(District) %>%
  summarise(
    d_log_irrigation = log(Net_Irrigated_Area[Year == 2019] + 1) -
                       log(Net_Irrigated_Area[Year == 1996] + 1),
    d_log_intensity  = log(Intensity_of_Cropping[Year == 2019] + 1) -
                       log(Intensity_of_Cropping[Year == 1996] + 1),
    d_log_holdings   = log(Number_of_Holdings[Year == 2019] + 1) -
                       log(Number_of_Holdings[Year == 1996] + 1),
    d_log_gca        = log(Gross_Cropped_Area[Year == 2019] + 1) -
                       log(Gross_Cropped_Area[Year == 1996] + 1),
    .groups = "drop"
  ) %>%
  mutate(across(starts_with("d_log_"), ~ as.numeric(scale(.)), .names = "z_{.col}"))

# Merge structural into merged panel
df_merged_struct <- df_merged %>%
  left_join(df_struct_2019,    by = "District") %>%
  left_join(df_struct_change,  by = "District")

# Z-score the 2019 levels too
df_merged_struct <- df_merged_struct %>%
  mutate(
    z_irrigation_2019 = as.numeric(scale(irrigation_2019)),
    z_intensity_2019  = as.numeric(scale(intensity_2019)),
    z_holdings_2019   = as.numeric(scale(holdings_2019)),
    z_gca_2019        = as.numeric(scale(gca_2019))
  )

cat(sprintf("Merged+struct rows: %d  |  struct districts matched: %d\n",
            nrow(df_merged_struct),
            sum(!is.na(df_merged_struct$irrigation_2019))))

## ── 1e. Save ──────────────────────────────────────────────────────────────── ##
out_dir <- file.path(ROOT, "data/Regression_data")
write_csv(df_merged,        file.path(out_dir, "df_2_merged_levels.csv"))
write_csv(df_merged_struct, file.path(out_dir, "df_2_merged_struct.csv"))
# Canonical alias used by downstream scripts
write_csv(df_merged,        file.path(out_dir, "df_2_merged_v2.csv"))
cat("Saved df_2_merged_levels.csv and df_2_merged_v2.csv\n")

## ── 1f. Report ────────────────────────────────────────────────────────────── ##
report_dir <- file.path(ROOT, "output/stage2/summary")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

report_lines <- c(
  "# Stage 2 Data Prep Report",
  sprintf("Generated: %s", Sys.time()),
  "",
  "## df_2_merged",
  sprintf("- N rows: %d", nrow(df_merged)),
  sprintf("- Years: %s", paste(sort(unique(df_merged$year)), collapse = ", ")),
  sprintf("- Districts: %d", length(unique(df_merged$District))),
  sprintf("- Seasons: %s", paste(sort(unique(df_merged$growing_season)), collapse = ", ")),
  "",
  "## N per gender x meal_type",
  capture.output(print(df_merged %>% count(gender, meal_type))),
  "",
  "## District mismatches (wage not in yield)",
  if (length(unmatched) == 0) "None" else paste("-", unmatched),
  "",
  "## Year overlap",
  sprintf("Wage years with yield data: %s",
          paste(sort(intersect(unique(df_merged$year), unique(df_yield$year))), collapse = ", ")),
  "",
  "## df_2_merged_struct",
  sprintf("- N rows: %d", nrow(df_merged_struct)),
  sprintf("- Struct districts matched: %d / %d",
          sum(!is.na(df_merged_struct$irrigation_2019)),
          nrow(df_merged_struct))
)

writeLines(report_lines, file.path(report_dir, "data_prep_report.md"))
cat("Saved data_prep_report.md\n")
cat("=== DATA PREP COMPLETE ===\n")

## ── 2. New controls merge (v2) ────────────────────────────────────────────── ##
cat("\n=== STAGE 2 DATA PREP v2: NEW CONTROLS ===\n")

n_base <- nrow(df_merged)
cat(sprintf("Base panel (df_2_merged): %d obs\n", n_base))

## Load new control files
land_ratio   <- read_csv(file.path(ROOT, "data/Regression_data/land_ratio_clean.csv"),
                         show_col_types = FALSE)
land_util    <- read_csv(file.path(ROOT, "data/Regression_data/land_utilization_clean.csv"),
                         show_col_types = FALSE) %>%
                select(District, year, cropping_intensity, Gross_cropped, Net_cropped)
crop_shares  <- read_csv(file.path(ROOT, "data/Regression_data/crop_area_shares_clean.csv"),
                         show_col_types = FALSE)
pop_panel    <- read_csv(file.path(ROOT, "data/Regression_data/population_panel_2013_2023.csv"),
                         show_col_types = FALSE) %>%
                select(District, year, pop_density)
pop_static   <- read_csv(file.path(ROOT, "data/Regression_data/pop_density_district_2020.csv"),
                         show_col_types = FALSE) %>%
                select(District, mean_pop_density_2020)
irrigation   <- if (file.exists(file.path(ROOT, "data/Regression_data/irrigation_panel_2017_2025.csv"))) {
  read_csv(file.path(ROOT, "data/Regression_data/irrigation_panel_2017_2025.csv"),
           show_col_types = FALSE) %>%
    select(District, year, irrigated_total_000_acres, irrigated_boro_000_acres,
           irrigated_aman_000_acres)
} else {
  tibble(District = character(), year = integer(),
         irrigated_total_000_acres = double(),
         irrigated_boro_000_acres = double(),
         irrigated_aman_000_acres = double())
}

## Sequential left joins on District + year (or District only for static)
df_v2 <- df_merged %>%
  left_join(land_ratio,  by = c("District", "year")) %>%
  left_join(land_util,   by = c("District", "year")) %>%
  left_join(crop_shares, by = c("District", "year")) %>%
  left_join(pop_panel,   by = c("District", "year")) %>%
  left_join(irrigation,  by = c("District", "year")) %>%
  left_join(pop_static,  by = "District") %>%
  mutate(
    log_irrigation_total = log(irrigated_total_000_acres + 1),
    log_irrigation_boro  = log(irrigated_boro_000_acres + 1),
    log_irrigation_aman  = log(irrigated_aman_000_acres + 1)
  )

cat(sprintf("After all left joins: %d obs (expected %d)\n", nrow(df_v2), n_base))

## Compute first-differences within District × growing_season × gender × meal_type
df_v2 <- df_v2 %>%
  arrange(District, growing_season, gender, meal_type, year) %>%
  group_by(District, growing_season, gender, meal_type) %>%
  mutate(
    diff_ratio_double       = ratio_double_cropped - lag(ratio_double_cropped),
    diff_cropping_intensity = cropping_intensity   - lag(cropping_intensity),
    diff_pop_density        = pop_density          - lag(pop_density)
  ) %>%
  ungroup()

## Missing-rate report for new variables
new_vars <- c("ratio_double_cropped", "cropping_intensity", "Gross_cropped",
              "Net_cropped", "share_Boro", "share_Aus", "share_Aman",
              "pop_density", "mean_pop_density_2020",
              "irrigated_total_000_acres", "log_irrigation_total",
              "diff_ratio_double", "diff_cropping_intensity", "diff_pop_density")

miss_report <- sapply(new_vars, function(v) {
  if (!v %in% names(df_v2)) return(NA_real_)
  mean(is.na(df_v2[[v]])) * 100
})

cat("\nMissing rate (%) per new variable:\n")
for (v in names(miss_report)) {
  flag <- if (!is.na(miss_report[v]) && miss_report[v] > 20) "  *** FLAG >20% — EXCLUDED FROM MAIN SPEC ***" else ""
  cat(sprintf("  %-30s %5.1f%%%s\n", v, miss_report[v], flag))
}

flagged_vars <- names(miss_report)[!is.na(miss_report) & miss_report > 20]
if (length(flagged_vars) > 0) {
  cat(sprintf("\nFLAGGED (>20%% missing): %s\n", paste(flagged_vars, collapse = ", ")))
} else {
  cat("\nNo variables exceed 20%% missing threshold.\n")
}

## Save v2 panel (canonical levels output)
write_csv(df_v2, file.path(out_dir, "df_2_merged_levels.csv"))
write_csv(df_v2, file.path(out_dir, "df_2_merged_v2.csv"))
cat(sprintf("Saved df_2_merged_levels.csv / df_2_merged_v2.csv  (%d obs, %d cols)\n",
            nrow(df_v2), ncol(df_v2)))

## Save v2 merge report
v2_report <- c(
  "# Stage 2 Panel v2 — New Controls Merge Report",
  sprintf("Generated: %s", Sys.time()),
  "",
  "## Base panel",
  sprintf("- Source: data/Regression_data/df_2_merged.csv"),
  sprintf("- N obs: %d", n_base),
  "",
  "## After merge",
  sprintf("- N obs: %d  (left joins preserve all base rows)", nrow(df_v2)),
  "",
  "## New variables merged",
  "| Variable | Source file | Type | Missing % | Flag |",
  "|----------|-------------|------|-----------|------|",
  paste0("| ratio_double_cropped | land_ratio_clean.csv | level | ",
         sprintf("%.1f%%", miss_report["ratio_double_cropped"]), " | ",
         if ("ratio_double_cropped" %in% flagged_vars) ">20% EXCLUDED" else "OK", " |"),
  paste0("| cropping_intensity | land_utilization_clean.csv | level | ",
         sprintf("%.1f%%", miss_report["cropping_intensity"]), " | ",
         if ("cropping_intensity" %in% flagged_vars) ">20% EXCLUDED" else "OK", " |"),
  paste0("| Gross_cropped | land_utilization_clean.csv | level | ",
         sprintf("%.1f%%", miss_report["Gross_cropped"]), " | ",
         if ("Gross_cropped" %in% flagged_vars) ">20% EXCLUDED" else "OK", " |"),
  paste0("| Net_cropped | land_utilization_clean.csv | level | ",
         sprintf("%.1f%%", miss_report["Net_cropped"]), " | ",
         if ("Net_cropped" %in% flagged_vars) ">20% EXCLUDED" else "OK", " |"),
  paste0("| share_Boro | crop_area_shares_clean.csv | level | ",
         sprintf("%.1f%%", miss_report["share_Boro"]), " | ",
         if ("share_Boro" %in% flagged_vars) ">20% EXCLUDED" else "OK", " |"),
  paste0("| share_Aus | crop_area_shares_clean.csv | level | ",
         sprintf("%.1f%%", miss_report["share_Aus"]), " | ",
         if ("share_Aus" %in% flagged_vars) ">20% EXCLUDED" else "OK", " |"),
  paste0("| share_Aman | crop_area_shares_clean.csv | level | ",
         sprintf("%.1f%%", miss_report["share_Aman"]), " | ",
         if ("share_Aman" %in% flagged_vars) ">20% EXCLUDED" else "OK", " |"),
  paste0("| pop_density | population_panel_2013_2023.csv | level | ",
         sprintf("%.1f%%", miss_report["pop_density"]), " | ",
         if ("pop_density" %in% flagged_vars) ">20% EXCLUDED" else "OK", " |"),
  paste0("| mean_pop_density_2020 | pop_density_district_2020.csv | static | ",
         sprintf("%.1f%%", miss_report["mean_pop_density_2020"]), " | ",
         if ("mean_pop_density_2020" %in% flagged_vars) ">20% EXCLUDED" else "OK", " |"),
  paste0("| diff_ratio_double | computed | FD | ",
         sprintf("%.1f%%", miss_report["diff_ratio_double"]), " | ",
         if ("diff_ratio_double" %in% flagged_vars) ">20% EXCLUDED" else "OK", " |"),
  paste0("| diff_cropping_intensity | computed | FD | ",
         sprintf("%.1f%%", miss_report["diff_cropping_intensity"]), " | ",
         if ("diff_cropping_intensity" %in% flagged_vars) ">20% EXCLUDED" else "OK", " |"),
  paste0("| diff_pop_density | computed | FD | ",
         sprintf("%.1f%%", miss_report["diff_pop_density"]), " | ",
         if ("diff_pop_density" %in% flagged_vars) ">20% EXCLUDED" else "OK", " |"),
  "",
  "## Flagged variables (>20% missing — excluded from main spec)",
  if (length(flagged_vars) == 0) "None" else paste("-", flagged_vars, collapse = "\n"),
  "",
  "## FD variables computed within",
  "  District × growing_season × gender × meal_type",
  "",
  "## Output",
  sprintf("- data/Regression_data/df_2_merged_v2.csv  (%d obs, %d cols)", nrow(df_v2), ncol(df_v2))
)

writeLines(v2_report, file.path(report_dir, "merge_v2_report.md"))
cat("Saved merge_v2_report.md\n")
cat("=== v2 MERGE COMPLETE ===\n")
