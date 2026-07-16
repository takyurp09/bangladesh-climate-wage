## 06_summary_stats.R
## Stage 1 summary statistics: climate panel + yield hat
## Table 1: Climate & yield by season (GDD/EDD/precip/log_yield/yield_hat)
## Table 2: Climate + yield_hat panel balance
## Output: output/stage1/tables/

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(knitr)
  library(kableExtra)
})

ROOT    <- here::here()
out_tbl <- file.path(ROOT, "output/stage1/tables")
out_sum <- file.path(ROOT, "output/stage1/summary")
dir.create(out_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(out_sum, recursive = TRUE, showWarnings = FALSE)

r2 <- function(x) round(x, 2)

cat("=== STAGE 1 SUMMARY STATISTICS ===\n")

## ── Load data ────────────────────────────────────────────────────────────── ##
clim <- read_csv(file.path(ROOT, "data/Regression_data/bangladesh_rice_regression_panel.csv"),
                 show_col_types = FALSE) %>%
  rename(growing_season = season, District = district)

yield_hat <- read_csv(file.path(ROOT, "output/stage1/fitted/yield_hat_2017_2023.csv"),
                      show_col_types = FALSE) %>%
  rename(District = district, growing_season = season)

cat(sprintf("Climate panel: N=%d rows, years %d-%d\n",
            nrow(clim), min(clim$year), max(clim$year)))
cat(sprintf("Yield hat:     N=%d rows, years %d-%d\n",
            nrow(yield_hat), min(yield_hat$year), max(yield_hat$year)))

## ════════════════════════════════════════════════════════════════════════════ ##
## TABLE 1 — CLIMATE & YIELD SUMMARY BY SEASON (2013–2023)                    ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("\n--- TABLE 1: Climate & yield by season ---\n")

clim2 <- clim %>%
  select(District, growing_season, year, gdd_10_30, edd_30, pr1, log_yield)

clim2_by_season <- clim2 %>%
  group_by(growing_season) %>%
  summarise(
    N_obs         = n(),
    GDD_mean      = r2(mean(gdd_10_30, na.rm = TRUE)),
    GDD_sd        = r2(sd(gdd_10_30,   na.rm = TRUE)),
    GDD_min       = r2(min(gdd_10_30,  na.rm = TRUE)),
    GDD_max       = r2(max(gdd_10_30,  na.rm = TRUE)),
    EDD_mean      = r2(mean(edd_30,    na.rm = TRUE)),
    EDD_sd        = r2(sd(edd_30,      na.rm = TRUE)),
    EDD_min       = r2(min(edd_30,     na.rm = TRUE)),
    EDD_max       = r2(max(edd_30,     na.rm = TRUE)),
    Precip_mean   = r2(mean(pr1,       na.rm = TRUE)),
    Precip_sd     = r2(sd(pr1,         na.rm = TRUE)),
    LogYield_mean = r2(mean(log_yield, na.rm = TRUE)),
    LogYield_sd   = r2(sd(log_yield,   na.rm = TRUE)),
    .groups = "drop"
  )

yhat_by_season <- yield_hat %>%
  group_by(growing_season) %>%
  summarise(
    YieldHat_mean = r2(mean(yield_hat, na.rm = TRUE)),
    YieldHat_sd   = r2(sd(yield_hat,   na.rm = TRUE)),
    .groups = "drop"
  )

tbl1 <- clim2_by_season %>%
  left_join(yhat_by_season, by = "growing_season")

col_labs <- c("Season", "N",
              "Mean", "SD", "Min", "Max",
              "Mean", "SD", "Min", "Max",
              "Mean", "SD",
              "Mean", "SD",
              "Mean (2017-23)", "SD (2017-23)")
names(tbl1) <- col_labs

kt1_html <- kable(tbl1, format = "html",
                  caption = "Climate and Yield Summary by Season (2013-2023)",
                  booktabs = TRUE) %>%
  add_header_above(c(" " = 2,
                     "GDD [10,30] (C-days)" = 4,
                     "EDD >30 (C-days)" = 4,
                     "Precipitation (mm)" = 2,
                     "log(Yield)" = 2,
                     "Yield Hat (Stage 1)" = 2)) %>%
  kable_styling(bootstrap_options = c("striped", "condensed"),
                full_width = FALSE, font_size = 12) %>%
  footnote(general = "Unit: district x season x year. Climate from ERA5 via sinusoidal interpolation (Schlenker-Roberts). Yield hat from Stage 1 season-specific FD regressions (2017-2023 only).",
           general_title = "Note: ", footnote_as_chunk = TRUE)

kt1_tex <- kable(tbl1, format = "latex",
                 caption = "Climate and Yield Summary by Season (2013--2023)",
                 booktabs = TRUE) %>%
  add_header_above(c(" " = 2,
                     "GDD [10,30] ($^\\circ$C-days)" = 4,
                     "EDD $>$30 ($^\\circ$C-days)" = 4,
                     "Precipitation (mm)" = 2,
                     "log(Yield)" = 2,
                     "Yield Hat (Stage 1)" = 2),
                   escape = FALSE) %>%
  kable_styling(latex_options = c("hold_position", "scale_down")) %>%
  footnote(general = "Unit: district $\\times$ season $\\times$ year. Climate from ERA5, sinusoidal interpolation (Schlenker-Roberts). Yield hat from Stage 1 season-specific FD regressions (2017--2023 only).",
           escape = FALSE, general_title = "Note: ", footnote_as_chunk = TRUE)

writeLines(as.character(kt1_html), file.path(out_tbl, "summary_stats_climate.html"))
kableExtra::save_kable(kt1_tex, file.path(out_tbl, "summary_stats_climate.tex"))
cat("Saved summary_stats_climate.tex/.html\n")

## ════════════════════════════════════════════════════════════════════════════ ##
## TABLE 2 — CLIMATE + YIELD PANEL BALANCE                                    ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("\n--- TABLE 2: Panel balance ---\n")

clim_districts <- length(unique(clim$District))
clim_years     <- sort(unique(clim$year))
yhat_districts <- length(unique(yield_hat$District))
yhat_years     <- sort(unique(yield_hat$year))

clim_miss_gdd <- mean(is.na(clim$gdd_10_30)) * 100
clim_miss_edd <- mean(is.na(clim$edd_30))    * 100
clim_miss_pr  <- mean(is.na(clim$pr1))       * 100

tbl2 <- data.frame(
  Statistic = c(
    "Districts in climate panel",
    "Districts with yield hat",
    "Years in climate panel",
    "Years with yield hat",
    "Total obs (climate panel)",
    "Total obs (yield hat)",
    "Missing: GDD [10,30] (%)",
    "Missing: EDD >30 (%)",
    "Missing: Precipitation (%)"
  ),
  Value = c(
    as.character(clim_districts),
    as.character(yhat_districts),
    paste(range(clim_years), collapse = "-"),
    paste(range(yhat_years), collapse = "-"),
    as.character(nrow(clim)),
    as.character(nrow(yield_hat)),
    sprintf("%.1f%%", clim_miss_gdd),
    sprintf("%.1f%%", clim_miss_edd),
    sprintf("%.1f%%", clim_miss_pr)
  ),
  stringsAsFactors = FALSE
)

kt2_html <- kable(tbl2, format = "html",
                  caption = "Climate and Yield Panel Balance",
                  booktabs = TRUE,
                  col.names = c("Statistic", "Value")) %>%
  kable_styling(bootstrap_options = c("striped", "condensed"),
                full_width = FALSE, font_size = 13) %>%
  footnote(general = "Climate panel: ERA5 2013-2023, district x season x year (64 districts). Yield hat: Stage 1 FD regressions, 2017-2023 only.",
           general_title = "Note: ", footnote_as_chunk = TRUE)

kt2_tex <- kable(tbl2, format = "latex",
                 caption = "Climate and Yield Panel Balance",
                 booktabs = TRUE,
                 col.names = c("Statistic", "Value")) %>%
  kable_styling(latex_options = "hold_position") %>%
  footnote(general = "Climate panel: ERA5 2013--2023, district $\\times$ season $\\times$ year. Yield hat: Stage 1 FD regressions, 2017--2023 only.",
           escape = FALSE, general_title = "Note: ", footnote_as_chunk = TRUE)

writeLines(as.character(kt2_html), file.path(out_tbl, "summary_stats_climate_balance.html"))
kableExtra::save_kable(kt2_tex, file.path(out_tbl, "summary_stats_climate_balance.tex"))
cat("Saved summary_stats_climate_balance.tex/.html\n")

## ── Console report ────────────────────────────────────────────────────────── ##
edd_season <- clim %>%
  group_by(growing_season) %>%
  summarise(mean_edd = r2(mean(edd_30,    na.rm = TRUE)),
            mean_gdd = r2(mean(gdd_10_30, na.rm = TRUE)),
            .groups = "drop")

cat("\nMean climate by season:\n")
print(edd_season)
cat(sprintf("\nClimate panel: %d districts x %d years (%d-%d) = %d obs\n",
            clim_districts, length(clim_years),
            min(clim_years), max(clim_years), nrow(clim)))
cat(sprintf("Yield hat: %d districts x %d years (%d-%d)\n",
            yhat_districts, length(yhat_years),
            min(yhat_years), max(yhat_years)))

cat("=== STAGE 1 SUMMARY STATS COMPLETE ===\n")
