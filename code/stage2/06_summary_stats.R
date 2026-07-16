## 06_summary_stats.R
## Stage 2 summary statistics: wage panel only
## Table 1: Wage panel by gender x meal_type x season
## Table 2: Wage panel balance
## Output: output/stage2/tables/summary_stats_*.tex/.html

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(knitr)
  library(kableExtra)
})

ROOT    <- here::here()
out_tbl <- file.path(ROOT, "output/stage2/tables")
out_sum <- file.path(ROOT, "output/stage2/summary")
dir.create(out_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(out_sum, recursive = TRUE, showWarnings = FALSE)

## ── helper: save tex + html ───────────────────────────────────────────────── ##
save_kable <- function(kt, stem, caption, footnote_txt) {
  # HTML
  kt_html <- kt %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                  full_width = FALSE, font_size = 13) %>%
    footnote(general = footnote_txt, general_title = "Note: ",
             footnote_as_chunk = TRUE)
  save_kable_html <- function(x, file) {
    writeLines(as.character(x), file)
  }
  save_kable_html(kt_html, file.path(out_tbl, paste0(stem, ".html")))

  # LaTeX
  kt_tex <- kt %>%
    kable_styling(latex_options = c("hold_position", "striped")) %>%
    footnote(general = footnote_txt, general_title = "Note: ",
             footnote_as_chunk = TRUE, escape = FALSE)
  save_kable_tex <- kableExtra::save_kable
  save_kable_tex(kt_tex, file.path(out_tbl, paste0(stem, ".tex")))

  cat(sprintf("Saved %s.tex/.html\n", stem))
}

## rounding helper
r2 <- function(x) round(x, 2)

cat("=== SUMMARY STATISTICS ===\n")

## ════════════════════════════════════════════════════════════════════════════ ##
## TABLE 1 — WAGE PANEL SUMMARY                                                ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("\n--- TABLE 1: Wage panel ---\n")

wage <- read_csv(file.path(ROOT, "data/Regression_data/df_2_merged_v2.csv"),
                 show_col_types = FALSE)

cat(sprintf("Wage panel: N=%d rows, years %s-%s\n",
            nrow(wage), min(wage$year), max(wage$year)))

## 1a: by gender
by_gender <- wage %>%
  group_by(gender) %>%
  summarise(
    N              = n(),
    real_wage_mean = r2(mean(real_wage,      na.rm = TRUE)),
    real_wage_sd   = r2(sd(real_wage,        na.rm = TRUE)),
    real_wage_min  = r2(min(real_wage,       na.rm = TRUE)),
    real_wage_max  = r2(max(real_wage,       na.rm = TRUE)),
    diff_wage_mean = r2(mean(diff_real_wage, na.rm = TRUE)),
    diff_wage_sd   = r2(sd(diff_real_wage,   na.rm = TRUE)),
    diff_lyhat_mean= r2(mean(diff_log_yield_hat, na.rm = TRUE)),
    diff_lyhat_sd  = r2(sd(diff_log_yield_hat,   na.rm = TRUE)),
    .groups = "drop"
  )

## 1b: by meal_type
by_meal <- wage %>%
  group_by(meal_type) %>%
  summarise(
    N              = n(),
    real_wage_mean = r2(mean(real_wage,      na.rm = TRUE)),
    real_wage_sd   = r2(sd(real_wage,        na.rm = TRUE)),
    real_wage_min  = r2(min(real_wage,       na.rm = TRUE)),
    real_wage_max  = r2(max(real_wage,       na.rm = TRUE)),
    diff_wage_mean = r2(mean(diff_real_wage, na.rm = TRUE)),
    diff_wage_sd   = r2(sd(diff_real_wage,   na.rm = TRUE)),
    diff_lyhat_mean= r2(mean(diff_log_yield_hat, na.rm = TRUE)),
    diff_lyhat_sd  = r2(sd(diff_log_yield_hat,   na.rm = TRUE)),
    .groups = "drop"
  )

## 1c: by season
by_season <- wage %>%
  group_by(growing_season) %>%
  summarise(
    N              = n(),
    real_wage_mean = r2(mean(real_wage,      na.rm = TRUE)),
    real_wage_sd   = r2(sd(real_wage,        na.rm = TRUE)),
    real_wage_min  = r2(min(real_wage,       na.rm = TRUE)),
    real_wage_max  = r2(max(real_wage,       na.rm = TRUE)),
    diff_wage_mean = r2(mean(diff_real_wage, na.rm = TRUE)),
    diff_wage_sd   = r2(sd(diff_real_wage,   na.rm = TRUE)),
    diff_lyhat_mean= r2(mean(diff_log_yield_hat, na.rm = TRUE)),
    diff_lyhat_sd  = r2(sd(diff_log_yield_hat,   na.rm = TRUE)),
    .groups = "drop"
  )

## 1d: overall
overall <- wage %>%
  summarise(
    gender         = "All",
    N              = n(),
    real_wage_mean = r2(mean(real_wage,      na.rm = TRUE)),
    real_wage_sd   = r2(sd(real_wage,        na.rm = TRUE)),
    real_wage_min  = r2(min(real_wage,       na.rm = TRUE)),
    real_wage_max  = r2(max(real_wage,       na.rm = TRUE)),
    diff_wage_mean = r2(mean(diff_real_wage, na.rm = TRUE)),
    diff_wage_sd   = r2(sd(diff_real_wage,   na.rm = TRUE)),
    diff_lyhat_mean= r2(mean(diff_log_yield_hat, na.rm = TRUE)),
    diff_lyhat_sd  = r2(sd(diff_log_yield_hat,   na.rm = TRUE))
  )

## Combine into one display table: Panel A=overall+gender, B=meal, C=season
col_labs <- c("Group", "N", "Mean", "SD", "Min", "Max",
              "Mean Δ", "SD Δ", "Mean ΔlnYhat", "SD ΔlnYhat")

panel_A <- bind_rows(
  overall %>% rename(group = gender),
  by_gender %>% rename(group = gender)
) %>% mutate(group = as.character(group))

panel_B <- by_meal %>% rename(group = meal_type) %>%
  mutate(group = paste0("  Meal: ", group))

panel_C <- by_season %>% rename(group = growing_season) %>%
  mutate(group = paste0("  ", group))

tbl1_combined <- bind_rows(
  data.frame(group = "Panel A: By gender", N = NA_integer_,
             real_wage_mean=NA, real_wage_sd=NA, real_wage_min=NA, real_wage_max=NA,
             diff_wage_mean=NA, diff_wage_sd=NA, diff_lyhat_mean=NA, diff_lyhat_sd=NA),
  panel_A,
  data.frame(group = "Panel B: By meal provision", N = NA_integer_,
             real_wage_mean=NA, real_wage_sd=NA, real_wage_min=NA, real_wage_max=NA,
             diff_wage_mean=NA, diff_wage_sd=NA, diff_lyhat_mean=NA, diff_lyhat_sd=NA),
  panel_B,
  data.frame(group = "Panel C: By growing season", N = NA_integer_,
             real_wage_mean=NA, real_wage_sd=NA, real_wage_min=NA, real_wage_max=NA,
             diff_wage_mean=NA, diff_wage_sd=NA, diff_lyhat_mean=NA, diff_lyhat_sd=NA),
  panel_C
)

names(tbl1_combined) <- col_labs

kt1 <- kable(tbl1_combined, format = "html",
             caption = "Table 1: Wage Panel Summary Statistics (2017–2023)",
             booktabs = TRUE, linesep = "") %>%
  add_header_above(c(" " = 2,
                     "Real Wage (BDT/day)" = 4,
                     "Δ Real Wage" = 2,
                     "Δ log(Yield Hat)" = 2))

kt1_tex <- kable(tbl1_combined, format = "latex",
                 caption = "Wage Panel Summary Statistics (2017--2023)",
                 booktabs = TRUE, linesep = "") %>%
  add_header_above(c(" " = 2,
                     "Real Wage (BDT/day)" = 4,
                     "$\\\\Delta$ Real Wage" = 2,
                     "$\\\\Delta$ log(Yield Hat)" = 2),
                   escape = FALSE) %>%
  kable_styling(latex_options = c("hold_position", "scale_down")) %>%
  footnote(general = "Unit: district $\\\\times$ season $\\\\times$ gender $\\\\times$ meal type $\\\\times$ year. Real wages deflated by district CPI. Yield hat from Stage 1 season-specific FD regressions.",
           escape = FALSE, general_title = "Note: ", footnote_as_chunk = TRUE)

kt1_html <- kt1 %>%
  kable_styling(bootstrap_options = c("striped", "condensed"),
                full_width = FALSE, font_size = 13) %>%
  footnote(general = "Unit: district x season x gender x meal type x year. Real wages deflated by district CPI. Yield hat from Stage 1 season-specific FD regressions.",
           general_title = "Note: ", footnote_as_chunk = TRUE)

writeLines(as.character(kt1_html), file.path(out_tbl, "summary_stats_wage.html"))
kableExtra::save_kable(kt1_tex, file.path(out_tbl, "summary_stats_wage.tex"))
cat("Saved summary_stats_wage.tex/.html\n")


## ════════════════════════════════════════════════════════════════════════════ ##
## TABLE 3 — PANEL BALANCE                                                     ##
## ════════════════════════════════════════════════════════════════════════════ ##
cat("\n--- TABLE 3: Balance ---\n")

wage_districts  <- length(unique(wage$District))
wage_years      <- sort(unique(wage$year))

# Complete cases in wage panel (main regression variables)
pct_complete <- mean(complete.cases(
  wage[, c("diff_real_wage", "diff_log_yield_hat", "gender", "meal_type")]
)) * 100

obs_per_district <- wage %>%
  group_by(District) %>%
  summarise(n = n(), .groups = "drop") %>%
  summarise(mean_n = r2(mean(n)), min_n = min(n), max_n = max(n))

tbl3 <- data.frame(
  Statistic = c(
    "Districts in wage panel",
    "Years in wage panel",
    "Total obs (wage panel)",
    "% obs with all regression vars non-missing",
    "Mean obs per district (wage)",
    "Min obs per district (wage)",
    "Max obs per district (wage)",
    "Missing: diff_log_yield_hat (%)",
    "Missing: real_wage (%)",
    "Missing: ratio_double_cropped (%)"
  ),
  Value = c(
    as.character(wage_districts),
    paste(range(wage_years), collapse = "–"),
    as.character(nrow(wage)),
    sprintf("%.1f%%", pct_complete),
    as.character(obs_per_district$mean_n),
    as.character(obs_per_district$min_n),
    as.character(obs_per_district$max_n),
    sprintf("%.1f%%", mean(is.na(wage$diff_log_yield_hat)) * 100),
    sprintf("%.1f%%", mean(is.na(wage$real_wage)) * 100),
    sprintf("%.1f%%", mean(is.na(wage$ratio_double_cropped)) * 100)
  ),
  stringsAsFactors = FALSE
)

kt3_html <- kable(tbl3, format = "html",
                  caption = "Table 2: Wage Panel Balance",
                  booktabs = TRUE,
                  col.names = c("Statistic", "Value")) %>%
  kable_styling(bootstrap_options = c("striped", "condensed"),
                full_width = FALSE, font_size = 13) %>%
  footnote(general = "Wage panel: BBS wage surveys, district x season x gender x meal type x year. Climate summary statistics in output/stage1/tables/summary_stats_climate.tex.",
           general_title = "Note: ", footnote_as_chunk = TRUE)

kt3_tex <- kable(tbl3, format = "latex",
                 caption = "Wage Panel Balance",
                 booktabs = TRUE,
                 col.names = c("Statistic", "Value")) %>%
  kable_styling(latex_options = "hold_position") %>%
  footnote(general = "Wage panel: BBS wage surveys, district $\\\\times$ season $\\\\times$ gender $\\\\times$ meal type $\\\\times$ year. Climate summary in Stage 1 output.",
           escape = FALSE, general_title = "Note: ", footnote_as_chunk = TRUE)

writeLines(as.character(kt3_html), file.path(out_tbl, "summary_stats_balance.html"))
kableExtra::save_kable(kt3_tex, file.path(out_tbl, "summary_stats_balance.tex"))
cat("Saved summary_stats_balance.tex/.html\n")

## ════════════════════════════════════════════════════════════════════════════ ##

## ── Console report ────────────────────────────────────────────────────── ##
cat("\n--- Summary ---\n")
wage_gm <- wage %>%
  group_by(gender, meal_type) %>%
  summarise(mean_rw = r2(mean(real_wage, na.rm = TRUE)), N = n(), .groups = "drop")
print(wage_gm)
cat(sprintf("\nWage panel: %d districts x %d years (%d-%d) = %d obs\n",
            wage_districts, length(wage_years),
            min(wage_years), max(wage_years), nrow(wage)))
cat(sprintf("Complete cases (regression vars): %.1f%%\n", pct_complete))
cat("=== SUMMARY STATS COMPLETE ===\n")
