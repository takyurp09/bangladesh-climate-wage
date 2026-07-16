# ============================================================
# rob1_thresholds_byseason.R
# New Table A6: Season-Specific GDD/EDD Threshold Sensitivity
#
# Runs threshold sensitivity SEPARATELY for Boro and Aus.
# Boro: 6 specs varying GDD lower/upper bounds [B1-B6]
# Aus:  4 specs varying EDD threshold [A1-A4]
# Aman excluded (no significant predictors — placebo only).
#
# FE:       year FE only (season-specific FD panel, matches 01_main_regressions.R)
# SE:       clustered by district
# Output:   output/stage1/tables/rob1_thresholds_byseason.tex
# ============================================================

suppressPackageStartupMessages({
  library(fixest)
  library(dplyr)
  library(readr)
  library(here)
})

ROOT <- here::here()
dir.create(file.path(ROOT, "output/stage1/tables"), recursive = TRUE, showWarnings = FALSE)

# ── 1. Data preparation ───────────────────────────────────────────────────────

df_raw <- read_csv(
  file.path(ROOT, "data/Regression_data/bangladesh_rice_regression_panel.csv"),
  show_col_types = FALSE
) %>%
  arrange(district, season, year) %>%
  group_by(district, season) %>%
  mutate(
    log_yield      = log(yield_per_ha),
    diff_log_yield = log_yield - lag(log_yield),
    # ── Boro: GDD lower-bound variants ────────────────────────
    diff_gdd_8_30  = gdd_8_30  - lag(gdd_8_30),   # B1
    diff_gdd_10_30 = gdd_10_30 - lag(gdd_10_30),  # B2 MAIN
    diff_gdd_12_30 = gdd_12_30 - lag(gdd_12_30),  # B3
    diff_gdd_15_30 = gdd_15_30 - lag(gdd_15_30),  # B4
    # ── Boro: GDD upper-bound variants ────────────────────────
    diff_gdd_10_32 = gdd_10_32 - lag(gdd_10_32),  # B5
    diff_gdd_10_35 = gdd_10_35 - lag(gdd_10_35),  # B6
    # ── Aus: EDD threshold variants ───────────────────────────
    diff_edd_28    = edd_28 - lag(edd_28),         # A1
    diff_edd_30    = edd_30 - lag(edd_30),         # A2 MAIN
    diff_edd_32    = edd_32 - lag(edd_32),         # A3
    diff_edd_35    = edd_35 - lag(edd_35)          # A4
  ) %>%
  ungroup() %>%
  filter(year >= 2013, !is.na(diff_log_yield))

df_boro <- df_raw %>% filter(season == "Boro", !is.na(diff_gdd_10_30))
df_aus  <- df_raw %>% filter(season == "Aus",  !is.na(diff_edd_30))

cat("Boro N:", nrow(df_boro), "| Aus N:", nrow(df_aus), "\n")

# ── 2. Boro regressions: vary GDD threshold ───────────────────────────────────

mb1 <- feols(diff_log_yield ~ diff_gdd_8_30  | year, data = df_boro, cluster = ~district)
mb2 <- feols(diff_log_yield ~ diff_gdd_10_30 | year, data = df_boro, cluster = ~district)
mb3 <- feols(diff_log_yield ~ diff_gdd_12_30 | year, data = df_boro, cluster = ~district)
mb4 <- feols(diff_log_yield ~ diff_gdd_15_30 | year, data = df_boro, cluster = ~district)
mb5 <- feols(diff_log_yield ~ diff_gdd_10_32 | year, data = df_boro, cluster = ~district)
mb6 <- feols(diff_log_yield ~ diff_gdd_10_35 | year, data = df_boro, cluster = ~district)

boro_models <- list(mb1, mb2, mb3, mb4, mb5, mb6)
cat("Boro regressions complete\n")

# ── 3. Aus regressions: vary EDD threshold ────────────────────────────────────

ma1 <- feols(diff_log_yield ~ diff_edd_28 | year, data = df_aus, cluster = ~district)
ma2 <- feols(diff_log_yield ~ diff_edd_30 | year, data = df_aus, cluster = ~district)
ma3 <- feols(diff_log_yield ~ diff_edd_32 | year, data = df_aus, cluster = ~district)
ma4 <- feols(diff_log_yield ~ diff_edd_35 | year, data = df_aus, cluster = ~district)

aus_models <- list(ma1, ma2, ma3, ma4)
cat("Aus regressions complete\n")

# ── 4. Extract stats ──────────────────────────────────────────────────────────

extract_stats <- function(m) {
  ct   <- coeftable(m)
  coef <- ct[1, "Estimate"]
  se   <- ct[1, "Std. Error"]
  pval <- ct[1, "Pr(>|t|)"]
  stars <- dplyr::case_when(
    pval < 0.01 ~ "***",
    pval < 0.05 ~ "**",
    pval < 0.10 ~ "*",
    TRUE        ~ ""
  )
  list(
    coef   = coef,
    se     = se,
    pval   = pval,
    stars  = stars,
    n      = nobs(m),
    within = unname(r2(m, "wr2"))
  )
}

bs  <- lapply(boro_models, extract_stats)
as_ <- lapply(aus_models,  extract_stats)

# Print results to console
cat("\n--- Boro results ---\n")
for (i in seq_along(bs)) {
  cat(sprintf("B%d: coef=%.6f%s, SE=%.6f, p=%.3f, N=%d, Within-R2=%.3f\n",
              i, bs[[i]]$coef, bs[[i]]$stars, bs[[i]]$se,
              bs[[i]]$pval, bs[[i]]$n, bs[[i]]$within))
}
cat("\n--- Aus results ---\n")
for (i in seq_along(as_)) {
  cat(sprintf("A%d: coef=%.6f%s, SE=%.6f, p=%.3f, N=%d, Within-R2=%.3f\n",
              i, as_[[i]]$coef, as_[[i]]$stars, as_[[i]]$se,
              as_[[i]]$pval, as_[[i]]$n, as_[[i]]$within))
}

# ── 5. LaTeX helpers ──────────────────────────────────────────────────────────

fmt_coef <- function(x)       sprintf("%.6f", x)
fmt_se   <- function(x)       sprintf("(%.6f)", x)
fmt_r2   <- function(x)       sprintf("%.3f", x)
fmt_n    <- function(x)       formatC(x, format = "d", big.mark = ",")
fmt_pval <- function(x) {
  if (is.na(x) || x < 0.001) "$<$0.001" else sprintf("%.3f", x)
}

# Row builder: label & N cells separated by " & ", closed by " \\"
row_tex <- function(label, cells, extra = "") {
  paste0(label, " & ", paste(cells, collapse = " & "), " \\\\", extra)
}

# ── 6. Assemble cell vectors ──────────────────────────────────────────────────

# Boro (6 columns, Panel A)
b_coef <- sapply(1:6, function(i) paste0(fmt_coef(bs[[i]]$coef), bs[[i]]$stars))
b_se   <- sapply(1:6, function(i) fmt_se(bs[[i]]$se))
b_n    <- sapply(1:6, function(i) fmt_n(bs[[i]]$n))
b_r2   <- sapply(1:6, function(i) fmt_r2(bs[[i]]$within))
b_pv   <- sapply(1:6, function(i) fmt_pval(bs[[i]]$pval))

# Aus (4 columns, padded to 6 for alignment in same tabular)
a_coef <- c(sapply(1:4, function(i) paste0(fmt_coef(as_[[i]]$coef), as_[[i]]$stars)), "", "")
a_se   <- c(sapply(1:4, function(i) fmt_se(as_[[i]]$se)), "", "")
a_n    <- c(sapply(1:4, function(i) fmt_n(as_[[i]]$n)), "", "")
a_r2   <- c(sapply(1:4, function(i) fmt_r2(as_[[i]]$within)), "", "")
a_pv   <- c(sapply(1:4, function(i) fmt_pval(as_[[i]]$pval)), "", "")

# ── 7. Build LaTeX table ──────────────────────────────────────────────────────

lines <- c(
  "% rob1_thresholds_byseason.tex",
  "% Auto-generated by code/stage1/rob1_thresholds_byseason.R",
  "%",
  "\\begin{tabular}{lcccccc}",
  "\\toprule",
  # ── Panel A header ──────────────────────────────────────────
  paste0("\\multicolumn{7}{c}{\\textit{Panel A: Boro Season",
         " \\textemdash{} GDD Threshold Sensitivity",
         " ($\\Delta\\log(\\text{yield}_{\\text{Boro}})$)}} \\\\"),
  "\\addlinespace[2pt]",
  "\\cmidrule(lr){1-7}",
  paste0(" & \\textbf{B1} & \\textbf{B2} & \\textbf{B3}",
         " & \\textbf{B4} & \\textbf{B5} & \\textbf{B6} \\\\"),
  paste0(" & $\\Delta$GDD & $\\Delta$GDD & $\\Delta$GDD",
         " & $\\Delta$GDD & $\\Delta$GDD & $\\Delta$GDD \\\\"),
  paste0(" & $[8,30]$ & $[10,30]$ & $[12,30]$",
         " & $[15,30]$ & $[10,32]$ & $[10,35]$ \\\\"),
  paste0(" & & \\multicolumn{1}{c}{\\footnotesize\\textit{(Main)}}",
         " & & & & \\\\"),
  "\\midrule",
  # ── Panel A data ────────────────────────────────────────────
  row_tex("$\\Delta$GDD", b_coef),
  row_tex("",             b_se,   "[4pt]"),
  row_tex("$N$",          b_n),
  row_tex("Within $R^{2}$", b_r2),
  row_tex("$p$-value",    b_pv),
  row_tex("Year FE",      rep("Yes", 6)),
  "\\midrule",
  # ── Panel B header ──────────────────────────────────────────
  paste0("\\multicolumn{7}{c}{\\textit{Panel B: Aus Season",
         " \\textemdash{} EDD Threshold Sensitivity",
         " ($\\Delta\\log(\\text{yield}_{\\text{Aus}})$)}} \\\\"),
  "\\addlinespace[2pt]",
  "\\cmidrule(lr){1-5}",
  paste0(" & \\textbf{A1} & \\textbf{A2} & \\textbf{A3}",
         " & \\textbf{A4} & & \\\\"),
  paste0(" & $\\Delta$EDD${>}28$ & $\\Delta$EDD${>}30$",
         " & $\\Delta$EDD${>}32$ & $\\Delta$EDD${>}35$ & & \\\\"),
  paste0(" & & \\multicolumn{1}{c}{\\footnotesize\\textit{(Main)}}",
         " & & & & \\\\"),
  "\\midrule",
  # ── Panel B data ────────────────────────────────────────────
  row_tex("$\\Delta$EDD", a_coef),
  row_tex("",             a_se,   "[4pt]"),
  row_tex("$N$",          a_n),
  row_tex("Within $R^{2}$", a_r2),
  row_tex("$p$-value",    a_pv),
  row_tex("Year FE",      c(rep("Yes", 4), "", "")),
  "\\bottomrule",
  "\\end{tabular}",
  # ── Notes ───────────────────────────────────────────────────
  "\\begin{minipage}{\\linewidth}\\smallskip",
  paste0("\\footnotesize\\textit{Notes:} ",
         "Each column is a separate season-specific FD regression. ",
         "Boro regressions use only GDD with the specified threshold. ",
         "Aus regressions use only EDD with the specified threshold. ",
         "Main specification column corresponds to the threshold ",
         "used to generate fitted values entering Stage~2. ",
         "Year FE included in all specifications. ",
         "SE clustered by district. ",
         "$^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$."),
  "\\end{minipage}"
)

out_path <- file.path(ROOT, "output/stage1/tables/rob1_thresholds_byseason.tex")
writeLines(lines, out_path)
cat("\nTable A6 saved to:", out_path, "\n")
