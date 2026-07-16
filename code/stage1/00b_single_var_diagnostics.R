suppressPackageStartupMessages({
  library(dplyr); library(readr); library(fixest); library(here)
})

# Manual VIF: 1 / (1 - R^2 of regressing x1 on x2 + controls)
manual_vif2 <- function(d) {
  r2_gdd <- summary(lm(diff_gdd_10_30 ~ diff_edd_30 + factor(year), data=d))$r.squared
  r2_edd <- summary(lm(diff_edd_30 ~ diff_gdd_10_30 + factor(year), data=d))$r.squared
  c(vif_gdd = 1/(1-r2_gdd), vif_edd = 1/(1-r2_edd))
}

ROOT <- here::here()
DATA <- file.path(ROOT, "data/Regression_data/bangladesh_rice_regression_panel.csv")

df <- read_csv(DATA, show_col_types = FALSE) %>%
  arrange(district, season, year) %>%
  group_by(district, season) %>%
  mutate(
    log_yield      = log(yield_per_ha),
    diff_gdd_10_30 = gdd_10_30 - lag(gdd_10_30)
  ) %>%
  ungroup() %>%
  filter(!is.na(diff_log_yield), !is.na(diff_gdd_10_30), !is.na(diff_edd_30))

seasons <- c("Boro", "Aus", "Aman")

results <- list()

for (s in seasons) {
  d <- filter(df, season == s)
  cat("\n====", s, "(N =", nrow(d), ")====\n")

  m_gdd   <- feols(diff_log_yield ~ diff_gdd_10_30              | year, data=d, cluster=~district, warn=FALSE, notes=FALSE)
  m_edd   <- feols(diff_log_yield ~ diff_edd_30                 | year, data=d, cluster=~district, warn=FALSE, notes=FALSE)
  m_joint <- feols(diff_log_yield ~ diff_gdd_10_30 + diff_edd_30 | year, data=d, cluster=~district, warn=FALSE, notes=FALSE)

  # Extract coef/SE/p
  extract_row <- function(m, var) {
    co <- coef(m); se <- se(m); pv <- pvalue(m)
    data.frame(
      coef = co[var], se = se[var], pval = pv[var],
      sig  = ifelse(pv[var]<0.01, "***", ifelse(pv[var]<0.05,"**", ifelse(pv[var]<0.10,"*","")))
    )
  }

  r_gdd_only  <- extract_row(m_gdd,   "diff_gdd_10_30")
  r_edd_only  <- extract_row(m_edd,   "diff_edd_30")
  r_joint_gdd <- extract_row(m_joint, "diff_gdd_10_30")
  r_joint_edd <- extract_row(m_joint, "diff_edd_30")

  # VIF (manual)
  vif_vals <- manual_vif2(d)
  vif_gdd  <- vif_vals["vif_gdd"]
  vif_edd  <- vif_vals["vif_edd"]

  # Sign flip?
  flip_gdd <- sign(r_gdd_only$coef) != sign(r_joint_gdd$coef)
  flip_edd <- sign(r_edd_only$coef) != sign(r_joint_edd$coef)

  cat(sprintf("  GDD-only : coef=%+.3e  SE=%.3e  p=%.3f %s\n",
              r_gdd_only$coef, r_gdd_only$se, r_gdd_only$pval, r_gdd_only$sig))
  cat(sprintf("  EDD-only : coef=%+.3e  SE=%.3e  p=%.3f %s\n",
              r_edd_only$coef, r_edd_only$se, r_edd_only$pval, r_edd_only$sig))
  cat(sprintf("  Joint GDD: coef=%+.3e  SE=%.3e  p=%.3f %s  VIF=%.2f\n",
              r_joint_gdd$coef, r_joint_gdd$se, r_joint_gdd$pval, r_joint_gdd$sig, vif_gdd))
  cat(sprintf("  Joint EDD: coef=%+.3e  SE=%.3e  p=%.3f %s  VIF=%.2f\n",
              r_joint_edd$coef, r_joint_edd$se, r_joint_edd$pval, r_joint_edd$sig, vif_edd))
  cat(sprintf("  GDD sign flip (single->joint): %s\n", ifelse(flip_gdd, "YES ⚠", "no")))
  cat(sprintf("  EDD sign flip (single->joint): %s\n", ifelse(flip_edd, "YES ⚠", "no")))

  results[[s]] <- list(
    n=nrow(d),
    gdd_only=r_gdd_only, edd_only=r_edd_only,
    joint_gdd=r_joint_gdd, joint_edd=r_joint_edd,
    vif_gdd=vif_gdd, vif_edd=vif_edd,
    flip_gdd=flip_gdd, flip_edd=flip_edd
  )
}

# ── Write markdown ────────────────────────────────────────────────────────────
OUT <- file.path(ROOT, "output/stage1/summary/single_var_models.md")

lines <- c(
  "# Single-Variable & Season-Split Diagnostics",
  "",
  "Spec: FD, year FE, cluster=district. VIF from equivalent `lm()` with factor(year).",
  ""
)

for (s in seasons) {
  r <- results[[s]]
  lines <- c(lines,
    paste0("## ", s, "  (N = ", r$n, ")"),
    "",
    "| Model | Variable | Coef | SE | p-value | Sig | VIF |",
    "|---|---|---|---|---|---|---|",
    sprintf("| GDD-only  | diff_GDD_10_30 | %+.4e | %.4e | %.3f | %s | — |",
            r$gdd_only$coef, r$gdd_only$se, r$gdd_only$pval, r$gdd_only$sig),
    sprintf("| EDD-only  | diff_EDD_30    | %+.4e | %.4e | %.3f | %s | — |",
            r$edd_only$coef, r$edd_only$se, r$edd_only$pval, r$edd_only$sig),
    sprintf("| Joint     | diff_GDD_10_30 | %+.4e | %.4e | %.3f | %s | %.2f |",
            r$joint_gdd$coef, r$joint_gdd$se, r$joint_gdd$pval, r$joint_gdd$sig, r$vif_gdd),
    sprintf("| Joint     | diff_EDD_30    | %+.4e | %.4e | %.3f | %s | %.2f |",
            r$joint_edd$coef, r$joint_edd$se, r$joint_edd$pval, r$joint_edd$sig, r$vif_edd),
    "",
    paste0("**GDD sign flip (single → joint):** ", ifelse(r$flip_gdd, "**YES ⚠**", "no")),
    paste0("  **EDD sign flip (single → joint):** ", ifelse(r$flip_edd, "**YES ⚠**", "no")),
    ""
  )
}

# Summary flag table
lines <- c(lines,
  "---",
  "",
  "## Summary: Sign Flips & VIF",
  "",
  "| Season | GDD-only sign | GDD-joint sign | GDD flip? | EDD-only sign | EDD-joint sign | EDD flip? | VIF(GDD) | VIF(EDD) |",
  "|---|---|---|---|---|---|---|---|---|"
)
for (s in seasons) {
  r <- results[[s]]
  lines <- c(lines, sprintf(
    "| %s | %s | %s | %s | %s | %s | %s | %.2f | %.2f |",
    s,
    ifelse(r$gdd_only$coef > 0, "+", "−"),
    ifelse(r$joint_gdd$coef > 0, "+", "−"),
    ifelse(r$flip_gdd, "**YES**", "no"),
    ifelse(r$edd_only$coef > 0, "+", "−"),
    ifelse(r$joint_edd$coef > 0, "+", "−"),
    ifelse(r$flip_edd, "**YES**", "no"),
    r$vif_gdd, r$vif_edd
  ))
}

writeLines(lines, OUT)
cat("\n✓ Saved to", OUT, "\n")
