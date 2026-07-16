## rob_build_master_table.R
## Builds rob_m5_threemeal.tex with TWO-COLUMN structure:
##   Col A: Differential effect (β_Three, tested H₀: β_Three = β_None)
##   Col B: Total effect (β_None + β_Three, delta-method SE)
## All FD specs use +: form. R3 (levels) uses * form.
## Run after rob_r1_r3_r6_r11.R has confirmed the +: numbers.

suppressPackageStartupMessages({
  library(here); library(dplyr); library(readr); library(fixest)
  library(knitr); library(kableExtra)
})
ROOT    <- here::here()
out_tbl <- file.path(ROOT, "output/stage2/tables")
dir.create(out_tbl, recursive=TRUE, showWarnings=FALSE)

## ── Delta-method helper ───────────────────────────────────────────────────── ##
dm_total <- function(model, base_nm, int_nm, G=NULL) {
  null <- list(dm_coef=NA_real_, dm_se=NA_real_, dm_p=NA_real_)
  if (is.null(model)) return(null)
  cf <- coef(model); V <- vcov(model)
  if (!base_nm %in% names(cf) || !int_nm %in% names(cf)) return(null)
  est  <- cf[base_nm] + cf[int_nm]
  var_e <- V[base_nm,base_nm] + V[int_nm,int_nm] + 2*V[base_nm,int_nm]
  se_e  <- sqrt(max(var_e, 0))
  if (is.null(G)) G <- n_distinct(names(fixef(model)[[1]]))
  t_e  <- est / se_e
  p_e  <- 2 * pt(-abs(t_e), df=G-1L)
  list(dm_coef=as.numeric(est), dm_se=as.numeric(se_e), dm_p=as.numeric(p_e))
}

## ── Build row from model ──────────────────────────────────────────────────── ##
make_row <- function(label, model,
                     base_nm = "diff_log_yield_hat",
                     int_nm  = "diff_log_yield_hat:meal_typeThree",
                     G       = NULL,
                     n_obs   = NULL,
                     fe_label = "yr+D\\^{}s",
                     outcome  = "$\\Delta$wage",
                     ## For R8/R9: supply pre-computed SE/p directly
                     override_diff_se = NULL, override_diff_p = NULL,
                     ## For R8/R9: no delta-method
                     no_dm = FALSE) {
  stars <- function(p) ifelse(is.na(p),"",
    ifelse(p<0.01,"***",ifelse(p<0.05,"**",ifelse(p<0.10,"*",""))))
  fmt2 <- function(x) ifelse(is.na(x), "---", sprintf("%.2f", x))
  fmt3 <- function(x) ifelse(is.na(x), "---", sprintf("%.3f", x))

  if (!is.null(model)) {
    cf <- coef(model); sv <- se(model); pv <- pvalue(model)
    if (int_nm %in% names(cf)) {
      diff_coef <- as.numeric(cf[int_nm])
      diff_se   <- if (!is.null(override_diff_se)) override_diff_se else as.numeric(sv[int_nm])
      diff_p    <- if (!is.null(override_diff_p))  override_diff_p  else as.numeric(pv[int_nm])
    } else {
      diff_coef <- NA; diff_se <- NA; diff_p <- NA
    }
    if (is.null(G)) G <- n_distinct(df$District)
    if (is.null(n_obs)) n_obs <- nobs(model)
    dm <- if (no_dm) list(dm_coef=NA, dm_se=NA, dm_p=NA) else
      dm_total(model, base_nm, int_nm, G)
  } else {
    diff_coef <- NA; diff_se <- NA; diff_p <- NA
    dm <- list(dm_coef=NA, dm_se=NA, dm_p=NA)
    if (is.null(n_obs)) n_obs <- NA
  }

  data.frame(
    Specification = label,
    `Diff. coef`  = paste0(fmt2(diff_coef), stars(diff_p)),
    `Diff. SE`    = fmt2(diff_se),
    `Diff. p`     = fmt3(diff_p),
    `Total coef`  = fmt2(dm$dm_coef),
    `Total SE`    = fmt2(dm$dm_se),
    `Total p`     = fmt3(dm$dm_p),
    `N`           = ifelse(is.na(n_obs),"---",format(as.integer(n_obs),big.mark=",")),
    check.names   = FALSE, stringsAsFactors = FALSE
  )
}

## ── Load data ─────────────────────────────────────────────────────────────── ##
df <- read_csv(file.path(ROOT,"data/Regression_data/df_2_merged_v2.csv"),
               show_col_types=FALSE) %>%
  mutate(meal_type = relevel(factor(meal_type),ref="None"),
         gender    = relevel(factor(gender),    ref="Female"),
         growing_season = factor(growing_season,levels=c("Boro","Aus","Aman")))

G_all <- n_distinct(df$District)  ## 63

## ── Fit all models with +: form (FD specs) ───────────────────────────────── ##
sf <- function(fml, data, clust=~District)
  tryCatch(feols(fml,data=data,cluster=clust,warn=FALSE,notes=FALSE),
           error=function(e){cat("FAILED:",conditionMessage(e),"\n");NULL})

## Baseline
m_base <- sf(diff_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal_type + gender |
               year + District^growing_season, df)

## R1: controls (census merge)
census <- read_csv(file.path(ROOT,"data/agricultural_census_2019.csv"),show_col_types=FALSE) %>%
  mutate(irrigation_share = Net_Irrigated_Area/Net_Cultivated_Area,
         avg_holdings     = Number_of_Holdings/Net_Cultivated_Area*1000,
         crop_intensity   = Intensity_of_Cropping) %>%
  select(District,irrigation_share,avg_holdings,crop_intensity)
df_r1 <- df %>% left_join(census,by="District") %>%
  filter(!is.na(irrigation_share),!is.na(avg_holdings),!is.na(crop_intensity))
G_r1 <- n_distinct(df_r1$District)
m_r1 <- sf(diff_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal_type + gender +
              irrigation_share + avg_holdings + crop_intensity |
              year + District^growing_season, df_r1)

## R2: District FE only
m_r2 <- sf(diff_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal_type + gender |
             year + District, df)

## R3: Levels — * form required (meal_type main effects needed for levels)
df_r3 <- df %>% arrange(District,growing_season,year) %>%
  group_by(District,growing_season) %>%
  mutate(log_yield_hat=cumsum(diff_log_yield_hat)) %>% ungroup()
m_r3 <- sf(real_wage ~ log_yield_hat * meal_type + gender |
             year + District^growing_season, df_r3)

## R4: Boro only
df_boro <- df %>% filter(growing_season=="Boro")
G_r4 <- n_distinct(df_boro$District)
m_r4 <- sf(diff_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal_type + gender |
             year + District, df_boro)

## R5: Drop top-5 variance
top5 <- df %>% group_by(District) %>%
  summarise(v=var(diff_log_yield_hat,na.rm=TRUE)) %>%
  arrange(desc(v)) %>% slice_head(n=5) %>% pull(District)
df_r5 <- df %>% filter(!District %in% top5)
G_r5 <- n_distinct(df_r5$District)
m_r5 <- sf(diff_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal_type + gender |
             year + District^growing_season, df_r5)

## R6a: Early 2017-2020
df_r6a <- df %>% filter(year<=2020)
G_r6a <- n_distinct(df_r6a$District)
m_r6a <- sf(diff_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal_type + gender |
              year + District^growing_season, df_r6a)

## R6b: Late 2021-2023
df_r6b <- df %>% filter(year>=2021)
G_r6b <- n_distinct(df_r6b$District)
m_r6b <- sf(diff_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal_type + gender |
              year + District^growing_season, df_r6b)

## R10: Log wage
m_r10 <- sf(diff_log_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal_type + gender |
              year + District^growing_season, df)

## R11: Nominal wage
wage_fd <- read_csv(file.path(ROOT,"data/Regression_data/wage_by_growing_season.csv"),
                    show_col_types=FALSE) %>%
  filter(!is.na(wage)) %>%
  group_by(District,growing_season,gender,meal_type,year) %>%
  summarise(nominal_wage=mean(wage,na.rm=TRUE),.groups="drop") %>%
  arrange(District,growing_season,gender,meal_type,year) %>%
  group_by(District,growing_season,gender,meal_type) %>%
  mutate(diff_nominal_wage=nominal_wage-lag(nominal_wage)) %>%
  ungroup() %>% filter(!is.na(diff_nominal_wage))
df_r11 <- df %>%
  select(District,growing_season,gender,meal_type,year,diff_log_yield_hat) %>%
  inner_join(wage_fd %>% select(District,growing_season,gender,meal_type,year,diff_nominal_wage),
             by=c("District","growing_season","gender","meal_type","year")) %>%
  mutate(meal_type=relevel(factor(meal_type),ref="None"),
         gender=relevel(factor(gender),ref="Female"))
G_r11 <- n_distinct(df_r11$District)
m_r11 <- sf(diff_nominal_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal_type + gender |
              year + District^growing_season, df_r11)

## Female / Male
df_f <- df %>% filter(gender=="Female")
df_m <- df %>% filter(gender=="Male")
G_f  <- n_distinct(df_f$District)
G_m  <- n_distinct(df_m$District)
m_f  <- sf(diff_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal_type |
             year + District^growing_season, df_f)
m_m  <- sf(diff_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal_type |
             year + District^growing_season, df_m)

## ── R9: Two-way cluster — check df ────────────────────────────────────────── ##
m_r9 <- tryCatch(
  feols(diff_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal_type + gender |
          year + District^growing_season,
        data=df, cluster=~District+year, warn=FALSE, notes=FALSE),
  error=function(e){cat("R9 error:",conditionMessage(e),"\n");NULL})
if (!is.null(m_r9)) {
  int_nm_r9 <- "diff_log_yield_hat:meal_typeThree"
  r9_int_coef <- coef(m_r9)[int_nm_r9]
  r9_int_se   <- se(m_r9)[int_nm_r9]
  r9_int_p_default <- pvalue(m_r9)[int_nm_r9]
  ## Manual p at df=62 (G-1)
  r9_t_df62 <- r9_int_coef / r9_int_se
  r9_p_df62 <- 2*pt(-abs(r9_t_df62), df=62)
  cat(sprintf("\nR9 Two-way cluster:\n"))
  cat(sprintf("  interaction coef=%.4f  SE=%.4f\n", r9_int_coef, r9_int_se))
  cat(sprintf("  p (fixest default, df=min_clust-1=6): %.4f\n", r9_int_p_default))
  cat(sprintf("  p (df=62, G_district-1):              %.4f\n", r9_p_df62))
  cat(sprintf("  [dk_se_comparison.tex showed p=0.028 using normal approx]\n"))
}

## Build r9_row before assembling full list
if (!is.null(m_r9)) {
  r9_int_coef2 <- as.numeric(coef(m_r9)["diff_log_yield_hat:meal_typeThree"])
  r9_int_se2   <- as.numeric(se(m_r9)["diff_log_yield_hat:meal_typeThree"])
  r9_int_p2    <- as.numeric(pvalue(m_r9)["diff_log_yield_hat:meal_typeThree"])
  r9_row <- data.frame(
    Specification = "R9: Two-way cluster$^{\\S}$",
    `Diff. coef`  = paste0(sprintf("%.2f", r9_int_coef2),
                           ifelse(r9_int_p2<0.01,"***",ifelse(r9_int_p2<0.05,"**",
                             ifelse(r9_int_p2<0.10,"*","")))),
    `Diff. SE`    = sprintf("%.2f", r9_int_se2),
    `Diff. p`     = sprintf("%.3f", r9_int_p2),
    `Total coef`="---",`Total SE`="---",`Total p`="---",
    `N`="5,946",
    check.names=FALSE, stringsAsFactors=FALSE)
} else {
  r9_row <- data.frame(
    Specification="R9: Two-way cluster$^{\\S}$",
    `Diff. coef`="---",`Diff. SE`="---",`Diff. p`="---",
    `Total coef`="---",`Total SE`="---",`Total p`="---",
    `N`="---",check.names=FALSE,stringsAsFactors=FALSE)
}

## ── Assemble rows ─────────────────────────────────────────────────────────── ##
rows <- list(
  make_row("\\textbf{Baseline M5}", m_base, G=G_all),
  make_row("R1: Add controls",      m_r1,   G=G_r1),
  make_row("R2: District FE only",  m_r2,   G=G_all),
  make_row("R3: Levels$^{\\dagger}$", m_r3, G=G_all,
           base_nm="log_yield_hat",
           int_nm ="log_yield_hat:meal_typeThree",
           outcome="wage"),
  make_row("R4: Boro only",         m_r4,   G=G_r4),
  make_row("R5: Drop top-5 var$^{\\ddagger}$", m_r5, G=G_r5),
  make_row("R6a: Early (2017--2020)", m_r6a, G=G_r6a),
  make_row("R6b: Late (2021--2023)", m_r6b,  G=G_r6b),
  ## R8: DK SE — from dk_se_comparison.tex Panel B; no delta-method for non-cluster VCV
  data.frame(
    Specification = "R8: DK SE (lag=2)$^{\\S}$",
    `Diff. coef`  = "$-222.60$**",
    `Diff. SE`    = "95.75",
    `Diff. p`     = "0.020",
    `Total coef`  = "---",`Total SE`  = "---",`Total p`   = "---",
    `N`           = "5,946",
    check.names=FALSE, stringsAsFactors=FALSE
  ),
  r9_row,
  make_row("R10: Log wage",   m_r10, G=G_all,
           outcome="$\\Delta\\log$(wage)"),
  make_row("R11: Nominal wage", m_r11, G=G_r11,
           outcome="$\\Delta$nom."),
  make_row("Male only",   m_m, G=G_m),
  make_row("Female only", m_f, G=G_f)
)

tbl <- bind_rows(rows)
print(tbl, row.names=FALSE)

## ── Count significance ────────────────────────────────────────────────────── ##
diff_p_vec <- as.numeric(sub("[*]+$","",tbl[["Diff. p"]]))
n_sig <- sum(!is.na(diff_p_vec) & diff_p_vec < 0.10, na.rm=TRUE)
n_total <- sum(!is.na(diff_p_vec))
cat(sprintf("\nDifferential significant at p<0.10: %d of %d rows\n", n_sig, n_total))

## ── Build LaTeX ──────────────────────────────────────────────────────────────  ##
header_note <- paste(
  "\\textit{Differential effect}: tests H$_0$: $\\beta_{\\text{Three}} = \\beta_{\\text{None}}$",
  "(implicit contract prediction). SE = $\\text{SE}(\\hat{\\beta}_{\\text{Three}})$.",
  "\\textit{Total Three-meal}: $\\hat{\\beta}_{\\text{None}} + \\hat{\\beta}_{\\text{Three}}$;",
  "SE via delta method including $\\text{Cov}(\\hat{\\beta}_{\\text{None}}, \\hat{\\beta}_{\\text{Three}}) = +6{,}554$.",
  "$^{*}$p$<$0.1, $^{**}$p$<$0.05, $^{***}$p$<$0.01.",
  "$^{\\dagger}$ R3 uses levels regression with $^*$ form (meal-type main effects required).",
  "$^{\\ddagger}$ Differential remains negative ($-189.6$, $p=0.056$) excluding five",
  "highest-variance districts (Netrokona, Sunamganj, Sherpur, Jamalpur, Gaibandha);",
  "the total effect is positive because the no-meal baseline is imprecisely estimated ($+270.2$,",
  "SE$=281.8$) in this subsample.",
  "$^{\\S}$ R8/R9 report interaction-term SE under alternative VCV (DK lag=2 / two-way cluster);",
  "delta-method total not computed for non-district-clustered VCV."
)

kt_tex <- kable(tbl, format="latex", booktabs=TRUE,
                escape=FALSE,
                caption="ROB-M5: Three-meal Pass-Through --- Differential vs.~Total Effect",
                col.names=c("Specification",
                            "Coef.","SE","$p$",
                            "Coef.","SE","$p$",
                            "$N$")) %>%
  add_header_above(c(" "=1,
                     "Differential vs.\\ No-meal (interaction)"=3,
                     "Total Three-meal (delta-method)"=3,
                     " "=1),
                   escape=FALSE) %>%
  kable_styling(latex_options=c("hold_position","scale_down")) %>%
  row_spec(1, bold=TRUE) %>%
  footnote(general=header_note, escape=FALSE,
           general_title="\\textit{Note:} ", footnote_as_chunk=TRUE)

kableExtra::save_kable(kt_tex, file.path(out_tbl,"rob_m5_threemeal.tex"))
cat(sprintf("\nSaved rob_m5_threemeal.tex\n"))
cat("Task 1 complete. All specs unified to +: form.\n")
