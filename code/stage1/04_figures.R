# stage1_publication_figures_fixed.R
# FIGURE 1: 4-panel summary | FIGURE 2: Spatial heterogeneity
# Loads from: output/stage1/models/stage1_main_models.RData
#             output/stage1/fitted/yield_hat_2017_2023.csv

Sys.setenv(PROJ_LIB = "")

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(fixest)
  library(ggplot2); library(patchwork); library(here)
})

ROOT <- here::here()
PLOTDIR <- file.path(ROOT, "output/stage1/plots")
dir.create(PLOTDIR, showWarnings = FALSE, recursive = TRUE)

# Okabe-Ito palette (colorblind-safe)
OI <- c(Boro = "#E69F00", Aus = "#56B4E9", Aman = "#009E73")

# ── Load models ───────────────────────────────────────────────────────────────
load(file.path(ROOT, "output/stage1/models/stage1_main_models.RData"))
# Objects: m_main, m_boro, m_aus, m_aman, m_levels, m_fe_a, m_fe_b, m_fe_c

# ── Load panel (for per-season joint models used in coefplots) ────────────────
df_est <- read_csv(
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
  filter(!is.na(diff_log_yield), !is.na(diff_gdd_10_30), !is.na(diff_edd_30))

# Joint joint per-season models (GDD + EDD, year FE) for coefplots
m_j_boro <- feols(diff_log_yield ~ diff_gdd_10_30 + diff_edd_30 | year,
                  data = df_est %>% filter(season == "Boro"), cluster = ~district, warn = FALSE, notes = FALSE)
m_j_aus  <- feols(diff_log_yield ~ diff_gdd_10_30 + diff_edd_30 | year,
                  data = df_est %>% filter(season == "Aus"),  cluster = ~district, warn = FALSE, notes = FALSE)
m_j_aman <- feols(diff_log_yield ~ diff_gdd_10_30 + diff_edd_30 | year,
                  data = df_est %>% filter(season == "Aman"), cluster = ~district, warn = FALSE, notes = FALSE)

# Extract coef + 95% CI for a given variable from list of fixest models
coef_ci <- function(model_list, var) {
  lapply(seq_along(model_list), function(i) {
    m  <- model_list[[i]]
    co <- coef(m); se <- se(m); pv <- pvalue(m)
    if (!var %in% names(co)) return(NULL)
    data.frame(label = names(model_list)[i],
               coef  = co[var], se = se[var], pval = pv[var],
               lo95  = co[var] - 1.96 * se[var],
               hi95  = co[var] + 1.96 * se[var])
  }) %>% bind_rows()
}

seas_models <- list(Boro = m_j_boro, Aus = m_j_aus, Aman = m_j_aman)

edd_coefs <- coef_ci(seas_models, "diff_edd_30") %>%
  mutate(label = factor(label, levels = c("Boro","Aus","Aman")),
         sig   = ifelse(pval < 0.05, "p<0.05", ifelse(pval < 0.1, "p<0.10", "ns")))
gdd_coefs <- coef_ci(seas_models, "diff_gdd_10_30") %>%
  mutate(label = factor(label, levels = c("Boro","Aus","Aman")),
         sig   = ifelse(pval < 0.05, "p<0.05", ifelse(pval < 0.1, "p<0.10", "ns")))

# ── Panel 1: EDD coef by season ───────────────────────────────────────────────
p1 <- ggplot(edd_coefs, aes(x = label, y = coef, color = label)) +
  geom_hline(yintercept = 0, color = "grey50", linetype = "dashed") +
  geom_errorbar(aes(ymin = lo95, ymax = hi95), width = 0.18, linewidth = 0.9) +
  geom_point(aes(shape = sig), size = 3.5) +
  scale_color_manual(values = OI, guide = "none") +
  scale_shape_manual(values = c("p<0.05"=16, "p<0.10"=17, "ns"=1), name = "") +
  labs(title = "A. DeltaEDD (>30°C) coefficient by season",
       x = NULL, y = "Coefficient (95% CI)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11),
        legend.position = "bottom")

# ── Panel 2: GDD coef by season ───────────────────────────────────────────────
p2 <- ggplot(gdd_coefs, aes(x = label, y = coef, color = label)) +
  geom_hline(yintercept = 0, color = "grey50", linetype = "dashed") +
  geom_errorbar(aes(ymin = lo95, ymax = hi95), width = 0.18, linewidth = 0.9) +
  geom_point(aes(shape = sig), size = 3.5) +
  scale_color_manual(values = OI, guide = "none") +
  scale_shape_manual(values = c("p<0.05"=16, "p<0.10"=17, "ns"=1), name = "") +
  labs(title = "B. DeltaGDD (10-30°C) coefficient by season",
       x = NULL, y = "Coefficient (95% CI)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11),
        legend.position = "bottom")

# ── Panel 3: FE comparison within R² ─────────────────────────────────────────
fe_r2 <- data.frame(
  spec = factor(c("(A) Year+District","(B) Year+Season","(C) Year+Dist×Season"),
                levels = c("(A) Year+District","(B) Year+Season","(C) Year+Dist×Season")),
  overall = c(r2(m_fe_a, "r2"), r2(m_fe_b, "r2"), r2(m_fe_c, "r2")),
  within  = c(r2(m_fe_a, "wr2"), r2(m_fe_b, "wr2"), r2(m_fe_c, "wr2"))
)

p3 <- ggplot(fe_r2, aes(x = spec, y = within, fill = spec)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = c("#457B9D","#A8DADC","#E63946")) +
  labs(title = "C. Within-R² by FE specification",
       x = NULL, y = "Within R²") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11),
        axis.text.x = element_text(size = 8))

# ── Panel 4: Fitted vs actual (2017-2023) ─────────────────────────────────────
df_hat <- read_csv(
  file.path(ROOT, "output/stage1/fitted/yield_hat_2017_2023.csv"),
  show_col_types = FALSE
) %>%
  mutate(diff_log_yield = log_yield - (log_yield - residual - yield_hat))
# yield_hat = fitted Deltalog_yield; residual = actual Deltalog_yield - yield_hat
# so actual Deltalog_yield = yield_hat + residual

df_hat$actual_fd <- df_hat$yield_hat + df_hat$residual

p4 <- ggplot(df_hat, aes(x = actual_fd, y = yield_hat, color = season)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, color = "grey30", linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  scale_color_manual(values = OI, name = "Season") +
  labs(title = "D. Fitted vs actual Deltalog(yield)  (2017-2023)",
       x = "Actual Deltalog(yield)", y = "Fitted Deltalog(yield)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11),
        legend.position = "bottom")

# ── Combine FIGURE 1 ──────────────────────────────────────────────────────────
fig1 <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title   = "First Stage: Climate Effects on Bangladesh Rice Yield",
    caption = "FD spec. Clustered SE by district. Estimation sample: 2013-2023.",
    theme   = theme(plot.title   = element_text(size = 13, face = "bold"),
                    plot.caption = element_text(size = 9))
  )

ggsave(file.path(PLOTDIR, "FIGURE1_summary_panel.png"), fig1, width = 12, height = 9, dpi = 300)
ggsave(file.path(PLOTDIR, "FIGURE1_summary_panel.pdf"), fig1, width = 12, height = 9)
cat("✓ FIGURE1 saved (PNG + PDF)\n")

# ── FIGURE 2: Spatial heterogeneity bar chart ──────────────────────────────────
cat("\n=== FIGURE 2: Spatial heterogeneity ===\n")

# Mean EDD per district for exposure ranking
edd_exposure <- df_est %>%
  group_by(district) %>%
  summarise(mean_edd = mean(edd_30, na.rm = TRUE), .groups = "drop") %>%
  slice_max(mean_edd, n = 20)

# District-level EDD coef (FD, year + season FE, no cluster since N small)
dist_coefs <- df_est %>%
  group_by(district) %>%
  group_modify(function(d, g) {
    if (nrow(d) < 15 || sd(d$diff_edd_30, na.rm=TRUE) == 0)
      return(data.frame(edd_coef=NA_real_, edd_se=NA_real_, edd_pval=NA_real_))
    tryCatch({
      m  <- feols(diff_log_yield ~ diff_edd_30 | year + season, data = d, warn=FALSE, notes=FALSE)
      co <- coef(m); s <- se(m); pv <- pvalue(m)
      data.frame(edd_coef=co["diff_edd_30"], edd_se=s["diff_edd_30"], edd_pval=pv["diff_edd_30"])
    }, error = function(e) data.frame(edd_coef=NA_real_, edd_se=NA_real_, edd_pval=NA_real_))
  }) %>%
  ungroup()

fig2_data <- edd_exposure %>%
  left_join(dist_coefs, by = "district") %>%
  filter(!is.na(edd_coef)) %>%
  mutate(
    sig_level = case_when(
      edd_pval < 0.01  ~ "p<0.01",
      edd_pval < 0.05  ~ "p<0.05",
      edd_pval < 0.10  ~ "p<0.10",
      TRUE             ~ "ns"
    ),
    sig_level = factor(sig_level, levels = c("p<0.01","p<0.05","p<0.10","ns")),
    district  = reorder(district, edd_coef)
  )

p_fig2 <- ggplot(fig2_data, aes(x = district, y = edd_coef, fill = sig_level)) +
  geom_col() +
  geom_errorbar(aes(ymin = edd_coef - 1.96*edd_se,
                    ymax = edd_coef + 1.96*edd_se), width = 0.3, linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "grey30") +
  scale_fill_manual(
    values = c("p<0.01"="#D73027","p<0.05"="#FC8D59","p<0.10"="#FEE08B","ns"="grey70"),
    name   = "Significance"
  ) +
  coord_flip() +
  labs(
    title    = "Spatial Heterogeneity in Heat Stress Effects\n(Top 20 Most EDD-Exposed Districts)",
    subtitle = "FD spec: Deltalog(yield) ~ DeltaEDD | year + season",
    x        = NULL, y = "DeltaEDD Coefficient (95% CI)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "right")

ggsave(file.path(PLOTDIR, "FIGURE2_spatial_heterogeneity.png"), p_fig2, width = 10, height = 7, dpi = 300)
ggsave(file.path(PLOTDIR, "FIGURE2_spatial_heterogeneity.pdf"), p_fig2, width = 10, height = 7)
cat("✓ FIGURE2 saved (PNG + PDF)\n")

cat("\n=== PUBLICATION FIGURES COMPLETE ===\n")

