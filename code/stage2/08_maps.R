## code/stage2/08_maps.R
## District-level regression-based maps (Stage 2)
## All maps derived from district-level regressions only
## Outputs → output/stage2/maps/

suppressPackageStartupMessages({
  library(sf)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(ggspatial)
  library(biscale)
  library(dplyr)
  library(fixest)
  library(here)
  library(readr)
  library(stringr)
})

ROOT   <- here::here()
MAPDIR <- file.path(ROOT, "output/stage2/maps")
SHP    <- file.path(ROOT, "data/climate_data/gadm41_BGD_shp/gadm41_BGD_2.shp")
WAGE   <- file.path(ROOT, "data/Regression_data/df_2_merged_v2.csv")
CLIM   <- file.path(ROOT, "data/Regression_data/bangladesh_rice_regression_panel.csv")
dir.create(MAPDIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(ROOT, "code/utils/district_names.R"))

cat("=== STAGE 2 MAPS ===\n")

# ── Load shapefile ────────────────────────────────────────────────────────────
shp <- st_read(SHP, quiet = TRUE) %>%
  st_transform(4326) %>%
  mutate(district_key = recode_shp(NAME_2))

bgd_outline <- shp %>% st_union()

# Verify join coverage
df_wage <- read_csv(WAGE, show_col_types = FALSE) %>%
  mutate(
    district_key = district_key(District),
    gender       = relevel(factor(gender),    ref = "Female"),
    meal_type    = relevel(factor(meal_type), ref = "None")
  )

shp_keys  <- unique(shp$district_key)
wage_keys <- unique(df_wage$district_key)
cat(sprintf("Shapefile districts: %d\n", length(shp_keys)))
cat(sprintf("Wage data districts: %d\n", length(wage_keys)))
unmatched <- setdiff(wage_keys, shp_keys)
if (length(unmatched) > 0) cat("Unmatched (wage→shp):", paste(unmatched, collapse=", "), "\n")

# ── Shared map theme ──────────────────────────────────────────────────────────
map_theme <- function() {
  theme_void(base_size = 9) +
    theme(
      legend.position      = "bottom",
      legend.title         = element_text(size = 7, face = "bold"),
      legend.text          = element_text(size = 6),
      legend.key.width     = unit(1.2, "cm"),
      legend.key.height    = unit(0.35, "cm"),
      legend.margin        = margin(t = 2, b = 2),
      plot.title           = element_text(size = 9, face = "bold", hjust = 0.5,
                                          margin = margin(b = 6)),
      plot.subtitle        = element_text(size = 7, hjust = 0.5,
                                          margin = margin(b = 2)),
      plot.caption         = element_text(size = 7, hjust = 0, colour = "grey40",
                                          margin = margin(t = 6)),
      plot.margin          = margin(t = 10, r = 20, b = 10, l = 10, unit = "mm")
    )
}

border_layer <- function() {
  list(
    geom_sf(data = shp, fill = NA, colour = "grey40", linewidth = 0.15),
    geom_sf(data = bgd_outline, fill = NA, colour = "black", linewidth = 0.5)
  )
}

save_map <- function(p, fname, w = 183, h = 200) {
  path <- file.path(MAPDIR, fname)
  ggsave(path, p, width = w, height = h, units = "mm", dpi = 300, bg = "white")
  cat(sprintf("  saved \u2713 %s  [%dx%dmm]\n", fname, w, h))
}

# ─────────────────────────────────────────────────────────────────────────────
# MAP 1 — District-Level Yield-Wage Pass-Through
# For each district: diff_real_wage ~ diff_log_yield_hat | year
# ─────────────────────────────────────────────────────────────────────────────
cat("\nMAP 1: District pass-through coefficients\n")

districts <- unique(df_wage$district_key)
N_THRESH  <- 15
P_THRESH  <- 0.2

run_district_m1 <- function(df) {
  tryCatch({
    # cluster = ~district_key invalid within a single district (1 cluster)
    # Use HC1 heteroskedasticity-robust SE instead
    m <- feols(diff_real_wage ~ diff_log_yield_hat | year,
               data = df, vcov = "HC1", warn = FALSE, notes = FALSE)
    cf <- coef(m)["diff_log_yield_hat"]
    se <- sqrt(diag(vcov(m)))["diff_log_yield_hat"]
    pv <- 2 * pnorm(-abs(cf / se))
    data.frame(coef = cf, se = se, pval = pv, n_obs = nobs(m))
  }, error = function(e) NULL)
}

passthrough_res <- lapply(districts, function(d) {
  sub <- df_wage %>% filter(district_key == d, !is.na(diff_real_wage), !is.na(diff_log_yield_hat))
  res <- run_district_m1(sub)
  if (is.null(res)) return(data.frame(district_key = d, coef = NA_real_, pval = NA_real_, n_obs = 0L))
  res$district_key <- d
  res
}) %>% bind_rows()

passthrough_res <- passthrough_res %>%
  mutate(coef_plot = ifelse(is.na(pval) | pval > P_THRESH | n_obs < N_THRESH, NA_real_, coef))

n_valid <- sum(!is.na(passthrough_res$coef_plot))
n_na    <- sum(is.na(passthrough_res$coef_plot))
cat(sprintf("  Districts with valid coef: %d | set to NA: %d\n", n_valid, n_na))

map1_df <- shp %>%
  left_join(passthrough_res, by = "district_key")

# Symmetric colour limits
lim1 <- max(abs(map1_df$coef_plot), na.rm = TRUE)
lim1 <- ceiling(lim1 / 10) * 10

p_map1 <- ggplot(map1_df) +
  geom_sf(aes(fill = coef_plot), colour = "grey40", linewidth = 0.15) +
  geom_sf(data = bgd_outline, fill = NA, colour = "black", linewidth = 0.5) +
  scale_fill_distiller(
    palette  = "RdBu",
    direction = -1,
    limits   = c(-lim1, lim1),
    na.value = "grey80",
    name     = "Coef (BDT/unit Δlog yield)",
    guide    = guide_colorbar(title.position = "top", title.hjust = 0.5,
                              barwidth = 8, barheight = 0.4)
  ) +
  annotation_scale(location = "br", width_hint = 0.25, text_cex = 0.6) +
  annotation_north_arrow(location = "tl", height = unit(0.7, "cm"),
                         width = unit(0.5, "cm"),
                         style = north_arrow_fancy_orienteering(text_size = 7)) +
  labs(
    title   = "Yield-Wage Pass-Through Coefficient by District",
    caption = str_wrap("Grey = insufficient data or p>0.2. FD spec, year FE. HC1 SE (within-district).", width = 80)
  ) +
  map_theme()

save_map(p_map1, "map1_passthrough_coef.png")

# ─────────────────────────────────────────────────────────────────────────────
# MAP 2 — Three-Meal Interaction Coefficient by District
# For each district: M5 spec → extract diff_log_yield_hat:meal_typeThree
# ─────────────────────────────────────────────────────────────────────────────
cat("\nMAP 2: Three-meal interaction coefficients\n")

run_district_m5 <- function(df) {
  tryCatch({
    # Need at least some Three-meal obs
    if (sum(df$meal_type == "Three", na.rm = TRUE) < 5) return(NULL)
    m <- feols(diff_real_wage ~ diff_log_yield_hat + diff_log_yield_hat:meal_type + gender | year,
               data = df, vcov = "HC1", warn = FALSE, notes = FALSE)
    cn  <- "diff_log_yield_hat:meal_typeThree"
    if (!cn %in% names(coef(m))) return(NULL)
    cf  <- coef(m)[cn]
    se  <- sqrt(diag(vcov(m)))[cn]
    pv  <- 2 * pnorm(-abs(cf / se))
    data.frame(coef = cf, se = se, pval = pv, n_obs = nobs(m))
  }, error = function(e) NULL)
}

meal_res <- lapply(districts, function(d) {
  sub <- df_wage %>% filter(district_key == d, !is.na(diff_real_wage), !is.na(diff_log_yield_hat))
  res <- run_district_m5(sub)
  if (is.null(res)) return(data.frame(district_key = d, coef = NA_real_, pval = NA_real_, n_obs = 0L))
  res$district_key <- d
  res
}) %>% bind_rows()

meal_res <- meal_res %>%
  mutate(coef_plot = ifelse(is.na(pval) | pval > P_THRESH | n_obs < N_THRESH, NA_real_, coef))

n_valid2 <- sum(!is.na(meal_res$coef_plot))
n_na2    <- sum(is.na(meal_res$coef_plot))
cat(sprintf("  Districts with valid coef: %d | set to NA: %d\n", n_valid2, n_na2))

map2_df <- shp %>%
  left_join(meal_res, by = "district_key")

lim2 <- max(abs(map2_df$coef_plot), na.rm = TRUE)
if (is.infinite(lim2) || lim2 == 0) lim2 <- 50
lim2 <- ceiling(lim2 / 10) * 10

p_map2 <- ggplot(map2_df) +
  geom_sf(aes(fill = coef_plot), colour = "grey40", linewidth = 0.15) +
  geom_sf(data = bgd_outline, fill = NA, colour = "black", linewidth = 0.5) +
  scale_fill_distiller(
    palette   = "RdBu",
    direction = -1,
    limits    = c(-lim2, lim2),
    na.value  = "grey80",
    name      = "Three-meal differential (BDT/day)",
    guide     = guide_colorbar(title.position = "top", title.hjust = 0.5,
                               barwidth = 8, barheight = 0.4)
  ) +
  annotation_scale(location = "br", width_hint = 0.25, text_cex = 0.6) +
  annotation_north_arrow(location = "tl", height = unit(0.7, "cm"),
                         width = unit(0.5, "cm"),
                         style = north_arrow_fancy_orienteering(text_size = 7)) +
  labs(
    title   = "Three-Meal Wage Differential vs. No-Meal Workers, by District",
    caption = str_wrap("Three-meal differential = interaction coefficient (three-meal vs. no-meal). Grey = insufficient data or p>0.2. FD spec M5, year FE. HC1 SE (within-district).", width = 80)
  ) +
  map_theme()

save_map(p_map2, "map2_threemeal_interaction.png")

# ─────────────────────────────────────────────────────────────────────────────
# MAP 3 — Bivariate: EDD Exposure × Pass-Through
# ─────────────────────────────────────────────────────────────────────────────
cat("\nMAP 3: Bivariate EDD × pass-through\n")

df_clim <- read_csv(CLIM, show_col_types = FALSE) %>%
  mutate(district_key = clean_name(district)) %>%
  filter(season == "Aus", year >= 2017, year <= 2023) %>%
  group_by(district_key) %>%
  summarise(mean_edd30 = mean(edd_30, na.rm = TRUE), .groups = "drop")

biv_df <- passthrough_res %>%
  select(district_key, pt_coef = coef_plot) %>%
  left_join(df_clim, by = "district_key") %>%
  filter(!is.na(pt_coef), !is.na(mean_edd30))

cat(sprintf("  Districts in bivariate: %d\n", nrow(biv_df)))

map3_df <- shp %>%
  left_join(biv_df, by = "district_key") %>%
  filter(!is.na(pt_coef), !is.na(mean_edd30))

if (nrow(map3_df) >= 6) {
  map3_df <- bi_class(map3_df, x = mean_edd30, y = pt_coef, style = "quantile", dim = 3)

  p_map3_main <- ggplot(map3_df) +
    geom_sf(data = shp, fill = "grey90", colour = "grey40", linewidth = 0.15) +
    geom_sf(aes(fill = bi_class), colour = "grey40", linewidth = 0.15, show.legend = FALSE) +
    geom_sf(data = bgd_outline, fill = NA, colour = "black", linewidth = 0.5) +
    bi_scale_fill(pal = "GrPink", dim = 3) +
    annotation_scale(location = "br", width_hint = 0.25, text_cex = 0.6) +
    annotation_north_arrow(location = "tl", height = unit(0.7, "cm"),
                           width = unit(0.5, "cm"),
                           style = north_arrow_fancy_orienteering(text_size = 7)) +
    labs(
      title   = "Heat Exposure vs Yield-Wage Pass-Through by District",
      caption = str_wrap("EDD = Extreme Degree Days >30\u00b0C, Aus season 2017\u20132023. Bivariate 3\u00d73 grid: low/med/high EDD \u00d7 negative/zero/positive pass-through.", width = 80)
    ) +
    map_theme()

  legend3 <- bi_legend(
    pal    = "GrPink",
    dim    = 3,
    xlab   = "EDD \u2192",
    ylab   = "Pass-Through \u2192",
    size   = 7
  )

  p_map3 <- p_map3_main + inset_element(legend3, left = 0.65, bottom = 0.02,
                                         right = 1.0, top = 0.35)
  save_map(p_map3, "map3_bivariate_edd_passthrough.png")
} else {
  cat("  WARNING: insufficient districts for bivariate map — skipping\n")
}

cat("\n=== MAP SUMMARY ===\n")
cat(sprintf("MAP 1 — Districts valid/NA: %d / %d\n", n_valid, n_na))
cat(sprintf("MAP 2 — Districts valid/NA: %d / %d\n", n_valid2, n_na2))
cat(sprintf("MAP 3 — Districts in bivariate: %d\n", if (exists("biv_df")) nrow(biv_df) else 0))
cat(sprintf("Output: %s\n", MAPDIR))
