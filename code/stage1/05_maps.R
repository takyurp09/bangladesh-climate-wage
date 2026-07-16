suppressPackageStartupMessages({
  library(sf)
  library(ggplot2)
  library(viridis)
  library(patchwork)
  library(ggspatial)
  library(dplyr)
  library(readr)
  library(fixest)
  library(here)
})

ROOT    <- here::here()
MAPDIR  <- file.path(ROOT, "output/stage1/maps")
SHP     <- file.path(ROOT, "data/climate_data/gadm41_BGD_shp/gadm41_BGD_2.shp")
CSV     <- file.path(ROOT, "data/Regression_data/bangladesh_rice_regression_panel.csv")

dir.create(MAPDIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(ROOT, "code/utils/district_names.R"))

# ── Load shapefile ────────────────────────────────────────────────────────────
shp <- st_read(SHP, quiet = TRUE) %>%
  mutate(district_key = recode_shp(NAME_2))

# Bangladesh outer boundary for thick border overlay
bgd_outline <- shp %>% st_union()

# ── Load & prep panel ─────────────────────────────────────────────────────────
df <- read_csv(CSV, show_col_types = FALSE) %>%
  mutate(
    district_key = district_key(district),
    log_yield    = log(yield_per_ha)
  ) %>%
  arrange(district_key, season, year) %>%
  group_by(district_key, season) %>%
  mutate(
    diff_gdd_10_30 = gdd_10_30 - lag(gdd_10_30),
    diff_log_yield_fd = log_yield - lag(log_yield)
  ) %>%
  ungroup()

# Report join coverage
shp_keys   <- unique(shp$district_key)
panel_keys <- unique(df$district_key)
unmatched_shp   <- setdiff(shp_keys, panel_keys)
unmatched_panel <- setdiff(panel_keys, shp_keys)
cat("\n=== JOIN DIAGNOSTICS ===\n")
cat("Shapefile districts:", length(shp_keys), "\n")
cat("Panel districts:    ", length(panel_keys), "\n")
cat("Unmatched in shapefile (→ grey on maps):", paste(unmatched_shp, collapse=", "), "\n")
cat("Unmatched in panel (→ no geometry):", paste(unmatched_panel, collapse=", "), "\n\n")

# ── Shared style helpers ──────────────────────────────────────────────────────
map_theme <- function() {
  theme_void(base_size = 10) +
    theme(
      legend.position    = "bottom",
      legend.key.width   = unit(0.6, "npc") * 0.3,
      legend.key.height  = unit(0.3, "cm"),
      legend.title       = element_text(size = 9),
      legend.text        = element_text(size = 8),
      plot.title         = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b=6)),
      strip.text         = element_text(size = 10, face = "bold"),
      panel.spacing      = unit(0.5, "lines")
    )
}

add_spatial_aids <- function(p) {
  p +
    annotation_scale(location = "bl", width_hint = 0.25, text_cex = 0.7) +
    annotation_north_arrow(
      location = "tr", which_north = "true",
      style = north_arrow_orienteering(text_size = 7),
      height = unit(0.8, "cm"), width = unit(0.8, "cm")
    )
}

border_geom <- geom_sf(data = bgd_outline, fill = NA, color = "black", linewidth = 0.6)

district_geom <- function(data_sf) {
  geom_sf(data = data_sf, aes(fill = fill_var),
          color = "grey30", linewidth = 0.2)
}

# ── MAP 1 — Mean EDD Exposure by District × Season ───────────────────────────
cat("Rendering Map 1 — Mean EDD exposure...\n")

m1_data <- df %>%
  group_by(district_key, season) %>%
  summarise(mean_edd = mean(edd_30, na.rm = TRUE), .groups = "drop")

m1_sf <- shp %>%
  left_join(m1_data, by = "district_key") %>%
  filter(!is.na(season)) %>%
  mutate(fill_var = mean_edd,
         season   = factor(season, levels = c("Boro","Aus","Aman")))

p1 <- ggplot() +
  geom_sf(data = m1_sf, aes(fill = fill_var), color = "grey30", linewidth = 0.2) +
  geom_sf(data = bgd_outline, fill = NA, color = "black", linewidth = 0.6) +
  facet_wrap(~season, ncol = 3) +
  scale_fill_gradientn(
    colors = c("#FFF7BC","#FEC44F","#D95F0E","#7F0000"),
    name = "Mean EDD (days)", na.value = "grey80",
    guide = guide_colorbar(barwidth = 8, barheight = 0.5, title.position = "top")
  ) +
  labs(title = "Mean Extreme Degree Days (>30°C) by District and Season") +
  map_theme()

p1 <- add_spatial_aids(p1)
ggsave(file.path(MAPDIR, "map1_edd_exposure.png"), p1, width = 12, height = 5.5, dpi = 300)
cat("✓ map1_edd_exposure.png\n")

# ── MAP 2 — Yield Trend by District ──────────────────────────────────────────
cat("Rendering Map 2 — Yield trend...\n")

m2_data <- df %>%
  mutate(period = case_when(
    year %in% 2013:2015 ~ "early",
    year %in% 2021:2023 ~ "late",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(period)) %>%
  group_by(district_key, period) %>%
  summarise(mean_yield = mean(log_yield, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = period, values_from = mean_yield) %>%
  mutate(fill_var = late - early)

m2_sf <- shp %>% left_join(m2_data, by = "district_key")

max_abs <- max(abs(m2_sf$fill_var), na.rm = TRUE)

p2 <- ggplot() +
  geom_sf(data = m2_sf, aes(fill = fill_var), color = "grey30", linewidth = 0.2) +
  border_geom +
  scale_fill_gradient2(
    low = "#D73027", mid = "white", high = "#2166AC",
    midpoint = 0, limits = c(-max_abs, max_abs),
    name = "Δ log(yield)", na.value = "grey80",
    guide = guide_colorbar(barwidth = 8, barheight = 0.5, title.position = "top")
  ) +
  labs(title = "Change in Log Yield: 2021–2023 vs 2013–2015") +
  map_theme()

p2 <- add_spatial_aids(p2)
ggsave(file.path(MAPDIR, "map2_yield_trend.png"), p2, width = 6, height = 6, dpi = 300)
cat("✓ map2_yield_trend.png\n")

# ── MAP 3 — District-Level EDD Coefficient ───────────────────────────────────
cat("Rendering Map 3 — District EDD coefficients...\n")

fd_df <- df %>%
  filter(!is.na(diff_edd_30), !is.na(diff_log_yield_fd))

m3_data <- fd_df %>%
  group_by(district_key) %>%
  group_modify(function(d, g) {
    if (nrow(d) < 20) return(data.frame(edd_coef = NA_real_, edd_p = NA_real_))
    if (sd(d$diff_edd_30, na.rm = TRUE) == 0) return(data.frame(edd_coef = NA_real_, edd_p = NA_real_))
    tryCatch({
      m <- feols(diff_log_yield_fd ~ diff_edd_30 | year, data = d,
                 warn = FALSE, notes = FALSE)
      co <- coef(m)["diff_edd_30"]
      pv <- pvalue(m)["diff_edd_30"]
      data.frame(edd_coef = co, edd_p = pv)
    }, error = function(e) data.frame(edd_coef = NA_real_, edd_p = NA_real_))
  }) %>%
  ungroup() %>%
  mutate(fill_var = ifelse(!is.na(edd_p) & edd_p <= 0.2, edd_coef, NA_real_))

m3_sf <- shp %>% left_join(m3_data, by = "district_key")

max3 <- max(abs(m3_sf$fill_var), na.rm = TRUE)
if (is.infinite(max3) || max3 == 0) max3 <- 1

p3 <- ggplot() +
  geom_sf(data = m3_sf, aes(fill = fill_var), color = "grey30", linewidth = 0.2) +
  border_geom +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#D73027",
    midpoint = 0, limits = c(-max3, max3),
    name = "EDD coef (p≤0.2)", na.value = "grey75",
    guide = guide_colorbar(barwidth = 8, barheight = 0.5, title.position = "top")
  ) +
  labs(title = "District-Level EDD Coefficient (FD, Year FE)\nGrey = p>0.2 or N<20") +
  map_theme()

p3 <- add_spatial_aids(p3)
ggsave(file.path(MAPDIR, "map3_edd_coef.png"), p3, width = 6, height = 6, dpi = 300)
cat("✓ map3_edd_coef.png\n")

# ── MAP 4 — Year of Maximum EDD Shock ────────────────────────────────────────
cat("Rendering Map 4 — Year of peak EDD...\n")

m4_data <- df %>%
  group_by(district_key, year) %>%
  summarise(total_edd = sum(edd_30, na.rm = TRUE), .groups = "drop") %>%
  group_by(district_key) %>%
  slice_max(total_edd, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(fill_var = factor(year))

m4_sf <- shp %>% left_join(m4_data, by = "district_key")

year_levels <- as.character(2013:2023)
pal_yr <- setNames(
  colorRampPalette(c("#313695","#74ADD1","#FEE090","#F46D43","#A50026"))(11),
  year_levels
)

p4 <- ggplot() +
  geom_sf(data = m4_sf, aes(fill = fill_var), color = "grey30", linewidth = 0.2) +
  border_geom +
  scale_fill_manual(
    values = pal_yr, na.value = "grey80",
    name = "Year", drop = FALSE,
    guide = guide_legend(nrow = 1, title.position = "top",
                         keywidth = unit(0.5,"cm"), keyheight = unit(0.3,"cm"))
  ) +
  labs(title = "Year of Peak Heat Stress by District") +
  map_theme() +
  theme(legend.key.width = unit(0.5, "cm"))

p4 <- add_spatial_aids(p4)
ggsave(file.path(MAPDIR, "map4_peak_edd_year.png"), p4, width = 7, height = 6.5, dpi = 300)
cat("✓ map4_peak_edd_year.png\n")

# ── MAP 5 — Mean GDD Exposure by District × Season ───────────────────────────
cat("Rendering Map 5 — Mean GDD exposure...\n")

m5_data <- df %>%
  group_by(district_key, season) %>%
  summarise(mean_gdd = mean(gdd_10_30, na.rm = TRUE), .groups = "drop")

m5_sf <- shp %>%
  left_join(m5_data, by = "district_key") %>%
  filter(!is.na(season)) %>%
  mutate(fill_var = mean_gdd,
         season   = factor(season, levels = c("Boro","Aus","Aman")))

p5 <- ggplot() +
  geom_sf(data = m5_sf, aes(fill = fill_var), color = "grey30", linewidth = 0.2) +
  geom_sf(data = bgd_outline, fill = NA, color = "black", linewidth = 0.6) +
  facet_wrap(~season, ncol = 3) +
  scale_fill_viridis_c(
    option = "mako", direction = 1,
    name = "Mean GDD (degree-days)", na.value = "grey80",
    guide = guide_colorbar(barwidth = 8, barheight = 0.5, title.position = "top")
  ) +
  labs(title = "Mean Growing Degree Days (10–30°C) by District and Season") +
  map_theme()

p5 <- add_spatial_aids(p5)
ggsave(file.path(MAPDIR, "map5_gdd_exposure.png"), p5, width = 12, height = 5.5, dpi = 300)
cat("✓ map5_gdd_exposure.png\n")

# ── MAP 6 — Within-District EDD-Yield Correlation ────────────────────────────
cat("Rendering Map 6 — EDD-yield correlation...\n")

m6_data <- df %>%
  group_by(district_key, season) %>%
  summarise(
    fill_var = cor(edd_30, log_yield, use = "complete.obs"),
    .groups = "drop"
  )

m6_sf <- shp %>%
  left_join(m6_data, by = "district_key") %>%
  filter(!is.na(season)) %>%
  mutate(season = factor(season, levels = c("Boro","Aus","Aman")))

p6 <- ggplot() +
  geom_sf(data = m6_sf, aes(fill = fill_var), color = "grey30", linewidth = 0.2) +
  geom_sf(data = bgd_outline, fill = NA, color = "black", linewidth = 0.6) +
  facet_wrap(~season, ncol = 3) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#D73027",
    midpoint = 0, limits = c(-1, 1),
    name = "Pearson r(EDD, log yield)", na.value = "grey80",
    guide = guide_colorbar(barwidth = 8, barheight = 0.5, title.position = "top")
  ) +
  labs(title = "Within-District Correlation: EDD and Log Yield") +
  map_theme()

p6 <- add_spatial_aids(p6)
ggsave(file.path(MAPDIR, "map6_edd_yield_corr.png"), p6, width = 12, height = 5.5, dpi = 300)
cat("✓ map6_edd_yield_corr.png\n")

cat("\n=== ALL 6 MAPS SAVED to", MAPDIR, "===\n")
