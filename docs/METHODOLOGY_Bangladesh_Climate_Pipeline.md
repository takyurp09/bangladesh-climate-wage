# Methodology: Bangladesh Rice Climate Indices Pipeline

**Authors:** Taky Tahmid  
**Date:** April 5, 2026  
**Institution:** University of Delaware  
**Pipeline Version:** 1.0

---

## Overview

We developed a spatially-explicit, season-specific climate indices pipeline for Bangladesh rice (Boro, Aus, Aman) following the methodology of Schlenker & Roberts (2009) and the CI2026 global aggregation framework. The pipeline converts hourly ERA5 reanalysis data to district-level seasonal climate exposures using Bangladesh-specific rice crop calendars, harvested-area weighted spatial aggregation, and sinusoidal degree-day calculations.

**Key Innovation:** This is the first application of the CI2026 aggregation methodology to Bangladesh using empirically-derived, season-specific rice calendars at 5-arcminute resolution.

---

## Pipeline Architecture

The pipeline consists of five sequential scripts that process raw climate data into regression-ready panels:

```
Script 1: ERA5 hourly → daily climate data
    ↓
Script 2: Daily climate → daily degree-day indices
    ↓
Script 3: Rice calendars → temporal aggregation windows
    ↓
Script 4: Daily indices → seasonal district aggregates (xagg)
    ↓
Script 5: Climate panel + yield data → regression dataset
```

**Total Runtime:** ~30 minutes (2013-2023, 64 districts, 3 seasons)  
**Output:** 2,065 observations × 66 variables (regression-ready panel)

---

## Script 1: ERA5 Hourly to Daily Conversion

### Purpose
Convert ERA5 hourly reanalysis data to daily temperature and precipitation for Bangladesh.

### Input Data
- **ERA5 hourly temperature** (`t2m_YYYY.nc`): 2-meter air temperature (Kelvin)
- **ERA5 hourly precipitation** (`tp_YYYY.nc`): Total precipitation (meters)
- **Spatial coverage:** Bangladesh bbox (20.5-26.7°N, 88.0-92.7°E)
- **Temporal coverage:** 2013-2023 (11 years)
- **Resolution:** 0.25° × 0.25° (~28 km)

### Processing Steps

1. **Temporal Subset:** Extract hourly data for Bangladesh bounding box
2. **Temperature Aggregation (hourly → daily):**
   - `tas` (daily mean) = mean of 24 hourly values
   - `tasmin` (daily minimum) = min of 24 hourly values
   - `tasmax` (daily maximum) = max of 24 hourly values
3. **Unit Conversion:** Kelvin → Celsius (subtract 273.15)
4. **Precipitation Aggregation (hourly → daily):**
   - Convert meters → millimeters (multiply by 1000)
   - Sum 24 hourly accumulations → daily total (mm/day)
5. **Output:** 4 files per year (tas, tasmin, tasmax, pr) × 11 years = **44 NetCDF files**

### Output Format
- **Files:** `tas_bangladesh_ERA5_YYYY.nc`, `tasmin_bangladesh_ERA5_YYYY.nc`, `tasmax_bangladesh_ERA5_YYYY.nc`, `pr_bangladesh_ERA5_YYYY.nc`
- **Dimensions:** time (365 days), lat (49), lon (38)
- **Variables:** tas/tasmin/tasmax (°C), pr (mm/day)

### Data Citation
**Hersbach et al.** (2023). ERA5 hourly data on single levels from 1940 to present. Copernicus Climate Change Service (C3S) Climate Data Store (CDS). https://doi.org/10.24381/cds.adbb2d47

---

## Script 2: Daily Climate Indices Computation

### Purpose
Calculate daily degree-day indices (EDD, GDD, HDD) and precipitation features using sinusoidal approximation of within-day temperature distribution.

### Methodology: Sinusoidal Degree-Days

Following **Schlenker & Roberts (2009)** and **Snyder (1985)**, we approximate the within-day temperature distribution as a sinusoid between daily minimum (`tasmin`) and maximum (`tasmax`). This allows analytical integration to compute time-weighted temperature exposure above/below thresholds.

#### Temperature Distribution Parameters
- **Daily mean:** M = (tasmax + tasmin) / 2
- **Daily half-range:** W = (tasmax - tasmin) / 2
- **Temperature at time θ:** T(θ) = M + W · sin(θ), where θ ∈ [-π, π]

### Indices Computed

#### 1. Exceedance Degree-Days (EDD)
**Definition:** Time-weighted integral of temperature ABOVE threshold.

**Formula:**
```
EDD(T_threshold) = (1/π) ∫[θ₁, π] [T(θ) - T_threshold] dθ
```
where θ₁ = arccos((T_threshold - M) / W)

**Thresholds:** 0, 4, 8, 12, 28, 30, 32, 35°C

**Interpretation:** 
- EDD_28: Heat stress during grain filling (Boro: March-April)
- EDD_30: Optimal temperature threshold for rice (Van Oort & Zwart, 2018)
- EDD_32, EDD_35: Extreme heat damage

**Case Handling:**
- If tasmax ≤ threshold: EDD = 0
- If tasmin ≥ threshold: EDD = M - threshold (entire day above)
- If threshold crosses day: Sinusoidal integral (6 cases handled)

#### 2. Growing Degree-Days (GDD)
**Definition:** Time-weighted integral of temperature between base and cap.

**Formula:**
```
GDD(T_base, T_cap) = (1/π) ∫[θ₁, θ₂] [T(θ) - T_base] dθ + (T_cap - T_base) · (θ₂ - π)/π
```
where:
- θ₁ = arccos((T_base - M) / W)
- θ₂ = arccos((T_cap - M) / W)

**Specifications (15 variants):**
- **Bases:** 8, 10, 12, 15, 20°C (rice range: 8-15°C typical)
- **Caps:** 30, 32, 35°C
- **Primary:** GDD_10_30 (standard for rice)

**Case Handling:**
- If tasmax ≤ base: GDD = 0
- If tasmin ≥ cap: GDD = cap - base
- If base < tasmin < tasmax < cap: GDD = M - base (simple linear)
- If day crosses base and/or cap: Sinusoidal integral (6 cases)

#### 3. Harmful Degree-Days (HDD)
**Definition:** Time-weighted integral of temperature BELOW threshold (cold damage).

**Formula:**
```
HDD(T_threshold) = (1/π) ∫[-π, θ₁] [T_threshold - T(θ)] dθ
```

**Thresholds:** 10, 15°C (Boro susceptible to cold in December-January)

#### 4. Precipitation Features
- **pr1:** Daily precipitation (mm/day) - direct from Script 1
- **pr2:** Daily precipitation squared (mm/day)² - for nonlinear precipitation effects

### Output
- **Files:** `indices_bangladesh_ERA5_YYYY.nc` (11 files, 2013-2023)
- **Dimensions:** time (365 days), lat (49), lon (38)
- **Variables:** 27 indices per day
  - 8 EDD variants
  - 15 GDD variants
  - 2 HDD variants
  - 2 precipitation features

### Methodological References
**Schlenker, W., & Roberts, M.J.** (2009). Nonlinear temperature effects indicate severe damages to US crop yields under climate change. *PNAS*, 106(37), 15594-15598. https://doi.org/10.1073/pnas.0906865106

**Snyder, R.L.** (1985). Hand calculating degree days. *Agricultural and Forest Meteorology*, 35, 353-358. https://doi.org/10.1016/0168-1923(85)90095-4

---

## Script 3: Bangladesh Rice Calendar Windows

### Purpose
Convert Bangladesh rice crop calendars (Boro, Aus, Aman NetCDFs) into temporal aggregation windows for each district.

### Input Data
- **Rice Calendars:** `Rice.Boro.crop.calendar.bgd.nc`, `Rice.Aus.crop.calendar.bgd.nc`, `Rice.Aman.crop.calendar.bgd.nc`
  - Resolution: 5 arcminutes (~10 km)
  - Variables: plant.start, plant.end, harvest.start, harvest.end (DOY 1-365)
  - Coverage: ~35-40% of Bangladesh land area (rice-growing regions)
- **Shapefile:** GADM v4.1 Bangladesh ADM2 (64 districts)

### Calendar Window Methodology (from CI2026)

#### Step 1: Extract Calendar DOYs per District
For each district polygon:
1. Extract polygon centroid coordinates
2. Find nearest calendar pixel (nearest-neighbor)
3. Read DOY values: plant_start, plant_end, harvest_start, harvest_end

#### Step 2: Convert to 2-Year Concatenated Indices
To handle year-wrapping seasons (e.g., Boro: Dec → May), we use a 2-year time axis (0-729 indices):

**Algorithm:**
```
T = 365
offset = 365  # Start second year at index 365

# Reference point: harvest end
harvest_end_day = harvest_end_doy + offset

# Ensure harvest_start ≤ harvest_end
harvest_start_day = harvest_start_doy + offset
if harvest_start_day > harvest_end_day:
    harvest_start_day -= T

# Ensure plant_end ≤ harvest_end
plant_end_day = plant_end_doy + offset
if plant_end_day > harvest_end_day:
    plant_end_day -= T

# Ensure plant_start ≤ plant_end
plant_start_day = plant_start_doy + offset
if plant_start_day > plant_end_day:
    plant_start_day -= T

# Between period
between_start_day = plant_end_day + 1
between_end_day = harvest_start_day - 1
```

#### Step 3: Define Three Periods
- **Plant:** Planting to end of establishment
- **Between:** Vegetative growth and flowering
- **Harvest:** Grain filling to harvest

**Output Format (per season):**
```python
windows = {
    "plant": {
        "start": array(64),  # 0-based start indices
        "end": array(64),    # 0-based end indices
        "days": array(64)    # day counts
    },
    "between": {...},
    "harvest": {...}
}
```

### Output
- **Files:** 3 pickle files
  - `calendar_windows_BGD_Rice.Boro.pkl`
  - `calendar_windows_BGD_Rice.Aus.pkl`
  - `calendar_windows_BGD_Rice.Aman.pkl`
- **Contents:** Temporal window arrays for 63-64 districts (1 district may have no rice)

### Calendar Validation
- **Boro:** Median plant_start = DOY 335 (Dec 1, nationally uniform, literature-based)
- **Aus:** Median plant_start = DOY 94 (Apr 4, ±31 days spatial variation)
- **Aman:** Median plant_start = DOY 207 (Jul 26, ±36 days spatial variation)

---

## Script 4: Seasonal Aggregation with xagg

### Purpose
Aggregate daily climate indices to seasonal district-level values using harvested-area weighted spatial aggregation and calendar-window temporal aggregation.

### Two-Stage Aggregation

#### Stage 1: Spatial Aggregation (Grid → District)

**Method:** `xagg` (xarray aggregation) with harvested-area weights

**Steps:**

1. **Load Harvested Area Rasters:**
   - **Irrigated rice (IRC):** `ANNUAL_AREA_HARVESTED_IRC_CROP3_HA.ASC` (for Boro)
   - **Rainfed rice (RFC):** `ANNUAL_AREA_HARVESTED_RFC_CROP3_HA.ASC` (for Aus, Aman)
   - Source: Monfreda et al. (2008), 5-arcminute resolution

2. **Create Pixel Overlap Weights:**
   ```python
   # Regrid harvested area to ERA5 grid (0.25°)
   regridder = xesmf.Regridder(harvest_grid, climate_grid, "bilinear")
   harvest_regridded = regridder(harvest_area)
   
   # Create xagg weights (polygon × pixel overlaps weighted by harvest area)
   weights = xagg.pixel_overlaps(climate_grid, district_polygons, 
                                  weights=harvest_regridded)
   ```

3. **Aggregate Climate Indices:**
   ```python
   # Area-weighted mean across pixels → district value
   district_daily = xagg.aggregate(climate_indices, weights)
   ```

**Result:** Daily climate indices (365 days) for each district, weighted by rice harvested area

**Harvested Area Purpose:**
- **Boro (IRC):** Pixels with more irrigated rice area contribute more to district aggregate
- **Aus/Aman (RFC):** Pixels with more rainfed rice area contribute more

#### Stage 2: Temporal Aggregation (Daily → Seasonal)

**Method:** Vectorized calendar-window masking

**Steps:**

1. **Concatenate Time Axis:** Create 2-year daily dataset (730 days = 365 × 2)
2. **Apply Calendar Windows:** For each period (plant, between, harvest):
   ```python
   # Create time × district mask
   mask = (time_index >= window_start) & (time_index < window_end)
   
   # Sum daily indices over masked days
   seasonal_value = daily_indices.where(mask, 0).sum(dim="time")
   ```

3. **Output:** Seasonal sums for each district × period × climate index

**Example:**
- District X, Boro season, plant period (DOY 335-365 + 1-30 = 61 days):
  - EDD_30 = sum of 61 daily EDD_30 values
  - GDD_10_30 = sum of 61 daily GDD_10_30 values
  - pr1 = sum of 61 daily precipitation values

### Season-Irrigation Mapping
| Season | Irrigation | Harvested Area Grid |
|--------|-----------|---------------------|
| Boro | Irrigated | IRC (fully irrigated dry season) |
| Aus | Rainfed | RFC (pre-monsoon) |
| Aman | Rainfed | RFC (monsoon) |

### Output
- **Per-year files:** `panel_bangladesh_rice_YYYY.parquet` (11 files)
- **Full panel:** `panel_bangladesh_rice_full.parquet`
- **Dimensions:** 6,336 rows = 64 districts × 3 seasons × 3 periods × 11 years
- **Variables:** 36 columns (27 climate indices + metadata)

### Methodological References
**Doelman, J.C. et al.** (CI2026, in preparation). Global climate-crop exposure aggregation framework.

**Monfreda, C., Ramankutty, N., & Foley, J.A.** (2008). Farming the planet: 2. Geographic distribution of crop areas, yields, physiological types, and net primary production in the year 2000. *Global Biogeochemical Cycles*, 22(1), GB1022. https://doi.org/10.1029/2007GB002947

---

## Script 5: Merge with Yield Data

### Purpose
Create final regression-ready dataset by merging climate panel with BBS rice yield data and adding robustness variables.

### Input Data

#### Climate Panel (from Script 4)
- 6,336 rows (districts × years × seasons × periods)
- 36 variables (climate indices + metadata)

#### Yield Data
- **Source:** Bangladesh Bureau of Statistics (BBS) Agricultural Yearbook
- **File:** `df_rice.csv`
- **Coverage:** 1970-2023, 64 districts, 3 seasons
- **Variables:**
  - District (name)
  - Year
  - Crop_type (Aman/Aus/Boro)
  - Area (hectares)
  - Production (metric tons)

### Processing Steps

#### 1. Harmonize District Names
**Issue:** Shapefile and BBS use different transliterations

**Mapping (13 corrections):**
```python
{
    'Bandarban': 'Banderban',
    'Bogra': 'Bogura',
    'Brahamanbaria': 'Brahmmanbaria',
    'Chittagong': 'Chattogram',
    'Comilla': 'Cumilla',
    "Cox'S Bazar": "Cox's bazar",
    'Jessore': 'Jashore',
    'Jhalokati': 'Jhalokathi',
    'Khagrachhari': 'Khagrachari',
    'Maulvibazar': 'Maulavi Bazar',
    'Nawabganj': 'Chapai Nawabganj',
    'Netrakona': 'Netrokona',
    'Pirojpur': 'Perojpur',
}
```

#### 2. Harmonize Season Names
```python
{
    'Rice.Boro': 'Boro',
    'Rice.Aus': 'Aus',
    'Rice.Aman': 'Aman',
}
```

#### 3. Aggregate Periods
Sum climate indices across plant + between + harvest → full growing season

```python
seasonal = climate.groupby(['district', 'year', 'season']).agg({
    'edd_*': 'sum',    # Sum degree-days across entire season
    'gdd_*': 'sum',
    'hdd_*': 'sum',
    'pr1': 'sum',      # Total seasonal precipitation
    'pr2': 'sum',
    'days': 'sum'      # Total growing days
})
```

#### 4. Calculate Yield Metrics
```python
yield_per_ha = production / area  # Metric tons per hectare
log_yield = log(yield_per_ha)     # Log transformation for regression
```

#### 5. Merge on (District, Year, Season)
```python
final = climate.merge(yield, on=['district', 'year', 'season'], how='inner')
```

**Result:** 2,065 matched observations (2013-2023)

#### 6. Create Robustness Variables

**First Differences (for fixed effects models):**
```python
diff_log_yield = log_yield[t] - log_yield[t-1]
diff_gdd_10_35 = gdd_10_35[t] - gdd_10_35[t-1]
diff_edd_* = edd_*[t] - edd_*[t-1]
diff_pr1 = pr1[t] - pr1[t-1]
```

**Interactions (nonlinear effects):**
```python
gdd_edd_interaction = gdd_10_35 × edd_30
gdd_pr1_interaction = gdd_10_35 × pr1
edd_pr1_interaction = edd_30 × pr1
```

**Polynomials (flexible functional forms):**
```python
gdd_10_35_sq = gdd_10_35²
edd_28_sq = edd_28²
edd_30_sq = edd_30²
edd_32_sq = edd_32²
edd_35_sq = edd_35²
```

**Quartile Indicators (distributional analysis):**
```python
gdd_quartile = quartile_rank(gdd_10_35)  # 1-4
edd_*_quartile = quartile_rank(edd_*)
pr1_quartile = quartile_rank(pr1)
```

**Binary Indicators (threshold effects):**
```python
high_gdd = (gdd_10_35 > median(gdd_10_35))  # 0/1
high_edd_* = (edd_* > median(edd_*))
high_pr1 = (pr1 > median(pr1))
```

### Final Dataset Structure

**Dimensions:** 2,065 observations × 66 variables

**Observation Unit:** District × Year × Season (2013-2023)

**Coverage:**
- 64 districts
- 11 years
- 3 seasons (Boro: 689 obs, Aman: 688 obs, Aus: 688 obs)

**Variable Categories:**
1. **Identifiers (5):** district, year, season, adm_code, country
2. **Yield Variables (4):** area, production, yield_per_ha, log_yield
3. **Climate Indices (27):** edd_*, gdd_*, hdd_*, pr1, pr2
4. **Robustness Variables (~30):**
   - First differences (8)
   - Interactions (3)
   - Polynomials (5)
   - Quartiles (6)
   - Binary indicators (6)

**Missing Values:**
- Yield: 9 observations (division by zero: area = 0)
- Lags: 189-201 observations (first year per district × season)

### Output Files
- **CSV:** `bangladesh_rice_regression_panel.csv` (human-readable)
- **Parquet:** `bangladesh_rice_regression_panel.parquet` (efficient storage)
- **Location:** `data/Regression_data/`

---

## Summary Statistics

### Climate Index Distributions (Full Season, 2013-2023)

| Variable | Mean | Std | Min | Median | Max | Unit |
|----------|------|-----|-----|--------|-----|------|
| **edd_28** | 193.4 | 127.8 | 0.0 | 165.1 | 610.8 | degree-days |
| **edd_30** | 85.9 | 60.7 | 0.0 | 73.8 | 334.2 | degree-days |
| **edd_32** | 28.1 | 23.9 | 0.0 | 22.3 | 161.1 | degree-days |
| **gdd_10_30** | 2853.3 | 1334.6 | 0.0 | 2526.7 | 5690.0 | degree-days |
| **pr1** | 909.8 | 766.7 | 0.0 | 727.4 | 4564.1 | mm |
| **log_yield** | 0.103 | 0.337 | -1.14 | 0.090 | 0.771 | log(MT/ha) |

### Season-Specific Patterns

**Boro (Dry Season, Irrigated):**
- Low precipitation (mean: ~200-300 mm)
- Moderate heat stress (edd_30: ~50-80 degree-days)
- High GDD accumulation (controlled irrigation → optimal growth)

**Aus (Pre-Monsoon, Rainfed):**
- Moderate precipitation (mean: ~500-800 mm)
- High heat stress (edd_30: ~100-120 degree-days, hottest season)
- Variable GDD (rainfall-dependent)

**Aman (Monsoon, Rainfed):**
- High precipitation (mean: ~1200-1600 mm)
- Low heat stress (edd_30: ~60-80 degree-days)
- High GDD accumulation (adequate rainfall)

---

## Validation & Quality Checks

### 1. Temperature Index Validation
- **Sinusoidal GDD caps correctly:** GDD_10_30 max = 20.0 (= cap 30 - base 10) ✓
- **EDD increases with threshold:** edd_35 < edd_32 < edd_30 < edd_28 ✓
- **Physical plausibility:** All degree-days ≥ 0 ✓

### 2. Precipitation Validation
- **Units correct:** pr1 max = 318 mm/day (extreme rainfall event) ✓
- **Squared term:** pr2 = pr1² (318² ≈ 101,124) ✓
- **Seasonal totals reasonable:** Boro: 200-500 mm, Aman: 1000-2000 mm ✓

### 3. Spatial Aggregation Validation
- **District count:** 64 (matches GADM shapefile) ✓
- **Coverage:** All 64 districts have climate data ✓
- **Harvest area weights:** IRC for Boro, RFC for Aus/Aman ✓

### 4. Temporal Aggregation Validation
- **Calendar windows:** 63-64 districts (1 may lack rice) ✓
- **Growing season days:** Boro: 120-150 days, Aus: 90-120 days, Aman: 120-150 days ✓
- **Year-wrapping:** Boro (Dec-May) handled correctly ✓

### 5. Merge Validation
- **Match rate:** 2,065 / 2,098 BBS observations (98.4%) ✓
- **Unmatched:** 33 observations (likely districts with no rice calendar coverage) ✓
- **Missing yields:** 9 / 2,065 (0.4%, division by zero) ✓

---

## Advantages Over Previous Approaches

### 1. Spatial Precision
- **Previous (R pipeline):** Unweighted mean across district pixels
- **Current (CI2026 method):** Harvested-area weighted aggregation
- **Benefit:** Districts with concentrated rice production properly weighted

### 2. Temporal Precision
- **Previous:** Hardcoded DOY windows (uniform across districts)
- **Current:** Spatially-varying calendars from MODIS (8-year median)
- **Benefit:** Captures regional variation in planting dates (±30 days)

### 3. Degree-Day Calculation
- **Previous:** Simple linear GDD = max(0, min(tas - base, cap - base))
- **Current:** Sinusoidal integration (Schlenker & Roberts 2009)
- **Benefit:** Accounts for within-day temperature distribution, more accurate heat stress

### 4. Season Differentiation
- **Previous:** Single rice season or manual DOY splits
- **Current:** Three distinct calendars (Boro/Aus/Aman) with irrigation types
- **Benefit:** Captures season-specific climate vulnerabilities

---

## Limitations & Future Work

### Current Limitations

1. **Calendar Spatial Resolution:** 5-arcminute (~10 km) may miss sub-district variation
2. **Aus Coverage Uncertainty:** Low production share (9%) → some districts may have incomplete calendars
3. **Irrigation Assumption:** Boro = 100% IRC, Aus/Aman = 100% RFC (reality: some irrigated Aman)
4. **Climate Grid Resolution:** ERA5 0.25° (~28 km) limits capture of microclimatic variation
5. **Static Harvested Area:** Monfreda (2008) weights represent year 2000, not 2013-2023 dynamics

### Potential Extensions

1. **Add More Crops:** Extend to wheat, potato, jute, pulses (requires new calendars)
2. **Finer Irrigation Mapping:** Use MIRCA2000 monthly irrigated/rainfed fractions
3. **Climate Model Ensemble:** Replicate with CMIP6 models for future projections
4. **Sub-District Analysis:** Use upazila (sub-district) boundaries for finer spatial detail
5. **Phenology Tracking:** Annual calendars instead of 8-year median (requires extensive validation)

---

## Software & Dependencies

### Core Libraries
- **Python 3.10+**
- **xarray 2023.1.0+** (multi-dimensional arrays)
- **pandas 2.0+** (dataframes)
- **geopandas 0.14+** (spatial operations)
- **xagg 0.3+** (spatial aggregation)
- **xesmf 0.8+** (regridding)
- **rasterio 1.3+** (raster I/O)
- **netCDF4 1.6+** (NetCDF I/O)

### Computational Requirements
- **RAM:** 16 GB minimum (32 GB recommended)
- **Storage:** ~50 GB (raw ERA5 + intermediate outputs)
- **Runtime:** ~30 minutes (full pipeline, 2013-2023)

### Code Availability
All scripts (01-05) available at:
```
wage_model/code/
├── 01_era5_hourly_to_daily_bangladesh.py
├── 02_compute_climate_indices_bangladesh.py
├── 03_bangladesh_calendar_windows.py
├── 04_aggregate_seasonal_bangladesh.py
└── 05_merge_yield_bangladesh.py
```

---

## Data Availability

### Input Data
- **ERA5:** Copernicus Climate Data Store (requires free account)
- **GADM v4.1:** https://gadm.org (public)
- **Monfreda et al. (2008):** Earthstat (public)
- **BBS Rice Yield:** Bangladesh Bureau of Statistics (public)
- **Bangladesh Rice Calendars:** Custom (available upon request)

### Output Data
- **Final Panel:** `bangladesh_rice_regression_panel.csv` (2,065 obs × 66 vars)
- **Intermediate Files:** Available upon request

---

## Citation

If using this methodology or dataset, please cite:

**Tahmid, T.** (2026). Bangladesh Rice Climate Indices Pipeline: Seasonal aggregation of ERA5 climate data using empirical crop calendars and harvested-area weighted spatial aggregation. University of Delaware. *[Dataset/Code Repository]*

**And core methodological papers:**

**Schlenker, W., & Roberts, M.J.** (2009). Nonlinear temperature effects indicate severe damages to US crop yields under climate change. *PNAS*, 106(37), 15594-15598.

**Monfreda, C., Ramankutty, N., & Foley, J.A.** (2008). Farming the planet: 2. Geographic distribution of crop areas, yields, physiological types, and net primary production in the year 2000. *Global Biogeochemical Cycles*, 22(1), GB1022.

---

## Author

**Taky Tahmid**  
University of Delaware

---

**Document Version:** 1.0  
**Last Updated:** April 5, 2026  
**Word Count:** ~4,200 words
