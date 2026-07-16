# Pipeline — ERA5 Climate Processing

Run once. Data already processed. Do not re-run unless raw data changes.

## Run order
1. `01_era5_download.py`        → ERA5 hourly → daily tas/pr
2. `02_climate_indices.py`      → GDD/EDD/HDD indices (Schlenker-Roberts sinusoidal)
3. `03_crop_calendar.py`        → rice calendars → district crop windows
4. `04_spatial_aggregation.py`  → area-weighted aggregation to district level
5. `05_merge_panel.py`          → merge climate + BBS yield → regression panel
6. `06_temperature_bins.py`     → 1°C temperature bin exposure

## Output
- `data/Regression_data/bangladesh_rice_regression_panel.csv` ← Stage 1 input
- `data/Regression_data/bangladesh_rice_temperature_bins.parquet`

## Notes
- GDD/EDD computed using Schlenker-Roberts sinusoidal interpolation
- Crop calendars from IRRI + local expert knowledge (see `03_crop_calendar.py`)
- Spatial aggregation weighted by harvested area (GAEZ raster)
