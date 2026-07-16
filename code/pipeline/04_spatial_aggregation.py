#!/usr/bin/env python3
"""
================================================================================
SCRIPT 4: Seasonal Aggregation with xagg (Bangladesh Rice)
================================================================================

Purpose: Aggregate daily climate indices to seasonal by ADM2 district

Input:
  - data/climate_data/era5_indices/netcdf/indices_bangladesh_ERA5_YYYY.nc
  - data/calendar_windows/calendar_windows_BGD_Rice.*.pkl
  - data/climate_data/gadm41_BGD_shp/gadm41_BGD_2.shp
  - Harvested area grids (IRC/RFC for rice)

Output:
  - data/seasonal_panels/panel_bangladesh_rice_YYYY.parquet (per year)
  - data/seasonal_panels/panel_bangladesh_rice_full.parquet (merged)

Method (CI2026):
  1. Create xagg weights using harvested area rasters
  2. Aggregate daily indices (grid → ADM2) with area-weighted mean
  3. Apply calendar windows (daily → seasonal) with vectorized masking
  4. Output: districts × seasons × indices

Test mode: Set TEST_YEAR to process single year only

Author: Taky Tahmid
Date: 2026-04-05
================================================================================
"""

from pathlib import Path
import os
import pickle
import time
import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import rioxarray
import xagg as xa
import xarray as xr
import xesmf as xe

xa.set_options(nan_to_zero_regridding=False)

# =============================================================================
# CONFIGURATION
# =============================================================================

COUNTRY_NAME = "bangladesh"
COUNTRY_TAG = "BGD"
ADM_LEVEL = 2

# TEST MODE: Set to year (e.g., 2013) to test single year, or None for all years
TEST_YEAR = None  # Change to None after successful test

YEAR_START = 2013
YEAR_END = 2023

# Paths
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent  # code/pipeline/ -> code/ -> project root

SHAPEFILE = PROJECT_ROOT / "data/climate_data/gadm41_BGD_shp/gadm41_BGD_2.shp"
INDICES_DIR = PROJECT_ROOT / "data/climate_data/era5_indices/netcdf"
WINDOWS_DIR = PROJECT_ROOT / "data/calendar_windows"
OUTPUT_DIR = PROJECT_ROOT / "data/seasonal_panels"

# Harvested area grids are large and are not stored in this public repository.
# Set HARVEST_AREA_GRID_DIR to the local directory containing these grids.
HARVEST_BASE = Path(os.environ.get(
    "HARVEST_AREA_GRID_DIR",
    PROJECT_ROOT / "data/harvested_area_grids"
))

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Rice crop code (GAEZ/Monfreda)
RICE_CODE = 3

# Rice seasons and irrigation types
RICE_SEASONS = {
    "Rice.Boro": "irrigated",    # Dry season, fully irrigated
    "Rice.Aus": "rainfed",        # Pre-monsoon, mostly rainfed
    "Rice.Aman": "rainfed",       # Monsoon, rainfed
}

# Climate indices to aggregate
CLIMATE_VARS = [
    "edd_0", "edd_4", "edd_8", "edd_12", "edd_28", "edd_30", "edd_32", "edd_35",
    "gdd_8_30", "gdd_8_32", "gdd_8_35", "gdd_10_30", "gdd_10_32", "gdd_10_35",
    "gdd_12_30", "gdd_12_32", "gdd_12_35", "gdd_15_30", "gdd_15_32", "gdd_15_35",
    "gdd_20_30", "gdd_20_32", "gdd_20_35",
    "hdd_10", "hdd_15",
    "pr1", "pr2"
]

# =============================================================================
# XAGG WEIGHT CREATION (CI2026 methodology)
# =============================================================================

def create_weights(harvestpath: Path, gdf: gpd.GeoDataFrame, grid_da: xr.DataArray) :
    """
    Create xagg pixel overlap weights using harvested area raster.
    
    Args:
        harvestpath: Path to harvested area .ASC file
        gdf: GeoDataFrame with district polygons
        grid_da: Sample grid from climate data (for regridding target)
    
    Returns:
        xagg PixelOverlaps object
    """
    
    print(f"    Creating weights from {harvestpath.name}...")
    
    # Load harvested area raster
    with rasterio.open(harvestpath) as src:
        data_array = rioxarray.open_rasterio(src)
        weights = data_array.squeeze().rename({"y": "latitude", "x": "longitude"})
        
        latitudes = weights["latitude"].values
        longitudes = weights["longitude"].values
        
        # Create bounds for regridding
        delta_lat = np.abs(latitudes[1] - latitudes[0]) / 2
        delta_lon = np.abs(longitudes[1] - longitudes[0]) / 2
        
        lat_bounds = np.array([[lat - delta_lat, lat + delta_lat] for lat in latitudes])
        lon_bounds = np.array([[lon - delta_lon, lon + delta_lon] for lon in longitudes])
        
        lat_bounds_da = xr.DataArray(
            lat_bounds,
            dims=["latitude", "bnds"],
            coords={"latitude": latitudes, "bnds": [0, 1]},
        )
        lon_bounds_da = xr.DataArray(
            lon_bounds,
            dims=["longitude", "bnds"],
            coords={"longitude": longitudes, "bnds": [0, 1]},
        )
        
        weights2 = weights.to_dataset(name="weights")
        weights2["lat_bounds"] = lat_bounds_da
        weights2["lon_bounds"] = lon_bounds_da
        weights2["latitude"].attrs["bounds"] = "lat_bounds"
        weights2["longitude"].attrs["bounds"] = "lon_bounds"
        
        weights3 = xr.decode_cf(weights2)
    
    # Regrid harvested area to climate grid
    regridder = xe.Regridder(weights3, grid_da, "bilinear", reuse_weights=False)
    weights4 = regridder(weights3)
    
    # Create pixel overlaps
    return xa.pixel_overlaps(grid_da, gdf, weights=weights4.weights)


# =============================================================================
# CALENDAR WINDOW AGGREGATION (CI2026 vectorized method)
# =============================================================================

def apply_calendar_windows_vectorized(ds: xr.Dataset, windows: dict, varlist: list) -> xr.Dataset:
    """
    Apply calendar windows to daily data (vectorized).
    
    Sums daily indices over plant/between/harvest periods for each district.
    Uses 2-year concatenated time axis to handle year-wrapping seasons.
    
    Args:
        ds: Daily dataset with climate indices
        windows: Calendar windows dict from pickle
        varlist: List of variables to aggregate
    
    Returns:
        Dataset with seasonal aggregates (dims: poly_idx, period)
    """
    
    # Concatenate dataset with itself (2-year axis)
    ds2 = xr.concat([ds, ds], dim="time")
    
    period_order = ["plant", "between", "harvest"]
    
    # Create time index array
    t_idx = xr.DataArray(
        np.arange(ds2.sizes["time"]),
        dims=("time",),
        coords={"time": ds2["time"]},
    )
    
    all_periods = []
    
    for period in period_order:
        start_da = xr.DataArray(windows[period]["start"], dims=("poly_idx",))
        end_da = xr.DataArray(windows[period]["end"], dims=("poly_idx",))
        days_da = xr.DataArray(windows[period]["days"], dims=("poly_idx",))
        
        # Broadcast time × poly_idx
        tt, ss = xr.broadcast(t_idx, start_da)
        _, ee = xr.broadcast(t_idx, end_da)
        
        # Create mask (time in window for each poly)
        mask = (tt >= ss) & (tt < ee)
        
        # Sum over time dimension
        summed = ds2[varlist].where(mask, other=0.0).sum(dim="time", skipna=True)
        
        summed = summed.assign_coords(period=period)
        summed["days"] = days_da
        
        all_periods.append(summed)
    
    # Concatenate periods
    out = xr.concat(all_periods, dim="period")
    out = out.assign_coords(period=("period", period_order))
    
    return out


# =============================================================================
# PROCESS ONE YEAR
# =============================================================================

def process_year(year: int, gdf: gpd.GeoDataFrame) -> pd.DataFrame:
    """
    Process one year: aggregate daily indices to seasonal panel.
    
    Args:
        year: Year to process
        gdf: GeoDataFrame with district polygons
    
    Returns:
        DataFrame with seasonal panel for this year
    """
    
    print(f"\n[{year}]")
    
    # Check for existing output
    out_file = OUTPUT_DIR / f"panel_bangladesh_rice_{year}.parquet"
    if out_file.exists():
        print(f"  ✓ Already exists, skipping")
        return pd.read_parquet(out_file)
    
    # Load indices file
    indices_file = INDICES_DIR / f"indices_{COUNTRY_NAME}_ERA5_{year}.nc"
    if not indices_file.exists():
        print(f"  ✗ Indices file not found: {indices_file.name}")
        return None
    
    print(f"  Loading indices...")
    start = time.time()
    ds = xr.open_dataset(indices_file)
    print(f"    → {ds.sizes['time']} days, {len(CLIMATE_VARS)} variables ({time.time()-start:.1f}s)")
    
    # Get sample grid for weights
    grid_for_weights = ds[CLIMATE_VARS[0]].isel(time=0, drop=True)
    
    # Get ADM codes/names
    adm_code_col = [c for c in gdf.columns if c.startswith("GID_")][0]
    adm_name_col = "NAME_2"
    adm_codes = gdf[adm_code_col].values
    adm_names = gdf[adm_name_col].values
    
    panel_dfs = []
    
    # Process each season × irrigation combination
    for season, irrigation_label in RICE_SEASONS.items():
        
        print(f"\n  [{season} - {irrigation_label}]")
        
        # Determine IRC/RFC
        infix = "IRC" if irrigation_label == "irrigated" else "RFC"
        
        # Create or load weights
        weightpath = OUTPUT_DIR / f"weights_rice_{infix}_{COUNTRY_TAG}.wm"
        
        if weightpath.exists():
            print(f"    Loading existing weights...")
            with open(weightpath, "rb") as fp:
                weightmap = pickle.load(fp)
        else:
            harvestpath = HARVEST_BASE / f"ANNUAL_AREA_HARVESTED_{infix}_CROP{RICE_CODE}_HA.ASC"
            if not harvestpath.exists():
                print(f"    ✗ Harvested area file not found: {harvestpath.name}")
                continue
            weightmap = create_weights(harvestpath, gdf, grid_for_weights)
            with open(weightpath, "wb") as fp:
                pickle.dump(weightmap, fp)
            print(f"    → Saved weights to {weightpath.name}")
        
        # Spatial aggregation (grid → districts)
        print(f"    Aggregating spatially (xagg)...")
        start = time.time()
        aggregated = xa.aggregate(ds[CLIMATE_VARS], weightmap).to_dataset()
        
        # Rename feature → poly_idx if needed
        if "feature" in aggregated.dims and "poly_idx" not in aggregated.dims:
            aggregated = aggregated.rename({"feature": "poly_idx"})
        
        if "poly_idx" not in aggregated.dims:
            print(f"    ✗ No poly_idx dimension after aggregation")
            continue
        
        print(f"    → Done ({time.time()-start:.1f}s)")
        
        # Load calendar windows
        windows_path = WINDOWS_DIR / f"calendar_windows_{COUNTRY_TAG}_{season}.pkl"
        if not windows_path.exists():
            print(f"    ✗ Calendar windows not found: {windows_path.name}")
            continue
        
        with open(windows_path, "rb") as f:
            cal_windows = pickle.load(f)
        
        # Temporal aggregation (daily → seasonal)
        print(f"    Aggregating temporally (calendar windows)...")
        start = time.time()
        seasonal = apply_calendar_windows_vectorized(aggregated[CLIMATE_VARS], cal_windows, CLIMATE_VARS)
        print(f"    → Done ({time.time()-start:.1f}s)")
        
        # Convert to DataFrame
        seasonal = seasonal.reset_coords(drop=True)
        df = seasonal.to_dataframe().reset_index()
        
        if df.empty:
            continue
        
        # Add metadata
        df["adm_code"] = adm_codes[df["poly_idx"].values]
        df["adm_name"] = adm_names[df["poly_idx"].values]
        df["year"] = year
        df["season"] = season
        df["irrigation"] = irrigation_label
        df["country"] = COUNTRY_TAG
        
        panel_dfs.append(df)
        print(f"    → {len(df)} rows (districts × periods)")
    
    if not panel_dfs:
        print(f"\n  ✗ No data aggregated for {year}")
        return None
    
    # Combine all seasons
    panel = pd.concat(panel_dfs, ignore_index=True)
    
    # Save per-year file
    panel.to_parquet(out_file)
    print(f"\n  ✓ Saved {len(panel)} rows to {out_file.name}")
    
    return panel


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("\n" + "="*80)
    print(f" Seasonal Aggregation with xagg ({COUNTRY_TAG} Rice)")
    print("="*80)
    
    if TEST_YEAR is not None:
        print(f"⚠️  TEST MODE: Processing year {TEST_YEAR} only")
        print(f"   Set TEST_YEAR = None to process all years ({YEAR_START}-{YEAR_END})")
    else:
        print(f"Years: {YEAR_START}-{YEAR_END}")
    
    print(f"Shapefile: {SHAPEFILE}")
    print(f"Indices:   {INDICES_DIR}")
    print(f"Windows:   {WINDOWS_DIR}")
    print(f"Harvest:   {HARVEST_BASE}")
    print(f"Output:    {OUTPUT_DIR}")
    print()
    
    # Load shapefile
    print("Loading shapefile...")
    gdf = gpd.read_file(SHAPEFILE)
    gdf = gdf.to_crs("EPSG:4326")
    gdf = gdf.reset_index(drop=True)
    print(f"  → {len(gdf)} districts (ADM{ADM_LEVEL})")
    
    overall_start = time.time()
    
    # Determine years to process
    if TEST_YEAR is not None:
        years = [TEST_YEAR]
    else:
        years = range(YEAR_START, YEAR_END + 1)
    
    # Process years
    all_panels = []
    for year in years:
        panel = process_year(year, gdf)
        if panel is not None:
            all_panels.append(panel)
    
    # Merge all years
    if all_panels:
        print("\n" + "="*80)
        print(" Merging all years...")
        print("="*80)
        
        full_panel = pd.concat(all_panels, ignore_index=True)
        full_output = OUTPUT_DIR / "panel_bangladesh_rice_full.parquet"
        full_panel.to_parquet(full_output)
        
        print(f"\n✓ Full panel shape: {full_panel.shape}")
        print(f"✓ Saved: {full_output}")
        print(f"\nBreakdown:")
        print(f"  Districts: {full_panel['adm_code'].nunique()}")
        print(f"  Years: {sorted(full_panel['year'].unique())}")
        print(f"  Seasons: {sorted(full_panel['season'].unique())}")
        print(f"  Periods: {sorted(full_panel['period'].unique())}")
    
    total_time = time.time() - overall_start
    print(f"\nTotal runtime: {total_time/60:.1f} minutes")
    print()


if __name__ == "__main__":
    main()
