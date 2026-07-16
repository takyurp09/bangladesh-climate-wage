#!/usr/bin/env python3
"""
================================================================================
SCRIPT 3: Bangladesh Rice Calendar Windows (Boro/Aus/Aman)
================================================================================

Purpose: Convert Bangladesh rice crop calendars to calendar windows (.pkl)

Input:
  - bangladesh_rice_calendar/data/final/Rice.Boro.crop.calendar.bgd.nc
  - bangladesh_rice_calendar/data/final/Rice.Aus.crop.calendar.bgd.nc
  - bangladesh_rice_calendar/data/final/Rice.Aman.crop.calendar.bgd.nc
  - data/climate_data/gadm41_BGD_shp/gadm41_BGD_2.shp

Output:
  - data/calendar_windows/calendar_windows_BGD_Rice.Boro.pkl
  - data/calendar_windows/calendar_windows_BGD_Rice.Aus.pkl
  - data/calendar_windows/calendar_windows_BGD_Rice.Aman.pkl

Each .pkl contains windows dict with keys: 'plant', 'between', 'harvest'
Each window has: 'start' (0-based index), 'end' (0-based index), 'days' (count)

Method (from CI2026):
  1. Load ADM2 shapefile (64 districts)
  2. For each district, extract calendar DOYs from nearest pixel
  3. Convert DOYs to 0-based indices for 2-year concatenated time axis
  4. Save as pickle for fast loading in Script 4

Author: Taky Tahmid
Date: 2026-04-05
================================================================================
"""

from pathlib import Path
import pickle
import geopandas as gpd
import numpy as np
import xarray as xr

# =============================================================================
# CONFIGURATION
# =============================================================================

COUNTRY_NAME = "bangladesh"
COUNTRY_TAG = "BGD"
ADM_LEVEL = 2

# Paths
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent  # code/pipeline/ -> code/ -> project root

SHAPEFILE = PROJECT_ROOT / "data/climate_data/gadm41_BGD_shp/gadm41_BGD_2.shp"
CALENDAR_DIR = PROJECT_ROOT / "bangladesh_rice_calendar/data/final"
WINDOWS_DIR = PROJECT_ROOT / "data/calendar_windows"

WINDOWS_DIR.mkdir(parents=True, exist_ok=True)

# Rice calendars (3 seasons)
RICE_CALENDARS = ["Rice.Boro", "Rice.Aus", "Rice.Aman"]

# =============================================================================
# CALENDAR WINDOW COMPUTATION (from CI2026)
# =============================================================================

def compute_windows_for_calendar(gdf: gpd.GeoDataFrame, cal_name: str) -> dict:
    """
    Compute plant/between/harvest windows for each polygon.
    
    Logic from CI2026 ALTERNATIVE_edds_calendar.py:
    - Uses 2-year concatenated time axis (0-729 indices)
    - Handles year-wrapping seasons (e.g., Boro: Dec-May)
    - Returns 0-based start/end indices + day counts
    
    Args:
        gdf: GeoDataFrame with district polygons
        cal_name: Calendar name (e.g., 'Rice.Boro')
    
    Returns:
        dict with 'plant', 'between', 'harvest' windows
    """
    
    # Load calendar NetCDF
    cal_path = CALENDAR_DIR / f"{cal_name}.crop.calendar.bgd.nc"
    
    if not cal_path.exists():
        raise FileNotFoundError(f"Calendar not found: {cal_path}")
    
    print(f"  Loading {cal_name}...")
    dscal = xr.open_dataset(cal_path)
    
    # Get centroids for nearest-pixel extraction
    centroids = gdf.geometry.centroid
    
    n_poly = len(gdf)
    
    # Initialize arrays
    plant_start = np.full(n_poly, np.nan, dtype=float)
    plant_end = np.full(n_poly, np.nan, dtype=float)
    between_start = np.full(n_poly, np.nan, dtype=float)
    between_end = np.full(n_poly, np.nan, dtype=float)
    harvest_start = np.full(n_poly, np.nan, dtype=float)
    harvest_end = np.full(n_poly, np.nan, dtype=float)
    days_plant = np.full(n_poly, np.nan, dtype=float)
    days_between = np.full(n_poly, np.nan, dtype=float)
    days_harvest = np.full(n_poly, np.nan, dtype=float)
    
    T = 365
    offset = T  # Start at day 365 to handle year-wrapping
    
    # Extract calendar values for each polygon
    for ii, row in gdf.iterrows():
        try:
            # Get nearest pixel calendar values
            pix = dscal.sel(
                longitude=centroids.iloc[ii].x,
                latitude=centroids.iloc[ii].y,
                method="nearest",
            )
            
            # Extract DOYs (1-365)
            hs_doy = int(pix["harvest.start"].values)
            he_doy = int(pix["harvest.end"].values)
            ps_doy = int(pix["plant.start"].values)
            pe_doy = int(pix["plant.end"].values)
            
            # Skip if NaN (ocean/outside growing area)
            if np.isnan([hs_doy, he_doy, ps_doy, pe_doy]).any():
                continue
            
            # Convert DOYs to 0-based indices on 2-year axis
            # Following CI2026 logic exactly
            
            # Harvest end (reference point)
            harvest_end_day = he_doy + offset
            
            # Harvest start (must be <= harvest end)
            harvest_start_day = hs_doy + offset
            if harvest_start_day > harvest_end_day:
                harvest_start_day -= T
            
            # Plant end (must be <= harvest end)
            plant_end_day = pe_doy + offset
            if plant_end_day > harvest_end_day:
                plant_end_day -= T
            
            # Plant start (must be <= plant end)
            plant_start_day = ps_doy + offset
            if plant_start_day > plant_end_day:
                plant_start_day -= T
            
            # Between period (plant end + 1 to harvest start - 1)
            between_start_day = plant_end_day + 1
            between_end_day = harvest_start_day - 1
            if between_end_day < between_start_day:
                between_end_day = between_start_day - 1
            
            # Convert to 0-based indices (subtract 1 from start, keep end as-is)
            plant_start[ii] = plant_start_day - 1
            plant_end[ii] = plant_end_day
            between_start[ii] = between_start_day - 1
            between_end[ii] = between_end_day
            harvest_start[ii] = harvest_start_day - 1
            harvest_end[ii] = harvest_end_day
            
            # Day counts
            days_plant[ii] = plant_end_day - plant_start_day
            days_between[ii] = between_end_day - between_start_day
            days_harvest[ii] = harvest_end_day - harvest_start_day
            
        except Exception as e:
            # Skip polygons with errors (ocean, missing data)
            continue
    
    # Create windows dict (CI2026 format)
    windows = {
        "plant": {
            "start": plant_start,
            "end": plant_end,
            "days": days_plant,
        },
        "between": {
            "start": between_start,
            "end": between_end,
            "days": days_between,
        },
        "harvest": {
            "start": harvest_start,
            "end": harvest_end,
            "days": days_harvest,
        },
    }
    
    # Count valid polygons
    valid_count = np.sum(~np.isnan(plant_start))
    print(f"    → Valid districts: {valid_count}/{n_poly}")
    
    return windows


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("\n" + "="*80)
    print(f" Precompute Calendar Windows for {COUNTRY_TAG} Rice")
    print("="*80)
    print(f"Shapefile: {SHAPEFILE}")
    print(f"Calendars: {CALENDAR_DIR}")
    print(f"Output:    {WINDOWS_DIR}")
    print()
    
    # Load shapefile
    if not SHAPEFILE.exists():
        raise FileNotFoundError(f"Shapefile not found: {SHAPEFILE}")
    
    print(f"Loading ADM{ADM_LEVEL} shapefile...")
    gdf = gpd.read_file(SHAPEFILE)
    gdf = gdf.to_crs("EPSG:4326")
    gdf = gdf.reset_index(drop=True)
    
    print(f"  → {len(gdf)} districts (ADM{ADM_LEVEL})")
    print()
    
    # Process each rice calendar
    for cal in RICE_CALENDARS:
        out_path = WINDOWS_DIR / f"calendar_windows_{COUNTRY_TAG}_{cal}.pkl"
        
        if out_path.exists():
            print(f"[{cal}]")
            print(f"  ✓ Already exists, skipping")
            continue
        
        print(f"[{cal}]")
        windows = compute_windows_for_calendar(gdf, cal)
        
        # Save pickle
        with open(out_path, "wb") as f:
            pickle.dump(windows, f)
        
        print(f"  ✓ Saved: {out_path.name}")
        print()
    
    print("="*80)
    print(" ✓ CALENDAR WINDOWS COMPLETE")
    print("="*80)
    
    # Summary
    pkl_files = sorted(WINDOWS_DIR.glob("calendar_windows_*.pkl"))
    print(f"\nTotal .pkl files: {len(pkl_files)}")
    for f in pkl_files:
        print(f"  {f.name}")
    print()


if __name__ == "__main__":
    main()
