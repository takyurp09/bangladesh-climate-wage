#!/usr/bin/env python3
"""
================================================================================
SCRIPT 1: ERA5 Hourly → Daily Conversion (Bangladesh)
================================================================================

Purpose: Convert ERA5 hourly t2m and tp to daily tas/tasmin/tasmax/pr

Input:
  - data/climate_data/era5_hourly/netcdf/t2m_YYYY.nc  (hourly 2m temperature, Kelvin)
  - data/climate_data/era5_hourly/netcdf/tp_YYYY.nc   (hourly total precipitation, meters)

Output:
  - data/climate_data/era5_daily/netcdf/tas_bangladesh_ERA5_YYYY.nc    (daily mean temp, °C)
  - data/climate_data/era5_daily/netcdf/tasmin_bangladesh_ERA5_YYYY.nc (daily min temp, °C)
  - data/climate_data/era5_daily/netcdf/tasmax_bangladesh_ERA5_YYYY.nc (daily max temp, °C)
  - data/climate_data/era5_daily/netcdf/pr_bangladesh_ERA5_YYYY.nc     (daily precip, mm/day)

Years: 2013-2023 (default, configurable)

Method:
  1. Read hourly ERA5 NetCDF
  2. Subset to Bangladesh bbox (20.5-26.7°N, 88.0-92.7°E)
  3. Temperature: K → °C, hourly → daily (mean/min/max)
  4. Precipitation: m → mm, hourly → daily (sum)
  5. Save as daily NetCDF files

Runtime: ~2-3 minutes per year

Author: Taky Tahmid
Date: 2026-04-05
================================================================================
"""

from pathlib import Path
import xarray as xr
import numpy as np

# =============================================================================
# CONFIGURATION
# =============================================================================

COUNTRY_NAME = "bangladesh"
YEAR_START = 2013
YEAR_END = 2023

# Bangladesh bounding box
BBOX = {
    "lat": slice(26.7, 20.5),  # Note: slice is inclusive, reversed for descending coords
    "lon": slice(88.0, 92.7)
}

# Paths (go up one level from code/ to wage_model/)
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent  # code/pipeline/ -> code/ -> project root
HOURLY_DIR = PROJECT_ROOT / "data/climate_data/era5_hourly/netcdf"
DAILY_DIR = PROJECT_ROOT / "data/climate_data/era5_daily/netcdf"

# Create output directory
DAILY_DIR.mkdir(parents=True, exist_ok=True)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def process_temperature_year(year: int) -> None:
    """
    Process one year of temperature data: hourly t2m → daily tas/tasmin/tasmax
    
    Args:
        year: Year to process (e.g., 2013)
    """
    
    # Input/output paths
    in_file = HOURLY_DIR / f"t2m_{year}.nc"
    out_tas = DAILY_DIR / f"tas_{COUNTRY_NAME}_ERA5_{year}.nc"
    out_tasmin = DAILY_DIR / f"tasmin_{COUNTRY_NAME}_ERA5_{year}.nc"
    out_tasmax = DAILY_DIR / f"tasmax_{COUNTRY_NAME}_ERA5_{year}.nc"
    
    # Skip if all outputs exist
    if out_tas.exists() and out_tasmin.exists() and out_tasmax.exists():
        print(f"  ✓ {year} temperature: already exists, skipping")
        return
    
    if not in_file.exists():
        print(f"  ✗ {year} temperature: input file not found - {in_file}")
        return
    
    print(f"  Processing {year} temperature...")
    
    # Load hourly data
    ds = xr.open_dataset(in_file)
    
    # Normalize coordinate names
    rename = {}
    if "valid_time" in ds.coords:
        rename["valid_time"] = "time"
    if "latitude" in ds.coords:
        rename["latitude"] = "lat"
    if "longitude" in ds.coords:
        rename["longitude"] = "lon"
    if rename:
        ds = ds.rename(rename)
    
    # Subset to Bangladesh
    ds = ds.sel(BBOX)
    
    if ds.sizes.get("time", 0) == 0:
        print(f"  ✗ {year} temperature: no data in bbox")
        return
    
    # Get t2m variable
    if "t2m" not in ds.data_vars:
        print(f"  ✗ {year} temperature: 't2m' variable not found")
        return
    
    t2m = ds["t2m"]
    
    # Convert K → °C
    t2m_c = (t2m - 273.15).astype("float32")
    
    # Ensure we have complete days (24-hour chunks)
    n_time = t2m_c.sizes["time"]
    n_days = n_time // 24
    t2m_c = t2m_c.isel(time=slice(0, n_days * 24))
    
    # Reshape to (days, 24, lat, lon)
    daily_time = t2m_c["time"].isel(time=slice(0, n_days * 24, 24))
    arr = t2m_c.data.reshape((n_days, 24) + t2m_c.data.shape[1:])
    
    # Compute daily aggregates
    tas = xr.DataArray(
        arr.mean(axis=1),
        dims=("time", "lat", "lon"),
        coords={"time": daily_time, "lat": t2m_c["lat"], "lon": t2m_c["lon"]},
        name="tas",
    )
    tas.attrs["units"] = "degC"
    tas.attrs["long_name"] = "Daily mean 2m air temperature"
    
    tasmin = xr.DataArray(
        arr.min(axis=1),
        dims=("time", "lat", "lon"),
        coords={"time": daily_time, "lat": t2m_c["lat"], "lon": t2m_c["lon"]},
        name="tasmin",
    )
    tasmin.attrs["units"] = "degC"
    tasmin.attrs["long_name"] = "Daily minimum 2m air temperature"
    
    tasmax = xr.DataArray(
        arr.max(axis=1),
        dims=("time", "lat", "lon"),
        coords={"time": daily_time, "lat": t2m_c["lat"], "lon": t2m_c["lon"]},
        name="tasmax",
    )
    tasmax.attrs["units"] = "degC"
    tasmax.attrs["long_name"] = "Daily maximum 2m air temperature"
    
    # Save individual files
    enc = {"zlib": True, "complevel": 1}
    
    tas.to_dataset().to_netcdf(out_tas, encoding={"tas": enc})
    tasmin.to_dataset().to_netcdf(out_tasmin, encoding={"tasmin": enc})
    tasmax.to_dataset().to_netcdf(out_tasmax, encoding={"tasmax": enc})
    
    print(f"    → Saved {n_days} days: tas, tasmin, tasmax")


def process_precipitation_year(year: int) -> None:
    """
    Process one year of precipitation data: hourly tp → daily pr
    
    Args:
        year: Year to process (e.g., 2013)
    """
    
    # Input/output paths
    in_file = HOURLY_DIR / f"tp_{year}.nc"
    out_pr = DAILY_DIR / f"pr_{COUNTRY_NAME}_ERA5_{year}.nc"
    
    # Skip if output exists
    if out_pr.exists():
        print(f"  ✓ {year} precipitation: already exists, skipping")
        return
    
    if not in_file.exists():
        print(f"  ✗ {year} precipitation: input file not found - {in_file}")
        return
    
    print(f"  Processing {year} precipitation...")
    
    # Load hourly data
    ds = xr.open_dataset(in_file)
    
    # Normalize coordinate names
    rename = {}
    if "valid_time" in ds.coords:
        rename["valid_time"] = "time"
    if "latitude" in ds.coords:
        rename["latitude"] = "lat"
    if "longitude" in ds.coords:
        rename["longitude"] = "lon"
    if rename:
        ds = ds.rename(rename)
    
    # Subset to Bangladesh
    ds = ds.sel(BBOX)
    
    if ds.sizes.get("time", 0) == 0:
        print(f"  ✗ {year} precipitation: no data in bbox")
        return
    
    # Get tp variable
    if "tp" not in ds.data_vars:
        print(f"  ✗ {year} precipitation: 'tp' variable not found")
        return
    
    tp = ds["tp"]
    
    # Convert m → mm (hourly accumulation)
    tp_mm_hourly = (tp * 1000.0).astype("float32")
    
    # Resample to daily (sum of 24 hourly values)
    pr_daily = tp_mm_hourly.resample(time="1D").sum()
    
    pr_daily.attrs["units"] = "mm/day"
    pr_daily.attrs["long_name"] = "Daily total precipitation"
    pr_daily = pr_daily.rename("pr")
    
    # Save
    enc = {"pr": {"zlib": True, "complevel": 1}}
    pr_daily.to_dataset().to_netcdf(out_pr, encoding=enc)
    
    print(f"    → Saved {pr_daily.sizes['time']} days: pr")


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("\n" + "="*80)
    print(f" ERA5 Hourly → Daily Conversion ({COUNTRY_NAME.upper()})")
    print("="*80)
    print(f"Years: {YEAR_START}-{YEAR_END}")
    print(f"Bbox: {BBOX['lat'].start}°N-{BBOX['lat'].stop}°N, {BBOX['lon'].start}°E-{BBOX['lon'].stop}°E")
    print(f"Input:  {HOURLY_DIR}")
    print(f"Output: {DAILY_DIR}")
    print()
    
    # Process all years
    for year in range(YEAR_START, YEAR_END + 1):
        print(f"\n[{year}]")
        
        # Temperature (tas, tasmin, tasmax)
        process_temperature_year(year)
        
        # Precipitation (pr)
        process_precipitation_year(year)
    
    print("\n" + "="*80)
    print(" ✓ CONVERSION COMPLETE")
    print("="*80)
    
    # Summary
    daily_files = sorted(DAILY_DIR.glob("*.nc"))
    print(f"\nTotal daily files created: {len(daily_files)}")
    
    tas_count = len(list(DAILY_DIR.glob("tas_*.nc")))
    tasmin_count = len(list(DAILY_DIR.glob("tasmin_*.nc")))
    tasmax_count = len(list(DAILY_DIR.glob("tasmax_*.nc")))
    pr_count = len(list(DAILY_DIR.glob("pr_*.nc")))
    
    print(f"  tas:    {tas_count} files")
    print(f"  tasmin: {tasmin_count} files")
    print(f"  tasmax: {tasmax_count} files")
    print(f"  pr:     {pr_count} files")
    print()


if __name__ == "__main__":
    main()
