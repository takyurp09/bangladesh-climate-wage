#!/usr/bin/env python3
"""
================================================================================
SCRIPT 2: Compute Climate Indices (Bangladesh)
================================================================================

Purpose: Compute EDD, GDD, HDD, pr1, pr2 from daily ERA5 data

Methods:
  - EDD: Sinusoidal exceedance degree-days (Schlenker & Roberts 2009)
  - GDD: Sinusoidal growing degree-days (Schlenker & Roberts 2009, Snyder 1985)
  - HDD: Sinusoidal harmful degree-days (above upper threshold)
  - pr1: Daily precipitation (mm/day)
  - pr2: pr1²

Key References:
  - Schlenker & Roberts (2009). PNAS 106(37), 15594-15598.
  - Snyder (1985). Agricultural and Forest Meteorology, 35, 353-358.

Author: Taky Tahmid
Date: 2026-04-05
================================================================================
"""

from pathlib import Path
import xarray as xr
import numpy as np
import time

# =============================================================================
# CONFIGURATION
# =============================================================================

COUNTRY_NAME = "bangladesh"
YEAR_START = 2013
YEAR_END = 2023

# Thresholds
EDD_THRESHOLDS = [0, 4, 8, 12, 28, 30, 32, 35]

# GDD specifications: (base, cap)
# Following Schlenker & Roberts framework for different crops
GDD_SPECS = [
    (8, 30), (8, 32), (8, 35),
    (10, 30), (10, 32), (10, 35),
    (12, 30), (12, 32), (12, 35),
    (15, 30), (15, 32), (15, 35),
    (20, 30), (20, 32), (20, 35),
]

HDD_THRESHOLDS = [10, 15]

# Paths
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent  # code/pipeline/ -> code/ -> project root
DAILY_DIR = PROJECT_ROOT / "data/climate_data/era5_daily/netcdf"
INDICES_DIR = PROJECT_ROOT / "data/climate_data/era5_indices/netcdf"

INDICES_DIR.mkdir(parents=True, exist_ok=True)

# =============================================================================
# SINUSOIDAL CLIMATE INDICES (Schlenker & Roberts 2009, Snyder 1985)
# =============================================================================

def calculate_sinusoidal_gdd(tasmax: xr.DataArray, tasmin: xr.DataArray, 
                             base: float, cap: float) -> xr.DataArray:
    """
    Sinusoidal Growing Degree-Days (Schlenker & Roberts 2009, Snyder 1985)
    
    Calculates time-weighted integral of temperature above base, capped at upper threshold.
    Uses sinusoidal approximation of within-day temperature distribution.
    
    Args:
        tasmax: Daily maximum temperature (°C)
        tasmin: Daily minimum temperature (°C)
        base: Base temperature (°C)
        cap: Upper threshold/cap (°C)
    
    Returns:
        Daily GDD (degree-days)
    """
    
    M = (tasmax + tasmin) / 2.0  # Daily mean
    W = (tasmax - tasmin) / 2.0  # Daily half-range
    
    eps = 1e-6  # Guard against flat days
    
    # Case 1: Entire day below base
    gdd = xr.where(tasmax <= base, 0.0, np.nan)
    
    # Case 2: Entire day above cap
    gdd = xr.where((tasmin >= cap) & (tasmax > base), cap - base, gdd)
    
    # Case 3: Day spans base but not cap (base < tmin < tmax < cap)
    mask3 = (tasmin > base) & (tasmax < cap)
    gdd = xr.where(mask3, M - base, gdd)
    
    # Case 4: Day crosses base from below (tmin < base < tmax < cap)
    mask4 = (tasmin < base) & (tasmax > base) & (tasmax < cap)
    Wsafe = xr.where(np.abs(W) < eps, np.nan, W)
    theta_base = np.arccos((base - M) / Wsafe)
    gdd4 = (W * np.sin(theta_base) + (M - base) * theta_base) / np.pi
    gdd = xr.where(mask4, gdd4, gdd)
    
    # Case 5: Day crosses cap from below (base < tmin < cap < tmax)
    mask5 = (tasmin > base) & (tasmin < cap) & (tasmax > cap)
    theta_cap = np.arccos((cap - M) / Wsafe)
    gdd5 = (W * np.sin(theta_cap) + (M - base) * theta_cap) / np.pi + (cap - base) * (np.pi - theta_cap) / np.pi
    gdd = xr.where(mask5, gdd5, gdd)
    
    # Case 6: Day spans both base and cap (tmin < base < cap < tmax)
    mask6 = (tasmin < base) & (tasmax > cap)
    theta_base6 = np.arccos((base - M) / Wsafe)
    theta_cap6 = np.arccos((cap - M) / Wsafe)
    gdd6 = ((W * np.sin(theta_base6) + (M - base) * theta_base6) / np.pi + 
            (cap - base) * (theta_base6 - theta_cap6) / np.pi)
    gdd = xr.where(mask6, gdd6, gdd)
    
    # Handle flat days (W ≈ 0)
    gdd = xr.where(
        np.abs(W) < eps,
        xr.where((M > base) & (M < cap), M - base,
                xr.where(M >= cap, cap - base, 0.0)),
        gdd
    )
    
    # Numerical safety
    gdd = gdd.clip(min=0.0, max=cap - base)
    
    return gdd.astype("float32")


def calculate_sinusoidal_edd(tasmax: xr.DataArray, tasmin: xr.DataArray, 
                             threshold: float) -> xr.DataArray:
    """
    Sinusoidal Exceedance Degree-Days (Schlenker & Roberts 2009)
    
    Calculates time-weighted integral of temperature ABOVE threshold.
    This is the harmful temperature exposure.
    
    Args:
        tasmax: Daily maximum temperature (°C)
        tasmin: Daily minimum temperature (°C)
        threshold: Threshold temperature (°C)
    
    Returns:
        Daily EDD (degree-days)
    """
    
    M = (tasmax + tasmin) / 2.0
    W = (tasmax - tasmin) / 2.0
    
    eps = 1e-6
    
    # Safe arccos path
    Wsafe = xr.where(np.abs(W) < eps, np.nan, W)
    R = (threshold - M) / Wsafe
    R = R.clip(-1.0, 1.0)
    
    theta = np.arccos(R)
    
    # Sinusoidal exceedance integral
    edd_mid = (W * np.sin(theta) + (M - threshold) * theta) / np.pi
    
    edd = xr.where(
        threshold >= tasmax, 0.0,
        xr.where(
            threshold <= tasmin, M - threshold,
            edd_mid
        )
    )
    
    # Handle flat days safely
    edd = xr.where(np.abs(W) < eps, xr.where(M > threshold, M - threshold, 0.0), edd)
    
    # Numerical safety
    edd = edd.clip(min=0.0)
    
    return edd.astype("float32")


def calculate_sinusoidal_hdd(tasmax: xr.DataArray, tasmin: xr.DataArray, 
                             threshold: float) -> xr.DataArray:
    """
    Sinusoidal Harmful Degree-Days (temperature BELOW threshold)
    
    Calculates time-weighted integral of temperature below threshold.
    This is for cold damage.
    
    Args:
        tasmax: Daily maximum temperature (°C)
        tasmin: Daily minimum temperature (°C)
        threshold: Threshold temperature (°C)
    
    Returns:
        Daily HDD (degree-days)
    """
    
    M = (tasmax + tasmin) / 2.0
    W = (tasmax - tasmin) / 2.0
    
    eps = 1e-6
    
    # Safe arccos path (inverted for cold)
    Wsafe = xr.where(np.abs(W) < eps, np.nan, W)
    R = (M - threshold) / Wsafe
    R = R.clip(-1.0, 1.0)
    
    theta = np.arccos(R)
    
    # Sinusoidal integral below threshold
    hdd_mid = (W * np.sin(theta) + (threshold - M) * theta) / np.pi
    
    hdd = xr.where(
        threshold <= tasmin, 0.0,
        xr.where(
            threshold >= tasmax, threshold - M,
            hdd_mid
        )
    )
    
    # Handle flat days safely
    hdd = xr.where(np.abs(W) < eps, xr.where(M < threshold, threshold - M, 0.0), hdd)
    
    # Numerical safety
    hdd = hdd.clip(min=0.0)
    
    return hdd.astype("float32")


# =============================================================================
# PROCESS ONE YEAR
# =============================================================================

def process_year(year: int) -> None:
    """
    Process one year: load daily data, compute all indices, save
    
    Args:
        year: Year to process (e.g., 2013)
    """
    
    # Output path
    out_file = INDICES_DIR / f"indices_{COUNTRY_NAME}_ERA5_{year}.nc"
    
    if out_file.exists():
        print(f"  ✓ {year}: already exists, skipping")
        return
    
    print(f"\n[{year}]")
    print(f"  Loading daily data...")
    
    # Input files
    tas_file = DAILY_DIR / f"tas_{COUNTRY_NAME}_ERA5_{year}.nc"
    tasmin_file = DAILY_DIR / f"tasmin_{COUNTRY_NAME}_ERA5_{year}.nc"
    tasmax_file = DAILY_DIR / f"tasmax_{COUNTRY_NAME}_ERA5_{year}.nc"
    pr_file = DAILY_DIR / f"pr_{COUNTRY_NAME}_ERA5_{year}.nc"
    
    # Check all files exist
    for f in [tas_file, tasmin_file, tasmax_file, pr_file]:
        if not f.exists():
            print(f"    ✗ Missing input file: {f.name}")
            return
    
    # Load data
    start = time.time()
    tas = xr.open_dataset(tas_file)["tas"]
    tasmin = xr.open_dataset(tasmin_file)["tasmin"]
    tasmax = xr.open_dataset(tasmax_file)["tasmax"]
    pr = xr.open_dataset(pr_file)["pr"]
    print(f"    → Loaded {tas.sizes['time']} days ({time.time()-start:.1f}s)")
    
    # Initialize output dataset
    indices = {}
    
    # --- EDD (Sinusoidal) ---
    print(f"  Computing EDD ({len(EDD_THRESHOLDS)} thresholds)...")
    start = time.time()
    for t in EDD_THRESHOLDS:
        var_name = f"edd_{t}"
        indices[var_name] = calculate_sinusoidal_edd(tasmax, tasmin, t)
        indices[var_name].attrs["long_name"] = f"Exceedance degree-days (threshold {t}°C, sinusoidal)"
        indices[var_name].attrs["units"] = "degree_days"
        indices[var_name].attrs["method"] = "Schlenker & Roberts (2009)"
    print(f"    → Done ({time.time()-start:.1f}s)")
    
    # --- GDD (Sinusoidal) ---
    print(f"  Computing GDD ({len(GDD_SPECS)} variants)...")
    start = time.time()
    for base, cap in GDD_SPECS:
        var_name = f"gdd_{base}_{cap}"
        indices[var_name] = calculate_sinusoidal_gdd(tasmax, tasmin, base, cap)
        indices[var_name].attrs["long_name"] = f"Growing degree-days (base {base}°C, cap {cap}°C, sinusoidal)"
        indices[var_name].attrs["units"] = "degree_days"
        indices[var_name].attrs["method"] = "Schlenker & Roberts (2009), Snyder (1985)"
    print(f"    → Done ({time.time()-start:.1f}s)")
    
    # --- HDD (Sinusoidal) ---
    print(f"  Computing HDD ({len(HDD_THRESHOLDS)} thresholds)...")
    start = time.time()
    for t in HDD_THRESHOLDS:
        var_name = f"hdd_{t}"
        indices[var_name] = calculate_sinusoidal_hdd(tasmax, tasmin, t)
        indices[var_name].attrs["long_name"] = f"Harmful degree-days (threshold {t}°C, sinusoidal)"
        indices[var_name].attrs["units"] = "degree_days"
        indices[var_name].attrs["method"] = "Sinusoidal (cold damage)"
    print(f"    → Done ({time.time()-start:.1f}s)")
    
    # --- PR ---
    print(f"  Computing precipitation indices...")
    start = time.time()
    indices["pr1"] = pr.astype("float32")
    indices["pr1"].attrs["long_name"] = "Daily precipitation"
    indices["pr1"].attrs["units"] = "mm/day"
    
    indices["pr2"] = (pr ** 2).astype("float32")
    indices["pr2"].attrs["long_name"] = "Daily precipitation squared"
    indices["pr2"].attrs["units"] = "(mm/day)^2"
    print(f"    → Done ({time.time()-start:.1f}s)")
    
    # Create dataset
    print(f"  Saving indices file...")
    start = time.time()
    ds_out = xr.Dataset(indices)
    ds_out.attrs["title"] = f"Climate Indices for Bangladesh ({year})"
    ds_out.attrs["source"] = "ERA5 daily data"
    ds_out.attrs["institution"] = "University of Delaware"
    ds_out.attrs["method"] = "Sinusoidal degree-days (Schlenker & Roberts 2009, Snyder 1985)"
    ds_out.attrs["created"] = str(np.datetime64('today'))
    
    # Save with compression
    encoding = {var: {"zlib": True, "complevel": 1} for var in ds_out.data_vars}
    ds_out.to_netcdf(out_file, encoding=encoding)
    
    n_vars = len(ds_out.data_vars)
    n_days = ds_out.sizes["time"]
    file_size_mb = out_file.stat().st_size / 1024 / 1024
    print(f"    → Saved {n_vars} variables, {n_days} days ({file_size_mb:.1f} MB, {time.time()-start:.1f}s)")


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("\n" + "="*80)
    print(f" Compute Climate Indices ({COUNTRY_NAME.upper()})")
    print(" Method: Sinusoidal Degree-Days (Schlenker & Roberts 2009)")
    print("="*80)
    print(f"Years: {YEAR_START}-{YEAR_END}")
    print(f"Input:  {DAILY_DIR}")
    print(f"Output: {INDICES_DIR}")
    print()
    print(f"Variables to compute:")
    print(f"  EDD (sinusoidal): {len(EDD_THRESHOLDS)} thresholds")
    print(f"  GDD (sinusoidal): {len(GDD_SPECS)} variants (base/cap pairs)")
    print(f"  HDD (sinusoidal): {len(HDD_THRESHOLDS)} thresholds")
    print(f"  PR:  2 variables (pr1, pr2)")
    print(f"  Total: {len(EDD_THRESHOLDS) + len(GDD_SPECS) + len(HDD_THRESHOLDS) + 2} variables per day")
    
    overall_start = time.time()
    
    # Process all years
    for year in range(YEAR_START, YEAR_END + 1):
        process_year(year)
    
    print("\n" + "="*80)
    print(" ✓ INDICES COMPUTATION COMPLETE")
    print("="*80)
    
    # Summary
    indices_files = sorted(INDICES_DIR.glob("indices_*.nc"))
    print(f"\nTotal indices files created: {len(indices_files)}")
    
    total_time = time.time() - overall_start
    print(f"Total runtime: {total_time/60:.1f} minutes")
    print()


if __name__ == "__main__":
    main()
