#!/usr/bin/env python3
"""
================================================================================
SCRIPT 7: Temperature Bin Aggregation & Schlenker-Roberts Figure
================================================================================

Purpose: Create temperature bin exposure variables and replicate Schlenker & 
         Roberts (2009) temperature-yield relationship figure

Method:
  1. For each district-season-year observation:
     - Load daily tasmin/tasmax for growing season (using calendar windows)
     - Calculate time spent in each 1°C temperature bin (0-1°C, 1-2°C, ..., 39-40°C)
     - Use sinusoidal approximation for within-day temperature distribution
  2. Create dataset with 40 temperature bin variables
  3. Merge with yield data
  4. Run regression: log_yield ~ temp_bins + district_FE + year_FE
  5. Plot coefficients (Schlenker-Roberts style figure)

Author: Taky Tahmid
Date: 2026-04-06
================================================================================
"""

from pathlib import Path
import pickle
import sys
import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent / "utils"))
from district_names import harmonize_district  # noqa: E402
import xagg as xa
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent  # code/pipeline/ -> code/ -> project root

DAILY_DIR = PROJECT_ROOT / "data/climate_data/era5_daily/netcdf"
WINDOWS_DIR = PROJECT_ROOT / "data/calendar_windows"
SHAPEFILE = PROJECT_ROOT / "data/climate_data/gadm41_BGD_shp/gadm41_BGD_2.shp"
YIELD_DATA = PROJECT_ROOT / "data/agriculturedata_bangladesh/crop/df_rice.csv"
OUTPUT_DIR = PROJECT_ROOT / "data/Regression_data"
FIGURE_DIR = PROJECT_ROOT / "figures"

OUTPUT_DIR.mkdir(exist_ok=True)
FIGURE_DIR.mkdir(exist_ok=True)

# Temperature bins (1°C intervals, 0-40°C)
TEMP_BINS = np.arange(0, 41, 1)
N_BINS = len(TEMP_BINS) - 1

YEARS = range(2013, 2024)
SEASONS = {
    'Rice.Boro': 'Boro',
    'Rice.Aus': 'Aus',
    'Rice.Aman': 'Aman'
}

# District harmonization imported from code/utils/district_names.py

# =============================================================================
# TEMPERATURE BIN CALCULATION (Sinusoidal)
# =============================================================================

def calculate_temp_bin_exposure(tasmin, tasmax, temp_bins):
    """
    Calculate fraction of day spent in each temperature bin using sinusoidal approximation.
    
    Args:
        tasmin: Daily minimum temperature (°C)
        tasmax: Daily maximum temperature (°C)
        temp_bins: Temperature bin edges (e.g., [0, 1, 2, ..., 40])
    
    Returns:
        Array of fractions (length = len(temp_bins) - 1)
    """
    
    n_bins = len(temp_bins) - 1
    exposure = np.zeros(n_bins)
    
    if np.isnan(tasmin) or np.isnan(tasmax) or tasmin >= tasmax:
        return exposure
    
    # Sinusoidal parameters
    M = (tasmax + tasmin) / 2  # Mean
    W = (tasmax - tasmin) / 2  # Half-range
    
    if W == 0:
        # Constant temperature all day
        bin_idx = np.searchsorted(temp_bins, M) - 1
        if 0 <= bin_idx < n_bins:
            exposure[bin_idx] = 1.0
        return exposure
    
    # For each bin, calculate fraction of day with temperature in that bin
    for i in range(n_bins):
        T_low = temp_bins[i]
        T_high = temp_bins[i + 1]
        
        # Temperature function: T(θ) = M + W·sin(θ), θ ∈ [-π/2, π/2]
        # (shifted to have min at θ=-π/2, max at θ=π/2)
        
        # Find θ values where T(θ) = T_low and T(θ) = T_high
        # sin(θ) = (T - M) / W
        
        # Case 1: Entire day below bin
        if tasmax <= T_low:
            continue
        
        # Case 2: Entire day above bin
        if tasmin >= T_high:
            continue
        
        # Case 3: Entire day within bin
        if tasmin >= T_low and tasmax <= T_high:
            exposure[i] = 1.0
            continue
        
        # Case 4: Bin partially overlaps day
        # Calculate θ boundaries
        sin_low = np.clip((T_low - M) / W, -1, 1)
        sin_high = np.clip((T_high - M) / W, -1, 1)
        
        theta_low = np.arcsin(sin_low)
        theta_high = np.arcsin(sin_high)
        
        # Handle multiple crossings (temperature crosses bin boundaries)
        if tasmin < T_low <= tasmax:
            # Crosses lower boundary
            if tasmax <= T_high:
                # Stays below upper boundary
                fraction = (np.pi/2 - theta_low) / np.pi
            else:
                # Crosses both boundaries
                fraction = (theta_high - theta_low) / np.pi
        elif tasmin < T_high <= tasmax:
            # Crosses upper boundary only
            fraction = (theta_high + np.pi/2) / np.pi
        else:
            # Within bin for entire temperature range
            fraction = (theta_high - theta_low) / np.pi
        
        exposure[i] = np.clip(fraction, 0, 1)
    
    # Normalize to sum to 1 (entire day)
    total = exposure.sum()
    if total > 0:
        exposure = exposure / total
    
    return exposure


# =============================================================================
# AGGREGATE TEMPERATURE BINS FOR SEASON
# =============================================================================

def aggregate_temp_bins_for_season(year, season, district_idx, windows_dict, shapefile):
    """
    Aggregate temperature bin exposure for one district-season-year.
    
    Returns:
        Dictionary with temp_bin_0_1, temp_bin_1_2, ..., temp_bin_39_40 (total days in each bin)
    """
    
    # Load daily temperature for this year
    tasmin_file = DAILY_DIR / f"tasmin_bangladesh_ERA5_{year}.nc"
    tasmax_file = DAILY_DIR / f"tasmax_bangladesh_ERA5_{year}.nc"
    
    if not tasmin_file.exists() or not tasmax_file.exists():
        return None
    
    tasmin_ds = xr.open_dataset(tasmin_file)
    tasmax_ds = xr.open_dataset(tasmax_file)
    
    # Get district polygon centroid (simple spatial aggregation)
    district_geom = shapefile.iloc[district_idx].geometry
    centroid = district_geom.centroid
    
    # Find nearest grid cell
    lat_idx = np.abs(tasmin_ds.lat.values - centroid.y).argmin()
    lon_idx = np.abs(tasmin_ds.lon.values - centroid.x).argmin()
    
    # Extract time series for this grid cell
    tasmin = tasmin_ds.tasmin.isel(lat=lat_idx, lon=lon_idx).values
    tasmax = tasmax_ds.tasmax.isel(lat=lat_idx, lon=lon_idx).values
    
    # Get calendar window indices

    # Concatenate to handle year-wrapping (like Script 4)
    tasmin_full = np.concatenate([tasmin, tasmin])
    tasmax_full = np.concatenate([tasmax, tasmax])

    # Initialize total exposure across all periods
    total_exposure = np.zeros(N_BINS)
                # Get window dict for this district
    window = {
        "plant": {"start": windows_dict["plant"]["start"][district_idx], "end": windows_dict["plant"]["end"][district_idx]},
        "between": {"start": windows_dict["between"]["start"][district_idx], "end": windows_dict["between"]["end"][district_idx]},
        "harvest": {"start": windows_dict["harvest"]["start"][district_idx], "end": windows_dict["harvest"]["end"][district_idx]}
    }
    
    for period in ['plant', 'between', 'harvest']:
        start_idx = window[period]['start']
        end_idx = window[period]['end']
        
        if np.isnan(start_idx) or np.isnan(end_idx):
            continue
        
        start_idx = int(start_idx)
        end_idx = int(end_idx)
        
        # Loop through days in this period
        for day_idx in range(start_idx, end_idx):
            tmin = tasmin_full[day_idx]
            tmax = tasmax_full[day_idx]
            
            # Calculate exposure for this day
            day_exposure = calculate_temp_bin_exposure(tmin, tmax, TEMP_BINS)
            total_exposure += day_exposure
    
    # Create result dictionary
    result = {f'temp_bin_{i}_{i+1}': total_exposure[i] for i in range(N_BINS)}
    
    tasmin_ds.close()
    tasmax_ds.close()
    
    return result


# =============================================================================
# MAIN AGGREGATION
# =============================================================================

def main():
    print("\n" + "="*80)
    print(" TEMPERATURE BIN AGGREGATION (Schlenker-Roberts Method)")
    print("="*80)
    
    # Load shapefile
    print("\nLoading shapefile...")
    gdf = gpd.read_file(SHAPEFILE)
    gdf = gdf.to_crs("EPSG:4326").reset_index(drop=True)
    n_districts = len(gdf)
    print(f"  → {n_districts} districts")
    
    # Load yield data
    print("\nLoading yield data...")
    yield_df = pd.read_csv(YIELD_DATA)
    yield_df = yield_df[(yield_df['Year'] >= 2013) & (yield_df['Year'] <= 2023)]
    yield_df = yield_df.rename(columns={
        'District': 'district',
        'Year': 'year',
        'Crop_type': 'season',
        'Area': 'area',
        'Production': 'production'
    })
    yield_df['yield_per_ha'] = yield_df['production'] / yield_df['area']
    yield_df['log_yield'] = np.log(yield_df['yield_per_ha'])
    yield_df = yield_df.replace([np.inf, -np.inf], np.nan)
    print(f"  → {len(yield_df)} yield observations")
    
    # Aggregate temperature bins for each observation
    print("\nAggregating temperature bins...")
    print("(This will take ~10-15 minutes)")
    
    all_rows = []
    
    for year in YEARS:
        print(f"\n[{year}]")
        
        for season_cal, season_yield in SEASONS.items():
            
            # Load calendar windows
            windows_file = WINDOWS_DIR / f"calendar_windows_BGD_{season_cal}.pkl"
            if not windows_file.exists():
                continue
            
            with open(windows_file, 'rb') as f:
                windows_dict = pickle.load(f)
            
            print(f"  {season_yield}...")
            
            for district_idx in range(n_districts):
                district_name = gdf.iloc[district_idx]['NAME_2']
                
                # Skip if no calendar
                if district_idx >= len(windows_dict['plant']['start']):
                    continue
                
                # Create window for this district
                window = {
                    'plant': {
                        'start': windows_dict['plant']['start'][district_idx],
                        'end': windows_dict['plant']['end'][district_idx]
                    },
                    'between': {
                        'start': windows_dict['between']['start'][district_idx],
                        'end': windows_dict['between']['end'][district_idx]
                    },
                    'harvest': {
                        'start': windows_dict['harvest']['start'][district_idx],
                        'end': windows_dict['harvest']['end'][district_idx]
                    }
                }
                
                # Aggregate temperature bins
                temp_bins = aggregate_temp_bins_for_season(
                    year, season_cal, district_idx, windows_dict, gdf
                )
                
                if temp_bins is None:
                    continue
                
                # Create row
                row = {
                    'district': district_name,
                    'year': year,
                    'season': season_yield,
                    **temp_bins
                }
                
                all_rows.append(row)
            
            print(f"    → {len([r for r in all_rows if r['year']==year and r['season']==season_yield])} obs")
    
    # Create DataFrame
    temp_bin_df = pd.DataFrame(all_rows)
    
    # Harmonize district names
    temp_bin_df['district'] = temp_bin_df['district'].map(harmonize_district)
    
    print(f"\n✓ Temperature bin dataset: {temp_bin_df.shape}")
    
    # Merge with yield
    print("\nMerging with yield data...")
    merged = pd.merge(temp_bin_df, yield_df[['district', 'year', 'season', 'area', 'production', 'yield_per_ha', 'log_yield']],
                      on=['district', 'year', 'season'], how='inner')
    
    print(f"  → {len(merged)} observations matched")
    
    # Save
    output_file = OUTPUT_DIR / "bangladesh_rice_temperature_bins.parquet"
    merged.to_parquet(output_file, index=False)
    print(f"\n✓ Saved: {output_file}")
    
    # Run regression and plot
    print("\n" + "="*80)
    print(" SCHLENKER-ROBERTS REGRESSION & FIGURE")
    print("="*80)
    
    create_schlenker_roberts_figure(merged)
    
    print("\n✓ COMPLETE!\n")


# =============================================================================
# SCHLENKER-ROBERTS FIGURE
# =============================================================================

def create_schlenker_roberts_figure(df):
    """
    Create Schlenker-Roberts style figure with temperature bins on x-axis.
    """
    
    # Remove missing yields
    df = df.dropna(subset=['log_yield'])
    
    # Demean log_yield by district-season (fixed effects)
    df['log_yield_demeaned'] = df.groupby(['district', 'season'])['log_yield'].transform(lambda x: x - x.mean())
    
    # Temperature bin columns
    temp_cols = [f'temp_bin_{i}_{i+1}' for i in range(N_BINS)]
    
    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    seasons = ['Boro', 'Aus', 'Aman']
    colors = {'Boro': 'blue', 'Aus': 'red', 'Aman': 'orange'}
    
    for idx, season in enumerate(seasons):
        ax = axes[idx]
        season_data = df[df['season'] == season].copy()
        
        if len(season_data) < 50:
            print(f"  ⚠️  {season}: Only {len(season_data)} observations, skipping")
            continue
        
        # Prepare regression
        X = season_data[temp_cols].values
        y = season_data['log_yield_demeaned'].values
        
        # Run regression (no intercept - coefficients are marginal effects)
        reg = LinearRegression(fit_intercept=False)
        reg.fit(X, y)
        
        coefs = reg.coef_
        
        # Normalize: Set exposure-weighted mean to zero (like Schlenker & Roberts)
        # Omit 29-30°C bin as reference (optimal for rice)
        ref_bin_idx = 29  # 29-30°C bin
        avg_exposure = season_data[temp_cols].mean().values
        # Normalize: center so exposure-weighted mean = 0 (Schlenker & Roberts method)
        weighted_mean = np.sum(coefs * avg_exposure) / np.sum(avg_exposure)
        coefs = coefs - weighted_mean
        coefs = coefs - weighted_mean
        temp_centers = TEMP_BINS[:-1] + 0.5
        
        # Plot step function
        ax.step(temp_centers, coefs, where='mid', color=colors[season], 
                linewidth=2, label='Step Function', alpha=0.7)
        
        # Plot polynomial smooth (8th order)
        valid_idx = ~np.isnan(coefs)
        if valid_idx.sum() >= 9:
            z = np.polyfit(temp_centers[valid_idx], coefs[valid_idx], 8)
            p = np.poly1d(z)
            temp_smooth = np.linspace(0, 40, 300)
            ax.plot(temp_smooth, p(temp_smooth), 'k-', linewidth=2.5, 
                    label='Polynomial (8th-order)')
            
            # Confidence band (simplified)
            residuals = y - reg.predict(X)
            se = np.std(residuals) / np.sqrt(len(y))
            ax.fill_between(temp_smooth, p(temp_smooth) - 1.96*se, 
                           p(temp_smooth) + 1.96*se, alpha=0.2, color='gray')
        
        # Histogram at bottom
        ax2 = ax.twinx()
        avg_exposure = season_data[temp_cols].mean().values
        ax2.bar(temp_centers, avg_exposure, width=0.8, alpha=0.3, 
                color='green', edgecolor='none')
        ax2.set_ylabel('Exposure (Days)', fontsize=9, color='green')
        ax2.tick_params(axis='y', labelcolor='green')
        ax2.set_ylim(0, avg_exposure.max() * 3)
        
        # Labels
        ax.set_xlabel('Temperature (Celsius)', fontsize=11, fontweight='bold')
        if idx == 0:
            ax.set_ylabel('Log Yield (MT/ha)', fontsize=11)
        ax.set_title(season, fontsize=13, fontweight='bold')
        ax.axhline(0, color='gray', linestyle='--', linewidth=0.8)
        ax.set_xlim(0, 40)
        ax.grid(True, alpha=0.2)
        ax.legend(fontsize=9, loc='upper left')
        ax.spines['top'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        
        print(f"  ✓ {season}: Regression R² = {reg.score(X, y):.3f}")
    
    plt.tight_layout()
    output_fig = FIGURE_DIR / "schlenker_roberts_temperature_yield.png"
    plt.savefig(output_fig, dpi=300, bbox_inches='tight')
    print(f"\n✓ Figure saved: {output_fig}")


# =============================================================================
# RUN
# =============================================================================

if __name__ == "__main__":
    main()
