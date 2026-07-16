#!/usr/bin/env python3
"""
================================================================================
SCRIPT 5: Merge Climate Panel with Yield Data (Bangladesh Rice)
================================================================================

Purpose: Create final regression-ready dataset

Input:
  - data/seasonal_panels/panel_bangladesh_rice_full.parquet
  - data/agriculturedata_bangladesh/crop/df_rice.csv

Output:
  - data/Regression_data/bangladesh_rice_regression_panel.csv
  - data/Regression_data/bangladesh_rice_regression_panel.parquet

Method:
  1. Load climate panel (6,336 rows: districts × years × seasons × periods)
  2. Load BBS yield data (2013-2023)
  3. Harmonize district names (climate → yield format)
  4. Map season names ("Rice.Boro" → "Boro")
  5. Merge on (district, year, season)
  6. Calculate yield metrics (yield_per_ha, log_yield)
  7. Aggregate periods (plant/between/harvest) → full season
  8. Add robustness variables (lags, interactions, polynomials, quartiles)

Author: Taky Tahmid
Date: 2026-04-05
================================================================================
"""

from pathlib import Path
import sys
import pandas as pd
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR.parent / "utils"))
from district_names import harmonize_district, district_key  # noqa: E402

# =============================================================================
# CONFIGURATION
# =============================================================================

PROJECT_ROOT = SCRIPT_DIR.parent.parent  # code/pipeline/ -> code/ -> project root

CLIMATE_PANEL = PROJECT_ROOT / "data/seasonal_panels/panel_bangladesh_rice_full.parquet"
YIELD_DATA = PROJECT_ROOT / "data/agriculturedata_bangladesh/crop/df_rice.csv"
OUTPUT_DIR = PROJECT_ROOT / "data/Regression_data"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Season name mapping (climate → yield format)
SEASON_MAPPING = {
    'Rice.Boro': 'Boro',
    'Rice.Aus': 'Aus',
    'Rice.Aman': 'Aman',
}

# Climate variables to aggregate across periods
CLIMATE_VARS = [
    "edd_0", "edd_4", "edd_8", "edd_12", "edd_28", "edd_30", "edd_32", "edd_35",
    "gdd_8_30", "gdd_8_32", "gdd_8_35", "gdd_10_30", "gdd_10_32", "gdd_10_35",
    "gdd_12_30", "gdd_12_32", "gdd_12_35", "gdd_15_30", "gdd_15_32", "gdd_15_35",
    "gdd_20_30", "gdd_20_32", "gdd_20_35",
    "hdd_10", "hdd_15",
    "pr1", "pr2"
]

# =============================================================================
# LOAD AND PROCESS DATA
# =============================================================================

def load_and_prepare_climate() -> pd.DataFrame:
    """Load climate panel and prepare for merging."""
    
    print("Loading climate panel...")
    df = pd.read_parquet(CLIMATE_PANEL)
    print(f"  → {len(df)} rows loaded")
    
    # Harmonize district names
    print("  Harmonizing district names...")
    df['district'] = df['adm_name'].map(harmonize_district)
    
    # Map season names
    print("  Mapping season names...")
    df['season'] = df['season'].map(SEASON_MAPPING)
    
    # Aggregate across periods (plant + between + harvest = full season)
    print("  Aggregating periods to full season...")
    
    # Group by district, year, season and sum climate indices
    agg_dict = {var: 'sum' for var in CLIMATE_VARS}
    agg_dict['days'] = 'sum'  # Total growing days
    agg_dict['adm_code'] = 'first'
    agg_dict['irrigation'] = 'first'
    agg_dict['country'] = 'first'
    
    df_seasonal = df.groupby(['district', 'year', 'season']).agg(agg_dict).reset_index()
    
    print(f"  → {len(df_seasonal)} seasonal observations (after aggregating periods)")
    
    return df_seasonal


def load_and_prepare_yield() -> pd.DataFrame:
    """Load yield data and prepare for merging."""
    
    print("\nLoading yield data...")
    df = pd.read_csv(YIELD_DATA)
    print(f"  → {len(df)} rows loaded")
    
    # Filter to 2013-2025 (extend as more years become available)
    df = df[(df['Year'] >= 2013) & (df['Year'] <= 2025)]
    print(f"  → {len(df)} rows after filtering to 2013-2025")
    
    # Rename columns to match
    df = df.rename(columns={
        'District': 'district',
        'Year': 'year',
        'Crop_type': 'season',
        'Area': 'area',
        'Production': 'production',
    })
    
    # Harmonize district names to BBS labels (fixes Cox's Bazar apostrophe mismatch)
    df['district'] = df['district'].map(harmonize_district)

    # Calculate yield
    df['yield_per_ha'] = df['production'] / df['area']
    df['log_yield'] = np.log(df['yield_per_ha'])
    
    # Remove inf/-inf from log(0)
    df = df.replace([np.inf, -np.inf], np.nan)
    
    print(f"  → Yield calculated for {df['yield_per_ha'].notna().sum()} observations")
    
    return df[['district', 'year', 'season', 'area', 'production', 'yield_per_ha', 'log_yield']]


# =============================================================================
# MERGE AND CREATE ROBUSTNESS VARIABLES
# =============================================================================

def create_robustness_variables(df: pd.DataFrame) -> pd.DataFrame:
    """Add robustness variables for regressions."""
    
    print("\nCreating robustness variables...")
    
    # Sort for lag operations
    df = df.sort_values(['district', 'season', 'year']).reset_index(drop=True)
    
    # First differences (by district × season)
    df['diff_log_yield'] = df.groupby(['district', 'season'])['log_yield'].diff()
    df['diff_gdd_10_35'] = df.groupby(['district', 'season'])['gdd_10_35'].diff()
    df['diff_edd_28'] = df.groupby(['district', 'season'])['edd_28'].diff()
    df['diff_edd_30'] = df.groupby(['district', 'season'])['edd_30'].diff()
    df['diff_edd_32'] = df.groupby(['district', 'season'])['edd_32'].diff()
    df['diff_edd_35'] = df.groupby(['district', 'season'])['edd_35'].diff()
    df['diff_pr1'] = df.groupby(['district', 'season'])['pr1'].diff()
    df['diff_pr2'] = df.groupby(['district', 'season'])['pr2'].diff()
    
    # Interactions
    df['gdd_edd_interaction'] = df['gdd_10_35'] * df['edd_30']
    df['gdd_pr1_interaction'] = df['gdd_10_35'] * df['pr1']
    df['edd_pr1_interaction'] = df['edd_30'] * df['pr1']
    
    # Polynomials
    df['gdd_10_35_sq'] = df['gdd_10_35'] ** 2
    df['edd_28_sq'] = df['edd_28'] ** 2
    df['edd_30_sq'] = df['edd_30'] ** 2
    df['edd_32_sq'] = df['edd_32'] ** 2
    df['edd_35_sq'] = df['edd_35'] ** 2
    
    # Quartiles
    df['gdd_quartile'] = pd.qcut(df['gdd_10_35'], q=4, labels=False, duplicates='drop') + 1
    df['edd_28_quartile'] = pd.qcut(df['edd_28'], q=4, labels=False, duplicates='drop') + 1
    df['edd_30_quartile'] = pd.qcut(df['edd_30'], q=4, labels=False, duplicates='drop') + 1
    df['edd_32_quartile'] = pd.qcut(df['edd_32'], q=4, labels=False, duplicates='drop') + 1
    df['edd_35_quartile'] = pd.qcut(df['edd_35'], q=4, labels=False, duplicates='drop') + 1
    df['pr1_quartile'] = pd.qcut(df['pr1'], q=4, labels=False, duplicates='drop') + 1
    
    # Binary indicators (above median)
    df['high_gdd'] = (df['gdd_10_35'] > df.groupby(['district', 'season'])['gdd_10_35'].transform('median')).astype(int)
    df['high_edd_28'] = (df['edd_28'] > df.groupby(['district', 'season'])['edd_28'].transform('median')).astype(int)
    df['high_edd_30'] = (df['edd_30'] > df.groupby(['district', 'season'])['edd_30'].transform('median')).astype(int)
    df['high_edd_32'] = (df['edd_32'] > df.groupby(['district', 'season'])['edd_32'].transform('median')).astype(int)
    df['high_edd_35'] = (df['edd_35'] > df.groupby(['district', 'season'])['edd_35'].transform('median')).astype(int)
    df['high_pr1'] = (df['pr1'] > df.groupby(['district', 'season'])['pr1'].transform('median')).astype(int)
    
    print(f"  → {df.shape[1]} total variables (after adding robustness vars)")
    
    return df


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("\n" + "="*80)
    print(" MERGE CLIMATE PANEL WITH YIELD DATA")
    print("="*80)
    print()
    
    # Load data
    climate = load_and_prepare_climate()
    yield_df = load_and_prepare_yield()
    
    # Check merge keys
    print("\n" + "="*80)
    print(" MERGE DIAGNOSTICS")
    print("="*80)
    
    climate_keys = set(zip(climate['district'], climate['year'], climate['season']))
    yield_keys = set(zip(yield_df['district'], yield_df['year'], yield_df['season']))
    
    print(f"\nClimate observations: {len(climate_keys)}")
    print(f"Yield observations:   {len(yield_keys)}")
    print(f"Intersection:         {len(climate_keys & yield_keys)}")
    print(f"Climate only:         {len(climate_keys - yield_keys)}")
    print(f"Yield only:           {len(yield_keys - climate_keys)}")
    
    # Merge
    print("\n" + "="*80)
    print(" MERGING")
    print("="*80)
    
    merged = pd.merge(
        climate,
        yield_df,
        on=['district', 'year', 'season'],
        how='inner'
    )
    
    print(f"\n✓ Merged dataset: {merged.shape}")
    print(f"  Districts: {merged['district'].nunique()}")
    print(f"  Years: {sorted(merged['year'].unique())}")
    print(f"  Seasons: {sorted(merged['season'].unique())}")
    
    # Add robustness variables
    final = create_robustness_variables(merged)
    
    # Save outputs
    print("\n" + "="*80)
    print(" SAVING OUTPUTS")
    print("="*80)
    
    csv_output = OUTPUT_DIR / "bangladesh_rice_regression_panel.csv"
    parquet_output = OUTPUT_DIR / "bangladesh_rice_regression_panel.parquet"
    
    final.to_csv(csv_output, index=False)
    final.to_parquet(parquet_output, index=False)
    
    print(f"\n✓ CSV saved:     {csv_output}")
    print(f"✓ Parquet saved: {parquet_output}")
    
    # Summary statistics
    print("\n" + "="*80)
    print(" FINAL DATASET SUMMARY")
    print("="*80)
    
    print(f"\nShape: {final.shape}")
    print(f"Columns: {final.shape[1]}")
    print(f"\nKey variables summary:")
    print(final[['log_yield', 'edd_28', 'edd_30', 'edd_32', 'gdd_10_35', 'pr1', 'area']].describe())
    
    print(f"\nMissing values:")
    missing = final.isnull().sum()
    print(missing[missing > 0].sort_values(ascending=False))
    
    print(f"\nObservations by season:")
    print(final.groupby('season').size())
    
    print(f"\nObservations by year:")
    print(final.groupby('year').size())
    
    print("\n" + "="*80)
    print(" ✓ PIPELINE COMPLETE!")
    print("="*80)
    print()


if __name__ == "__main__":
    main()
