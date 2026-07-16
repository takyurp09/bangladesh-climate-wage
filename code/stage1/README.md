# Stage 1: Climate to Rice Yields

This stage estimates season-specific climate-yield relationships for Bangladesh rice production.

## Purpose

The Stage 1 scripts construct fitted rice-yield shocks that are used as the main exposure variable in the wage pass-through analysis.

## Main Files

1. `00_diagnostics.R`  
   Data quality checks and variable diagnostics.

2. `00b_single_var_diagnostics.R`  
   Single-variable checks for growing degree days and extreme degree days.

3. `01_main_regressions.R`  
   Main season-specific climate-yield regressions and fitted-value construction.

4. `02_robustness.R`  
   Alternative thresholds and robustness specifications.

5. `03_temperature_bins.R`  
   Temperature-bin robustness following the climate-econometrics literature.

6. `04_figures.R`, `05_maps.R`, `06_summary_stats.R`  
   Publication figures, maps, and summary statistics.

## Canonical Specification

The current public version follows the canonical July 2026 project specification:

- Outcome: log rice yield
- Exposure: season-specific GDD/EDD measures
- Fixed effects: year and district
- Clustering: district
- Unit: district by rice growing season

First-difference variants are retained as robustness checks where relevant.

## Output

The key output is fitted log yield, used by Stage 2. Generated outputs are intentionally excluded from the public repository.
