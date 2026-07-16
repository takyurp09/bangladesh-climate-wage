# Stage 2: Rice Yield Shocks to Agricultural Wages

This stage estimates whether fitted rice-yield shocks pass through to agricultural wages in Bangladesh.

## Purpose

The Stage 2 scripts test whether wage pass-through differs by labor arrangement. Meal provision is used as a proxy for worker attachment:

- no meal: more casual labor
- one/two meals: intermediate arrangements
- three meals: more attached or resident labor

## Main Files

1. `00_dataprep.R`  
   Builds the wage-yield-control panel.

2. `01_main_regressions.R`  
   Main levels regressions and meal-type interaction models.

3. `02_by_season.R`  
   Season-specific results for Boro, Aus, and Aman.

4. `03_heterogeneity.R`  
   District and structural heterogeneity analyses.

5. `04_robustness.R`  
   Alternative specifications and inference checks.

6. `05_figures.R`, `06_summary_stats.R`, `08_maps.R`, `09_plots.R`  
   Figures, tables, and spatial summaries.

7. `10_meal_monotonicity.R` to `17_identification_diagnostics.R`  
   Additional mechanism, robustness, and identification diagnostics.

## Canonical Specification

The current public version uses:

```text
log_real_wage ~ log_yield_hat + gender + meal_type | year + District^growing_season
```

The main heterogeneity model interacts fitted log yield with meal type:

```text
log_real_wage ~ log_yield_hat * meal_type + gender | year + District^growing_season
```

Standard errors are clustered by district.

## Output

Generated tables, model objects, and figures are intentionally excluded from the public repository. Users with the source data can recreate them by running the scripts in order.
