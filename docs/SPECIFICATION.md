# Econometric Specification (Canonical)

> **PERMANENT PROJECT DECISION (July 2026):** This repository uses the
> **LEVELS** specification for all production code, tables, and figures.
> First-differences is retained only for appendix robustness checks explicitly
> labelled "FD". Do not revert to FD for new outputs.

**Effective: July 2026.** All production code, tables, and figures use the **levels**
specification unless explicitly labelled as a first-differences robustness check.

## Stage 1 (climate → yield)

Season-specific panel fixed effects in **levels**:

| Season | Regression | Fixed effects |
|--------|------------|---------------|
| Boro | `log_yield ~ gdd_10_30 \| year + district` | year + district |
| Aus | `log_yield ~ edd_30 \| year + district` | year + district |
| Aman | `log_yield ~ gdd_10_30 + edd_30 \| year + district` | year + district |

Output: `yield_hat` = fitted log yield (levels), saved to
`output/stage1/fitted/yield_hat_levels_2017_*.csv`.

First-differences Stage 1 is retained only in `output/stage1/tables/appendix_fd.*`.

## Stage 2 (yield → wages)

| Model | Regression | Fixed effects |
|-------|------------|---------------|
| M3 | `log_real_wage ~ log_yield_hat + gender + meal_type` | `year + District^growing_season` |
| M5 | `log_real_wage ~ log_yield_hat * meal_type + gender` | `year + District^growing_season` |

- **Main regressor:** `log_yield_hat` (fitted log yield, levels)
- **Dependent variable:** `log_real_wage`
- **Clustering:** district
- **Sample:** 2017–2025 (64 districts)
- **Headline estimand (M5):** three-meal minus zero-meal interaction on `log_yield_hat`

## Deprecated (do not use for new outputs)

- `diff_real_wage ~ diff_log_yield_hat` (first-differences Stage 2)
- BDT/day pass-through as primary estimand
- Sample end year 2023 as default

## Data files

| File | Purpose |
|------|---------|
| `data/Regression_data/df_2_merged_levels.csv` | Canonical Stage 2 panel |
| `data/Regression_data/df_2_merged_v2.csv` | Alias (same content) |
| `data/Regression_data/irrigation_panel_2017_2025.csv` | District irrigation controls |

## Environment

**Always activate `crop_env` before any R or Python command:**

```bash
conda activate crop_env
```

All pipeline scripts, figure builds, and table generation assume this environment.
Do not run project R code outside `crop_env` unless debugging environment issues.

## Pipeline

```bash
conda activate crop_env
bash code/run_levels_pipeline.sh
Rscript code/paper/run_all_figures.R
Rscript code/paper/build_all_appendix_tables.R
```
