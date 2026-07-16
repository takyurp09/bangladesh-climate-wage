# Reproducibility and Version-Control Workflow

This repository is organized as a public code portfolio for a research pipeline. It is designed to show how the project is structured, version-controlled, and documented, while keeping raw data and generated outputs outside Git.

## Workflow

1. Clone the repository.

```bash
git clone https://github.com/takyurp09/bangladesh-climate-wage.git
cd bangladesh-climate-wage
```

2. Create the software environment.

```bash
conda env create -f environment.yml
conda activate bangladesh-climate-wage
```

3. Obtain the required input data from the original providers.

See `data/README.md` for data-source notes. Raw data are not committed to this repository.

4. Place data in the expected local directories.

The code assumes the following broad structure:

```text
data/
├── climate_data/
├── calendar_windows/
├── harvested_area_grids/
└── Regression_data/
```

For harvested-area grids, users may set:

```bash
export HARVEST_AREA_GRID_DIR=/path/to/harvested_area_grids
```

5. Run the pipeline in stages.

```bash
# Climate exposure construction
python code/pipeline/01_era5_download.py
python code/pipeline/02_climate_indices.py
python code/pipeline/03_crop_calendar.py
python code/pipeline/04_spatial_aggregation.py
python code/pipeline/05_merge_panel.py
python code/pipeline/06_temperature_bins.py

# Climate-to-yield analysis
Rscript code/stage1/00_diagnostics.R
Rscript code/stage1/01_main_regressions.R
Rscript code/stage1/02_robustness.R
Rscript code/stage1/03_temperature_bins.R
Rscript code/stage1/04_figures.R
Rscript code/stage1/05_maps.R
Rscript code/stage1/06_summary_stats.R

# Yield-to-wage analysis
Rscript code/stage2/00_dataprep.R
Rscript code/stage2/01_main_regressions.R
Rscript code/stage2/02_by_season.R
Rscript code/stage2/03_heterogeneity.R
Rscript code/stage2/04_robustness.R
Rscript code/stage2/05_figures.R
Rscript code/stage2/06_summary_stats.R
```

## Version-Control Rules

Commit code, documentation, configuration files, and small public-facing notes.

Do not commit:

- raw data
- derived datasets
- generated regression outputs
- generated figures and tables
- manuscript drafts
- local paths
- API keys or credentials
- cache folders and local environment files

The `.gitignore` file enforces these rules for common data, output, and cache formats.

## Recommended Commit Practice

Use small, descriptive commits:

```bash
git status
git add README.md docs/REPRODUCIBILITY.md
git commit -m "Document reproducibility workflow"
git push
```

Good commit messages describe the research workflow change, for example:

- `Add climate exposure aggregation script`
- `Document Stage 2 wage pass-through models`
- `Update data availability notes`
- `Add robustness diagnostics`

## Public Portfolio Scope

This public version is intended to demonstrate:

- reproducible project structure
- Git/GitHub version control
- R and Python workflow design
- transparent data availability practices
- separation of code from raw data and generated outputs

It is not intended to be a complete public release of all thesis materials.
