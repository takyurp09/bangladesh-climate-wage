# Data Availability

This repository is a public code portfolio version of an applied research pipeline. Raw and derived data are not committed because several inputs are large, licensed, or obtained from external data providers.

## Main Data Inputs

- ERA5 hourly temperature and precipitation from the Copernicus Climate Data Store
- Bangladesh rice yield and agricultural statistics from public statistical sources
- Bangladesh agricultural wage data from official statistical publications
- District boundaries and geospatial layers used for spatial aggregation
- Crop calendars and harvested-area weights used to aggregate climate exposure by rice growing season

## Reproducibility Notes

The code is organized so that users can inspect the full workflow:

1. `code/pipeline/` builds district-season climate exposure measures.
2. `code/stage1/` estimates climate-to-yield relationships.
3. `code/stage2/` estimates yield-shock pass-through to agricultural wages.

To fully reproduce the analysis, users must obtain the underlying data from the original providers and place them in the expected local directories. The public repository intentionally excludes raw data, generated regression outputs, and manuscript drafts.

## Privacy and Licensing

This public version includes only code and documentation intended for open review.
