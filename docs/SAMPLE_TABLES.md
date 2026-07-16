# Sample Tables

The following examples are small real-data extracts from processed outputs. They document data structure and variable definitions without releasing raw data, full panels, or complete result tables.

## Climate Exposure

| year | District | growing_season | gdd | edd | precip |
|---:|---|---|---:|---:|---:|
| 1970 | Bagerhat | Aman | 1964.4762 | 0.0000 | 898.3837 |
| 1970 | Bagerhat | Aus | 2166.7418 | 0.0000 | 1394.8040 |
| 1970 | Bagerhat | Boro | 1240.1576 | 0.0000 | 133.7011 |

## Rice Yield Panel

| District | year | growing_season | yield |
|---|---:|---|---:|
| Bagerhat | 2007 | Aman | 0.7286 |
| Bagerhat | 2007 | Aus | 0.7191 |
| Bagerhat | 2007 | Boro | 1.1668 |

## Wage Panel

| District | year | Month | growing_season | gender | meal_type | real_wage |
|---|---:|---|---|---|---|---:|
| Bagerhat | 2017 |  | Aman | Female | None | 325.8479 |
| Bagerhat | 2018 |  | Aman | Female | None | 310.9945 |
| Bagerhat | 2019 |  | Aman | Female | None | 377.8527 |

## Analysis Panel Extract

| District | year | growing_season | gender | meal_type | real_wage | log_yield_hat | ratio_double_cropped |
|---|---:|---|---|---|---:|---:|---:|
| Bagerhat | 2017 | Aman | Female | None | 325.8479 | -0.2173 | 0.2601 |
| Bagerhat | 2018 | Aman | Female | None | 310.9945 | -0.2027 | 0.2609 |
| Bagerhat | 2019 | Aman | Female | None | 377.8527 | -0.2051 | 0.2231 |

## Contract-Slope Diagnostic

| contract | estimate_bdt | se_bdt | p | n |
|---|---:|---:|---:|---:|
| Zero meals | -31.9623 | 5.6890 | 2.128e-08 | 2674 |
| Three meals | -46.7335 | 7.0084 | 3.501e-11 | 1684 |

## Interpretation

- `gdd`: growing degree days during the rice growing-season window.
- `edd`: extreme degree days during the rice growing-season window.
- `log_yield_hat`: fitted log yield from the climate-to-yield first stage.
- `meal_type`: labor-contract proxy used for wage pass-through heterogeneity.
- `estimate_bdt`: diagnostic slope in Bangladeshi taka units from a selected visualization table.
