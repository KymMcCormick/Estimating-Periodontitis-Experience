# MCE periodontal measurement project

This project scaffold separates the workflow into two parts:

1. Create respondent-level MCE and comparator periodontal measures from NHANES periodontal site-level data.
2. Run the main analyses comparing MCE with mean CAL and single-threshold extent measures.

## Expected derived input

The measure-building script expects a site-level file at:

`data/derived/nhanes_core_site_level.rds`

with at least these variables, or equivalent names that you map in `scripts/02_build_mce_measure.R`:

- `seqn`: participant identifier
- `age`: age in years
- `tooth`: tooth position identifier
- `cal`: clinical attachment loss value
- `wtmec6yr`: NHANES 6-year MEC examination weight
- `sdmvstra`: NHANES strata
- `sdmvpsu`: NHANES PSU

Optional but used in later scripts:

- `hba1c`: glycated haemoglobin
- `diabetes`: diabetes status, coded 0/1

## Run

From the project root:

```r
source("run_all.R")
```

## Definition used here

MCE is implemented as the average of respondent-level CAL extent measures across ordered CAL thresholds. With thresholds 3, 4, 5, and 6 mm:

`MCE = mean(extent_CAL3, extent_CAL4, extent_CAL5, extent_CAL6)`

Equivalently, it is the mean proportion of ordered CAL thresholds crossed across observed periodontal units. The default unit in this scaffold is the retained tooth, using maximum CAL per tooth. You can switch to site-level calculation in `scripts/02_build_mce_measure.R`.
