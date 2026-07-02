# perioepe

**perioepe** is an early-stage repository for estimating **Expected Periodontitis Experience (EPE)** from observed periodontal examination data when tooth loss is informative.

The first iteration uses:

1. observed tooth-level maximum clinical attachment loss (CAL) among retained teeth,
2. missing tooth counts,
3. age, and
4. an age-specific periodontal-attribution model for extraction reasons.

The current model does **not** reconstruct individual disease histories or identify the true cause of loss for any given missing tooth. It estimates expected periodontal experience under transparent assumptions about informative missingness.

## Repository status

This is a research-code repository in transition toward an R package. The current scripts are designed to reproduce the first EPE analytic dataset from NHANES 2009--2014 files. The `R/` folder contains reusable functions that can later become package functions.

## Recommended run order

Place the required NHANES `.xpt` files in `data/raw/`, then run:

```r
source("scripts/00_setup.R")
source("scripts/01_prepare_nhanes_core.R")
source("scripts/02_build_epe_iteration01.R")
```

Or from a terminal:

```bash
Rscript scripts/01_prepare_nhanes_core.R
Rscript scripts/02_build_epe_iteration01.R
```

## Required NHANES input files

The first iteration expects the following files in `data/raw/`:

| File | NHANES cycle | Purpose |
|---|---:|---|
| `DEMO_F.xpt` | 2009--2010 | Demographics, age, survey design variables |
| `DEMO_G.xpt` | 2011--2012 | Demographics, age, survey design variables |
| `DEMO_H.xpt` | 2013--2014 | Demographics, age, survey design variables |
| `GHB_F.xpt` | 2009--2010 | HbA1c / glycohemoglobin |
| `GHB_G.xpt` | 2011--2012 | HbA1c / glycohemoglobin |
| `GHB_H.xpt` | 2013--2014 | HbA1c / glycohemoglobin |
| `OHXPER_F.xpt` | 2009--2010 | Periodontal examination |
| `OHXPER_G.xpt` | 2011--2012 | Periodontal examination |
| `OHXPER_H.xpt` | 2013--2014 | Periodontal examination |

NHANES raw data are not stored in this repository.

## Main outputs

The primary script writes:

| Output | Description |
|---|---|
| `data/processed/nhanes_demographics_hba1c.rds` | Participant-level age, HbA1c, cycle and survey-design variables |
| `data/processed/incident_periodontal_attribution_lookup.csv` | Age-specific incident periodontal-attribution curve, `q_incident` |
| `data/processed/cumulative_periodontal_attribution_lookup.csv` | Age-specific cumulative missing-tooth attribution curve, `pi_cumulative` |
| `data/processed/analytic_dataset_epe.csv` | Main analytic dataset |
| `data/processed/analytic_dataset_epe.rds` | Main analytic dataset, R format |
| `outputs/figures/` | Diagnostic figures for the attribution model and missing-tooth curve |

## Naming convention

This repository uses the following rule:

- **NHANES source variables** retain their original uppercase names, for example `SEQN`, `RIDAGEYR`, `WTMEC2YR`, `SDMVPSU`, `SDMVSTRA`, and `LBXGH`.
- **Derived variables** use `snake_case`, for example `n_missing_teeth`, `pi_cumulative`, `expected_missing_burden`, and `mean_extent_3to6mm`.
- **Participant-level datasets** use the prefix `nhanes_` or `analytic_`.
- **Model lookup tables** use the suffix `_lookup`.
- **Scripts** are numbered by pipeline order.

See `docs/naming-conventions.md` for details.

## Conceptual model

The first iteration distinguishes two quantities:

- `q_incident(a)`: probability that an extraction at age `a` is attributed to periodontitis.
- `pi_cumulative(A)`: expected periodontal-attribution weight for a tooth already missing at examination age `A`.

EPE then combines observed retained-tooth CAL with the expected contribution from missing teeth:

```text
epe_total = observed_retained_burden + expected_missing_burden
```

where:

```text
expected_missing_burden = n_missing_teeth * pi_cumulative(age) * assigned_missing_cal
```

The default analytic variable `epe` is currently set to `epe_mean_tooth`, which standardises `epe_total` by the number of modelled teeth.

## Project structure

```text
perioepe/
├── R/                         # reusable functions, future package functions
├── scripts/                   # executable research pipeline
│   └── legacy/                # uploaded scripts preserved for provenance
├── data/
│   ├── raw/                   # local NHANES XPT files; not tracked by Git
│   └── processed/             # generated analytic data; not tracked by Git
├── outputs/
│   ├── figures/               # generated figures; not tracked by Git
│   └── tables/                # generated tables; not tracked by Git
├── docs/                      # model notes, naming conventions, roadmap
└── tests/                     # early unit tests for package transition
```

## Roadmap

Planned extensions include:

1. validation against diabetes / HbA1c patterning,
2. age-by-tooth-position imputation models,
3. neighbourhood-informed tooth-level models,
4. sensitivity analyses for the assigned burden of missing teeth,
5. generalised functions for non-NHANES datasets, and
6. conversion to a documented R package.

See `docs/roadmap-to-r-package.md`.

## Important caution

EPE is an expectation-based measurement framework. It is intended to make assumptions explicit and testable; it should not be interpreted as a definitive reconstruction of the true reason for any individual tooth loss event.
