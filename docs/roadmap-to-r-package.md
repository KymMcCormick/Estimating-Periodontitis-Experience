# Roadmap to an R package

The repository is currently a reproducibility pipeline. The aim is to evolve it into an R package that other researchers can apply to their own datasets.

## Phase 1: Reproducible research repository

- [x] Preserve uploaded scripts in `scripts/legacy/`.
- [x] Create stable project structure.
- [x] Add helper functions for tooth-level maximum CAL and observed-site measures.
- [x] Document input files, outputs, and naming conventions.
- [ ] Run the full pipeline locally with NHANES raw files.
- [ ] Add model validation scripts.
- [ ] Add manuscript-ready tables and figures.

## Phase 2: Function extraction

Move repeated script sections into functions:

- [ ] `fit_incident_attribution_model()`
- [ ] `build_cumulative_attribution_lookup()`
- [ ] `predict_pi_cumulative()`
- [ ] `calc_epe_cumulative()`
- [ ] `classify_diabetes_status()`
- [ ] `prepare_nhanes_periodontal_data()`

## Phase 3: Package interface

Aim for a public-facing interface such as:

```r
library(perioepe)

epe <- estimate_epe(
  tooth_cal = tooth_cal_data,
  age = age,
  model = "iteration01",
  advanced_cal = 10
)
```

## Phase 4: Generalisation beyond NHANES

Future functions should accept user data with a standard minimal structure:

| Required field | Meaning |
|---|---|
| `id` | Participant identifier |
| `age` | Examination age |
| `CAL_02` ... `CAL_31` | Tooth-level observed maximum CAL; missing teeth coded `NA` |

Dataset-specific scripts can then map raw column names to this standard structure.

## Phase 5: Documentation and testing

- [ ] Add roxygen2 documentation.
- [ ] Add examples using simulated data.
- [ ] Add unit tests for EPE calculation and helper functions.
- [ ] Add package website with `pkgdown`.
- [ ] Add citation and release notes.
