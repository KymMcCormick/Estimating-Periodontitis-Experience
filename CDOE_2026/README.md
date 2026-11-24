# Estimating Periodontitis Experience Under Informative Tooth Loss

Author: McCormick, K. M. 2026

## Overview

This directory contains the R analytical code accompanying:

McCormick, K. M., Guzzo, G., & Amarasena, N. *Estimating Periodontitis Experience Under Informative Tooth Loss*.  Community Dentistry and Oral Epidemiology.

Expected Periodontitis Experience (EPE) is a model-based framework for estimating population-level cumulative periodontal experience when tooth loss has removed evidence of prior periodontal destruction. EPE combines observed clinical attachment loss (CAL) on retained teeth with an expected periodontal contribution from missing tooth positions under explicit attribution and severity assumptions.

The primary implementation uses the maximum CAL across four interproximal sites per retained tooth, age-specific cumulative missing-tooth attribution weights calibrated using published extraction-cause studies and NHANES tooth-loss accumulation, and a fixed assigned CAL severity of 10 mm. Contributions are summed across 28 eligible tooth positions (excluding third molars) and divided by 28. In this analytical implementation, tooth positions without a valid interproximal CAL value are treated as missing.

EPE is intended for population-level research. It is not an individual diagnostic instrument, a validated attribution of individual tooth loss, or a reconstruction of individual disease histories.

## Version and scope

The scripts in `CDOE_2026/` implement the analytical specification accompanying the CDOE manuscript. Earlier implementations elsewhere in this repository may use different calibration or measurement definitions. Use this directory to reproduce the manuscript analyses.

## Data requirements

The analyses use publicly available National Health and Nutrition Examination Survey (NHANES) data for 2009–2014. Download the following SAS transport (XPT) files from the [NHANES website](https://wwwn.cdc.gov/nchs/nhanes/) and place them in `CDOE_2026/data/raw/`.

| Dataset | 2009–2010 | 2011–2012 | 2013–2014 |
|---|---|---|---|
| Demographics | DEMO_F.XPT | DEMO_G.XPT | DEMO_H.XPT |
| Glycohemoglobin | GHB_F.XPT | GHB_G.XPT | GHB_H.XPT |
| Periodontal examination | OHXPER_F.XPT | OHXPER_G.XPT | OHXPER_H.XPT |

The raw datasets are not distributed with this code. The empirical extraction-attribution schedules used for calibration are defined in `00_setup.R`.

## Software requirements

R 4.1.0 or later is required for the native pipe syntax. Install the required packages before running the analysis:

```r
install.packages(c(
  "dplyr", "tidyr", "purrr", "readr", "stringr", "tibble",
  "haven", "survey", "ggplot2", "patchwork", "scales"
))
```

## Reproducing the analyses

1. Clone or download this repository.
2. Open R or RStudio and set the working directory to `CDOE_2026/` (the directory containing `99_run_all.R`).
3. Create `data/raw/` if necessary and place the nine NHANES files in it.
4. Install the required packages and start a clean R session.
5. Run:

```r
source("99_run_all.R")
```

The master script executes the analysis in sequence, rebuilds the derived datasets and validates the expected reporting files. Output directories are created automatically. **The master script deletes and recreates `outputs/reporting/` at the beginning of each run.** Do not store unrelated files there.

Alternatively, from a terminal with R available on the system path, run `Rscript 99_run_all.R` while in `CDOE_2026/`.

## Script organisation

| Script | Purpose |
|---|---|
| `00_setup.R` | Shared dependencies, paths, study constants, attribution schedules and analytical functions |
| `01_prepare_nhanes_core.R` | Prepare participant demographics, survey variables and HbA1c |
| `02_build_epe_analytic_data_interproximal_corrected.R` | Construct the analytical dataset and primary calibration objects |
| `03_produce_main_outputs.R` | Generate main manuscript tables and figures |
| `04_produce_supp_calibration_outputs.R` | Generate supplementary calibration outputs |
| `05_produce_supp_table_s3.R` | Summarise age-specific EPE components |
| `06_produce_supp_tables_s4_s5.R` | Estimate late-life curvature, rank correlations and exploratory joint HbA1c models |
| `07_produce_supp_sensitivity_outputs.R` | Evaluate sensitivity to assigned CAL severity |
| `10_sensitivity_leave_one_source_out_downstream.R` | Evaluate sensitivity to extraction-attribution source selection |
| `99_run_all.R` | Execute and validate the complete analytical pipeline |

## Outputs and reproducibility checks

Intermediate analytical objects are saved in `data/derived/`. The pipeline produces 21 reporting files in `outputs/reporting/`, including:

- Main manuscript Tables 1–3 and Figures 1–2.
- Underlying numerical data for Figure 2.
- Supplementary Tables S1–S8 (Table S2 is exported as three separate panels).
- Supplementary Figures S1–S3 and S5–S6.

Supplementary Figure S4 is a methodological schematic included in the manuscript's supporting information and is not generated by the analytical pipeline.

Selected reference results for checking reproducibility:

| Quantity | Manuscript value |
|---|---:|
| Broader periodontal analytic sample | 10,714 |
| Primary HbA1c analytic sample | 10,324 |
| Incident periodontal attribution at age 30 | 12.0% |
| Incident periodontal attribution at age 50 | 50.6% |
| Incident periodontal attribution at age 80 | 61.7% |
| Cumulative missing-tooth attribution at age 80 | 45.3% |
| Age-adjusted EPE–HbA1c coefficient (percentage points per 1-SD increase) | 0.144 |
| Design-weighted EPE–HbA1c rank correlation | 0.291 |

Values are rounded to manuscript precision. Confidence intervals involving EPE are conditional on the fitted attribution model and assigned CAL severity and do not incorporate the full uncertainty in the external calibration evidence.

## Citation and related materials

If using this code or the EPE framework in research, please cite the accompanying CDOE manuscript. The analytical materials associated with the manuscript are also referenced at [OSF](https://osf.io/x2yma/).

## Licence

See the `LICENSE` file in the root of this repository for licensing terms.

## Contact

Kym M. McCormick  
Australian Research Centre for Population Oral Health (ARCPOH)  
Adelaide University, Adelaide, South Australia, Australia  
kym.mccormick@adelaide.edu.au
