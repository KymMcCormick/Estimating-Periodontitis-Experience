# Downstream validation and reporting analyses

This document describes scripts added after the first EPE construction pipeline.
The scripts assume that `scripts/02_build_epe_iteration01.R` has already created
`data/processed/analytic_dataset_epe.rds`.

## Run order

```r
source("scripts/00_setup.R")
source("scripts/01_prepare_nhanes_core.R")
source("scripts/02_build_epe_iteration01.R")
source("scripts/03_diagnose_incident_attribution_model.R")
source("scripts/04_summarise_epe_distribution.R")
source("scripts/05_assess_hba1c_associational_coherence.R")
source("scripts/06_assess_diabetes_discrimination_auc.R")
```

## Scripts

### `03_diagnose_incident_attribution_model.R`

Fits the source-adjusted incident periodontal-attribution curve `q_incident(a)`
from empirical extraction-reason schedules and writes diagnostic plots/tables.
This script is useful for inspecting the attribution model independently of the
full EPE construction script.

### `04_summarise_epe_distribution.R`

Creates a survey-weighted age-group table for periodontal observability, mean
CAL, extent CAL >=4 mm, and EPE.

### `05_assess_hba1c_associational_coherence.R`

Assesses whether EPE shows coherent external association with HbA1c relative to
observed-site periodontal measures. It estimates survey-weighted rank-correlation
approximations and survey-weighted regression models.

### `06_assess_diabetes_discrimination_auc.R`

Compares EPE and observed-site periodontal measures for discrimination of
diabetes status, defined as HbA1c >= 6.5%, with prediabetes excluded. AUC
confidence intervals use a stratified PSU bootstrap approximation.

The bootstrap size is set inside the script. For quick testing, reduce `B`; for
final analyses, use a larger value such as 1000 or 2000.

## Naming conventions

These scripts use the repository-standard variable name `epe`. Earlier local
scripts used `EPE`. A short compatibility block is included in downstream scripts
so that old local analytic datasets can still be read, but new code should use
`epe`.

## Data policy

Generated CSV tables and figures are written to `outputs/`. Depending on project
policy, these may or may not be committed. Raw NHANES data and processed RDS/CSV
analytic datasets should remain local and should not be committed.
