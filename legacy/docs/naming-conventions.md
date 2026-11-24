# Naming conventions

This repository should stay consistent as the imputation model develops.

## Files

Use lowercase names with underscores:

```text
01_prepare_nhanes_core.R
02_build_epe_iteration01.R
incident_periodontal_attribution_lookup.csv
cumulative_periodontal_attribution_lookup.csv
analytic_dataset_epe.csv
```

## Variables

### Keep NHANES source variables unchanged

Do not rename source variables until they are transformed into derived variables. Examples:

| Source variable | Meaning |
|---|---|
| `SEQN` | NHANES participant identifier |
| `RIDAGEYR` | Age in years |
| `WTMEC2YR` | MEC examination weight for a two-year cycle |
| `SDMVPSU` | Primary sampling unit |
| `SDMVSTRA` | Sampling stratum |
| `LBXGH` | HbA1c / glycohemoglobin |

### Use snake_case for derived variables

Examples:

| Derived variable | Meaning |
|---|---|
| `source_cycle` | NHANES cycle code: `F`, `G`, or `H` |
| `n_present_teeth` | Number of modelled teeth with observed CAL |
| `n_missing_teeth` | Number of modelled teeth without observed CAL |
| `q_incident` | Age-specific incident periodontal-attribution probability |
| `pi_cumulative` | Age-specific cumulative missing-tooth attribution weight |
| `assigned_missing_cal` | CAL burden assigned to a missing tooth before weighting |
| `expected_missing_burden` | Expected CAL contribution from missing teeth |
| `epe_total` | Observed retained burden plus expected missing burden |
| `epe_mean_tooth` | `epe_total` divided by number of modelled teeth |

## Tooth-level CAL columns

Tooth-level maximum CAL variables use:

```text
CAL_02, CAL_03, ..., CAL_31
```

Third molars are excluded by default: `01`, `16`, `17`, and `32`.

## Model iterations

Use explicit iteration labels rather than vague terms such as `final` or `new`:

```text
02_build_epe_iteration01.R
02_build_epe_iteration02_age_position.R
02_build_epe_iteration03_neighbourhood.R
```

## Output datasets

Use one row per participant unless otherwise stated. Suggested names:

```text
analytic_dataset_epe.csv
analytic_dataset_epe_iteration02.csv
tooth_level_prediction_dataset.csv
```
