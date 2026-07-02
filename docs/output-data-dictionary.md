# Output data dictionary

Main output: `data/processed/analytic_dataset_epe.csv` and `.rds`.

| Variable | Meaning |
|---|---|
| `SEQN` | NHANES participant identifier |
| `source_cycle` | NHANES cycle code |
| `observed_retained_burden` | Sum of observed tooth-level maximum CAL over retained modelled teeth |
| `n_present_teeth` | Number of modelled teeth with observed tooth-level CAL |
| `n_missing_teeth` | Number of modelled teeth without observed tooth-level CAL |
| `max_observed_cal` | Maximum tooth-level CAL observed for the participant |
| `assigned_missing_cal` | CAL value assigned to missing teeth before probabilistic weighting |
| `pi_cumulative` | Cumulative periodontal-attribution weight at participant age |
| `expected_missing_burden` | Expected periodontal burden contributed by missing teeth |
| `epe_total` | Total expected periodontitis experience burden |
| `epe_mean_tooth` | Mean expected burden per modelled tooth |
| `epe` | Default EPE variable; currently `epe_mean_tooth` |
| `measure_model` | Model label |
| `RIDAGEYR` | Age in years |
| `LBXGH` | HbA1c / glycohemoglobin |
| `WTMEC2YR` | NHANES two-year MEC weight |
| `SDMVPSU` | NHANES PSU |
| `SDMVSTRA` | NHANES stratum |
| `Age` | Copy of `RIDAGEYR` for plotting/reporting |
| `DiabetesStatus` | HbA1c-derived category: Normal, Prediabetes, Diabetes |
| `AgeGroup_Table2` | Age group used for tabulation |
| `mean_CAL_site` | Mean observed site-level CAL |
| `extent_ge_3mm` | Percent observed sites with CAL >= 3 mm |
| `extent_ge_4mm` | Percent observed sites with CAL >= 4 mm |
| `extent_ge_5mm` | Percent observed sites with CAL >= 5 mm |
| `extent_ge_6mm` | Percent observed sites with CAL >= 6 mm |
| `mean_extent_3to6mm` | Average of extent thresholds 3--6 mm |
| `n_sites_observed` | Number of observed periodontal CAL sites |
| `mean_CAL` | Copy of `mean_CAL_site` for convenience |

Note: some variables retain uppercase because they originate in the current analysis script. Future refactoring should migrate all derived variables to snake_case while preserving NHANES source names.
