# Model specification: iteration 01

## Purpose

Estimate Expected Periodontitis Experience (EPE) from observed periodontal data when missing teeth are treated as informatively missing.

## Conceptual distinction

The current model distinguishes incident extraction attribution from cumulative missing-tooth attribution.

| Quantity | Interpretation |
|---|---|
| `q_incident(a)` | Probability that an extraction occurring at age `a` was attributed to periodontitis |
| `pi_cumulative(A)` | Expected periodontal-attribution weight for a tooth already missing at examination age `A` |

NHANES observes teeth that are already missing at the examination age. Therefore, the model uses `pi_cumulative(A)` rather than directly assigning `q_incident(A)` to each missing tooth.

## Current inputs

The first iteration uses:

- age,
- number of missing modelled teeth,
- observed tooth-level maximum CAL among retained teeth, and
- age-patterned extraction-reason evidence from external studies.

## EPE calculation

For each participant:

```text
observed_retained_burden = sum(observed tooth-level max CAL)
expected_missing_burden = n_missing_teeth * pi_cumulative(age) * assigned_missing_cal
epe_total = observed_retained_burden + expected_missing_burden
epe_mean_tooth = epe_total / n_teeth_modelled
```

The default `epe` variable currently equals `epe_mean_tooth`.

## Current missing-tooth burden rule

The first iteration assigns each missing tooth a potential CAL burden of:

```text
assigned_missing_cal = max(advanced_CAL, max_observed_cal)
```

where `advanced_CAL = 10` by default.

This is intentionally simple. Future iterations should evaluate sensitivity to this value and replace it with tooth-position-specific or distributional assumptions.

## Interpretation

EPE is an expected-burden measure, not an individual-level causal diagnosis. It estimates periodontal experience under stated assumptions about informative missingness.

## Planned extensions

1. Age-by-position missing-tooth attribution.
2. Tooth-level disease prediction using observed neighbours and homologous teeth.
3. Sensitivity analysis for assigned missing CAL.
4. External validation in other datasets.
5. Package functions for applying EPE to non-NHANES datasets.
