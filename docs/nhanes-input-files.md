# NHANES input files

The first iteration uses NHANES 2009--2014 adults aged 30 years and older.

## Required files

| File | Cycle | Role |
|---|---|---|
| `DEMO_F.xpt` | 2009--2010 | age, weights, strata, PSU |
| `DEMO_G.xpt` | 2011--2012 | age, weights, strata, PSU |
| `DEMO_H.xpt` | 2013--2014 | age, weights, strata, PSU |
| `GHB_F.xpt` | 2009--2010 | HbA1c |
| `GHB_G.xpt` | 2011--2012 | HbA1c |
| `GHB_H.xpt` | 2013--2014 | HbA1c |
| `OHXPER_F.xpt` | 2009--2010 | periodontal CAL site measures |
| `OHXPER_G.xpt` | 2011--2012 | periodontal CAL site measures |
| `OHXPER_H.xpt` | 2013--2014 | periodontal CAL site measures |

## Expected directory

```text
data/raw/
```

Raw data are ignored by Git.

## Pooled survey weights

For 2009--2014 analyses combining three two-year cycles, the scripts create:

```r
WTMEC6YR = WTMEC2YR / 3
```

This is used in the missing-tooth accumulation curve for design-based uncertainty.
