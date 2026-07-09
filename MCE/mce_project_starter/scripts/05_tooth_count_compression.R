# scripts/05_tooth_count_compression.R
# ============================================================
# Compare measures by observed tooth/unit count
# ============================================================

nhanes_mce <- readRDS("data/derived/nhanes_mce_measures.rds")
nhanes_design <- make_nhanes_design(nhanes_mce)

# In this scaffold, n_observed_units is retained teeth when MCE is built using unit = "tooth".
nhanes_mce <- nhanes_mce |>
  dplyr::mutate(
    observed_unit_group = cut(
      n_observed_units,
      breaks = c(0, 9, 19, 27, Inf),
      labels = c("1-9", "10-19", "20-27", "28+"),
      right = TRUE
    )
  )

nhanes_design <- make_nhanes_design(nhanes_mce)

outcomes <- c("mean_cal", "mce", "extent_cal_ge_3", "extent_cal_ge_4", "extent_cal_ge_5", "extent_cal_ge_6")
outcomes <- intersect(outcomes, names(nhanes_mce))

by_tooth_count <- lapply(outcomes, function(outcome) {
  formula <- stats::as.formula(paste0("~", outcome))
  survey::svyby(
    formula,
    ~ observed_unit_group,
    nhanes_design,
    survey::svymean,
    na.rm = TRUE,
    keep.names = FALSE
  ) |>
    dplyr::mutate(measure = outcome)
}) |>
  dplyr::bind_rows()

save_csv(by_tooth_count, "outputs/tables/measures_by_observed_unit_count.csv")
