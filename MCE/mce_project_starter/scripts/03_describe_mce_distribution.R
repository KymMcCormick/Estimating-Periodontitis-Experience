# scripts/03_describe_mce_distribution.R
# ============================================================
# Describe MCE distribution
# ============================================================

nhanes_mce <- readRDS("data/derived/nhanes_mce_measures.rds")
nhanes_design <- make_nhanes_design(nhanes_mce)

measure_vars <- c("mean_cal", "mce", "mce_percent", "extent_cal_ge_3", "extent_cal_ge_4", "extent_cal_ge_5", "extent_cal_ge_6")
measure_vars <- intersect(measure_vars, names(nhanes_mce))

summary_table <- lapply(measure_vars, function(var) {
  formula <- stats::as.formula(paste0("~", var))
  mean_est <- survey::svymean(formula, nhanes_design, na.rm = TRUE)
  quant_est <- survey::svyquantile(formula, nhanes_design, quantiles = c(0.25, 0.5, 0.75), na.rm = TRUE)

  data.frame(
    measure = var,
    mean = as.numeric(mean_est),
    se = as.numeric(SE(mean_est)),
    q25 = as.numeric(quant_est[[1]][1]),
    median = as.numeric(quant_est[[1]][2]),
    q75 = as.numeric(quant_est[[1]][3])
  )
}) |>
  dplyr::bind_rows()

save_csv(summary_table, "outputs/tables/mce_distribution_summary.csv")
