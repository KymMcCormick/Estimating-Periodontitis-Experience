# scripts/06_hba1c_alignment.R
# ============================================================
# Associational coherence with HbA1c
# ============================================================

nhanes_mce <- readRDS("data/derived/nhanes_mce_measures.rds") |>
  dplyr::mutate(age_c30 = age - 30)

if (!"hba1c" %in% names(nhanes_mce)) {
  stop("HbA1c column not found. Add hba1c to the core data or skip this script.", call. = FALSE)
}

nhanes_design <- make_nhanes_design(nhanes_mce)

measures <- c("mean_cal", "mce", "extent_cal_ge_3", "extent_cal_ge_4", "extent_cal_ge_5", "extent_cal_ge_6")
measures <- intersect(measures, names(nhanes_mce))

hba1c_models <- lapply(measures, function(measure) {
  formula <- stats::as.formula(paste0("hba1c ~ ", measure, " + age_c30 + I(age_c30^2)"))
  model <- survey::svyglm(formula, design = nhanes_design, family = gaussian())
  coef_table <- summary(model)$coefficients

  data.frame(
    measure = measure,
    term = rownames(coef_table),
    estimate = coef_table[, "Estimate"],
    se = coef_table[, "Std. Error"],
    p_value = coef_table[, "Pr(>|t|)"],
    row.names = NULL
  )
}) |>
  dplyr::bind_rows()

save_csv(hba1c_models, "outputs/tables/hba1c_alignment_models.csv")
