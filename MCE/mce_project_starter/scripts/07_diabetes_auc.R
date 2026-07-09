# scripts/07_diabetes_auc.R
# ============================================================
# Diabetes discrimination analysis
# ============================================================
# Note: survey-weighted AUC is not implemented here. This script
# provides a transparent unweighted screening AUC. Replace or extend
# with the exact AUC approach used in the manuscript.
# ============================================================

if (!requireNamespace("pROC", quietly = TRUE)) {
  stop("Install package 'pROC' before running diabetes AUC analysis.", call. = FALSE)
}

nhanes_mce <- readRDS("data/derived/nhanes_mce_measures.rds") |>
  dplyr::filter(!is.na(diabetes))

if (!"diabetes" %in% names(nhanes_mce)) {
  stop("Diabetes column not found. Add diabetes to the core data or skip this script.", call. = FALSE)
}

measures <- c("mean_cal", "mce", "extent_cal_ge_3", "extent_cal_ge_4", "extent_cal_ge_5", "extent_cal_ge_6")
measures <- intersect(measures, names(nhanes_mce))

auc_table <- lapply(measures, function(measure) {
  complete_data <- nhanes_mce |>
    dplyr::filter(!is.na(.data[[measure]]), !is.na(diabetes))

  roc_obj <- pROC::roc(
    response = complete_data$diabetes,
    predictor = complete_data[[measure]],
    quiet = TRUE
  )

  data.frame(
    measure = measure,
    n = nrow(complete_data),
    auc = as.numeric(pROC::auc(roc_obj))
  )
}) |>
  dplyr::bind_rows()

save_csv(auc_table, "outputs/tables/diabetes_auc_unweighted.csv")
