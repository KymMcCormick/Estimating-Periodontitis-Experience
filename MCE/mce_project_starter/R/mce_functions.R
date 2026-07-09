# R/mce_functions.R
# ============================================================
# Functions for constructing MCE and comparator measures
# ============================================================

check_required_columns <- function(data, required) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(
      "Missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

make_mce_measures <- function(data,
                              id_col = "seqn",
                              tooth_col = "tooth",
                              cal_col = "cal",
                              thresholds = c(3, 4, 5, 6),
                              unit = c("tooth", "site"),
                              covariates = c("age", "wtmec6yr", "sdmvstra", "sdmvpsu", "hba1c", "diabetes")) {
  unit <- match.arg(unit)

  check_required_columns(data, c(id_col, tooth_col, cal_col))

  available_covariates <- intersect(covariates, names(data))

  id_sym <- rlang::sym(id_col)
  tooth_sym <- rlang::sym(tooth_col)
  cal_sym <- rlang::sym(cal_col)

  if (unit == "tooth") {
    analysis_units <- data |>
      dplyr::filter(!is.na(!!cal_sym)) |>
      dplyr::group_by(!!id_sym, !!tooth_sym) |>
      dplyr::summarise(
        cal_unit = max(!!cal_sym, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    analysis_units <- data |>
      dplyr::filter(!is.na(!!cal_sym)) |>
      dplyr::mutate(cal_unit = !!cal_sym) |>
      dplyr::select(!!id_sym, cal_unit)
  }

  extent_list <- lapply(thresholds, function(threshold) {
    analysis_units |>
      dplyr::group_by(!!id_sym) |>
      dplyr::summarise(
        !!paste0("extent_cal_ge_", threshold) := mean(cal_unit >= threshold, na.rm = TRUE),
        .groups = "drop"
      )
  })

  extent_data <- Reduce(
    function(x, y) dplyr::left_join(x, y, by = id_col),
    extent_list
  )

  extent_cols <- paste0("extent_cal_ge_", thresholds)

  respondent_measures <- analysis_units |>
    dplyr::group_by(!!id_sym) |>
    dplyr::summarise(
      n_observed_units = dplyr::n(),
      mean_cal = mean(cal_unit, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::left_join(extent_data, by = id_col) |>
    dplyr::mutate(
      mce = rowMeans(dplyr::across(dplyr::all_of(extent_cols)), na.rm = TRUE),
      mce_percent = 100 * mce
    )

  covariate_data <- data |>
    dplyr::select(dplyr::all_of(c(id_col, available_covariates))) |>
    dplyr::distinct(!!id_sym, .keep_all = TRUE)

  respondent_measures |>
    dplyr::left_join(covariate_data, by = id_col)
}
