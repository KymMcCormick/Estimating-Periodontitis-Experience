# R/calc_epe_cumulative.R
# Expected Periodontitis Experience under cumulative missing-tooth attribution.

calc_epe_cumulative <- function(
    cal_df,
    age_vec,
    pi_fun,
    advanced_cal = 10
) {
  cal_mat <- as.matrix(cal_df)

  observed_retained_burden <- rowSums(cal_mat, na.rm = TRUE)
  n_present_teeth <- rowSums(!is.na(cal_mat))
  n_missing_teeth <- rowSums(is.na(cal_mat))
  n_teeth_modelled <- ncol(cal_mat)

  max_observed_cal <- apply(
    cal_mat,
    1,
    function(x) {
      if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
    }
  )

  assigned_missing_cal <- pmax(
    advanced_cal,
    max_observed_cal,
    na.rm = TRUE
  )

  assigned_missing_cal[is.na(assigned_missing_cal)] <- advanced_cal

  pi_cumulative <- pi_fun(age_vec)

  expected_missing_burden <-
    n_missing_teeth * pi_cumulative * assigned_missing_cal

  epe_total <- observed_retained_burden + expected_missing_burden
  epe_mean_tooth <- epe_total / n_teeth_modelled

  tibble::tibble(
    observed_retained_burden = observed_retained_burden,
    n_present_teeth = n_present_teeth,
    n_missing_teeth = n_missing_teeth,
    max_observed_cal = max_observed_cal,
    assigned_missing_cal = assigned_missing_cal,
    pi_cumulative = pi_cumulative,
    expected_missing_burden = expected_missing_burden,
    epe_total = epe_total,
    epe_mean_tooth = epe_mean_tooth,
    epe = epe_mean_tooth
  )
}
