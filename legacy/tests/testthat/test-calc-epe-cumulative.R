test_that("calc_epe_cumulative returns expected simple values", {
  source("R/calc_epe_cumulative.R")

  cal_df <- data.frame(
    CAL_02 = c(2, NA),
    CAL_03 = c(4, NA),
    CAL_04 = c(NA, 6)
  )

  pi_fun <- function(age) rep(0.5, length(age))

  out <- calc_epe_cumulative(
    cal_df = cal_df,
    age_vec = c(40, 70),
    pi_fun = pi_fun,
    advanced_cal = 10
  )

  expect_equal(out$n_missing_teeth, c(1, 2))
  expect_equal(out$observed_retained_burden, c(6, 6))
  expect_equal(out$expected_missing_burden, c(5, 10))
  expect_equal(out$epe_total, c(11, 16))
  expect_equal(out$epe_mean_tooth, c(11 / 3, 16 / 3))
})
