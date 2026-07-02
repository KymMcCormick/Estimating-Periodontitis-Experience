test_that("weighted_pava returns non-decreasing fitted values", {
  source("R/weighted_pava.R")

  out <- weighted_pava(c(1, 3, 2, 4), w = c(1, 1, 1, 1))

  expect_true(all(diff(out) >= 0))
})
