# scripts/04_age_trajectories.R
# ============================================================
# Survey-weighted quadratic age trajectories
# ============================================================

nhanes_mce <- readRDS("data/derived/nhanes_mce_measures.rds") |>
  dplyr::mutate(age_c30 = age - 30)

nhanes_design <- make_nhanes_design(nhanes_mce)

outcomes <- c("mean_cal", "mce", "extent_cal_ge_3", "extent_cal_ge_4", "extent_cal_ge_5", "extent_cal_ge_6")
outcomes <- intersect(outcomes, names(nhanes_mce))

fit_age_model <- function(outcome) {
  formula <- stats::as.formula(paste0(outcome, " ~ age_c30 + I(age_c30^2)"))
  survey::svyglm(formula, design = nhanes_design, family = gaussian())
}

age_models <- setNames(lapply(outcomes, fit_age_model), outcomes)

age_grid <- data.frame(age = seq(30, 80, by = 1)) |>
  dplyr::mutate(age_c30 = age - 30)

prediction_data <- lapply(names(age_models), function(outcome) {
  pred <- predict(age_models[[outcome]], newdata = age_grid, se.fit = TRUE)
  data.frame(
    outcome = outcome,
    age = age_grid$age,
    estimate = as.numeric(pred$fit),
    se = as.numeric(pred$se.fit)
  ) |>
    dplyr::mutate(
      lower = estimate - 1.96 * se,
      upper = estimate + 1.96 * se
    )
}) |>
  dplyr::bind_rows()

save_csv(prediction_data, "outputs/tables/age_trajectory_predictions.csv")

age_plot <- ggplot2::ggplot(prediction_data, ggplot2::aes(x = age, y = estimate)) +
  ggplot2::geom_line() +
  ggplot2::facet_wrap(~ outcome, scales = "free_y") +
  ggplot2::labs(
    x = "Age in years",
    y = "Predicted value",
    title = "Survey-weighted quadratic age trajectories"
  ) +
  ggplot2::theme_minimal()

save_plot(age_plot, "outputs/figures/age_trajectories.png", width = 9, height = 6)
