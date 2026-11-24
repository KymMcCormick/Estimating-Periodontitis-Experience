# scripts/04_age_trajectories.R
# Survey-weighted quadratic age trajectories for MCE and comparison periodontal measures.
#
# Purpose:
#   Compare how MCE, mean CAL, and single-threshold extent measures change with age.
#   This is the main age-patterning analysis after constructing the MCE analytic dataset.
#
# Inputs:
#   data/processed/analytic_dataset_mce.rds
#
# Outputs:
#   outputs/tables/mce_age_model_coefficients.csv
#   outputs/tables/mce_age_quadratic_tests.csv
#   outputs/tables/mce_age_trajectory_predictions.csv
#   outputs/tables/mce_age_specific_slopes.csv
#   outputs/figures/mce_age_trajectories_extent_scale.png
#   outputs/figures/mce_age_trajectories_all_faceted.png
#   outputs/figures/mce_age_specific_slopes.png

source("scripts/00_setup.R")

# ------------------------------------------------------------
# 1. Load analytic dataset
# ------------------------------------------------------------

mce_path <- "data/processed/analytic_dataset_mce.rds"

if (!file.exists(mce_path)) {
  stop(
    "Missing analytic dataset: ", mce_path, "\n",
    "Run scripts/02_build_mce_measure.R before this script."
  )
}

analytic_mce <- readRDS(mce_path)

required_vars <- c(
  "SEQN", "Age", "mce", "mean_CAL",
  "extent_ge_3mm", "extent_ge_4mm", "extent_ge_5mm", "extent_ge_6mm",
  "WTMEC6YR", "SDMVPSU", "SDMVSTRA"
)

missing_vars <- setdiff(required_vars, names(analytic_mce))

if (length(missing_vars) > 0) {
  stop(
    "The analytic MCE dataset is missing required variable(s):\n",
    paste0("  - ", missing_vars, collapse = "\n")
  )
}

analysis_data <- analytic_mce |>
  dplyr::mutate(
    Age = as.numeric(Age),
    age_c = Age - 30,
    age_c_sq = age_c^2
  ) |>
  dplyr::filter(
    !is.na(Age),
    Age >= 30,
    !is.na(WTMEC6YR),
    !is.na(SDMVPSU),
    !is.na(SDMVSTRA)
  )

message("Rows in MCE analytic dataset: ", nrow(analytic_mce))
message("Rows available for age trajectory models: ", nrow(analysis_data))

if (nrow(analysis_data) == 0) {
  stop("No complete rows available for survey-weighted age trajectory models.")
}

# ------------------------------------------------------------
# 2. Survey design
# ------------------------------------------------------------

nhanes_design <- survey::svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC6YR,
  nest = TRUE,
  data = analysis_data
)

# ------------------------------------------------------------
# 3. Measures to model
# ------------------------------------------------------------

measure_lookup <- tibble::tribble(
  ~measure,          ~measure_label,          ~measure_family,
  "mce",            "MCE",                  "Extent scale",
  "mean_CAL",       "Mean CAL",             "Mean CAL scale",
  "extent_ge_3mm",  "Extent CAL >=3 mm",    "Extent scale",
  "extent_ge_4mm",  "Extent CAL >=4 mm",    "Extent scale",
  "extent_ge_5mm",  "Extent CAL >=5 mm",    "Extent scale",
  "extent_ge_6mm",  "Extent CAL >=6 mm",    "Extent scale"
)

trajectory_measures <- measure_lookup$measure

# ------------------------------------------------------------
# 4. Helper functions
# ------------------------------------------------------------

fit_age_model <- function(outcome, design) {
  f <- stats::as.formula(paste0(outcome, " ~ age_c + age_c_sq"))
  survey::svyglm(f, design = design, family = stats::gaussian())
}

fit_linear_age_model <- function(outcome, design) {
  f <- stats::as.formula(paste0(outcome, " ~ age_c"))
  survey::svyglm(f, design = design, family = stats::gaussian())
}

predict_svyglm_identity <- function(model, newdata) {
  # Manual prediction avoids differences in predict.svyglm output across
  # survey package versions. The models use an identity link.
  tt <- stats::delete.response(stats::terms(model))
  X <- stats::model.matrix(tt, newdata)

  beta <- stats::coef(model)
  V <- stats::vcov(model)

  # Keep columns aligned even if model.matrix has attributes/ordering differences.
  X <- X[, names(beta), drop = FALSE]

  fit <- as.numeric(X %*% beta)
  se <- sqrt(diag(X %*% V %*% t(X)))

  tibble::tibble(
    estimate = fit,
    se = se,
    ci_low = estimate - 1.96 * se,
    ci_high = estimate + 1.96 * se
  )
}

calc_age_slope <- function(model, ages) {
  beta <- stats::coef(model)
  V <- stats::vcov(model)

  needed <- c("age_c", "age_c_sq")
  if (!all(needed %in% names(beta))) {
    stop("Model does not contain age_c and age_c_sq coefficients.")
  }

  purrr::map_dfr(ages, function(age_value) {
    age_c_value <- age_value - 30

    # slope = b_age_c + 2 * b_age_c_sq * age_c
    L <- rep(0, length(beta))
    names(L) <- names(beta)
    L["age_c"] <- 1
    L["age_c_sq"] <- 2 * age_c_value

    estimate <- as.numeric(sum(L * beta))
    se <- sqrt(as.numeric(t(L) %*% V %*% L))

    tibble::tibble(
      age = age_value,
      slope = estimate,
      se = se,
      ci_low = slope - 1.96 * se,
      ci_high = slope + 1.96 * se
    )
  })
}

extract_coefficients <- function(model, outcome) {
  broom::tidy(model) |>
    dplyr::mutate(measure = outcome, .before = 1)
}

extract_quadratic_test <- function(model_quad, model_linear, outcome) {
  # Primary test uses design-based Wald test for the quadratic term.
  wald <- tryCatch(
    survey::regTermTest(model_quad, ~age_c_sq),
    error = function(e) e
  )

  if (inherits(wald, "error")) {
    return(tibble::tibble(
      measure = outcome,
      test = "Quadratic age term",
      F = NA_real_,
      df_num = NA_real_,
      df_denom = NA_real_,
      p_value = NA_real_,
      note = conditionMessage(wald)
    ))
  }

  tibble::tibble(
    measure = outcome,
    test = "Quadratic age term",
    F = as.numeric(wald$Ftest),
    df_num = as.numeric(wald$df),
    df_denom = as.numeric(wald$ddf),
    p_value = as.numeric(wald$p),
    note = NA_character_
  )
}

# ------------------------------------------------------------
# 5. Fit age models
# ------------------------------------------------------------

age_models_quad <- purrr::set_names(trajectory_measures) |>
  purrr::map(fit_age_model, design = nhanes_design)

age_models_linear <- purrr::set_names(trajectory_measures) |>
  purrr::map(fit_linear_age_model, design = nhanes_design)

model_coefficients <- purrr::imap_dfr(
  age_models_quad,
  extract_coefficients
) |>
  dplyr::left_join(measure_lookup, by = "measure") |>
  dplyr::select(measure, measure_label, measure_family, dplyr::everything())

readr::write_csv(
  model_coefficients,
  "outputs/tables/mce_age_model_coefficients.csv"
)

message("Saved: outputs/tables/mce_age_model_coefficients.csv")

quadratic_tests <- purrr::map_dfr(
  trajectory_measures,
  function(m) extract_quadratic_test(
    model_quad = age_models_quad[[m]],
    model_linear = age_models_linear[[m]],
    outcome = m
  )
) |>
  dplyr::left_join(measure_lookup, by = "measure") |>
  dplyr::select(measure, measure_label, measure_family, dplyr::everything())

readr::write_csv(
  quadratic_tests,
  "outputs/tables/mce_age_quadratic_tests.csv"
)

message("Saved: outputs/tables/mce_age_quadratic_tests.csv")

# ------------------------------------------------------------
# 6. Predicted age trajectories
# ------------------------------------------------------------

age_grid <- tibble::tibble(
  Age = seq(30, 80, by = 1)
) |>
  dplyr::mutate(
    age_c = Age - 30,
    age_c_sq = age_c^2
  )

age_predictions <- purrr::imap_dfr(
  age_models_quad,
  function(model, outcome) {
    predict_svyglm_identity(model, age_grid) |>
      dplyr::bind_cols(age_grid |> dplyr::select(Age)) |>
      dplyr::mutate(measure = outcome, .before = 1)
  }
) |>
  dplyr::left_join(measure_lookup, by = "measure") |>
  dplyr::select(measure, measure_label, measure_family, Age, estimate, se, ci_low, ci_high)

readr::write_csv(
  age_predictions,
  "outputs/tables/mce_age_trajectory_predictions.csv"
)

message("Saved: outputs/tables/mce_age_trajectory_predictions.csv")

# Predictions at selected ages for easy reporting.
selected_age_predictions <- age_predictions |>
  dplyr::filter(Age %in% c(30, 40, 50, 60, 70, 80))

readr::write_csv(
  selected_age_predictions,
  "outputs/tables/mce_age_trajectory_predictions_selected_ages.csv"
)

message("Saved: outputs/tables/mce_age_trajectory_predictions_selected_ages.csv")

# ------------------------------------------------------------
# 7. Age-specific slopes
# ------------------------------------------------------------

slope_ages <- c(40, 50, 60, 70, 80)

age_slopes <- purrr::imap_dfr(
  age_models_quad,
  function(model, outcome) {
    calc_age_slope(model, slope_ages) |>
      dplyr::mutate(measure = outcome, .before = 1)
  }
) |>
  dplyr::left_join(measure_lookup, by = "measure") |>
  dplyr::select(measure, measure_label, measure_family, age, slope, se, ci_low, ci_high)

readr::write_csv(
  age_slopes,
  "outputs/tables/mce_age_specific_slopes.csv"
)

message("Saved: outputs/tables/mce_age_specific_slopes.csv")

# ------------------------------------------------------------
# 8. Figures
# ------------------------------------------------------------

extent_prediction_data <- age_predictions |>
  dplyr::filter(measure_family == "Extent scale") |>
  dplyr::mutate(
    measure_label = factor(
      measure_label,
      levels = c(
        "Extent CAL >=3 mm",
        "MCE",
        "Extent CAL >=4 mm",
        "Extent CAL >=5 mm",
        "Extent CAL >=6 mm"
      )
    )
  )

extent_trajectory_plot <- ggplot2::ggplot(
  extent_prediction_data,
  ggplot2::aes(x = Age, y = estimate, linetype = measure_label)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = ci_low, ymax = ci_high, group = measure_label),
    alpha = 0.10,
    colour = NA
  ) +
  ggplot2::geom_line(linewidth = 1.05) +
  ggplot2::scale_x_continuous(breaks = seq(30, 80, 10)) +
  ggplot2::coord_cartesian(ylim = c(0, NA)) +
  ggplot2::labs(
    title = "Survey-weighted age trajectories for MCE and extent measures",
    subtitle = "Quadratic age models, NHANES 2009-2014 adults aged >=30 years",
    x = "Age",
    y = "Predicted measure value",
    linetype = "Measure"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "bottom",
    plot.title = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  filename = "outputs/figures/mce_age_trajectories_extent_scale.png",
  plot = extent_trajectory_plot,
  width = 10,
  height = 7,
  dpi = 300
)

message("Saved: outputs/figures/mce_age_trajectories_extent_scale.png")

faceted_trajectory_plot <- ggplot2::ggplot(
  age_predictions,
  ggplot2::aes(x = Age, y = estimate)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = ci_low, ymax = ci_high),
    alpha = 0.12,
    colour = NA
  ) +
  ggplot2::geom_line(linewidth = 1.0) +
  ggplot2::facet_wrap(~measure_label, scales = "free_y", ncol = 3) +
  ggplot2::scale_x_continuous(breaks = seq(30, 80, 10)) +
  ggplot2::labs(
    title = "Survey-weighted quadratic age trajectories",
    subtitle = "Facets use free y-axes because mean CAL and extent measures are on different scales",
    x = "Age",
    y = "Predicted measure value"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  filename = "outputs/figures/mce_age_trajectories_all_faceted.png",
  plot = faceted_trajectory_plot,
  width = 11,
  height = 7.5,
  dpi = 300
)

message("Saved: outputs/figures/mce_age_trajectories_all_faceted.png")

slope_plot_data <- age_slopes |>
  dplyr::filter(measure_family == "Extent scale") |>
  dplyr::mutate(
    measure_label = factor(
      measure_label,
      levels = c(
        "Extent CAL >=3 mm",
        "MCE",
        "Extent CAL >=4 mm",
        "Extent CAL >=5 mm",
        "Extent CAL >=6 mm"
      )
    )
  )

slope_plot <- ggplot2::ggplot(
  slope_plot_data,
  ggplot2::aes(x = age, y = slope, linetype = measure_label)
) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.35, alpha = 0.6) +
  ggplot2::geom_line(linewidth = 1.0) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = ci_low, ymax = ci_high),
    width = 0.8,
    alpha = 0.75
  ) +
  ggplot2::scale_x_continuous(breaks = slope_ages) +
  ggplot2::labs(
    title = "Model-implied annual age slopes",
    subtitle = "Slope is the predicted annual change in measure value at selected ages",
    x = "Age",
    y = "Annual age slope",
    linetype = "Measure"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "bottom",
    plot.title = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  filename = "outputs/figures/mce_age_specific_slopes.png",
  plot = slope_plot,
  width = 10,
  height = 7,
  dpi = 300
)

message("Saved: outputs/figures/mce_age_specific_slopes.png")

# ------------------------------------------------------------
# 9. Console summary
# ------------------------------------------------------------

message("")
message("Selected age predictions:")
print(
  selected_age_predictions |>
    dplyr::select(measure_label, Age, estimate, se) |>
    dplyr::arrange(measure_label, Age),
  n = Inf
)

message("")
message("Age-specific slopes:")
print(
  age_slopes |>
    dplyr::select(measure_label, age, slope, se, ci_low, ci_high) |>
    dplyr::arrange(measure_label, age),
  n = Inf
)

message("")
message("Completed MCE age trajectory analysis.")
