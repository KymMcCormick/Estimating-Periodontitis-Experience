# ============================================================
# 04_produce_supp_calibration_outputs.R
# Produce supplementary calibration tables and figures.
# ============================================================
# Reporting outputs:
#   outputs/reporting/supp_table_s1.csv
#   outputs/reporting/supp_table_s2_panel_a.csv
#   outputs/reporting/supp_table_s2_panel_b.csv
#   outputs/reporting/supp_table_s2_panel_c.csv
#   outputs/reporting/supp_figure_s1.png
#   outputs/reporting/supp_figure_s2.png
#   outputs/reporting/supp_figure_s3.png
#
# Supplementary Figure S4 is a conceptual schematic and is not generated
# from the analytic data, so it is deliberately outside this pipeline.

#source("00_setup.R")

if (!file.exists(ANALYTIC_DATA_PATH) || !file.exists(CALIBRATION_DATA_PATH)) {
  stop("Run 01_prepare_nhanes_core.R and 02_build_epe_analytic_data.R first.")
}

analytic <- readRDS(ANALYTIC_DATA_PATH) |> add_reporting_variables()
calibration <- readRDS(CALIBRATION_DATA_PATH)

model_order <- c(
  "All schedules",
  "Without Schedule A",
  "Without Schedule B",
  "Without Schedule C"
)

model_specs <- tibble::tribble(
  ~model, ~excluded_schedule,
  "All schedules", NA_character_,
  "Without Schedule A", "Schedule A",
  "Without Schedule B", "Schedule B",
  "Without Schedule C", "Schedule C"
)

incident_models <- purrr::map(
  model_specs$excluded_schedule,
  ~fit_incident_attribution(if (is.na(.x)) NULL else .x)
) |>
  purrr::set_names(model_specs$model)

incident_curves <- purrr::imap_dfr(incident_models, function(result, model) {
  result$curve |>
    dplyr::mutate(model = factor(model, levels = model_order), .before = 1)
})

cumulative_curves <- purrr::imap_dfr(incident_models, function(result, model) {
  make_cumulative_attribution(
    incident_curve = result$curve,
    missing_curve = calibration$missing_accumulation$curve
  ) |>
    dplyr::mutate(model = factor(model, levels = model_order), .before = 1)
})

model_datasets <- purrr::map(model_order, function(model) {
  cumulative_curve <- cumulative_curves |>
    dplyr::filter(as.character(.data$model) == .env$model)

  recalculate_epe_from_components(
    analytic,
    cumulative_curve = cumulative_curve,
    floor = PRIMARY_MISSING_CAL_FLOOR
  ) |>
    dplyr::mutate(model = factor(model, levels = model_order))
}) |>
  purrr::set_names(model_order)

primary_check <- model_datasets[["All schedules"]] |>
  dplyr::summarise(
    maximum_absolute_difference = max(abs(epe_model - epe), na.rm = TRUE),
    mean_absolute_difference = mean(abs(epe_model - epe), na.rm = TRUE)
  )

print(primary_check)
if (primary_check$maximum_absolute_difference > 1e-8) {
  stop(
    "The all-schedules recalculation does not reproduce the primary EPE. ",
    "Do not use the supplementary sensitivity outputs until the upstream ",
    "calibration objects and analytic dataset are rebuilt together."
  )
}


# ============================================================
# Supplementary Table S1. Selected incident estimates
# ============================================================

s1_ages <- c(18, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80)

supp_table_s1 <- incident_models[["All schedules"]]$curve |>
  dplyr::filter(age %in% s1_ages) |>
  dplyr::transmute(
    Age = age,
    `Source-adjusted` = format_percent_value(q_incident, 1),
    `Minimum source-specific` = format_percent_value(q_min_source, 1),
    `Maximum source-specific` = format_percent_value(q_max_source, 1)
  )

save_report_csv(supp_table_s1, "supp_table_s1.csv")


# ============================================================
# Supplementary Table S2, Panel A. Incident probabilities
# ============================================================

s2_ages <- c(30, 40, 50, 60, 70, 80)

supp_table_s2_panel_a <- incident_curves |>
  dplyr::filter(age %in% s2_ages) |>
  dplyr::select(model, age, q_incident) |>
  tidyr::pivot_wider(names_from = age, values_from = q_incident) |>
  dplyr::transmute(
    Model = as.character(model),
    `Age 30 years` = format_percent_value(`30`, 1),
    `Age 40 years` = format_percent_value(`40`, 1),
    `Age 50 years` = format_percent_value(`50`, 1),
    `Age 60 years` = format_percent_value(`60`, 1),
    `Age 70 years` = format_percent_value(`70`, 1),
    `Age 80 years` = format_percent_value(`80`, 1)
  )

save_report_csv(supp_table_s2_panel_a, "supp_table_s2_panel_a.csv")


# ============================================================
# Supplementary Table S2, Panel B. Downstream EPE components
# ============================================================

summarise_model_components <- function(data, model) {
  analysis_data <- data |>
    dplyr::filter(
      is.finite(epe_model),
      is.finite(expected_missing_component_model),
      is.finite(WTMEC6YR),
      WTMEC6YR > 0,
      !is.na(SDMVPSU),
      !is.na(SDMVSTRA)
    )
  design <- make_survey_design(analysis_data)

  mean_epe_result <- survey_mean_ci(design, "epe_model")
  mean_missing_result <- survey_mean_ci(design, "expected_missing_component_model")

  selected_pi <- cumulative_curves |>
    dplyr::filter(as.character(.data$model) == .env$model, age %in% c(60, 70, 80)) |>
    dplyr::select(age, pi_cumulative) |>
    tidyr::pivot_wider(names_from = age, values_from = pi_cumulative)

  tibble::tibble(
    model = factor(model, levels = model_order),
    pi_60 = selected_pi$`60`,
    pi_70 = selected_pi$`70`,
    pi_80 = selected_pi$`80`,
    mean_epe = mean_epe_result$estimate,
    mean_epe_low = mean_epe_result$conf_low,
    mean_epe_high = mean_epe_result$conf_high,
    mean_missing = mean_missing_result$estimate,
    mean_missing_low = mean_missing_result$conf_low,
    mean_missing_high = mean_missing_result$conf_high
  )
}

panel_b_raw <- purrr::imap_dfr(model_datasets, summarise_model_components)

supp_table_s2_panel_b <- panel_b_raw |>
  dplyr::transmute(
    Model = as.character(model),
    `Cumulative attribution weight, age 60` = format_percent_value(pi_60, 1),
    `Cumulative attribution weight, age 70` = format_percent_value(pi_70, 1),
    `Cumulative attribution weight, age 80` = format_percent_value(pi_80, 1),
    `Mean EPE (95% CI)` = format_ci(mean_epe, mean_epe_low, mean_epe_high, 2),
    `Expected missing-tooth component (95% CI)` = format_ci(
      mean_missing,
      mean_missing_low,
      mean_missing_high,
      2
    )
  )

save_report_csv(supp_table_s2_panel_b, "supp_table_s2_panel_b.csv")


# ============================================================
# Supplementary Table S2, Panel C. Curvature and HbA1c alignment
# ============================================================

summarise_model_curvature <- function(data, model) {
  purrr::map_dfr(GLYCAEMIC_STATUS_LEVELS, function(status) {
    fit_late_life_contrasts(data, outcome = "epe_model", status = status) |>
      dplyr::filter(quantity == "late_life_curvature") |>
      dplyr::transmute(
        model = factor(model, levels = model_order),
        DiabetesStatus,
        estimate,
        conf_low,
        conf_high
      )
  })
}

curvature_raw <- purrr::imap_dfr(model_datasets, summarise_model_curvature)

curvature_wide <- curvature_raw |>
  dplyr::mutate(status = tolower(DiabetesStatus)) |>
  dplyr::select(model, status, estimate, conf_low, conf_high) |>
  tidyr::pivot_wider(
    names_from = status,
    values_from = c(estimate, conf_low, conf_high),
    names_glue = "{status}_{.value}"
  )

summarise_model_hba1c <- function(data, model) {
  analysis_data <- data |>
    dplyr::filter(
      is.finite(LBXGH),
      is.finite(Age),
      is.finite(epe_model),
      is.finite(mean_CAL),
      is.finite(extent_ge_4mm),
      is.finite(pct_sites_observed),
      is.finite(WTMEC6YR),
      WTMEC6YR > 0,
      !is.na(SDMVPSU),
      !is.na(SDMVSTRA)
    ) |>
    dplyr::mutate(
      age_c = Age - 55,
      z_epe_model = z_weighted(epe_model, WTMEC6YR)
    )

  design <- make_survey_design(analysis_data)
  fit <- survey::svyglm(
    LBXGH ~ z_epe_model + age_c + I(age_c^2),
    design = design
  )
  term <- tidy_svy_term(fit, "z_epe_model")

  tibble::tibble(
    model = factor(model, levels = model_order),
    beta = term$estimate,
    beta_low = term$conf_low,
    beta_high = term$conf_high,
    rho = survey_rank_correlation(analysis_data, "epe_model", "LBXGH")
  )
}

hba1c_raw <- purrr::imap_dfr(model_datasets, summarise_model_hba1c)

panel_c_raw <- curvature_wide |>
  dplyr::left_join(hba1c_raw, by = "model") |>
  dplyr::arrange(model)

supp_table_s2_panel_c <- panel_c_raw |>
  dplyr::transmute(
    Model = as.character(model),
    `Normal glycaemia curvature (95% CI)` = format_ci(
      normal_estimate, normal_conf_low, normal_conf_high, 3
    ),
    `Prediabetes curvature (95% CI)` = format_ci(
      prediabetes_estimate, prediabetes_conf_low, prediabetes_conf_high, 3
    ),
    `Diabetes curvature (95% CI)` = format_ci(
      diabetes_estimate, diabetes_conf_low, diabetes_conf_high, 3
    ),
    `Age-adjusted HbA1c β (95% CI)` = format_ci(beta, beta_low, beta_high, 3),
    `Rank correlation with HbA1c` = format_number(rho, 3)
  )

save_report_csv(supp_table_s2_panel_c, "supp_table_s2_panel_c.csv")


# ============================================================
# Supplementary Figure S1. Incident model and source schedules
# ============================================================

primary_model <- incident_models[["All schedules"]]

supp_figure_s1 <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = ATTRIBUTION_SCHEDULES,
    ggplot2::aes(
      x = age_low,
      xend = age_high_plot,
      y = p,
      yend = p,
      colour = schedule
    ),
    linewidth = 1,
    alpha = 0.60
  ) +
  ggplot2::geom_point(
    data = ATTRIBUTION_SCHEDULES,
    ggplot2::aes(
      x = age_mid,
      y = p,
      colour = schedule,
      shape = schedule,
      size = n
    ),
    alpha = 0.85
  ) +
  ggplot2::geom_line(
    data = primary_model$prediction_grid,
    ggplot2::aes(
      x = age_mid,
      y = q_source,
      colour = schedule,
      group = schedule
    ),
    linewidth = 0.7,
    alpha = 0.35
  ) +
  ggplot2::geom_ribbon(
    data = primary_model$curve,
    ggplot2::aes(x = age, ymin = q_min_source, ymax = q_max_source),
    fill = "grey65",
    alpha = 0.25
  ) +
  ggplot2::geom_line(
    data = primary_model$curve,
    ggplot2::aes(x = age, y = q_incident, linetype = "Source-adjusted incident model"),
    colour = "black",
    linewidth = 1.2
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0, 0.8, 0.2),
    limits = c(0, 0.85)
  ) +
  ggplot2::scale_x_continuous(breaks = seq(10, 80, 10), limits = c(10, 80)) +
  ggplot2::scale_size_continuous(range = c(2, 6)) +
  ggplot2::labs(
   # title = "Incident periodontal-attribution model",
   # subtitle = "Empirical age-band probabilities with source-adjusted fitted curve",
    x = "Age",
    y = "Probability extraction was attributed to periodontitis",
    colour = "Empirical schedule",
    shape = "Empirical schedule",
    size = "Age-band n",
    linetype = "Model"
  ) +
  study_theme(base_size = 12) +
  ggplot2::theme(legend.box = "vertical")

save_report_plot(supp_figure_s1, "supp_figure_s1.png", width = 10, height = 7)


# ============================================================
# Supplementary Figure S2. Leave-one-source-out incident curves
# ============================================================

supp_figure_s2 <- ggplot2::ggplot(
  incident_curves,
  ggplot2::aes(x = age, y = q_incident, linetype = model)
) +
  ggplot2::geom_line(colour = "black", linewidth = 0.95) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0, 0.8, 0.2),
    limits = c(0, 0.85)
  ) +
  ggplot2::scale_x_continuous(breaks = seq(10, 80, 10), limits = c(10, 80)) +
  ggplot2::labs(
   # title = "Leave-one-source-out sensitivity of incident attribution model",
    x = "Age",
    y = "Incident periodontal-attribution probability",
    linetype = "Model"
  ) +
  study_theme(base_size = 12)

save_report_plot(supp_figure_s2, "supp_figure_s2.png", width = 9, height = 6)


# ============================================================
# Supplementary Figure S3. Missing-tooth accumulation by age
# ============================================================

missing_points <- calibration$missing_accumulation$points
missing_curve <- calibration$missing_accumulation$curve

supp_figure_s3 <- ggplot2::ggplot() +
  ggplot2::geom_errorbar(
    data = missing_points,
    ggplot2::aes(
      x = age,
      ymin = pmax(mean_missing - se_missing, 0),
      ymax = mean_missing + se_missing
    ),
    width = 0.30,
    linewidth = 0.35,
    alpha = 0.55
  ) +
  ggplot2::geom_point(
    data = missing_points,
    ggplot2::aes(x = age, y = mean_missing),
    size = 1.8,
    alpha = 0.75
  ) +
  ggplot2::geom_line(
    data = missing_curve,
    ggplot2::aes(x = age, y = mean_missing_ordered),
    linewidth = 1.2
  ) +
  ggplot2::scale_x_continuous(breaks = seq(30, 80, 10)) +
  ggplot2::scale_y_continuous(breaks = seq(0, 15, 5)) +
  ggplot2::coord_cartesian(ylim = c(0, 15)) +
  ggplot2::labs(
   # title = "Estimated accumulation of missing teeth by age",
    #subtitle = "Points show survey-weighted age-specific means with SE; line shows smoothed non-decreasing curve",
    x = "Age",
    y = "Mean number of missing teeth"
  ) +
  study_theme(base_size = 12)

save_report_plot(supp_figure_s3, "supp_figure_s3.png", width = 9, height = 6)
