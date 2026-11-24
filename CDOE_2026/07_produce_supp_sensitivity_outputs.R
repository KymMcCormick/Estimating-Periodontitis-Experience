# ============================================================
# 07_produce_supp_sensitivity_outputs.R
# Sensitivity to the minimum assigned CAL value for missing teeth.
# ============================================================
# Reporting outputs:
#   outputs/reporting/supp_table_s6.csv
#   outputs/reporting/supp_table_s7.csv
#   outputs/reporting/supp_table_s8.csv
#   outputs/reporting/supp_figure_s5.png
#   outputs/reporting/supp_figure_s6.png
#
# The discontinued diabetes-status AUC analysis is not calculated here and
# no AUC-related table, figure, model, helper, or package is retained.

source("00_setup.R")

if (!file.exists(ANALYTIC_DATA_PATH)) {
  stop("Run 01_prepare_nhanes_core.R and 02_build_epe_analytic_data.R first.")
}

calibration <- readRDS(CALIBRATION_DATA_PATH)
analytic <- readRDS(ANALYTIC_DATA_PATH) |> add_reporting_variables()

make_floor_data <- function(floor) {
  
  recalculate_epe_from_components(
    data = analytic,
    cumulative_curve = calibration$cumulative,
    floor = floor
  ) |>
    dplyr::mutate(
      cal_floor = floor,
      assigned_missing_cal_sens = assigned_missing_cal_model,
      expected_missing_component_sens = expected_missing_component_model,
      observed_retained_component_sens = observed_retained_component_model,
      epe_sens = epe_model
    )
}

summarise_one_floor <- function(floor) {
  floor_data <- make_floor_data(floor) |>
    dplyr::filter(
      is.finite(epe_sens),
      is.finite(WTMEC6YR),
      WTMEC6YR > 0,
      !is.na(SDMVPSU),
      !is.na(SDMVSTRA)
    )

  if (floor == PRIMARY_MISSING_CAL_FLOOR) {
    maximum_difference <- max(abs(floor_data$epe_sens - floor_data$epe), na.rm = TRUE)
    message("Maximum EPE reproduction difference at 10 mm: ", signif(maximum_difference, 5))
    if (maximum_difference > 1e-8) {
      stop("The 10-mm sensitivity specification does not reproduce primary EPE.")
    }
  }

  broad_design <- make_survey_design(floor_data)
  mean_assigned_result <- survey_mean_ci(broad_design, "assigned_missing_cal_sens")
  mean_missing_result <- survey_mean_ci(broad_design, "expected_missing_component_sens")
  mean_epe_result <- survey_mean_ci(broad_design, "epe_sens")

  # Use the same complete-case definition as main Table 3 so that the
  # 10-mm row reproduces the reported primary HbA1c coefficients.
  hba1c_data <- floor_data |>
    dplyr::filter(
      is.finite(LBXGH),
      is.finite(Age),
      is.finite(epe_sens),
      is.finite(mean_CAL),
      is.finite(extent_ge_4mm),
      is.finite(pct_sites_observed)
    ) |>
    dplyr::mutate(
      age_c = Age - 55,
      z_epe_sens = z_weighted(epe_sens, WTMEC6YR)
    )

  hba1c_design <- make_survey_design(hba1c_data)
  fit_unadjusted <- survey::svyglm(LBXGH ~ z_epe_sens, design = hba1c_design)
  fit_age_adjusted <- survey::svyglm(
    LBXGH ~ z_epe_sens + age_c + I(age_c^2),
    design = hba1c_design
  )

  beta_unadjusted_result <- tidy_svy_term(fit_unadjusted, "z_epe_sens")
  beta_age_adjusted_result <- tidy_svy_term(fit_age_adjusted, "z_epe_sens")
  rho <- survey_rank_correlation(hba1c_data, "epe_sens", "LBXGH")

  curvature <- purrr::map_dfr(GLYCAEMIC_STATUS_LEVELS, function(status) {
    fit_late_life_contrasts(
      floor_data,
      outcome = "epe_sens",
      status = status,
      z_critical = 1.96
    ) |>
      dplyr::mutate(cal_floor = floor, .before = 1)
  })

  summary <- tibble::tibble(
    cal_floor = floor,
    mean_assigned = mean_assigned_result$estimate,
    mean_assigned_low = mean_assigned_result$conf_low,
    mean_assigned_high = mean_assigned_result$conf_high,
    mean_missing = mean_missing_result$estimate,
    mean_missing_low = mean_missing_result$conf_low,
    mean_missing_high = mean_missing_result$conf_high,
    mean_epe = mean_epe_result$estimate,
    mean_epe_low = mean_epe_result$conf_low,
    mean_epe_high = mean_epe_result$conf_high,
    beta_unadjusted = beta_unadjusted_result$estimate,
    beta_unadjusted_low = beta_unadjusted_result$conf_low,
    beta_unadjusted_high = beta_unadjusted_result$conf_high,
    beta_age_adjusted = beta_age_adjusted_result$estimate,
    beta_age_adjusted_low = beta_age_adjusted_result$conf_low,
    beta_age_adjusted_high = beta_age_adjusted_result$conf_high,
    rho = rho
  )

  list(summary = summary, curvature = curvature)
}

floor_results <- purrr::map(MISSING_CAL_FLOORS, summarise_one_floor)
summary_raw <- purrr::map_dfr(floor_results, "summary")
curvature_raw <- purrr::map_dfr(floor_results, "curvature")


# ============================================================
# Supplementary Table S6. EPE components by CAL floor
# ============================================================

supp_table_s6 <- summary_raw |>
  dplyr::transmute(
    `CAL floor (mm)` = cal_floor,
    `Mean assigned missing CAL (mm)` = format_ci(
      mean_assigned, mean_assigned_low, mean_assigned_high, 2
    ),
    `Expected missing-tooth component` = format_ci(
      mean_missing, mean_missing_low, mean_missing_high, 2
    ),
    `EPE, mean-tooth scale` = format_ci(
      mean_epe, mean_epe_low, mean_epe_high, 2
    )
  )

save_report_csv(supp_table_s6, "supp_table_s6.csv")


# ============================================================
# Supplementary Table S7. HbA1c alignment by CAL floor
# ============================================================

supp_table_s7 <- summary_raw |>
  dplyr::transmute(
    `CAL floor (mm)` = cal_floor,
    `HbA1c β, unadjusted (95% CI)` = format_ci(
      beta_unadjusted, beta_unadjusted_low, beta_unadjusted_high, 3
    ),
    `HbA1c β, age-adjusted (95% CI)` = format_ci(
      beta_age_adjusted, beta_age_adjusted_low, beta_age_adjusted_high, 3
    ),
    `Rank correlation with HbA1c` = format_number(rho, 3)
  )

save_report_csv(supp_table_s7, "supp_table_s7.csv")


# ============================================================
# Supplementary Table S8. Late-life curvature by CAL floor
# ============================================================

s8_wide <- curvature_raw |>
  dplyr::select(
    cal_floor,
    DiabetesStatus,
    quantity,
    estimate,
    conf_low,
    conf_high
  ) |>
  tidyr::pivot_wider(
    names_from = quantity,
    values_from = c(estimate, conf_low, conf_high)
  ) |>
  dplyr::mutate(
    DiabetesStatus = factor(DiabetesStatus, levels = GLYCAEMIC_STATUS_LEVELS)
  ) |>
  dplyr::arrange(cal_floor, DiabetesStatus)

supp_table_s8 <- s8_wide |>
  dplyr::transmute(
    `CAL floor (mm)` = cal_floor,
    `Glycaemic status` = GLYCAEMIC_STATUS_LABELS[as.character(DiabetesStatus)],
    `Fitted change, 60–70 years (95% CI)` = format_ci(
      estimate_change_60_70,
      conf_low_change_60_70,
      conf_high_change_60_70,
      3
    ),
    `Fitted change, 70–80 years (95% CI)` = format_ci(
      estimate_change_70_80,
      conf_low_change_70_80,
      conf_high_change_70_80,
      3
    ),
    `Late-life curvature (95% CI)` = format_ci(
      estimate_late_life_curvature,
      conf_low_late_life_curvature,
      conf_high_late_life_curvature,
      3
    )
  )

save_report_csv(supp_table_s8, "supp_table_s8.csv")


# ============================================================
# Supplementary Figure S5. Age-adjusted HbA1c coefficient
# ============================================================

supp_figure_s5 <- ggplot2::ggplot(
  summary_raw,
  ggplot2::aes(x = cal_floor, y = beta_age_adjusted)
) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = beta_age_adjusted_low,
      ymax = beta_age_adjusted_high
    ),
    width = 0.15
  ) +
  ggplot2::geom_vline(
    xintercept = PRIMARY_MISSING_CAL_FLOOR,
    linetype = "dashed"
  ) +
  ggplot2::scale_x_continuous(breaks = MISSING_CAL_FLOORS) +
  ggplot2::labs(
   # title = "Sensitivity of EPE–HbA1c association to assigned missing-tooth CAL",
    x = "Assigned CAL severity for missing teeth (mm)",
    y = "Age-adjusted HbA1c coefficient for EPE"
  ) +
  study_theme(base_size = 12)

save_report_plot(supp_figure_s5, "supp_figure_s5.png", width = 8, height = 5)


# ============================================================
# Supplementary Figure S6. Late-life curvature by CAL floor
# ============================================================

curvature_plot_data <- curvature_raw |>
  dplyr::filter(quantity == "late_life_curvature") |>
  dplyr::mutate(
    status_label = factor(
      GLYCAEMIC_STATUS_LABELS[as.character(DiabetesStatus)],
      levels = unname(GLYCAEMIC_STATUS_LABELS)
    )
  )

supp_figure_s6 <- ggplot2::ggplot(
  curvature_plot_data,
  ggplot2::aes(
    x = cal_floor,
    y = estimate,
    group = status_label,
    linetype = status_label,
    shape = status_label
  )
) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
  ggplot2::geom_vline(
    xintercept = PRIMARY_MISSING_CAL_FLOOR,
    linetype = "dashed"
  ) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = conf_low, ymax = conf_high),
    width = 0.12,
    linewidth = 0.4
  ) +
  ggplot2::facet_wrap(
    ~status_label,
    ncol = 1,
    scales = "fixed"
  ) +
  ggplot2::scale_x_continuous(breaks = MISSING_CAL_FLOORS) +
  ggplot2::labs(
  #  title = "Sensitivity of EPE late-life curvature to assigned missing-tooth CAL",
    x = "Assigned CAL severity for missing teeth (mm)",
    y = "Late-life curvature",
    linetype = "Glycaemic status",
    shape = "Glycaemic status"
  ) +
  study_theme(base_size = 12)+
  ggplot2::theme(
    legend.position = "none",
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank()
  )

save_report_plot(supp_figure_s6, "supp_figure_s6.png", width = 8, height = 5)
