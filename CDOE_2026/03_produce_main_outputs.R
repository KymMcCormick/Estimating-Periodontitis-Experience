# ============================================================
# 03_produce_main_outputs.R
# Produce every numbered table and figure in the main manuscript.
# ============================================================
# Reporting outputs:
#   outputs/reporting/table_1.csv
#   outputs/reporting/table_2.csv
#   outputs/reporting/table_3.csv
#   outputs/reporting/figure_1.png
#   outputs/reporting/figure_2.png

#source("00_setup.R")

if (!file.exists(ANALYTIC_DATA_PATH) || !file.exists(CALIBRATION_DATA_PATH)) {
  stop("Run 01_prepare_nhanes_core.R and 02_build_epe_analytic_data.R first.")
}

analytic <- readRDS(ANALYTIC_DATA_PATH)
calibration <- readRDS(CALIBRATION_DATA_PATH)


# ============================================================
# Table 1. Empirical extraction-attribution schedules
# ============================================================

table_1 <- ATTRIBUTION_SCHEDULES |>
  dplyr::transmute(
    `Attribution schedule` = as.character(schedule),
    Author = author,
    Country = country,
    `Age range (years)` = age_band,
    n = n,
    `Conditional probability` = formatC(p, format = "f", digits = 3)
  )

save_report_csv(table_1, "table_1.csv")


# ============================================================
# Table 2. Survey-weighted periodontal measures by age group
# ============================================================

table_2_data <- analytic |>
  add_reporting_variables() |>
  dplyr::filter(
    !is.na(AgeGroup),
    is.finite(WTMEC6YR),
    WTMEC6YR > 0,
    !is.na(SDMVPSU),
    !is.na(SDMVSTRA),
    is.finite(pct_sites_observed),
    is.finite(mean_CAL),
    is.finite(extent_ge_4mm),
    is.finite(epe)
  )

table_2_design <- make_survey_design(table_2_data)

table_2_measures <- tibble::tribble(
  ~variable, ~column, ~digits,
  "pct_sites_observed", "Sites observed (%) (95% CI)", 1,
  "mean_CAL", "Mean CAL (mm) (95% CI)", 2,
  "extent_ge_4mm", "% CAL ≥ 4 mm (95% CI)", 1,
  "epe", "EPE (95% CI)", 2
)

make_table_2_row <- function(group = NULL) {
  if (is.null(group)) {
    design <- table_2_design
    data <- table_2_data
    label <- "Overall"
  } else {
    data <- table_2_data |> dplyr::filter(as.character(AgeGroup) == group)
    design <- make_survey_design(data)
    label <- group
  }

  values <- purrr::pmap_dfc(table_2_measures, function(variable, column, digits) {
    estimate <- survey_mean_ci(design, variable)
    tibble::tibble(!!column := format_ci(
      estimate$estimate,
      estimate$conf_low,
      estimate$conf_high,
      digits = digits
    ))
  })

  dplyr::bind_cols(
    tibble::tibble(
      `Age group (years)` = label,
      n = format(nrow(data), big.mark = ",", scientific = FALSE, trim = TRUE)
    ),
    values
  )
}

table_2 <- dplyr::bind_rows(
  make_table_2_row(),
  purrr::map_dfr(AGE_GROUP_LEVELS, make_table_2_row)
)

save_report_csv(table_2, "table_2.csv")


# ============================================================
# Table 3. Separate HbA1c regression models
# ============================================================

hba1c_data <- prepare_hba1c_complete_case(analytic)
if (nrow(hba1c_data) == 0) stop("No complete observations for the HbA1c analysis.")
hba1c_design <- make_survey_design(hba1c_data)

regression_specs <- tibble::tribble(
  ~measure, ~z_variable,
  "Expected Periodontitis Experience", "z_epe",
  "Mean CAL", "z_mean_CAL",
  "Extent CAL ≥4 mm", "z_extent_ge_4mm"
)

fit_table_3_measure <- function(measure, z_variable) {
  unadjusted <- survey::svyglm(
    stats::as.formula(paste0("LBXGH ~ ", z_variable)),
    design = hba1c_design
  )

  age_adjusted <- survey::svyglm(
    stats::as.formula(paste0(
      "LBXGH ~ ", z_variable, " + age_c + I(age_c^2)"
    )),
    design = hba1c_design
  )

  unadjusted_term <- tidy_svy_term(unadjusted, z_variable)
  adjusted_term <- tidy_svy_term(age_adjusted, z_variable)

  attenuation <- 100 * (
    unadjusted_term$estimate - adjusted_term$estimate
  ) / unadjusted_term$estimate

  tibble::tibble(
    `Periodontal measure` = measure,
    `Unadjusted β (95% CI)` = format_ci(
      unadjusted_term$estimate,
      unadjusted_term$conf_low,
      unadjusted_term$conf_high,
      digits = 3
    ),
    `Age-adjusted β (95% CI)` = format_ci(
      adjusted_term$estimate,
      adjusted_term$conf_low,
      adjusted_term$conf_high,
      digits = 3
    ),
    `Attenuation after age adjustment` = paste0(format_number(attenuation, 1), "%")
  )
}

table_3 <- purrr::pmap_dfr(regression_specs, fit_table_3_measure)
save_report_csv(table_3, "table_3.csv")


# ============================================================
# Figure 1. Incident and cumulative attribution curves
# ============================================================

figure_1_data <- calibration$cumulative |>
  dplyr::select(age, q_incident, pi_cumulative) |>
  tidyr::pivot_longer(
    cols = c(q_incident, pi_cumulative),
    names_to = "curve",
    values_to = "value"
  ) |>
  dplyr::mutate(
    curve = dplyr::recode(
      curve,
      pi_cumulative = "Cumulative missing-tooth attribution π(A)",
      q_incident = "Incident extraction attribution q(a)"
    ),
    curve = factor(
      curve,
      levels = c(
        "Cumulative missing-tooth attribution π(A)",
        "Incident extraction attribution q(a)"
      )
    )
  )

figure_1 <- ggplot2::ggplot(
  figure_1_data,
  ggplot2::aes(x = age, y = value, linetype = curve)
) +
  ggplot2::geom_line(linewidth = 1.2, colour = "black") +
  ggplot2::scale_linetype_manual(values = c("solid", "dashed")) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0, 0.8, 0.2),
    limits = c(0, 0.85)
  ) +
  ggplot2::scale_x_continuous(breaks = seq(30, 80, 10)) +
  ggplot2::labs(
    x = "Age",
    y = "Probability / attribution weight",
    linetype = "Curve"
  ) +
  study_theme(base_size = 12)

save_report_plot(figure_1, "figure_1.png", width = 9, height = 6)


# ============================================================
# Figure 2. Periodontal measures by age group and glycaemic status
# ============================================================

figure_2_data <- hba1c_data |>
  dplyr::mutate(
    status_label = factor(
      GLYCAEMIC_STATUS_LABELS[as.character(DiabetesStatus)],
      levels = unname(GLYCAEMIC_STATUS_LABELS)
    ),
    age_group_midpoint = dplyr::case_when(
      AgeGroup == "30–39" ~ 35,
      AgeGroup == "40–49" ~ 45,
      AgeGroup == "50–59" ~ 55,
      AgeGroup == "60–69" ~ 65,
      AgeGroup == "70–79" ~ 75,
      AgeGroup == "80+" ~ 80,
      TRUE ~ NA_real_
    )
  )

save_report_csv(figure_2_data, "figure_2_data.csv")

figure_2_measures <- tibble::tribble(
  ~variable, ~title,
  "epe", "A. Expected Periodontitis Experience",
  "mean_CAL", "B. Mean CAL (mm)",
  "extent_ge_4mm", "C. Sites with CAL ≥4 mm (%)",
  "pct_sites_observed", "D. Periodontal sites observed (%)"
)

status_offsets <- c(
  `Normal glycaemia` = -0.55,
  Prediabetes = 0,
  Diabetes = 0.55
)

estimate_age_group_points <- function(variable) {
  tidyr::expand_grid(
    status_label = levels(figure_2_data$status_label),
    AgeGroup = AGE_GROUP_LEVELS
  ) |>
    purrr::pmap_dfr(function(status_label, AgeGroup) {
      group_data <- figure_2_data |>
        dplyr::filter(
          as.character(.data$status_label) == .env$status_label,
          as.character(.data$AgeGroup) == .env$AgeGroup
        )

      if (nrow(group_data) == 0) return(tibble::tibble())

      design <- make_survey_design(group_data)
      estimate_result <- survey_mean_ci(design, variable)
      midpoint <- unique(group_data$age_group_midpoint)

      tibble::tibble(
        variable = variable,
        status_label = factor(status_label, levels = levels(figure_2_data$status_label)),
        AgeGroup = factor(AgeGroup, levels = AGE_GROUP_LEVELS),
        age_group_midpoint = midpoint[1],
        x = midpoint[1] + status_offsets[[status_label]],
        estimate = estimate_result$estimate,
        conf_low = estimate_result$conf_low,
        conf_high = estimate_result$conf_high
      )
    })
}

figure_2_points <- purrr::map_dfr(
  figure_2_measures$variable,
  estimate_age_group_points
)

fit_age_curve <- function(variable, status) {
  group_data <- figure_2_data |>
    dplyr::filter(
      as.character(status_label) == status,
      is.finite(.data[[variable]])
    )
  design <- make_survey_design(group_data)
  fit <- survey::svyglm(
    stats::as.formula(paste0(variable, " ~ Age + I(Age^2)")),
    design = design
  )
  prediction_data <- tibble::tibble(Age = seq(30, 80, by = 0.25))

  tibble::tibble(
    variable = variable,
    status_label = factor(status, levels = levels(figure_2_data$status_label)),
    Age = prediction_data$Age,
    fitted = as.numeric(stats::predict(fit, newdata = prediction_data))
  )
}

figure_2_curves <- tidyr::expand_grid(
  variable = figure_2_measures$variable,
  status = levels(figure_2_data$status_label)
) |>
  purrr::pmap_dfr(fit_age_curve)

shape_values <- c(`Normal glycaemia` = 16, Prediabetes = 17, Diabetes = 15)
linetype_values <- c(`Normal glycaemia` = "solid", Prediabetes = "dashed", Diabetes = "dotdash")
colour_values <- c(`Normal glycaemia` = "grey35", Prediabetes = "grey60", Diabetes = "black")

make_figure_2_panel <- function(variable, title) {
  points <- figure_2_points |> dplyr::filter(.data$variable == .env$variable)
  curves <- figure_2_curves |> dplyr::filter(.data$variable == .env$variable)
  x_label <- if (variable %in% c("extent_ge_4mm", "pct_sites_observed")) {
    "Age group, years"
  } else {
    NULL
  }
  y_label <- if (variable %in% c("epe", "extent_ge_4mm")) {
    "Survey-weighted mean"
  } else {
    NULL
  }

  ggplot2::ggplot() +
    ggplot2::geom_line(
      data = curves,
      ggplot2::aes(
        x = Age,
        y = fitted,
        linetype = status_label,
        colour = status_label
      ),
      linewidth = 0.9
    ) +
    ggplot2::geom_errorbar(
      data = points,
      ggplot2::aes(
        x = x,
        ymin = conf_low,
        ymax = conf_high,
        colour = status_label
      ),
      width = 0.35,
      linewidth = 0.45
    ) +
    ggplot2::geom_point(
      data = points,
      ggplot2::aes(
        x = x,
        y = estimate,
        shape = status_label,
        colour = status_label
      ),
      size = 2.5
    ) +
    ggplot2::scale_shape_manual(values = shape_values) +
    ggplot2::scale_linetype_manual(values = linetype_values) +
    ggplot2::scale_colour_manual(values = colour_values) +
    ggplot2::scale_x_continuous(
      breaks = c(35, 45, 55, 65, 75, 80),
      labels = AGE_GROUP_LEVELS
    ) +
    ggplot2::labs(
      title = title,
      x = x_label,
      y = y_label,
      shape = "Glycaemic status",
      linetype = "Glycaemic status",
      colour = "Glycaemic status"
    ) +
    study_theme(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 11, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
}

figure_2_panels <- purrr::pmap(figure_2_measures, make_figure_2_panel)

figure_2 <- patchwork::wrap_plots(plotlist = figure_2_panels, ncol = 2, guides = "collect") +
  patchwork::plot_annotation(
    theme = ggplot2::theme(
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5)
    )
  ) &
  ggplot2::theme(legend.position = "bottom")

save_report_plot(figure_2, "figure_2.png", width = 9, height = 9)
