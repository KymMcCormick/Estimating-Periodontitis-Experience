# ============================================================
# 05_produce_supp_table_s3.R
# Produce Supplementary Table S3: components contributing to EPE.
# ============================================================
# Reporting output:
#   outputs/reporting/supp_table_s3.csv

source("00_setup.R")

if (!file.exists(ANALYTIC_DATA_PATH)) {
  stop("Run 01_prepare_nhanes_core.R and 02_build_epe_analytic_data.R first.")
}

component_data <- readRDS(ANALYTIC_DATA_PATH) |>
  add_reporting_variables() |>
  dplyr::mutate(
    expected_contribution_per_missing_tooth = pi_cumulative * assigned_missing_cal,
    observed_retained_component = observed_retained_burden / N_ELIGIBLE_TEETH,
    expected_missing_component = expected_missing_burden / N_ELIGIBLE_TEETH
  ) |>
  dplyr::filter(
    !is.na(AgeGroup),
    is.finite(WTMEC6YR),
    WTMEC6YR > 0,
    !is.na(SDMVPSU),
    !is.na(SDMVSTRA),
    is.finite(n_present_teeth),
    is.finite(n_missing_teeth),
    is.finite(pi_cumulative),
    is.finite(assigned_missing_cal),
    is.finite(expected_contribution_per_missing_tooth),
    is.finite(observed_retained_component),
    is.finite(expected_missing_component)
  )

component_check <- component_data |>
  dplyr::summarise(
    maximum_absolute_difference = max(
      abs(epe - observed_retained_component - expected_missing_component),
      na.rm = TRUE
    )
  )

print(component_check)
if (component_check$maximum_absolute_difference > 1e-8) {
  stop("The EPE components do not reproduce participant-level EPE.")
}

component_design <- make_survey_design(component_data)

component_measures <- tibble::tribble(
  ~variable, ~column, ~digits,
  "n_present_teeth", "Retained teeth", 1,
  "n_missing_teeth", "Missing tooth positions", 1,
  "pi_cumulative", "Cumulative missing-tooth attribution weight (%)", 1,
  "assigned_missing_cal", "Assigned advanced CAL for missing teeth (mm)", 2,
  "expected_contribution_per_missing_tooth", "Expected contribution per missing tooth position (mm)", 2,
  "observed_retained_component", "Observed retained-tooth component", 2,
  "expected_missing_component", "Expected missing-tooth component", 2
)

make_component_row <- function(group = NULL) {
  if (is.null(group)) {
    data <- component_data
    design <- component_design
    label <- "Overall"
  } else {
    data <- component_data |> dplyr::filter(AgeGroup == group)
    design <- make_survey_design(data)
    label <- group
  }

  values <- purrr::pmap_dfc(component_measures, function(variable, column, digits) {
    estimate <- survey_mean_ci(design, variable)

    if (variable == "pi_cumulative") {
      value <- format_ci(
        100 * estimate$estimate,
        100 * estimate$conf_low,
        100 * estimate$conf_high,
        digits
      )
    } else {
      value <- format_ci(
        estimate$estimate,
        estimate$conf_low,
        estimate$conf_high,
        digits
      )
    }

    tibble::tibble(!!column := value)
  })

  dplyr::bind_cols(
    tibble::tibble(
      `Age group` = label,
      n = nrow(data)
    ),
    values
  )
}

supp_table_s3 <- dplyr::bind_rows(
  make_component_row(),
  purrr::map_dfr(AGE_GROUP_LEVELS, make_component_row)
)

save_report_csv(supp_table_s3, "supp_table_s3.csv")
