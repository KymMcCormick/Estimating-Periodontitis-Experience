# ============================================================
# 06_produce_supp_tables_s4_s5.R
# Produce Supplementary Table S4 (late-life curvature) and
# Supplementary Table S5 (HbA1c rank correlations and joint model).
# ============================================================
# Reporting outputs:
#   outputs/reporting/supp_table_s4.csv
#   outputs/reporting/supp_table_s5.csv

source("00_setup.R")

if (!file.exists(ANALYTIC_DATA_PATH)) {
  stop("Run 01_prepare_nhanes_core.R and 02_build_epe_analytic_data.R first.")
}

analytic <- readRDS(ANALYTIC_DATA_PATH) |> add_reporting_variables()


# ============================================================
# Supplementary Table S4. Late-life changes and curvature
# ============================================================

s4_data <- analytic |>
  dplyr::filter(
    !is.na(DiabetesStatus),
    is.finite(Age),
    is.finite(WTMEC6YR),
    WTMEC6YR > 0,
    !is.na(SDMVPSU),
    !is.na(SDMVSTRA)
  )

s4_measures <- tibble::tribble(
  ~outcome, ~measure,
  "epe", "EPE",
  "extent_ge_4mm", "% CAL ≥4 mm",
  "mean_CAL", "Mean CAL"
)

s4_raw <- tidyr::expand_grid(
  measure_row = seq_len(nrow(s4_measures)),
  DiabetesStatus = GLYCAEMIC_STATUS_LEVELS
) |>
  dplyr::left_join(
    s4_measures |> dplyr::mutate(measure_row = dplyr::row_number()),
    by = "measure_row"
  ) |>
  dplyr::select(-measure_row) |>
  purrr::pmap_dfr(function(DiabetesStatus, outcome, measure) {
    fit_late_life_contrasts(
      data = s4_data,
      outcome = outcome,
      status = DiabetesStatus,
      z_critical = 1.96
    ) |>
      dplyr::mutate(measure = measure)
  })

s4_wide <- s4_raw |>
  dplyr::select(
    DiabetesStatus,
    measure,
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
    DiabetesStatus = factor(DiabetesStatus, levels = GLYCAEMIC_STATUS_LEVELS),
    measure = factor(measure, levels = s4_measures$measure)
  ) |>
  dplyr::arrange(DiabetesStatus, measure)

supp_table_s4 <- s4_wide |>
  dplyr::transmute(
    `Glycaemic status` = GLYCAEMIC_STATUS_LABELS[as.character(DiabetesStatus)],
    Measure = as.character(measure),
    `Fitted change, 60–70 years` = format_number(estimate_change_60_70, 2),
    `Fitted change, 70–80 years` = format_number(estimate_change_70_80, 2),
    `Late-life curvature (95% CI)` = format_ci(
      estimate_late_life_curvature,
      conf_low_late_life_curvature,
      conf_high_late_life_curvature,
      2
    )
  )

save_report_csv(supp_table_s4, "supp_table_s4.csv")


# ============================================================
# Supplementary Table S5. Rank correlations and joint model
# ============================================================

hba1c_data <- prepare_hba1c_complete_case(analytic)
hba1c_design <- make_survey_design(hba1c_data)

rank_rows <- tibble::tribble(
  ~measure, ~variable,
  "Expected Periodontitis Experience", "epe",
  "Mean CAL", "mean_CAL",
  "Extent CAL ≥4 mm", "extent_ge_4mm",
  "Periodontal sites observed", "pct_sites_observed"
) |>
  dplyr::mutate(
    rank_correlation = purrr::map_dbl(
      variable,
      ~survey_rank_correlation(hba1c_data, .x, "LBXGH")
    )
  )

joint_fit <- survey::svyglm(
  LBXGH ~ z_epe + z_mean_CAL + z_extent_ge_4mm + age_c + I(age_c^2),
  design = hba1c_design
)

joint_terms <- dplyr::bind_rows(
  tidy_svy_term(joint_fit, "z_epe") |>
    dplyr::mutate(measure = "Expected Periodontitis Experience"),
  tidy_svy_term(joint_fit, "z_mean_CAL") |>
    dplyr::mutate(measure = "Mean CAL"),
  tidy_svy_term(joint_fit, "z_extent_ge_4mm") |>
    dplyr::mutate(measure = "Extent CAL ≥4 mm")
) |>
  dplyr::select(measure, estimate, conf_low, conf_high)

supp_table_s5 <- rank_rows |>
  dplyr::left_join(joint_terms, by = "measure") |>
  dplyr::transmute(
    Measure = measure,
    `Rank correlation with HbA1c` = format_number(rank_correlation, 3),
    `Joint age-adjusted β (95% CI)` = dplyr::if_else(
      is.na(estimate),
      "—",
      format_ci(estimate, conf_low, conf_high, 3)
    )
  )

save_report_csv(supp_table_s5, "supp_table_s5.csv")
