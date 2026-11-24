# ============================================================
# 02_build_epe_analytic_data.R
# Construct the primary EPE analytic dataset and calibration objects.
# ============================================================
# Internal outputs only:
#   data/derived/analytic_dataset_epe.rds
#   data/derived/calibration_primary.rds
#
# This script does not create manuscript tables or figures.

source("00_setup.R")

if (!file.exists(PARTICIPANT_CORE_PATH)) {
  stop("Run 01_prepare_nhanes_core.R first: ", PARTICIPANT_CORE_PATH)
}

participant_core <- readRDS(PARTICIPANT_CORE_PATH)
assert_has_columns(
  participant_core,
  c("SEQN", "source_cycle", "Age", "WTMEC2YR", "WTMEC6YR", "SDMVPSU", "SDMVSTRA", "LBXGH")
)


ohx_files <- paste0("OHXPER_", NHANES_CYCLES, ".xpt")

# EPE retained-tooth values: maximum CAL across four interproximal sites
# per tooth (mesiofacial, distofacial, mesiolingual, and distolingual).
# Mid-facial and mid-lingual sites are not used for this tooth-level maximum.
tooth_cal <- purrr::map_dfr(ohx_files, function(filename) {
  data <- read_nhanes_xpt(filename)
  make_epe_max_cal_by_tooth(data) |>
    dplyr::mutate(source_cycle = cycle_from_filename(filename), .after = SEQN)
})

# Conventional comparison measures and periodontal observability: all six
# recorded NHANES sites per eligible tooth.
site_measures <- purrr::map_dfr(ohx_files, function(filename) {
  data <- read_nhanes_xpt(filename)
  calc_observed_site_measures(data) |>
    dplyr::mutate(source_cycle = cycle_from_filename(filename), .after = SEQN)
})

periodontal <- tooth_cal |>
  dplyr::inner_join(
    site_measures,
    by = c("SEQN", "source_cycle")
  )

base_data <- participant_core |>
  dplyr::inner_join(
    periodontal,
    by = c("SEQN", "source_cycle")
  ) |>
  dplyr::arrange(source_cycle, SEQN)

cal_columns <- names(base_data)[grepl("^CAL_\\d{2}$", names(base_data))]
if (length(cal_columns) != N_ELIGIBLE_TEETH) {
  stop(
    "Expected ", N_ELIGIBLE_TEETH, " tooth-level CAL columns but found ",
    length(cal_columns), "."
  )
}

# The study is restricted to dentate adults with at least one observed
# tooth-level CAL contribution.
cal_matrix_preliminary <- as.matrix(base_data[, cal_columns, drop = FALSE])
n_present_preliminary <- rowSums(!is.na(cal_matrix_preliminary))
keep_dentate <- n_present_preliminary > 0

base_data <- base_data[keep_dentate, , drop = FALSE] |>
  dplyr::mutate(
    n_present_teeth_preliminary = n_present_preliminary[keep_dentate],
    n_missing_teeth_preliminary = N_ELIGIBLE_TEETH - n_present_teeth_preliminary
  )

# Primary incident model q(a), using the single canonical schedule table.
primary_incident <- fit_incident_attribution()

# NHANES missing-tooth accumulation, estimated in the broader periodontal
# sample irrespective of HbA1c availability.
missing_accumulation <- make_missing_accumulation(
  base_data |>
    dplyr::transmute(
      Age,
      n_missing_teeth = n_missing_teeth_preliminary,
      WTMEC6YR,
      SDMVPSU,
      SDMVSTRA
    )
)

primary_cumulative <- make_cumulative_attribution(
  incident_curve = primary_incident$curve,
  missing_curve = missing_accumulation$curve
)

# Participant-level EPE.
epe_components <- calculate_epe(
  tooth_cal_data = base_data |> dplyr::select(dplyr::all_of(cal_columns)),
  age = base_data$Age,
  cumulative_curve = primary_cumulative,
  floor = PRIMARY_MISSING_CAL_FLOOR
)

analytic <- dplyr::bind_cols(
  base_data |>
    dplyr::select(
      SEQN,
      source_cycle,
      Age,
      LBXGH,
      WTMEC2YR,
      WTMEC6YR,
      SDMVPSU,
      SDMVSTRA,
      mean_CAL,
      extent_ge_4mm,
      n_sites_observed
    ),
  epe_components
) |>
  add_reporting_variables() |>
  dplyr::arrange(source_cycle, SEQN)

assert_has_columns(
  analytic,
  c(
    "SEQN", "Age", "LBXGH", "WTMEC2YR", "WTMEC6YR", "SDMVPSU", "SDMVSTRA",
    "mean_CAL", "extent_ge_4mm", "n_sites_observed", "pct_sites_observed",
    "n_present_teeth", "n_missing_teeth", "observed_retained_burden",
    "max_observed_cal", "assigned_missing_cal", "pi_cumulative",
    "expected_missing_burden", "epe", "AgeGroup", "DiabetesStatus"
  )
)

if (any(analytic$n_present_teeth + analytic$n_missing_teeth != N_ELIGIBLE_TEETH)) {
  stop("The retained and missing tooth counts do not sum to 28 for every participant.")
}

saveRDS(analytic, ANALYTIC_DATA_PATH)

saveRDS(
  list(
    site_definitions = list(
      epe_retained_tooth = paste0(
        "Maximum CAL across four interproximal sites per retained tooth: ",
        paste(EPE_INTERPROXIMAL_CAL_SUFFIXES, collapse = ", "),
        ". Mid-facial and mid-lingual sites excluded."
      ),
      conventional_measures = paste0(
        "Mean CAL, extent CAL >=4 mm, and sites observed use all six ",
        "NHANES periodontal sites per eligible tooth."
      )
    ),
    attribution_schedules = ATTRIBUTION_SCHEDULES,
    incident = primary_incident,
    missing_accumulation = missing_accumulation,
    cumulative = primary_cumulative,
    primary_missing_cal_floor = PRIMARY_MISSING_CAL_FLOOR
  ),
  CALIBRATION_DATA_PATH
)

message("Saved: ", ANALYTIC_DATA_PATH)
message("Saved: ", CALIBRATION_DATA_PATH)
message(
  "EPE retained-tooth definition: maximum CAL across four interproximal ",
  "sites per tooth (112 possible sites; mid-facial and mid-lingual excluded)."
)
message(
  "Conventional measures: all six NHANES sites per tooth ",
  "(168 possible sites)."
)
message("Broader periodontal analytic sample: ", nrow(analytic))
message("HbA1c complete: ", sum(is.finite(analytic$LBXGH)))

if (nrow(analytic) != 10714) {
  warning(
    "The broader analytic sample is ", nrow(analytic),
    ", not the 10,714 reported in the current manuscript."
  )
}

if (sum(is.finite(analytic$LBXGH)) != 10324) {
  warning(
    "The HbA1c-complete sample is ", sum(is.finite(analytic$LBXGH)),
    ", not the 10,324 reported in the current manuscript."
  )
}
