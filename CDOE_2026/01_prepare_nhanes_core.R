# ============================================================
# 01_prepare_nhanes_core.R
# Prepare participant-level NHANES demographics, survey variables,
# and HbA1c for the 2009–2014 analysis cycles.
# ============================================================
# Internal output only:
#   data/derived/nhanes_participant_core.rds

source("00_setup.R")


demo_files <- paste0("DEMO_", NHANES_CYCLES, ".xpt")
ghb_files <- paste0("GHB_", NHANES_CYCLES, ".xpt")


demographics <- purrr::map_dfr(demo_files, function(filename) {
  read_nhanes_xpt(filename) |>
    dplyr::select(
      SEQN,
      RIDAGEYR,
      WTMEC2YR,
      SDMVPSU,
      SDMVSTRA
    ) |>
    dplyr::mutate(source_cycle = cycle_from_filename(filename))
}) |>
  dplyr::transmute(
    SEQN,
    source_cycle,
    Age = as.numeric(RIDAGEYR),
    WTMEC2YR = as.numeric(WTMEC2YR),
    WTMEC6YR = WTMEC2YR / 3,
    SDMVPSU,
    SDMVSTRA
  ) |>
  dplyr::filter(is.finite(Age), Age >= 30)


hba1c <- purrr::map_dfr(ghb_files, function(filename) {
  read_nhanes_xpt(filename) |>
    dplyr::select(SEQN, LBXGH) |>
    dplyr::mutate(source_cycle = cycle_from_filename(filename))
}) |>
  dplyr::mutate(LBXGH = as.numeric(LBXGH)) |>
  dplyr::group_by(SEQN, source_cycle) |>
  dplyr::summarise(
    LBXGH = if (all(is.na(LBXGH))) NA_real_ else mean(LBXGH, na.rm = TRUE),
    .groups = "drop"
  )


participant_core <- demographics |>
  dplyr::left_join(hba1c, by = c("SEQN", "source_cycle")) |>
  dplyr::arrange(source_cycle, SEQN)


if (anyDuplicated(participant_core[c("SEQN", "source_cycle")]) > 0) {
  stop("Participant core contains duplicated SEQN/cycle rows.")
}

saveRDS(participant_core, PARTICIPANT_CORE_PATH)

message("Saved: ", PARTICIPANT_CORE_PATH)
message("Participants aged ≥30 years in participant core: ", nrow(participant_core))
message("Participants with HbA1c: ", sum(is.finite(participant_core$LBXGH)))

