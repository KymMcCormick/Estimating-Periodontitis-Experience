# scripts/01_prepare_nhanes_core.R
# Prepare participant-level NHANES core data for periodontal analysis projects.
# Includes demographics, survey design variables, examination weights, cycle, and HbA1c.

source("scripts/00_setup.R")

# ------------------------------------------------------------
# 1. Helpers
# ------------------------------------------------------------

peek_header <- function(path, n = 60) {
  rawToChar(readBin(path, what = "raw", n = n))
}

stop_if_html <- function(path) {
  hdr <- peek_header(path, n = 60)
  if (grepl("^\\s*<!DOCTYPE html>|^\\s*<html", hdr, ignore.case = TRUE)) {
    stop("File is HTML, not XPT: ", path)
  }
  invisible(TRUE)
}

read_nhanes_xpt <- function(filename) {
  path <- file.path("data/raw", filename)
  if (!file.exists(path)) stop("Missing file: ", path)
  stop_if_html(path)
  haven::read_xpt(path)
}

cycle_from_filename <- function(filename) {
  out <- sub("^.*_([A-Z])\\.XPT$", "\\1", toupper(basename(filename)))
  
  if (identical(out, toupper(basename(filename)))) {
    stop("Could not extract NHANES cycle from filename: ", filename)
  }
  
  out
}

na_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}


# ------------------------------------------------------------
# 2. NHANES files
# ------------------------------------------------------------

demo_files <- c("DEMO_F.xpt", "DEMO_G.xpt", "DEMO_H.xpt")
ghb_files  <- c("GHB_F.xpt",  "GHB_G.xpt",  "GHB_H.xpt")

cols_demo <- c("SEQN", "RIDAGEYR", "WTMEC2YR", "SDMVPSU", "SDMVSTRA")
cols_ghb  <- c("SEQN", "LBXGH")


# ------------------------------------------------------------
# 3. Demographics and survey variables
# ------------------------------------------------------------

demo <- purrr::map_dfr(demo_files, function(f) {
  read_nhanes_xpt(f) %>%
    dplyr::select(dplyr::any_of(cols_demo)) %>%
    dplyr::mutate(source_cycle = cycle_from_filename(f))
}) %>%
  dplyr::arrange(SEQN, source_cycle) %>%
  dplyr::filter(!is.na(RIDAGEYR), RIDAGEYR >= 30) %>%
  dplyr::mutate(
    age = as.numeric(RIDAGEYR),
    wtmec2yr = as.numeric(WTMEC2YR),
    wtmec6yr = wtmec2yr / 3,
    sdmvpsu = SDMVPSU,
    sdmvstra = SDMVSTRA
  )


# ------------------------------------------------------------
# 4. HbA1c
# ------------------------------------------------------------

ghb <- purrr::map_dfr(ghb_files, function(f) {
  read_nhanes_xpt(f) %>%
    dplyr::select(dplyr::any_of(cols_ghb)) %>%
    dplyr::mutate(source_cycle = cycle_from_filename(f))
}) %>%
  dplyr::mutate(hba1c = as.numeric(LBXGH)) %>%
  dplyr::group_by(SEQN, source_cycle) %>%
  dplyr::summarise(
    hba1c = na_mean(hba1c),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 5. Merge and save
# ------------------------------------------------------------

nhanes_core <- demo %>%
  dplyr::left_join(ghb, by = c("SEQN", "source_cycle")) %>%
  dplyr::select(
    SEQN,
    source_cycle,
    age,
    wtmec2yr,
    wtmec6yr,
    sdmvpsu,
    sdmvstra,
    hba1c,
    dplyr::everything()
  )

saveRDS(nhanes_core, "data/processed/nhanes_participant_core.rds")

message("Saved: data/processed/nhanes_participant_core.rds")