# scripts/01_prepare_nhanes_core.R
# Clean and merge NHANES demographics (age) and HbA1c into a participant-level dataset.

source("scripts/00_setup.R")
source("R/utils_io.R")

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# ---- helpers ----
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
  # DEMO_F.xpt -> F
  sub("^.*_([A-Z])\\.xpt$", "\\1", toupper(filename))
}

# ---- file lists ----
demo_files <- c("DEMO_F.xpt","DEMO_G.xpt","DEMO_H.xpt")
ghb_files  <- c("GHB_F.xpt","GHB_G.xpt","GHB_H.xpt")

cols_demo <- c("SEQN","RIDAGEYR", "WTMEC2YR", "SDMVPSU", "SDMVSTRA")
cols_ghb  <- c("SEQN","LBXGH")

# ---- Demographics (age) ----
demo <- purrr::map_dfr(demo_files, function(f) {
  df <- read_nhanes_xpt(f)
  df %>%
    dplyr::select(dplyr::any_of(cols_demo)) %>%
    dplyr::mutate(source_cycle = cycle_from_filename(f))
}) %>%
  # optional: if any duplicates occur within a cycle, keep first non-missing age
  dplyr::arrange(SEQN, source_cycle) %>%
  dplyr::filter(!is.na(RIDAGEYR) & RIDAGEYR >= 30)  # manuscript: adults aged ≥30

# ---- HbA1c (glycohemoglobin) ----
ghb <- purrr::map_dfr(ghb_files, function(f) {
  df <- read_nhanes_xpt(f)
  df %>%
    dplyr::select(dplyr::any_of(cols_ghb)) %>%
    dplyr::mutate(source_cycle = cycle_from_filename(f))
}) %>%
  dplyr::mutate(LBXGH = as.numeric(LBXGH)) %>%
  # If there is ever >1 record per SEQN within a cycle, average within-cycle first.
  dplyr::group_by(SEQN, source_cycle) %>%
  dplyr::summarise(LBXGH = mean(LBXGH, na.rm = TRUE), .groups = "drop")

# ---- merge ----
demo_ghb <- demo %>%
  dplyr::left_join(ghb, by = c("SEQN", "source_cycle"))

saveRDS(demo_ghb, "data/processed/nhanes_demographics_hba1c.rds")
message("Saved: data/processed/nhanes_demographics_hba1c.rds")


