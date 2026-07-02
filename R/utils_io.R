# R/utils_io.R
# Shared input/output helpers for NHANES XPT files.

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

read_nhanes_xpt <- function(filename, raw_dir = "data/raw") {
  path <- file.path(raw_dir, filename)
  if (!file.exists(path)) stop("Missing file: ", path)
  stop_if_html(path)
  haven::read_xpt(path)
}

cycle_from_filename <- function(filename) {
  # Example: DEMO_F.xpt -> F
  sub("^.*_([A-Z])\\.xpt$", "\\1", toupper(filename))
}

standardise_cycle_code <- function(x) {
  x <- toupper(x)
  dplyr::case_when(
    grepl("_F\\.", x) | grepl("F$", x) ~ "F",
    grepl("_G\\.", x) | grepl("G$", x) ~ "G",
    grepl("_H\\.", x) | grepl("H$", x) ~ "H",
    TRUE ~ NA_character_
  )
}
