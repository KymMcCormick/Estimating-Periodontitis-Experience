# scripts/00_setup.R
# ============================================================
# Project setup
# ============================================================

options(survey.lonely.psu = "adjust")

required_packages <- c(
  "dplyr",
  "ggplot2",
  "readr",
  "rlang",
  "survey",
  "tidyr"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop(
    "Install missing package(s) before running the pipeline: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

required_dirs <- c(
  "data/raw",
  "data/derived",
  "outputs/tables",
  "outputs/figures"
)

invisible(lapply(required_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

source("R/mce_functions.R")
source("R/nhanes_helpers.R")

message("Setup complete.")
