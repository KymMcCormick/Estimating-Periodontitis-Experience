# scripts/00_setup.R
# Package loading and directory setup for reproducible NHANES periodontal analysis projects.

# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

required_pkgs <- c(
  "dplyr", "tidyr", "purrr", "readr", "stringr", "tibble",
  "haven", "httr", "survey",
  "ggplot2", "patchwork", "scales",
  "broom", "pROC"
)

to_install <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(to_install) > 0) {
  install.packages(to_install)
}

invisible(lapply(required_pkgs, library, character.only = TRUE))


# ------------------------------------------------------------
# 2. Reproducibility settings
# ------------------------------------------------------------

set.seed(123)

# Standard handling for lonely PSUs in survey analyses
options(survey.lonely.psu = "adjust")


# ------------------------------------------------------------
# 3. Project directories
# ------------------------------------------------------------

project_dirs <- c(
  "data/raw",
  "data/processed",
  "data/derived",
  "outputs/figures",
  "outputs/tables",
  "outputs/models"
)

invisible(lapply(project_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))