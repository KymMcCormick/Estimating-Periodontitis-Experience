# scripts/00_setup.R
# Package loading and directory setup for the EPE reproducibility pipeline.

required_pkgs <- c(
  "dplyr", "tidyr", "purrr", "readr", "stringr", "tibble",
  "haven", "httr", "survey",
  "ggplot2", "patchwork", "scales",
  "broom", "pROC"
)

to_install <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) {
  install.packages(to_install)
}

invisible(lapply(required_pkgs, library, character.only = TRUE))

# For deterministic reproduction
set.seed(123)

# Standard project directories
project_dirs <- c(
  "data/raw",
  "data/processed",
  "outputs/figures",
  "outputs/tables"
)

invisible(lapply(project_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
