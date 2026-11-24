# scripts/00_setup.R
# Package loading and directory setup for reproducible NHANES periodontal analysis projects.

# ============================================================
# 1. Required packages
# ============================================================

required_pkgs <- c(
  "dplyr", "tidyr", "purrr", "readr", "stringr", "tibble",
  "haven", "httr", "survey",
  "ggplot2", "patchwork", "scales",
  "broom", "pROC"
)

# ============================================================
# 2. Install and load packages
# ============================================================

to_install <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(to_install) > 0) {
  install.packages(to_install)
}

invisible(lapply(required_pkgs, library, character.only = TRUE))

# ============================================================
# 3. Reproducibility and survey settings
# ============================================================

set.seed(123)

# Standard handling for strata with a single PSU in survey analyses.
options(survey.lonely.psu = "adjust")

# ============================================================
# 4. Standard project directories
# ============================================================

project_dirs <- c(
  "data/raw",
  "data/processed",
  "data/derived",
  "outputs/figures",
  "outputs/tables",
  "outputs/models",
  "outputs/logs"
)

invisible(lapply(project_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# ============================================================
# 5. Lightweight project checks
# ============================================================

if (!dir.exists("data/raw")) {
  warning("data/raw directory was not created successfully.")
}

message("Setup complete.")
message("Working directory: ", getwd())
message("Packages loaded: ", paste(required_pkgs, collapse = ", "))
