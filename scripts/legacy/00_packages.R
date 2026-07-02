# code/00_packages.R
# Package loading for the reproducibility pipeline.

required_pkgs <- c(
  "dplyr","tidyr","purrr","readr","stringr",
  "haven","httr",
  "ggplot2","patchwork","scales",
  "broom","pROC"
)

to_install <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) {
  install.packages(to_install)
}

invisible(lapply(required_pkgs, library, character.only = TRUE))

# For deterministic reproduction
set.seed(123)
