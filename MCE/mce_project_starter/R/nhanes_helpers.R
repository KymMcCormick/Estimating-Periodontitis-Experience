# R/nhanes_helpers.R
# ============================================================
# Helper functions for NHANES survey analyses
# ============================================================

make_nhanes_design <- function(data,
                               weight_col = "wtmec6yr",
                               strata_col = "sdmvstra",
                               psu_col = "sdmvpsu") {
  check_required_columns(data, c(weight_col, strata_col, psu_col))

  survey::svydesign(
    ids = stats::as.formula(paste0("~", psu_col)),
    strata = stats::as.formula(paste0("~", strata_col)),
    weights = stats::as.formula(paste0("~", weight_col)),
    nest = TRUE,
    data = data
  )
}

save_csv <- function(data, path) {
  readr::write_csv(data, path)
  message("Saved: ", path)
  invisible(path)
}

save_plot <- function(plot, path, width = 7, height = 5, dpi = 300) {
  ggplot2::ggsave(path, plot = plot, width = width, height = height, dpi = dpi)
  message("Saved: ", path)
  invisible(path)
}
