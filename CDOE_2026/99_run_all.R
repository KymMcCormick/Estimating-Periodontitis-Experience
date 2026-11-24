# ============================================================
# 99_run_all.R
# Rebuild the EPE analysis and produce only the numbered reporting outputs.
# ============================================================

source("00_setup.R")

# Remove stale outputs from superseded analyses, including all former AUC
# products. Internal derived data are rebuilt but not deleted automatically.
unlink(DIR_REPORT, recursive = TRUE, force = TRUE)
dir.create(DIR_REPORT, recursive = TRUE, showWarnings = FALSE)

scripts <- c(
  "01_prepare_nhanes_core.R",
  "02_build_epe_analytic_data_interproximal_corrected.R",
  "03_produce_main_outputs.R",
  "04_produce_supp_calibration_outputs.R",
  "05_produce_supp_table_s3.R",
  "06_produce_supp_tables_s4_s5.R",
  "07_produce_supp_sensitivity_outputs.R",
)

for (script in scripts) {
  message("\n============================================================")
  message("Running ", script)
  message("============================================================")
  source(script, local = new.env(parent = globalenv()))
}

expected_outputs <- c(
  "table_1.csv",
  "table_2.csv",
  "table_3.csv",
  "figure_1.png",
  "figure_2.png",
  "figure_2_data.csv",
  "supp_table_s1.csv",
  "supp_table_s2_panel_a.csv",
  "supp_table_s2_panel_b.csv",
  "supp_table_s2_panel_c.csv",
  "supp_figure_s1.png",
  "supp_figure_s2.png",
  "supp_figure_s3.png",
  "supp_table_s3.csv",
  "supp_table_s4.csv",
  "supp_table_s5.csv",
  "supp_table_s6.csv",
  "supp_table_s7.csv",
  "supp_table_s8.csv",
  "supp_figure_s5.png",
  "supp_figure_s6.png"
)

actual_outputs <- sort(list.files(DIR_REPORT))
missing_outputs <- setdiff(expected_outputs, actual_outputs)
unexpected_outputs <- setdiff(actual_outputs, expected_outputs)

if (length(missing_outputs) > 0) {
  stop(
    "The pipeline completed without producing: ",
    paste(missing_outputs, collapse = ", ")
  )
}

if (length(unexpected_outputs) > 0) {
  stop(
    "Unexpected files were produced in the reporting directory: ",
    paste(unexpected_outputs, collapse = ", ")
  )
}

message(
  "\nPipeline complete. Exactly ",
  length(expected_outputs),
  " reporting files were produced in:"
)
message("  ", DIR_REPORT)