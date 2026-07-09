# run_all.R
# ============================================================
# Run complete MCE analysis pipeline
# ============================================================
# Purpose:
#   Source each numbered project script in order so the full
#   analysis can be reproduced from the project root.
#
# Usage:
#   Open this project in RStudio, confirm the working directory
#   is the project root, then run:
#      source("run_all.R")
#
# Notes:
#   - This script stops immediately if any step fails.
#   - Scripts use project-relative paths.
#   - Avoid setwd() inside analysis scripts.
# ============================================================

start_time <- Sys.time()

message("============================================================")
message("Starting full MCE analysis pipeline")
message("Started at: ", format(start_time))
message("Working directory: ", getwd())
message("============================================================")

analysis_scripts <- file.path(
  "scripts",
  c(
    "00_setup.R",
    "01_prepare_nhanes_core.R",
    "02_build_mce_measure.R",
    "03_describe_mce_distribution.R",
    "04_age_trajectories.R",
    "05_tooth_count_compression.R",
    "06_hba1c_alignment.R",
    "07_diabetes_auc.R"
  )
)

missing_scripts <- analysis_scripts[!file.exists(analysis_scripts)]

if (length(missing_scripts) > 0) {
  stop(
    "The following script(s) could not be found:\n",
    paste0("  - ", missing_scripts, collapse = "\n"),
    "\n\nCheck that run_all.R is in the project root."
  )
}

for (script in analysis_scripts) {
  step_start <- Sys.time()

  message("")
  message("------------------------------------------------------------")
  message("Running: ", script)
  message("Started: ", format(step_start))
  message("------------------------------------------------------------")

  tryCatch(
    {
      source(script, local = FALSE)

      step_end <- Sys.time()
      message("Completed: ", script)
      message(
        "Step duration: ",
        round(as.numeric(difftime(step_end, step_start, units = "mins")), 2),
        " minutes"
      )
    },
    error = function(e) {
      message("")
      message("ERROR in: ", script)
      message(conditionMessage(e))
      message("")
      stop("Pipeline stopped because a script failed.", call. = FALSE)
    }
  )
}

end_time <- Sys.time()

message("")
message("============================================================")
message("Full MCE analysis pipeline completed successfully")
message("Finished at: ", format(end_time))
message(
  "Total duration: ",
  round(as.numeric(difftime(end_time, start_time, units = "mins")), 2),
  " minutes"
)
message("============================================================")
