# run_all.R
# ============================================================
# Run complete EPE analysis pipeline
# ============================================================
#
# Purpose:
#   Source each numbered project script in order so the full
#   analysis can be reproduced from the project root.
#
# Usage:
#   1. Open the project root in RStudio.
#   2. Confirm your working directory is the GitHub project root.
#   3. Run:
#        source("run_all.R")
#
# Notes:
#   - This script stops immediately if any step fails.
#   - Individual scripts should use project-relative paths.
#   - Avoid using setwd() inside analysis scripts.
# ============================================================


# ------------------------------------------------------------
# 1. Basic setup
# ------------------------------------------------------------

start_time <- Sys.time()

message("============================================================")
message("Starting full EPE analysis pipeline")
message("Started at: ", format(start_time))
message("Working directory: ", getwd())
message("============================================================")


# ------------------------------------------------------------
# 2. List scripts in intended order
# ------------------------------------------------------------

analysis_scripts <- c(
  "00_setup.R",
  "01_prepare_nhanes_core.R",
  "02_build_epe_iteration01.R",
  "03_diagnose_incident_attribution_model.R",
  "04_summarise_epe_distribution.R",
  "05_assess_hba1c_associational_coherence.R",
  "06_assess_diabetes_discrimination_auc.R",
  "07_summarise_epe_components_by_age.R"
)


# If your scripts are stored inside a scripts/ folder, use this instead:
# analysis_scripts <- file.path("scripts", analysis_scripts)


# ------------------------------------------------------------
# 3. Check that all scripts exist before running
# ------------------------------------------------------------

missing_scripts <- analysis_scripts[!file.exists(analysis_scripts)]

if (length(missing_scripts) > 0) {
  stop(
    "The following script(s) could not be found:\n",
    paste0("  - ", missing_scripts, collapse = "\n"),
    "\n\nCheck whether run_all.R is in the project root and whether your ",
    "scripts are stored in the root folder or inside scripts/."
  )
}


# ------------------------------------------------------------
# 4. Run scripts sequentially
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# 5. Finish
# ------------------------------------------------------------

end_time <- Sys.time()

message("")
message("============================================================")
message("Full EPE analysis pipeline completed successfully")
message("Finished at: ", format(end_time))
message(
  "Total duration: ",
  round(as.numeric(difftime(end_time, start_time, units = "mins")), 2),
  " minutes"
)
message("============================================================")
