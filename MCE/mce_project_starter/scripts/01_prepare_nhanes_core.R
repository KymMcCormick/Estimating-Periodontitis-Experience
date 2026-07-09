# scripts/01_prepare_nhanes_core.R
# ============================================================
# Prepare NHANES core periodontal site-level data
# ============================================================
# This script is the bridge from your existing NHANES code into
# the MCE project.
#
# Recommended approach:
#   1. Copy your old preparation script into this project.
#   2. Edit it so that its final output is saved as:
#        data/derived/nhanes_core_site_level.rds
#
# Required columns for later scripts:
#   seqn, age, tooth, cal, wtmec6yr, sdmvstra, sdmvpsu
#
# Optional columns used later:
#   hba1c, diabetes
# ============================================================

core_path <- "data/derived/nhanes_core_site_level.rds"

if (!file.exists(core_path)) {
  stop(
    "Core site-level NHANES file not found at: ", core_path, "\n\n",
    "Next step: adapt your existing NHANES preparation code so that it writes this file.\n",
    "The file should contain one row per observed periodontal site, with participant ID, tooth position, CAL, survey design variables, and analysis covariates.",
    call. = FALSE
  )
}

nhanes_core <- readRDS(core_path)

check_required_columns(
  nhanes_core,
  c("seqn", "age", "tooth", "cal", "wtmec6yr", "sdmvstra", "sdmvpsu")
)

message("Loaded core NHANES site-level data: ", nrow(nhanes_core), " rows.")
