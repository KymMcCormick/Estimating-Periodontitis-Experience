# scripts/02_build_mce_measure.R
# ============================================================
# Build MCE and comparator periodontal measures
# ============================================================

core_path <- "data/derived/nhanes_core_site_level.rds"
measure_path <- "data/derived/nhanes_mce_measures.rds"

nhanes_core <- readRDS(core_path)

# Default: retained-tooth unit, using maximum CAL per tooth.
# Change unit = "site" if the manuscript defines MCE across sites instead.
nhanes_mce <- make_mce_measures(
  data = nhanes_core,
  id_col = "seqn",
  tooth_col = "tooth",
  cal_col = "cal",
  thresholds = c(3, 4, 5, 6),
  unit = "tooth",
  covariates = c("age", "wtmec6yr", "sdmvstra", "sdmvpsu", "hba1c", "diabetes")
)

saveRDS(nhanes_mce, measure_path)
readr::write_csv(nhanes_mce, "outputs/tables/nhanes_mce_measures_preview.csv")

message("Saved respondent-level MCE data: ", measure_path)
message("Rows: ", nrow(nhanes_mce))
