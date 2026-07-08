# ============================================================
# Supplementary Table S3:
# Component description of Expected Periodontitis Experience
# by age group
# ============================================================
#
# Purpose:
#   Describe how Expected Periodontitis Experience (EPE) is
#   assembled from observed retained-tooth burden and expected
#   missing-tooth burden across age groups.
#
# Inputs:
#   data/processed/analytic_dataset_epe.rds
#
# Outputs:
#   outputs/tables/supp_table_epe_components_by_age.csv
#   outputs/tables/supp_table_epe_components_by_age_long.csv
#   outputs/tables/epe_component_diagnostics.csv
#   outputs/tables/supp_table_epe_components_by_age.rds
#
# Notes:
#   - EPE is expressed on a mean-tooth scale.
#   - Component contributions are expressed on the same mean-tooth
#     scale, so that:
#
#       observed retained-tooth component
#       + expected missing-tooth component
#       = EPE
#
#   - Third molars are excluded upstream during EPE construction.
#   - This table is descriptive; it is intended to make the EPE
#     construction transparent rather than to test an association.
# ============================================================

source("scripts/00_setup.R")

library(tidyverse)
library(survey)

options(survey.lonely.psu = "adjust")

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------

analytic <- readRDS("data/processed/analytic_dataset_epe.rds")

# Compatibility with earlier local versions that used uppercase EPE.
if ("EPE" %in% names(analytic) && !"epe" %in% names(analytic)) {
  analytic <- analytic %>%
    dplyr::mutate(epe = EPE)
}


# ------------------------------------------------------------
# 2. Check required variables
# ------------------------------------------------------------

required_vars <- c(
  "AgeGroup_Table2",
  "WTMEC2YR",
  "SDMVPSU",
  "SDMVSTRA",
  "n_present_teeth",
  "n_missing_teeth",
  "observed_retained_burden",
  "expected_missing_burden",
  "pi_cumulative",
  "assigned_missing_cal",
  "epe"
)

missing_vars <- setdiff(required_vars, names(analytic))

if (length(missing_vars) > 0) {
  stop(
    "analytic_dataset_epe.rds is missing required variables: ",
    paste(missing_vars, collapse = ", "),
    ". Re-run scripts/02_build_epe_iteration01.R before this script."
  )
}


# ------------------------------------------------------------
# 3. Prepare data
# ------------------------------------------------------------

component_data <- analytic %>%
  mutate(
    AgeGroup_Table2 = as.character(AgeGroup_Table2),
    AgeGroup_Table2 = stringr::str_replace_all(AgeGroup_Table2, "–", "-"),
    AgeGroup_Table2 = factor(
      AgeGroup_Table2,
      levels = c("30-39", "40-49", "50-59", "60-69", "70-79", "80+")
    ),

    # NHANES 2009-2014 combines three 2-year cycles.
    WTMEC6YR = WTMEC2YR / 3,

    # EPE denominator: eligible tooth positions excluding third molars.
    n_eligible_teeth = n_present_teeth + n_missing_teeth,

    # Component contributions expressed on the same mean-tooth scale as EPE.
    observed_retained_component = observed_retained_burden / n_eligible_teeth,
    expected_missing_component = expected_missing_burden / n_eligible_teeth,
    epe_component_sum = observed_retained_component + expected_missing_component,

    # More interpretable descriptive versions.
    pi_cumulative_pct = pi_cumulative * 100,
    expected_missing_contribution_per_missing_tooth =
      pi_cumulative * assigned_missing_cal
  ) %>%
  filter(
    !is.na(AgeGroup_Table2),
    !is.na(WTMEC6YR),
    WTMEC6YR > 0,
    !is.na(SDMVPSU),
    !is.na(SDMVSTRA),
    is.finite(n_present_teeth),
    is.finite(n_missing_teeth),
    is.finite(n_eligible_teeth),
    n_eligible_teeth > 0,
    is.finite(observed_retained_burden),
    is.finite(expected_missing_burden),
    is.finite(pi_cumulative),
    is.finite(assigned_missing_cal),
    is.finite(epe)
  )

message("Rows in EPE component table data: ", nrow(component_data))

if (nrow(component_data) == 0) {
  stop(
    "No usable observations remain for the EPE component table. ",
    "Check analytic_dataset_epe.rds and required EPE component variables."
  )
}


# ------------------------------------------------------------
# 4. Diagnostics
# ------------------------------------------------------------

component_diagnostics <- component_data %>%
  summarise(
    n_unweighted = n(),
    min_n_eligible_teeth = min(n_eligible_teeth, na.rm = TRUE),
    max_n_eligible_teeth = max(n_eligible_teeth, na.rm = TRUE),
    n_not_28_eligible_teeth = sum(n_eligible_teeth != 28, na.rm = TRUE),
    max_abs_epe_component_difference =
      max(abs(epe - epe_component_sum), na.rm = TRUE),
    mean_abs_epe_component_difference =
      mean(abs(epe - epe_component_sum), na.rm = TRUE)
  )

print(component_diagnostics)

readr::write_csv(
  component_diagnostics,
  "outputs/tables/epe_component_diagnostics.csv"
)

if (component_diagnostics$max_abs_epe_component_difference > 1e-8) {
  warning(
    "EPE does not exactly equal observed_retained_component + ",
    "expected_missing_component for at least one participant. ",
    "This may reflect rounding, a changed EPE definition, or mismatched variables."
  )
}

if (component_diagnostics$n_not_28_eligible_teeth > 0) {
  warning(
    "Some participants do not have 28 eligible tooth positions. ",
    "Check whether upstream exclusions or tooth-position definitions changed."
  )
}


# ------------------------------------------------------------
# 5. Survey design
# ------------------------------------------------------------

des <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC6YR,
  nest = TRUE,
  data = component_data
)


# ------------------------------------------------------------
# 6. Measures to include
# ------------------------------------------------------------

measure_info <- tribble(
  ~var,                                           ~label,                                                   ~digits,
  "n_present_teeth",                              "Retained teeth, mean",                                  1,
  "n_missing_teeth",                              "Missing tooth positions, mean",                         1,
  "pi_cumulative_pct",                            "Cumulative missing-tooth attribution weight (%)",       1,
  "assigned_missing_cal",                         "Assigned advanced CAL for missing teeth (mm)",           2,
  "expected_missing_contribution_per_missing_tooth", "Expected contribution per missing tooth position (mm)", 2,
  "observed_retained_component",                  "Observed retained-tooth component",                      2,
  "expected_missing_component",                   "Expected missing-tooth component",                       2,
  "epe",                                          "EPE, mean-tooth scale",                                 2
)


# ------------------------------------------------------------
# 7. Helper functions
# ------------------------------------------------------------

format_mean_ci <- function(mean, lower, upper, digits = 2) {

  if (any(is.na(c(mean, lower, upper)))) {
    return(NA_character_)
  }

  paste0(
    formatC(mean, format = "f", digits = digits),
    " (",
    formatC(lower, format = "f", digits = digits),
    "-",
    formatC(upper, format = "f", digits = digits),
    ")"
  )
}

get_svy_mean_ci <- function(design_object, varname, digits = 2) {

  f_var <- as.formula(paste0("~", varname))

  est <- svymean(
    f_var,
    design = design_object,
    na.rm = TRUE
  )

  mean_val <- as.numeric(coef(est)[1])
  ci_vals <- as.numeric(confint(est)[1, ])

  format_mean_ci(
    mean = mean_val,
    lower = ci_vals[1],
    upper = ci_vals[2],
    digits = digits
  )
}


# ------------------------------------------------------------
# 8. Build table by age group
# ------------------------------------------------------------

age_groups <- levels(component_data$AgeGroup_Table2)

table_age_long <- map_dfr(age_groups, function(age_group) {

  group_design <- subset(
    des,
    AgeGroup_Table2 == age_group
  )

  group_n <- component_data %>%
    filter(AgeGroup_Table2 == age_group) %>%
    nrow()

  map_dfr(seq_len(nrow(measure_info)), function(i) {

    varname <- measure_info$var[i]
    label <- measure_info$label[i]
    digits <- measure_info$digits[i]

    tibble(
      AgeGroup_Table2 = age_group,
      n = group_n,
      measure = label,
      value = get_svy_mean_ci(
        design_object = group_design,
        varname = varname,
        digits = digits
      )
    )
  })
})

table_age_wide <- table_age_long %>%
  pivot_wider(
    names_from = measure,
    values_from = value
  )


# ------------------------------------------------------------
# 9. Overall row
# ------------------------------------------------------------

overall_row <- map_dfr(seq_len(nrow(measure_info)), function(i) {

  varname <- measure_info$var[i]
  label <- measure_info$label[i]
  digits <- measure_info$digits[i]

  tibble(
    AgeGroup_Table2 = "Overall",
    n = nrow(component_data),
    measure = label,
    value = get_svy_mean_ci(
      design_object = des,
      varname = varname,
      digits = digits
    )
  )
}) %>%
  pivot_wider(
    names_from = measure,
    values_from = value
  )


supp_table_epe_components_by_age <- bind_rows(
  overall_row,
  table_age_wide
)


# ------------------------------------------------------------
# 10. Print and save
# ------------------------------------------------------------

print(supp_table_epe_components_by_age, n = Inf)

readr::write_csv(
  supp_table_epe_components_by_age,
  "outputs/tables/supp_table_epe_components_by_age.csv"
)

readr::write_csv(
  table_age_long,
  "outputs/tables/supp_table_epe_components_by_age_long.csv"
)

saveRDS(
  supp_table_epe_components_by_age,
  "outputs/tables/supp_table_epe_components_by_age.rds"
)

message("Saved: outputs/tables/supp_table_epe_components_by_age.csv")
message("Saved: outputs/tables/supp_table_epe_components_by_age_long.csv")
message("Saved: outputs/tables/epe_component_diagnostics.csv")
message("Saved: outputs/tables/supp_table_epe_components_by_age.rds")


# ------------------------------------------------------------
# 11. Suggested table note for manuscript/supplement
# ------------------------------------------------------------

cat(
  "\nSuggested table note:\n",
  "Values are survey-weighted means with 95% confidence intervals in parentheses, ",
  "estimated using the NHANES complex survey design. Component contributions ",
  "are expressed on the same mean-tooth scale as EPE. The observed retained-tooth ",
  "component equals the sum of retained-tooth maximum CAL values divided by the ",
  "number of eligible tooth positions. The expected missing-tooth component equals ",
  "the expected missing-tooth burden divided by the number of eligible tooth ",
  "positions. EPE is the sum of these two mean-tooth components. Third molars ",
  "were excluded. CAL = clinical attachment loss; EPE = Expected Periodontitis ",
  "Experience.\n",
  sep = ""
)

