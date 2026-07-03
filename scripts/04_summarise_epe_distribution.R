# ============================================================
# Supplementary Table S1:
# Population survey-weighted means and 95% CIs by age group
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
# 2. Prepare data
# ------------------------------------------------------------

analytic_table <- analytic %>%
  mutate(
    AgeGroup_Table2 = as.character(AgeGroup_Table2),
    AgeGroup_Table2 = stringr::str_replace_all(AgeGroup_Table2, "–", "-"),
    
    AgeGroup_Table2 = factor(
      AgeGroup_Table2,
      levels = c("30-39", "40-49", "50-59", "60-69", "70-79", "80+")
    ),
    
    WTMEC6YR = WTMEC2YR / 3,
    
    # 28 teeth x 6 sites, third molars excluded
    pct_sites_observed = (n_sites_observed / 168) * 100
  ) %>%
  filter(
    !is.na(AgeGroup_Table2),
    !is.na(WTMEC6YR),
    WTMEC6YR > 0,
    !is.na(SDMVPSU),
    !is.na(SDMVSTRA)
  )


# ------------------------------------------------------------
# 3. Survey design
# ------------------------------------------------------------

des <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC6YR,
  nest = TRUE,
  data = analytic_table
)


# ------------------------------------------------------------
# 4. Measures to include
# ------------------------------------------------------------

measure_info <- tribble(
  ~var,                  ~label,                                  ~digits,
  "pct_sites_observed",  "Sites observed (%)",                     1,
  "mean_CAL",            "Mean CAL (mm)",                          2,
  "extent_ge_4mm",       "Extent CAL >=4 mm (%)",                  1,
  "epe",                 "EPE, mean-tooth scale",                  2
)


# ------------------------------------------------------------
# 5. Helper functions
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
# 6. Build table by age group
# ------------------------------------------------------------

age_groups <- levels(analytic_table$AgeGroup_Table2)

table_age <- map_dfr(age_groups, function(age_group) {
  
  group_design <- subset(
    des,
    AgeGroup_Table2 == age_group
  )
  
  group_n <- analytic_table %>%
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


table_age_wide <- table_age %>%
  pivot_wider(
    names_from = measure,
    values_from = value
  )


# ------------------------------------------------------------
# 7. Overall row
# ------------------------------------------------------------

overall_row <- map_dfr(seq_len(nrow(measure_info)), function(i) {
  
  varname <- measure_info$var[i]
  label <- measure_info$label[i]
  digits <- measure_info$digits[i]
  
  tibble(
    AgeGroup_Table2 = "Overall",
    n = nrow(analytic_table),
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


supp_table_mean_ci <- bind_rows(
  overall_row,
  table_age_wide
)


# ------------------------------------------------------------
# 8. Print and save
# ------------------------------------------------------------

print(supp_table_mean_ci, n = Inf)

readr::write_csv(
  supp_table_mean_ci,
  "outputs/tables/supp_table_population_means_ci_periodontal_measures_by_age.csv"
)