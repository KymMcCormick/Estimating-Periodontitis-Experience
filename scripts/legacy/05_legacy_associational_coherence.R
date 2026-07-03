# ============================================================
# Associational coherence with HbA1c
# ============================================================
#
# Purpose:
#   Test whether Expected Periodontitis Experience (EPE)
#   shows coherent external association with HbA1c relative to
#   conventional observed-site periodontal measures.
#
# Measures:
#   EPE
#   mean_CAL
#   extent_ge_4mm
#
# Outcome:
#   LBXGH = HbA1c (%)
# ============================================================

source("00_packages.R")

library(tidyverse)
library(survey)

options(survey.lonely.psu = "adjust")

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------

analytic <- readRDS("data/derived/analytic_dataset_epe.rds")


# ------------------------------------------------------------
# 2. Prepare analytic dataset
# ------------------------------------------------------------

assoc_data <- analytic %>%
  mutate(
    WTMEC6YR = WTMEC2YR / 3,
    age_c = Age - 55,
    pct_sites_observed = (n_sites_observed / 168) * 100
  ) %>%
  filter(
    !is.na(LBXGH),
    !is.na(Age),
    !is.na(EPE),
    !is.na(mean_CAL),
    !is.na(extent_ge_4mm),
    !is.na(pct_sites_observed),
    !is.na(WTMEC6YR),
    WTMEC6YR > 0,
    !is.na(SDMVPSU),
    !is.na(SDMVSTRA)
  )


# ------------------------------------------------------------
# 3. Weighted standardisation
# ------------------------------------------------------------
# Standardising allows coefficients to be compared across measures.
# Coefficients are interpreted as change in HbA1c percentage points
# per one weighted SD increase in the periodontal measure.

weighted_mean <- function(x, w) {
  keep <- !is.na(x) & !is.na(w) & w > 0
  sum(w[keep] * x[keep]) / sum(w[keep])
}

weighted_sd <- function(x, w) {
  keep <- !is.na(x) & !is.na(w) & w > 0
  m <- weighted_mean(x[keep], w[keep])
  sqrt(sum(w[keep] * (x[keep] - m)^2) / sum(w[keep]))
}

z_weighted <- function(x, w) {
  (x - weighted_mean(x, w)) / weighted_sd(x, w)
}


assoc_data <- assoc_data %>%
  mutate(
    z_EPE = z_weighted(EPE, WTMEC6YR),
    z_mean_CAL = z_weighted(mean_CAL, WTMEC6YR),
    z_extent_ge_4mm = z_weighted(extent_ge_4mm, WTMEC6YR),
    z_sites_observed = z_weighted(pct_sites_observed, WTMEC6YR)
  )


# ------------------------------------------------------------
# 4. Survey design
# ------------------------------------------------------------

des <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC6YR,
  nest = TRUE,
  data = assoc_data
)


# ------------------------------------------------------------
# 5. Measure labels
# ------------------------------------------------------------

measure_info <- tribble(
  ~raw_var,              ~z_var,              ~label,
  "EPE",                 "z_EPE",             "Expected Periodontitis Experience",
  "mean_CAL",            "z_mean_CAL",        "Mean CAL",
  "extent_ge_4mm",       "z_extent_ge_4mm",   "Extent CAL >=4 mm",
  "pct_sites_observed",  "z_sites_observed",  "Periodontal sites observed"
)


# ------------------------------------------------------------
# 6. Design-weighted rank correlations with HbA1c
# ------------------------------------------------------------
# This is a survey-weighted rank-correlation approximation:
# rank variables are created first, then their survey-weighted
# Pearson correlation is estimated using svyvar().

rank_data <- assoc_data %>%
  mutate(
    rank_HbA1c = rank(LBXGH, ties.method = "average"),
    rank_EPE = rank(EPE, ties.method = "average"),
    rank_mean_CAL = rank(mean_CAL, ties.method = "average"),
    rank_extent_ge_4mm = rank(extent_ge_4mm, ties.method = "average"),
    rank_sites_observed = rank(pct_sites_observed, ties.method = "average")
  )

rank_des <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC6YR,
  nest = TRUE,
  data = rank_data
)


get_rank_cor <- function(rank_measure_var, label, design_object) {
  
  f <- as.formula(
    paste0("~rank_HbA1c + ", rank_measure_var)
  )
  
  v <- svyvar(
    f,
    design = design_object,
    na.rm = TRUE
  )
  
  cor_mat <- cov2cor(as.matrix(v))
  
  tibble(
    measure = label,
    rank_correlation_with_HbA1c = cor_mat[1, 2]
  )
}


correlation_table <- tribble(
  ~rank_var,              ~label,
  "rank_EPE",             "Expected Periodontitis Experience",
  "rank_mean_CAL",        "Mean CAL",
  "rank_extent_ge_4mm",   "Extent CAL >=4 mm",
  "rank_sites_observed",  "Periodontal sites observed"
) %>%
  pmap_dfr(
    function(rank_var, label) {
      get_rank_cor(
        rank_measure_var = rank_var,
        label = label,
        design_object = rank_des
      )
    }
  ) %>%
  mutate(
    rank_correlation_with_HbA1c =
      round(rank_correlation_with_HbA1c, 3)
  )

print(correlation_table, n = Inf)


# ------------------------------------------------------------
# 7. Regression helpers
# ------------------------------------------------------------

tidy_svy_term <- function(fit, term, model, measure) {
  
  coef_table <- as.data.frame(summary(fit)$coefficients)
  coef_table$term <- rownames(coef_table)
  rownames(coef_table) <- NULL
  
  names(coef_table)[1:4] <- c(
    "estimate",
    "std_error",
    "t_value",
    "p_value"
  )
  
  ci <- as.data.frame(confint(fit))
  ci$term <- rownames(ci)
  rownames(ci) <- NULL
  names(ci)[1:2] <- c("conf_low", "conf_high")
  
  coef_table %>%
    left_join(ci, by = "term") %>%
    filter(term == !!term) %>%
    transmute(
      model = model,
      measure = measure,
      estimate,
      conf_low,
      conf_high,
      std_error,
      p_value,
      estimate_ci = paste0(
        formatC(estimate, format = "f", digits = 3),
        " (",
        formatC(conf_low, format = "f", digits = 3),
        ", ",
        formatC(conf_high, format = "f", digits = 3),
        ")"
      ),
      p_value_formatted = case_when(
        p_value < 0.001 ~ "<0.001",
        TRUE ~ formatC(p_value, format = "f", digits = 3)
      )
    )
}


fit_single_model <- function(z_var, label, design_object, age_adjusted = FALSE) {
  
  if (age_adjusted) {
    f <- as.formula(
      paste0("LBXGH ~ ", z_var, " + age_c + I(age_c^2)")
    )
    model_label <- "Age-adjusted"
  } else {
    f <- as.formula(
      paste0("LBXGH ~ ", z_var)
    )
    model_label <- "Unadjusted"
  }
  
  fit <- svyglm(
    f,
    design = design_object
  )
  
  tidy_svy_term(
    fit = fit,
    term = z_var,
    model = model_label,
    measure = label
  )
}


# ------------------------------------------------------------
# 8. Separate regression models
# ------------------------------------------------------------
# Main external alignment analysis.
# Coefficients are HbA1c percentage-point differences per 1 SD
# increase in each periodontal measure.

regression_results <- bind_rows(
  
  # Unadjusted models
  pmap_dfr(
    measure_info %>%
      filter(label != "Periodontal sites observed"),
    function(raw_var, z_var, label) {
      fit_single_model(
        z_var = z_var,
        label = label,
        design_object = des,
        age_adjusted = FALSE
      )
    }
  ),
  
  # Age-adjusted sensitivity models
  pmap_dfr(
    measure_info %>%
      filter(label != "Periodontal sites observed"),
    function(raw_var, z_var, label) {
      fit_single_model(
        z_var = z_var,
        label = label,
        design_object = des,
        age_adjusted = TRUE
      )
    }
  )
)

print(regression_results, n = Inf)


# ------------------------------------------------------------
# 9. Optional joint models
# ------------------------------------------------------------
# Interpret cautiously because EPE, mean CAL, and extent are
# overlapping periodontal summaries.

joint_unadjusted <- svyglm(
  LBXGH ~ z_EPE + z_mean_CAL + z_extent_ge_4mm,
  design = des
)

joint_age_adjusted <- svyglm(
  LBXGH ~ z_EPE + z_mean_CAL + z_extent_ge_4mm + age_c + I(age_c^2),
  design = des
)


joint_results <- bind_rows(
  tidy_svy_term(
    joint_unadjusted,
    term = "z_EPE",
    model = "Joint unadjusted",
    measure = "Expected Periodontitis Experience"
  ),
  tidy_svy_term(
    joint_unadjusted,
    term = "z_mean_CAL",
    model = "Joint unadjusted",
    measure = "Mean CAL"
  ),
  tidy_svy_term(
    joint_unadjusted,
    term = "z_extent_ge_4mm",
    model = "Joint unadjusted",
    measure = "Extent CAL >=4 mm"
  ),
  tidy_svy_term(
    joint_age_adjusted,
    term = "z_EPE",
    model = "Joint age-adjusted",
    measure = "Expected Periodontitis Experience"
  ),
  tidy_svy_term(
    joint_age_adjusted,
    term = "z_mean_CAL",
    model = "Joint age-adjusted",
    measure = "Mean CAL"
  ),
  tidy_svy_term(
    joint_age_adjusted,
    term = "z_extent_ge_4mm",
    model = "Joint age-adjusted",
    measure = "Extent CAL >=4 mm"
  )
)

print(joint_results, n = Inf)


# ------------------------------------------------------------
# 10. Save outputs
# ------------------------------------------------------------

readr::write_csv(
  correlation_table,
  "outputs/tables/hba1c_rank_correlations_periodontal_measures.csv"
)

readr::write_csv(
  regression_results,
  "outputs/tables/hba1c_single_measure_regressions.csv"
)

readr::write_csv(
  joint_results,
  "outputs/tables/hba1c_joint_measure_regressions.csv"
)