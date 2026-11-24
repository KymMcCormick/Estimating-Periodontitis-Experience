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
    !is.na(epe),
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
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  sum(w[keep] * x[keep]) / sum(w[keep])
}

weighted_sd <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  m <- weighted_mean(x[keep], w[keep])
  v <- sum(w[keep] * (x[keep] - m)^2) / sum(w[keep])
  if (!is.finite(v) || v <= 0) return(NA_real_)
  sqrt(v)
}

z_weighted <- function(x, w) {
  m <- weighted_mean(x, w)
  s <- weighted_sd(x, w)
  if (!is.finite(m) || !is.finite(s) || s <= 0) {
    return(rep(NA_real_, length(x)))
  }
  (x - m) / s
}


assoc_data <- assoc_data %>%
  mutate(
    z_epe = z_weighted(epe, WTMEC6YR),
    z_mean_CAL = z_weighted(mean_CAL, WTMEC6YR),
    z_extent_ge_4mm = z_weighted(extent_ge_4mm, WTMEC6YR),
    z_sites_observed = z_weighted(pct_sites_observed, WTMEC6YR)
  )

# Early diagnostic checks. These make failures much clearer than the
# downstream survey/glm error produced when a model has zero usable rows.
message("HbA1c coherence analytic rows after filtering: ", nrow(assoc_data))
if (nrow(assoc_data) == 0) {
  stop(
    "No usable observations remain for the HbA1c coherence analysis. ",
    "Check that analytic_dataset_epe.rds contains non-missing LBXGH, Age, epe, ",
    "mean_CAL, extent_ge_4mm, n_sites_observed, WTMEC2YR, SDMVPSU, and SDMVSTRA."
  )
}

standardisation_check <- assoc_data %>%
  summarise(
    n = n(),
    n_hba1c = sum(is.finite(LBXGH)),
    n_epe = sum(is.finite(epe)),
    sd_epe = weighted_sd(epe, WTMEC6YR),
    n_z_epe = sum(is.finite(z_epe)),
    n_z_mean_CAL = sum(is.finite(z_mean_CAL)),
    n_z_extent_ge_4mm = sum(is.finite(z_extent_ge_4mm)),
    n_z_sites_observed = sum(is.finite(z_sites_observed))
  )
print(standardisation_check)

bad_z <- c("z_epe", "z_mean_CAL", "z_extent_ge_4mm", "z_sites_observed")[
  vapply(
    c("z_epe", "z_mean_CAL", "z_extent_ge_4mm", "z_sites_observed"),
    function(v) sum(is.finite(assoc_data[[v]])) == 0,
    logical(1)
  )
]
if (length(bad_z) > 0) {
  stop(
    "The following standardised variables have no finite values: ",
    paste(bad_z, collapse = ", "),
    ". This usually means the source variable is entirely missing or has zero weighted variance."
  )
}


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
  "epe",                 "z_epe",             "Expected Periodontitis Experience",
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
    rank_epe = rank(epe, ties.method = "average"),
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
  "rank_epe",             "Expected Periodontitis Experience",
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
  
  # Prevent the unhelpful survey/glm error:
  # "no observations informative" / "object 'fit' not found".
  model_data <- design_object$variables
  keep <- is.finite(model_data$LBXGH) & is.finite(model_data[[z_var]])
  if (age_adjusted) {
    keep <- keep & is.finite(model_data$age_c)
  }
  
  if (sum(keep) == 0) {
    stop(
      "No finite observations available for ", label, " (", z_var, ") in the ",
      model_label, " model. Check missingness and weighted standardisation."
    )
  }
  
  design_model <- subset(
    update(design_object, .model_keep = keep),
    .model_keep
  )
  
  fit <- survey::svyglm(
    f,
    design = design_model
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

tibble::as_tibble(regression_results) |>
  print(n = Inf)


# ------------------------------------------------------------
# 9. Optional joint models
# ------------------------------------------------------------
# Interpret cautiously because EPE, mean CAL, and extent are
# overlapping periodontal summaries.

joint_keep <- with(
  des$variables,
  is.finite(LBXGH) &
    is.finite(z_epe) &
    is.finite(z_mean_CAL) &
    is.finite(z_extent_ge_4mm)
)

joint_age_keep <- with(
  des$variables,
  joint_keep & is.finite(age_c)
)

if (sum(joint_keep) == 0) {
  stop("No finite observations available for the joint unadjusted model.")
}
if (sum(joint_age_keep) == 0) {
  stop("No finite observations available for the joint age-adjusted model.")
}

joint_des <- subset(
  update(des, .joint_keep = joint_keep),
  .joint_keep
)

joint_age_des <- subset(
  update(des, .joint_age_keep = joint_age_keep),
  .joint_age_keep
)

joint_unadjusted <- survey::svyglm(
  LBXGH ~ z_epe + z_mean_CAL + z_extent_ge_4mm,
  design = joint_des
)

joint_age_adjusted <- survey::svyglm(
  LBXGH ~ z_epe + z_mean_CAL + z_extent_ge_4mm + age_c + I(age_c^2),
  design = joint_age_des
)


joint_results <- bind_rows(
  tidy_svy_term(
    joint_unadjusted,
    term = "z_epe",
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
    term = "z_epe",
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

tibble::as_tibble(joint_results) |>
  print(n = Inf)


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
