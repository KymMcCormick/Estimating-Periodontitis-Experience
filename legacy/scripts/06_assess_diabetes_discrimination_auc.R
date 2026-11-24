# ============================================================
# Discrimination of diabetes status using periodontal measures
# ============================================================
#
# Purpose:
#   Examine whether Expected Periodontitis Experience (EPE) and
#   conventional observed-site periodontal measures discriminate
#   diabetes status, defined using HbA1c >= 6.5%.
#
# Analysis set:
#   Normal glycaemia vs diabetes only; prediabetes excluded.
#
# Outputs:
#   1. diabetes_discrimination_auc_table.csv
#   2. diabetes_discrimination_auc_differences.csv
#   3. diabetes_discrimination_sample_summary.csv
#
# Notes:
#   - AUCs are survey-weighted.
#   - Confidence intervals are estimated using a stratified PSU
#     bootstrap approximation.
#   - Use B = 200 while testing and B = 1000 or 2000 for final runs.
# ============================================================

source("scripts/00_setup.R")

library(tidyverse)
library(survey)

options(survey.lonely.psu = "adjust")

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# User settings
# ------------------------------------------------------------

B <- 1000
set.seed(20260702)


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

auc_data <- analytic %>%
  mutate(
    WTMEC6YR = WTMEC2YR / 3,

    glycaemia_group = case_when(
      LBXGH < 5.7 ~ "Normal glycaemia",
      LBXGH >= 5.7 & LBXGH < 6.5 ~ "Prediabetes",
      LBXGH >= 6.5 ~ "Diabetes",
      TRUE ~ NA_character_
    ),

    diabetes_binary = case_when(
      glycaemia_group == "Diabetes" ~ 1,
      glycaemia_group == "Normal glycaemia" ~ 0,
      TRUE ~ NA_real_
    ),

    age_c = Age - 55,
    pct_sites_observed = (n_sites_observed / 168) * 100
  ) %>%
  filter(
    glycaemia_group %in% c("Normal glycaemia", "Diabetes"),
    !is.na(diabetes_binary),
    !is.na(Age),
    !is.na(epe),
    !is.na(mean_CAL),
    !is.na(extent_ge_4mm),
    !is.na(WTMEC6YR),
    WTMEC6YR > 0,
    !is.na(SDMVPSU),
    !is.na(SDMVSTRA)
  )


# ------------------------------------------------------------
# 3. Sample summary
# ------------------------------------------------------------

sample_summary <- auc_data %>%
  group_by(glycaemia_group) %>%
  summarise(
    n_unweighted = n(),
    weighted_n_approx = sum(WTMEC6YR, na.rm = TRUE),
    .groups = "drop"
  )

print(sample_summary, n = Inf)

readr::write_csv(
  sample_summary,
  "outputs/tables/diabetes_discrimination_sample_summary.csv"
)


# ------------------------------------------------------------
# 4. Weighted standardisation
# ------------------------------------------------------------
# Standardising allows coefficients or model terms to be compared
# across periodontal measures. For AUC, monotone transformations do
# not change discrimination in single-measure models, but using
# standardised variables keeps the modelling consistent with the
# HbA1c regression analysis.

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

auc_data <- auc_data %>%
  mutate(
    z_epe = z_weighted(epe, WTMEC6YR),
    z_mean_CAL = z_weighted(mean_CAL, WTMEC6YR),
    z_extent_ge_4mm = z_weighted(extent_ge_4mm, WTMEC6YR)
  )


# ------------------------------------------------------------
# 5. Survey design
# ------------------------------------------------------------

des <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC6YR,
  nest = TRUE,
  data = auc_data
)


# ------------------------------------------------------------
# 6. Weighted AUC function
# ------------------------------------------------------------
# Weighted probability that a randomly selected diabetes case has
# a higher predicted score than a randomly selected non-case.
# Ties are assigned half credit.

weighted_auc <- function(y, score, w) {

  keep <- !is.na(y) & !is.na(score) & !is.na(w) & w > 0
  y <- y[keep]
  score <- score[keep]
  w <- w[keep]

  if (!all(y %in% c(0, 1))) {
    stop("Outcome must be coded 0/1.")
  }

  total_case_w <- sum(w[y == 1])
  total_control_w <- sum(w[y == 0])

  if (total_case_w == 0 || total_control_w == 0) {
    return(NA_real_)
  }

  ord <- order(score)
  y <- y[ord]
  score <- score[ord]
  w <- w[ord]

  score_groups <- split(seq_along(score), score)

  control_less <- 0
  numerator <- 0

  for (idx in score_groups) {
    group_case_w <- sum(w[idx][y[idx] == 1])
    group_control_w <- sum(w[idx][y[idx] == 0])

    numerator <- numerator +
      group_case_w * (control_less + 0.5 * group_control_w)

    control_less <- control_less + group_control_w
  }

  numerator / (total_case_w * total_control_w)
}


# ------------------------------------------------------------
# 7. Model specifications
# ------------------------------------------------------------

model_specs <- tribble(
  ~model,          ~measure,             ~formula_string,
  "Unadjusted",    "EPE",                "diabetes_binary ~ z_epe",
  "Unadjusted",    "Mean CAL",           "diabetes_binary ~ z_mean_CAL",
  "Unadjusted",    "Extent CAL >=4 mm",  "diabetes_binary ~ z_extent_ge_4mm",
  "Age-adjusted",  "EPE",                "diabetes_binary ~ z_epe + age_c + I(age_c^2)",
  "Age-adjusted",  "Mean CAL",           "diabetes_binary ~ z_mean_CAL + age_c + I(age_c^2)",
  "Age-adjusted",  "Extent CAL >=4 mm",  "diabetes_binary ~ z_extent_ge_4mm + age_c + I(age_c^2)",
  "Age-only",      "Age only",           "diabetes_binary ~ age_c + I(age_c^2)"
)


# ------------------------------------------------------------
# 8. Point estimates using survey-weighted logistic models
# ------------------------------------------------------------

get_svy_auc <- function(formula_string, design_object, data_object) {

  fit <- svyglm(
    as.formula(formula_string),
    design = design_object,
    family = quasibinomial()
  )

  pred <- as.numeric(predict(fit, type = "response"))

  weighted_auc(
    y = data_object$diabetes_binary,
    score = pred,
    w = data_object$WTMEC6YR
  )
}

auc_point <- model_specs %>%
  mutate(
    auc = map_dbl(
      formula_string,
      ~ get_svy_auc(
        formula_string = .x,
        design_object = des,
        data_object = auc_data
      )
    )
  )

print(auc_point, n = Inf)


# ------------------------------------------------------------
# 9. Stratified PSU bootstrap for AUC confidence intervals
# ------------------------------------------------------------

psu_frame <- auc_data %>%
  distinct(SDMVSTRA, SDMVPSU)

get_boot_sample <- function(data, psu_frame) {

  boot_psus <- psu_frame %>%
    group_by(SDMVSTRA) %>%
    summarise(
      sampled_psu = list(sample(SDMVPSU, size = n(), replace = TRUE)),
      .groups = "drop"
    ) %>%
    unnest(sampled_psu) %>%
    count(SDMVSTRA, SDMVPSU = sampled_psu, name = "boot_mult")

  data %>%
    inner_join(boot_psus, by = c("SDMVSTRA", "SDMVPSU")) %>%
    mutate(boot_w = WTMEC6YR * boot_mult)
}

get_glm_auc <- function(formula_string, data_object, weight_var = "boot_w") {

  fit <- tryCatch(
    glm(
      as.formula(formula_string),
      data = data_object,
      weights = data_object[[weight_var]],
      family = quasibinomial()
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(NA_real_)
  }

  pred <- tryCatch(
    as.numeric(predict(fit, type = "response")),
    error = function(e) rep(NA_real_, nrow(data_object))
  )

  weighted_auc(
    y = data_object$diabetes_binary,
    score = pred,
    w = data_object[[weight_var]]
  )
}

boot_auc <- map_dfr(seq_len(B), function(b) {

  boot_data <- get_boot_sample(
    data = auc_data,
    psu_frame = psu_frame
  )

  model_specs %>%
    mutate(
      replicate = b,
      auc = map_dbl(
        formula_string,
        ~ get_glm_auc(
          formula_string = .x,
          data_object = boot_data,
          weight_var = "boot_w"
        )
      )
    )
})


# ------------------------------------------------------------
# 10. AUC table
# ------------------------------------------------------------

auc_ci <- boot_auc %>%
  group_by(model, measure) %>%
  summarise(
    auc_low = quantile(auc, 0.025, na.rm = TRUE),
    auc_high = quantile(auc, 0.975, na.rm = TRUE),
    n_boot = sum(!is.na(auc)),
    .groups = "drop"
  )

auc_table <- auc_point %>%
  left_join(auc_ci, by = c("model", "measure")) %>%
  mutate(
    auc_formatted = paste0(
      formatC(auc, format = "f", digits = 3),
      " (",
      formatC(auc_low, format = "f", digits = 3),
      "-",
      formatC(auc_high, format = "f", digits = 3),
      ")"
    )
  ) %>%
  select(
    model,
    measure,
    auc,
    auc_low,
    auc_high,
    auc_formatted,
    n_boot
  )

print(auc_table, n = Inf)

readr::write_csv(
  auc_table,
  "outputs/tables/diabetes_discrimination_auc_table.csv"
)


# ------------------------------------------------------------
# 11. Bootstrap differences in AUC
# ------------------------------------------------------------

make_point_diff <- function(model_name, comparison_name, reference_measure) {

  wide <- auc_point %>%
    filter(model == model_name) %>%
    select(measure, auc) %>%
    pivot_wider(names_from = measure, values_from = auc)

  tibble(
    model = model_name,
    comparison = comparison_name,
    delta_auc = wide$EPE - wide[[reference_measure]]
  )
}

point_diff <- bind_rows(
  make_point_diff("Unadjusted", "EPE - Mean CAL", "Mean CAL"),
  make_point_diff("Unadjusted", "EPE - Extent CAL >=4 mm", "Extent CAL >=4 mm"),
  make_point_diff("Age-adjusted", "EPE - Mean CAL", "Mean CAL"),
  make_point_diff("Age-adjusted", "EPE - Extent CAL >=4 mm", "Extent CAL >=4 mm")
)

boot_diff <- boot_auc %>%
  filter(model %in% c("Unadjusted", "Age-adjusted")) %>%
  select(replicate, model, measure, auc) %>%
  pivot_wider(names_from = measure, values_from = auc) %>%
  mutate(
    `EPE - Mean CAL` = EPE - `Mean CAL`,
    `EPE - Extent CAL >=4 mm` = EPE - `Extent CAL >=4 mm`
  ) %>%
  select(
    replicate,
    model,
    `EPE - Mean CAL`,
    `EPE - Extent CAL >=4 mm`
  ) %>%
  pivot_longer(
    cols = c(`EPE - Mean CAL`, `EPE - Extent CAL >=4 mm`),
    names_to = "comparison",
    values_to = "delta_auc"
  )

diff_ci <- boot_diff %>%
  group_by(model, comparison) %>%
  summarise(
    delta_low = quantile(delta_auc, 0.025, na.rm = TRUE),
    delta_high = quantile(delta_auc, 0.975, na.rm = TRUE),
    p_boot = {
      d <- delta_auc[!is.na(delta_auc)]
      2 * min(mean(d <= 0), mean(d >= 0))
    },
    n_boot = sum(!is.na(delta_auc)),
    .groups = "drop"
  )

diff_table <- point_diff %>%
  left_join(diff_ci, by = c("model", "comparison")) %>%
  mutate(
    delta_auc_formatted = paste0(
      formatC(delta_auc, format = "f", digits = 3),
      " (",
      formatC(delta_low, format = "f", digits = 3),
      "-",
      formatC(delta_high, format = "f", digits = 3),
      ")"
    ),
    p_boot_formatted = case_when(
      is.na(p_boot) ~ NA_character_,
      p_boot < 0.001 ~ "<0.001",
      TRUE ~ formatC(p_boot, format = "f", digits = 3)
    )
  ) %>%
  select(
    model,
    comparison,
    delta_auc,
    delta_low,
    delta_high,
    delta_auc_formatted,
    p_boot,
    p_boot_formatted,
    n_boot
  )

print(diff_table, n = Inf)

readr::write_csv(
  diff_table,
  "outputs/tables/diabetes_discrimination_auc_differences.csv"
)


# ------------------------------------------------------------
# 12. Compact publication table
# ------------------------------------------------------------

publication_auc_table <- auc_table %>%
  select(model, measure, auc_formatted) %>%
  pivot_wider(
    names_from = model,
    values_from = auc_formatted
  ) %>%
  select(
    measure,
    Unadjusted,
    `Age-adjusted`,
    `Age-only`
  )

print(publication_auc_table, n = Inf)

readr::write_csv(
  publication_auc_table,
  "outputs/tables/diabetes_discrimination_auc_publication_table.csv"
)

# End of file
