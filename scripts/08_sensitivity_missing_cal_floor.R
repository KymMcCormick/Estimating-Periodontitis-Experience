# ============================================================
# Sensitivity analysis:
# Assigned minimum CAL value for missing tooth positions
# ============================================================
#
# Purpose:
#   Test whether the main EPE findings depend on the minimum
#   assigned advanced CAL value used for missing tooth positions.
#
# Primary EPE specification:
#   assigned_missing_cal_i = max(10, max_observed_cal_i)
#
# Sensitivity specifications:
#   assigned_missing_cal_i(c) = max(c, max_observed_cal_i)
#
# Interpretation:
#   This is a robustness check, not an optimisation exercise.
#   Do not select the value of c that maximises the HbA1c association
#   or diabetes AUC. The primary value should be chosen a priori on
#   conceptual and clinical grounds.
#
# Inputs:
#   data/processed/analytic_dataset_epe.rds
#
# Outputs:
#   outputs/tables/supp_table_missing_cal_floor_sensitivity.csv
#   outputs/tables/supp_table_missing_cal_floor_sensitivity_delta_from_10mm.csv
#   outputs/tables/supp_table_missing_cal_floor_late_life_curvature.csv
#   outputs/figures/missing_cal_floor_sensitivity_hba1c_beta.png
#   outputs/figures/missing_cal_floor_sensitivity_auc.png
#   outputs/figures/missing_cal_floor_sensitivity_late_life_curvature.png
# ============================================================

source("scripts/00_setup.R")

library(tidyverse)
library(survey)

options(survey.lonely.psu = "adjust")

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 1. User settings
# ------------------------------------------------------------

# Primary value used in the manuscript.
primary_missing_cal_floor <- 10

# Use c(6, 8, 10, 12, 14) for a compact supplementary table.
# Use 6:14 for a fuller sensitivity curve.
missing_cal_floors <- c(6, 8, 10, 12, 14)
# missing_cal_floors <- 6:14


# ------------------------------------------------------------
# 2. Load analytic dataset
# ------------------------------------------------------------

analytic <- readRDS("data/processed/analytic_dataset_epe.rds")

# Compatibility with earlier local versions that used uppercase EPE.
if ("EPE" %in% names(analytic) && !"epe" %in% names(analytic)) {
  analytic <- analytic %>%
    dplyr::mutate(epe = EPE)
}

required_vars <- c(
  "SEQN",
  "Age",
  "DiabetesStatus",
  "LBXGH",
  "WTMEC2YR",
  "SDMVPSU",
  "SDMVSTRA",
  "n_present_teeth",
  "n_missing_teeth",
  "observed_retained_burden",
  "pi_cumulative",
  "max_observed_cal",
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
# 3. Helper functions
# ------------------------------------------------------------

format_ci <- function(est, lcl, ucl, digits = 3) {
  paste0(
    formatC(est, format = "f", digits = digits),
    " (",
    formatC(lcl, format = "f", digits = digits),
    "-",
    formatC(ucl, format = "f", digits = digits),
    ")"
  )
}

weighted_mean <- function(x, w) {
  keep <- complete.cases(x, w)
  weighted.mean(x[keep], w[keep], na.rm = TRUE)
}

weighted_sd <- function(x, w) {
  keep <- complete.cases(x, w)
  x <- x[keep]
  w <- w[keep]
  m <- weighted.mean(x, w, na.rm = TRUE)
  sqrt(weighted.mean((x - m)^2, w, na.rm = TRUE))
}

weighted_cor <- function(x, y, w) {
  keep <- complete.cases(x, y, w)

  x <- x[keep]
  y <- y[keep]
  w <- w[keep]

  mx <- weighted.mean(x, w)
  my <- weighted.mean(y, w)

  cov_xy <- weighted.mean((x - mx) * (y - my), w)
  var_x <- weighted.mean((x - mx)^2, w)
  var_y <- weighted.mean((y - my)^2, w)

  cov_xy / sqrt(var_x * var_y)
}

weighted_spearman <- function(x, y, w) {
  keep <- complete.cases(x, y, w)

  weighted_cor(
    rank(x[keep], ties.method = "average"),
    rank(y[keep], ties.method = "average"),
    w[keep]
  )
}

weighted_auc <- function(y, score, w) {
  keep <- complete.cases(y, score, w)

  y <- y[keep]
  score <- score[keep]
  w <- w[keep]

  if (!all(y %in% c(0, 1))) {
    stop("Outcome y must be coded 0/1.")
  }

  ord <- order(score)
  y <- y[ord]
  score <- score[ord]
  w <- w[ord]

  pos <- y == 1
  neg <- y == 0

  w_pos_total <- sum(w[pos])
  w_neg_total <- sum(w[neg])

  if (w_pos_total == 0 || w_neg_total == 0) {
    return(NA_real_)
  }

  auc_num <- 0
  neg_cum_before <- 0

  for (s in unique(score)) {
    in_tie <- score == s

    w_neg_tie <- sum(w[in_tie & neg])
    w_pos_tie <- sum(w[in_tie & pos])

    auc_num <- auc_num +
      w_pos_tie * (neg_cum_before + 0.5 * w_neg_tie)

    neg_cum_before <- neg_cum_before + w_neg_tie
  }

  auc_num / (w_pos_total * w_neg_total)
}

make_design <- function(data) {
  svydesign(
    ids = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~WTMEC6YR,
    nest = TRUE,
    data = data
  )
}

calc_epe_with_floor <- function(data, missing_cal_floor) {
  data %>%
    mutate(
      missing_cal_floor = missing_cal_floor,
      WTMEC6YR = WTMEC2YR / 3,

      # EPE denominator: eligible tooth positions excluding third molars.
      n_eligible_teeth = n_present_teeth + n_missing_teeth,

      assigned_missing_cal_sens = pmax(
        missing_cal_floor,
        max_observed_cal,
        na.rm = TRUE
      ),

      # This line is mainly protective. The current analytic sample should
      # include dentate participants, so max_observed_cal should not be NA.
      assigned_missing_cal_sens = if_else(
        is.na(assigned_missing_cal_sens),
        as.numeric(missing_cal_floor),
        assigned_missing_cal_sens
      ),

      expected_missing_burden_sens =
        n_missing_teeth * pi_cumulative * assigned_missing_cal_sens,

      observed_retained_component_sens =
        observed_retained_burden / n_eligible_teeth,

      expected_missing_component_sens =
        expected_missing_burden_sens / n_eligible_teeth,

      epe_sens =
        observed_retained_component_sens + expected_missing_component_sens,

      expected_missing_contribution_per_missing_tooth_sens =
        pi_cumulative * assigned_missing_cal_sens
    )
}

add_weighted_standardised_epe <- function(data) {
  m <- weighted_mean(data$epe_sens, data$WTMEC6YR)
  s <- weighted_sd(data$epe_sens, data$WTMEC6YR)

  data %>%
    mutate(epe_sens_z = (epe_sens - m) / s)
}

get_beta <- function(fit, term = "epe_sens_z") {
  est <- coef(fit)[term]
  ci <- confint(fit)[term, ]

  tibble(
    beta = as.numeric(est),
    lcl = as.numeric(ci[1]),
    ucl = as.numeric(ci[2])
  )
}

get_svy_mean <- function(design_object, varname) {
  est <- svymean(
    as.formula(paste0("~", varname)),
    design = design_object,
    na.rm = TRUE
  )

  as.numeric(coef(est)[1])
}

get_svy_mean_ci <- function(design_object, varname) {
  est <- svymean(
    as.formula(paste0("~", varname)),
    design = design_object,
    na.rm = TRUE
  )

  tibble(
    mean = as.numeric(coef(est)[1]),
    lcl = as.numeric(confint(est)[1, 1]),
    ucl = as.numeric(confint(est)[1, 2])
  )
}

late_life_curvature_one_group <- function(design_object, glycaemic_status) {

  group_design <- subset(
    design_object,
    DiabetesStatus == glycaemic_status
  )

  fit <- svyglm(
    epe_sens ~ Age + I(Age^2),
    design = group_design
  )

  newdat <- tibble(Age = c(60, 70, 80))
  X <- model.matrix(delete.response(terms(fit)), newdat)

  L_60_70 <- X[2, ] - X[1, ]
  L_70_80 <- X[3, ] - X[2, ]
  L_curvature <- X[3, ] - 2 * X[2, ] + X[1, ]

  contrast_est <- function(L) {
    est <- drop(L %*% coef(fit))
    se <- sqrt(drop(L %*% vcov(fit) %*% L))

    tibble(
      estimate = est,
      se = se,
      lcl = est - 1.96 * se,
      ucl = est + 1.96 * se
    )
  }

  bind_rows(
    contrast_est(L_60_70) %>% mutate(quantity = "change_60_70"),
    contrast_est(L_70_80) %>% mutate(quantity = "change_70_80"),
    contrast_est(L_curvature) %>% mutate(quantity = "late_life_curvature")
  ) %>%
    mutate(DiabetesStatus = as.character(glycaemic_status)) %>%
    select(DiabetesStatus, quantity, estimate, se, lcl, ucl)
}

analyse_one_floor <- function(missing_cal_floor, data) {

  sens_data <- data %>%
    calc_epe_with_floor(missing_cal_floor = missing_cal_floor) %>%
    filter(
      !is.na(WTMEC6YR),
      WTMEC6YR > 0,
      !is.na(SDMVPSU),
      !is.na(SDMVSTRA),
      is.finite(epe_sens),
      is.finite(Age)
    ) %>%
    add_weighted_standardised_epe()

  des <- make_design(sens_data)

  # Check that floor = 10 reproduces the primary EPE variable.
  if (missing_cal_floor == primary_missing_cal_floor) {
    max_abs_diff <- max(
      abs(sens_data$epe_sens - sens_data$epe),
      na.rm = TRUE
    )

    message(
      "Maximum absolute difference between recalculated EPE at ",
      primary_missing_cal_floor,
      " mm and saved epe variable: ",
      signif(max_abs_diff, 4)
    )

    if (max_abs_diff > 1e-8) {
      warning(
        "Recalculated EPE at the primary missing-CAL floor does not ",
        "match the saved epe variable exactly. Check whether the EPE ",
        "definition changed upstream."
      )
    }
  }

  # Descriptive means.
  mean_epe <- get_svy_mean_ci(des, "epe_sens") %>%
    rename(
      mean_epe = mean,
      mean_epe_lcl = lcl,
      mean_epe_ucl = ucl
    )

  mean_assigned_cal <- get_svy_mean_ci(des, "assigned_missing_cal_sens") %>%
    rename(
      mean_assigned_missing_cal = mean,
      mean_assigned_missing_cal_lcl = lcl,
      mean_assigned_missing_cal_ucl = ucl
    )

  mean_expected_component <- get_svy_mean_ci(des, "expected_missing_component_sens") %>%
    rename(
      mean_expected_missing_component = mean,
      mean_expected_missing_component_lcl = lcl,
      mean_expected_missing_component_ucl = ucl
    )

  # HbA1c models.
  hba1c_data <- sens_data %>%
    filter(
      !is.na(LBXGH),
      is.finite(epe_sens_z),
      is.finite(Age)
    )

  hba1c_des <- make_design(hba1c_data)

  fit_hba1c_unadj <- svyglm(
    LBXGH ~ epe_sens_z,
    design = hba1c_des
  )

  fit_hba1c_ageadj <- svyglm(
    LBXGH ~ epe_sens_z + Age + I(Age^2),
    design = hba1c_des
  )

  beta_unadj <- get_beta(fit_hba1c_unadj) %>%
    rename(
      beta_hba1c_unadj = beta,
      beta_hba1c_unadj_lcl = lcl,
      beta_hba1c_unadj_ucl = ucl
    )

  beta_ageadj <- get_beta(fit_hba1c_ageadj) %>%
    rename(
      beta_hba1c_ageadj = beta,
      beta_hba1c_ageadj_lcl = lcl,
      beta_hba1c_ageadj_ucl = ucl
    )

  rho_hba1c <- weighted_spearman(
    hba1c_data$epe_sens,
    hba1c_data$LBXGH,
    hba1c_data$WTMEC6YR
  )

  # Diabetes-status discrimination: Normal vs Diabetes, excluding Prediabetes.
  auc_data <- sens_data %>%
    filter(
      DiabetesStatus %in% c("Normal", "Diabetes"),
      !is.na(DiabetesStatus),
      is.finite(epe_sens_z),
      is.finite(Age)
    ) %>%
    mutate(
      diabetes_binary = if_else(DiabetesStatus == "Diabetes", 1, 0)
    )

  auc_des <- make_design(auc_data)

  fit_auc_unadj <- svyglm(
    diabetes_binary ~ epe_sens_z,
    design = auc_des,
    family = quasibinomial()
  )

  fit_auc_ageadj <- svyglm(
    diabetes_binary ~ epe_sens_z + Age + I(Age^2),
    design = auc_des,
    family = quasibinomial()
  )

  auc_unadj <- weighted_auc(
    y = auc_data$diabetes_binary,
    score = as.numeric(predict(fit_auc_unadj, type = "response")),
    w = auc_data$WTMEC6YR
  )

  auc_ageadj <- weighted_auc(
    y = auc_data$diabetes_binary,
    score = as.numeric(predict(fit_auc_ageadj, type = "response")),
    w = auc_data$WTMEC6YR
  )

  # Late-life curvature for EPE by glycaemic status.
  curvature <- map_dfr(
    levels(droplevels(sens_data$DiabetesStatus)),
    ~late_life_curvature_one_group(des, .x)
  ) %>%
    mutate(missing_cal_floor = missing_cal_floor) %>%
    relocate(missing_cal_floor)

  curvature_wide <- curvature %>%
    filter(quantity == "late_life_curvature") %>%
    select(missing_cal_floor, DiabetesStatus, estimate) %>%
    mutate(DiabetesStatus = paste0("curvature_", DiabetesStatus)) %>%
    pivot_wider(
      names_from = DiabetesStatus,
      values_from = estimate
    )

  summary_row <- tibble(
    missing_cal_floor = missing_cal_floor,
    n_unweighted = nrow(sens_data),
    n_hba1c = nrow(hba1c_data),
    n_auc = nrow(auc_data),
    rho_hba1c = rho_hba1c,
    auc_unadj = auc_unadj,
    auc_ageadj = auc_ageadj
  ) %>%
    bind_cols(
      mean_epe,
      mean_assigned_cal,
      mean_expected_component,
      beta_unadj,
      beta_ageadj
    ) %>%
    left_join(curvature_wide, by = "missing_cal_floor")

  list(
    summary = summary_row,
    curvature = curvature
  )
}


# ------------------------------------------------------------
# 4. Run sensitivity analysis
# ------------------------------------------------------------

sensitivity_results <- map(
  missing_cal_floors,
  ~analyse_one_floor(missing_cal_floor = .x, data = analytic)
)

missing_cal_floor_sensitivity <- map_dfr(
  sensitivity_results,
  "summary"
)

missing_cal_floor_curvature <- map_dfr(
  sensitivity_results,
  "curvature"
)


# ------------------------------------------------------------
# 5. Delta from primary 10 mm specification
# ------------------------------------------------------------

primary_row <- missing_cal_floor_sensitivity %>%
  filter(missing_cal_floor == primary_missing_cal_floor)

if (nrow(primary_row) != 1) {
  stop("Primary missing-CAL floor not found in sensitivity table.")
}

numeric_vars <- names(missing_cal_floor_sensitivity)[
  vapply(missing_cal_floor_sensitivity, is.numeric, logical(1))
]

numeric_vars <- setdiff(numeric_vars, "missing_cal_floor")

missing_cal_floor_sensitivity_delta <- missing_cal_floor_sensitivity

for (v in numeric_vars) {
  missing_cal_floor_sensitivity_delta[[paste0("delta_", v)]] <-
    missing_cal_floor_sensitivity_delta[[v]] - primary_row[[v]]
}


# ------------------------------------------------------------
# 6. Manuscript-friendly formatted table
# ------------------------------------------------------------

missing_cal_floor_sensitivity_formatted <- missing_cal_floor_sensitivity %>%
  transmute(
    `Minimum assigned CAL value for missing teeth (mm)` = missing_cal_floor,
    `Mean EPE` = format_ci(mean_epe, mean_epe_lcl, mean_epe_ucl, digits = 2),
    `Expected missing-tooth component` = format_ci(
      mean_expected_missing_component,
      mean_expected_missing_component_lcl,
      mean_expected_missing_component_ucl,
      digits = 2
    ),
    `HbA1c beta, unadjusted` = format_ci(
      beta_hba1c_unadj,
      beta_hba1c_unadj_lcl,
      beta_hba1c_unadj_ucl,
      digits = 3
    ),
    `HbA1c beta, age-adjusted` = format_ci(
      beta_hba1c_ageadj,
      beta_hba1c_ageadj_lcl,
      beta_hba1c_ageadj_ucl,
      digits = 3
    ),
    `Rank correlation with HbA1c` = formatC(
      rho_hba1c,
      format = "f",
      digits = 3
    ),
    `AUC, unadjusted` = formatC(
      auc_unadj,
      format = "f",
      digits = 3
    ),
    `AUC, age-adjusted` = formatC(
      auc_ageadj,
      format = "f",
      digits = 3
    ),
    `Late-life curvature, normal glycaemia` = formatC(
      curvature_Normal,
      format = "f",
      digits = 3
    ),
    `Late-life curvature, prediabetes` = formatC(
      curvature_Prediabetes,
      format = "f",
      digits = 3
    ),
    `Late-life curvature, diabetes` = formatC(
      curvature_Diabetes,
      format = "f",
      digits = 3
    )
  )


# ------------------------------------------------------------
# 7. Save tables
# ------------------------------------------------------------

readr::write_csv(
  missing_cal_floor_sensitivity,
  "outputs/tables/supp_table_missing_cal_floor_sensitivity.csv"
)

readr::write_csv(
  missing_cal_floor_sensitivity_delta,
  "outputs/tables/supp_table_missing_cal_floor_sensitivity_delta_from_10mm.csv"
)

readr::write_csv(
  missing_cal_floor_sensitivity_formatted,
  "outputs/tables/supp_table_missing_cal_floor_sensitivity_formatted.csv"
)

readr::write_csv(
  missing_cal_floor_curvature,
  "outputs/tables/supp_table_missing_cal_floor_late_life_curvature.csv"
)

print(missing_cal_floor_sensitivity_formatted, n = Inf)


# ------------------------------------------------------------
# 8. Figures
# ------------------------------------------------------------

hba1c_plot <- ggplot(
  missing_cal_floor_sensitivity,
  aes(x = missing_cal_floor, y = beta_hba1c_ageadj)
) +
  geom_line() +
  geom_point() +
  geom_errorbar(
    aes(
      ymin = beta_hba1c_ageadj_lcl,
      ymax = beta_hba1c_ageadj_ucl
    ),
    width = 0.15
  ) +
  geom_vline(
    xintercept = primary_missing_cal_floor,
    linetype = "dashed"
  ) +
  labs(
    title = "Sensitivity of EPE-HbA1c association to assigned missing-tooth CAL",
    x = "Minimum assigned CAL value for missing teeth (mm)",
    y = "Age-adjusted HbA1c coefficient for EPE"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  "outputs/figures/missing_cal_floor_sensitivity_hba1c_beta.png",
  plot = hba1c_plot,
  width = 8,
  height = 5,
  dpi = 300
)

auc_plot <- ggplot(
  missing_cal_floor_sensitivity,
  aes(x = missing_cal_floor, y = auc_ageadj)
) +
  geom_line() +
  geom_point() +
  geom_vline(
    xintercept = primary_missing_cal_floor,
    linetype = "dashed"
  ) +
  labs(
    title = "Sensitivity of diabetes-status discrimination to assigned missing-tooth CAL",
    x = "Minimum assigned CAL value for missing teeth (mm)",
    y = "Age-adjusted AUC"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  "outputs/figures/missing_cal_floor_sensitivity_auc.png",
  plot = auc_plot,
  width = 8,
  height = 5,
  dpi = 300
)

curvature_plot_data <- missing_cal_floor_curvature %>%
  filter(quantity == "late_life_curvature")

curvature_plot <- ggplot(
  curvature_plot_data,
  aes(
    x = missing_cal_floor,
    y = estimate,
    group = DiabetesStatus,
    linetype = DiabetesStatus,
    shape = DiabetesStatus
  )
) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_line() +
  geom_point() +
  geom_vline(
    xintercept = primary_missing_cal_floor,
    linetype = "dashed"
  ) +
  labs(
    title = "Sensitivity of EPE late-life curvature to assigned missing-tooth CAL",
    x = "Minimum assigned CAL value for missing teeth (mm)",
    y = "Late-life curvature",
    linetype = "Glycaemic status",
    shape = "Glycaemic status"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  "outputs/figures/missing_cal_floor_sensitivity_late_life_curvature.png",
  plot = curvature_plot,
  width = 8,
  height = 5,
  dpi = 300
)


# ------------------------------------------------------------
# 9. Suggested note
# ------------------------------------------------------------

cat(
  "\nSuggested supplementary note:\n",
  "Sensitivity analyses recalculated EPE after varying the minimum assigned ",
  "advanced CAL value for missing tooth positions. The primary specification ",
  "assigned each missing tooth position an advanced CAL value equal to the ",
  "greater of 10 mm or the participant's maximum observed tooth-level CAL. ",
  "Sensitivity specifications replaced the 10-mm floor with alternative values ",
  "and repeated the HbA1c association, diabetes-status discrimination, and ",
  "late-life curvature analyses. These analyses were used to assess robustness ",
  "to the missing-tooth burden parameter rather than to optimise the parameter ",
  "against an external outcome.\n",
  sep = ""
)

message("Saved: outputs/tables/supp_table_missing_cal_floor_sensitivity.csv")
message("Saved: outputs/tables/supp_table_missing_cal_floor_sensitivity_delta_from_10mm.csv")
message("Saved: outputs/tables/supp_table_missing_cal_floor_sensitivity_formatted.csv")
message("Saved: outputs/tables/supp_table_missing_cal_floor_late_life_curvature.csv")
message("Saved figures to outputs/figures/")
