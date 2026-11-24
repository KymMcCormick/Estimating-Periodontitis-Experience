# scripts/06_hba1c_alignment.R
# ============================================================
# Assess HbA1c alignment for MCE and comparison periodontal measures
# ============================================================
#
# Purpose:
#   Compare how MCE, mean CAL, and single-threshold extent measures align
#   with HbA1c as an external metabolic construct.
#
# Inputs:
#   data/processed/analytic_dataset_mce.rds
#
# Outputs:
#   outputs/tables/mce_hba1c_model_coefficients.csv
#   outputs/tables/mce_hba1c_weighted_correlations.csv
#   outputs/tables/mce_hba1c_prediction_curves_age60.csv
#   outputs/tables/mce_hba1c_by_mce_quintile.csv
#
#   outputs/figures/mce_hba1c_coefficients.png
#   outputs/figures/mce_hba1c_prediction_curves_age60.png
#   outputs/figures/mce_hba1c_by_mce_quintile.png
#
# Notes:
#   - This is an associational measurement-behaviour analysis, not a causal
#     model of diabetes or periodontal disease.
#   - Coefficients are expressed as the expected difference in HbA1c for a
#     one-survey-weighted-SD higher periodontal measure.
#   - Three models are fitted:
#       1. Unadjusted
#       2. Age-adjusted quadratic
#       3. Age + missing-teeth adjusted
# ============================================================

source("scripts/00_setup.R")

# ------------------------------------------------------------
# 1. Load MCE analytic dataset
# ------------------------------------------------------------

input_file <- "data/processed/analytic_dataset_mce.rds"

if (!file.exists(input_file)) {
  stop(
    "Missing input file: ", input_file, "\n",
    "Run scripts/02_build_mce_measure.R before this script."
  )
}

dat <- readRDS(input_file)

# ------------------------------------------------------------
# 2. Standardise required variable names robustly
# ------------------------------------------------------------

# HbA1c
if ("hba1c" %in% names(dat)) {
  dat$hba1c_analysis <- as.numeric(dat$hba1c)
} else if ("LBXGH" %in% names(dat)) {
  dat$hba1c_analysis <- as.numeric(dat$LBXGH)
} else {
  stop("Could not find HbA1c variable. Expected hba1c or LBXGH.")
}

# Age
if ("age" %in% names(dat)) {
  dat$age_analysis <- as.numeric(dat$age)
} else if ("Age" %in% names(dat)) {
  dat$age_analysis <- as.numeric(dat$Age)
} else if ("RIDAGEYR" %in% names(dat)) {
  dat$age_analysis <- as.numeric(dat$RIDAGEYR)
} else {
  stop("Could not find age variable. Expected age, Age, or RIDAGEYR.")
}

# Weights
if ("WTMEC6YR" %in% names(dat)) {
  dat$wtmec6yr_analysis <- as.numeric(dat$WTMEC6YR)
} else if ("wtmec6yr" %in% names(dat)) {
  dat$wtmec6yr_analysis <- as.numeric(dat$wtmec6yr)
} else if ("WTMEC2YR" %in% names(dat)) {
  dat$wtmec6yr_analysis <- as.numeric(dat$WTMEC2YR) / 3
} else if ("wtmec2yr" %in% names(dat)) {
  dat$wtmec6yr_analysis <- as.numeric(dat$wtmec2yr) / 3
} else {
  stop("Could not find NHANES MEC weight. Expected WTMEC6YR, wtmec6yr, WTMEC2YR, or wtmec2yr.")
}

# PSU and strata
if ("SDMVPSU" %in% names(dat)) {
  dat$sdmvpsu_analysis <- dat$SDMVPSU
} else if ("sdmvpsu" %in% names(dat)) {
  dat$sdmvpsu_analysis <- dat$sdmvpsu
} else {
  stop("Could not find PSU variable. Expected SDMVPSU or sdmvpsu.")
}

if ("SDMVSTRA" %in% names(dat)) {
  dat$sdmvstra_analysis <- dat$SDMVSTRA
} else if ("sdmvstra" %in% names(dat)) {
  dat$sdmvstra_analysis <- dat$sdmvstra
} else {
  stop("Could not find strata variable. Expected SDMVSTRA or sdmvstra.")
}

# Missing teeth
if ("n_missing_teeth" %in% names(dat)) {
  dat$n_missing_teeth_analysis <- as.numeric(dat$n_missing_teeth)
} else if ("n_present_teeth" %in% names(dat)) {
  dat$n_missing_teeth_analysis <- 28 - as.numeric(dat$n_present_teeth)
} else {
  stop("Could not find tooth-count variable. Expected n_missing_teeth or n_present_teeth.")
}

# Age terms centred at 30 years, matching earlier scripts.
dat <- dat |>
  dplyr::mutate(
    age_c = age_analysis - 30,
    age_c_sq = age_c^2
  )

# ------------------------------------------------------------
# 3. Measures to compare
# ------------------------------------------------------------

measure_lookup <- tibble::tribble(
  ~measure,          ~measure_label,        ~plot_set,
  "mce",             "MCE",                 "extent_scale",
  "extent_ge_3mm",   "Extent CAL >=3 mm",   "extent_scale",
  "extent_ge_4mm",   "Extent CAL >=4 mm",   "extent_scale",
  "extent_ge_5mm",   "Extent CAL >=5 mm",   "extent_scale",
  "extent_ge_6mm",   "Extent CAL >=6 mm",   "extent_scale",
  "mean_CAL",        "Mean CAL",            "mean_cal_scale"
)

missing_measures <- measure_lookup$measure[!measure_lookup$measure %in% names(dat)]

if (length(missing_measures) > 0) {
  stop(
    "The following measure variable(s) are missing from analytic_dataset_mce.rds:\n",
    paste0("  - ", missing_measures, collapse = "\n")
  )
}

# ------------------------------------------------------------
# 4. Helpers
# ------------------------------------------------------------

extract_svyvar_scalar <- function(x) {
  x_num <- as.numeric(x)
  if (length(x_num) < 1 || !is.finite(x_num[1])) {
    return(NA_real_)
  }
  x_num[1]
}

weighted_cor <- function(x, y, w) {
  ok <- is.finite(x) & is.finite(y) & is.finite(w) & w > 0
  x <- x[ok]
  y <- y[ok]
  w <- w[ok]

  if (length(x) < 3) return(NA_real_)

  wx <- sum(w * x) / sum(w)
  wy <- sum(w * y) / sum(w)

  cov_xy <- sum(w * (x - wx) * (y - wy)) / sum(w)
  var_x <- sum(w * (x - wx)^2) / sum(w)
  var_y <- sum(w * (y - wy)^2) / sum(w)

  if (var_x <= 0 || var_y <= 0) return(NA_real_)

  cov_xy / sqrt(var_x * var_y)
}

extract_measure_coefficient <- function(fit, measure, measure_label, model_label) {
  beta <- stats::coef(fit)
  vc <- stats::vcov(fit)

  if (!"measure_z" %in% names(beta)) {
    stop("Model does not contain measure_z term.")
  }

  est <- unname(beta["measure_z"])
  se <- sqrt(unname(vc["measure_z", "measure_z"]))

  tibble::tibble(
    measure = measure,
    measure_label = measure_label,
    model = model_label,
    term = "measure_z",
    estimate = est,
    std.error = se,
    statistic = est / se,
    p.value = 2 * stats::pt(abs(est / se), df = fit$df.residual, lower.tail = FALSE),
    conf.low = est - stats::qt(0.975, df = fit$df.residual) * se,
    conf.high = est + stats::qt(0.975, df = fit$df.residual) * se
  )
}

make_prediction_curve <- function(
    fit,
    model_label,
    measure_values,
    measure_mean,
    measure_sd,
    measure,
    measure_label,
    prediction_age = 60,
    n_missing_teeth_z_value = 0
) {
  nd <- tibble::tibble(
    measure_value = measure_values,
    measure_z = (measure_values - measure_mean) / measure_sd,
    age_c = prediction_age - 30,
    age_c_sq = (prediction_age - 30)^2,
    n_missing_teeth_z = n_missing_teeth_z_value
  )

  tt <- stats::delete.response(stats::terms(fit))
  X <- stats::model.matrix(tt, data = nd)
  beta <- stats::coef(fit)
  V <- stats::vcov(fit)

  X <- X[, names(beta), drop = FALSE]
  V <- V[names(beta), names(beta), drop = FALSE]

  pred <- as.vector(X %*% beta)
  se <- sqrt(diag(X %*% V %*% t(X)))

  tibble::tibble(
    measure = measure,
    measure_label = measure_label,
    model = model_label,
    prediction_age = prediction_age,
    n_missing_teeth_z = n_missing_teeth_z_value,
    measure_value = measure_values,
    predicted_hba1c = pred,
    se = se,
    ci_low = pred - 1.96 * se,
    ci_high = pred + 1.96 * se
  )
}

# ------------------------------------------------------------
# 5. Survey-weighted correlations and HbA1c regression models
# ------------------------------------------------------------

correlation_results <- list()
coefficient_results <- list()
prediction_results <- list()

for (i in seq_len(nrow(measure_lookup))) {

  measure <- measure_lookup$measure[i]
  measure_label <- measure_lookup$measure_label[i]

  df_m <- dat |>
    dplyr::filter(
      is.finite(.data[[measure]]),
      is.finite(hba1c_analysis),
      is.finite(age_analysis),
      is.finite(n_missing_teeth_analysis),
      is.finite(wtmec6yr_analysis),
      !is.na(sdmvpsu_analysis),
      !is.na(sdmvstra_analysis),
      wtmec6yr_analysis > 0
    )

  if (nrow(df_m) < 100) {
    warning("Skipping ", measure, ": fewer than 100 complete cases.")
    next
  }

  des_tmp <- survey::svydesign(
    ids = ~sdmvpsu_analysis,
    strata = ~sdmvstra_analysis,
    weights = ~wtmec6yr_analysis,
    nest = TRUE,
    data = df_m
  )

  measure_formula <- stats::as.formula(paste0("~", measure))

  measure_mean <- as.numeric(
    stats::coef(
      survey::svymean(measure_formula, des_tmp, na.rm = TRUE)
    )[1]
  )

  measure_var <- survey::svyvar(
    measure_formula,
    des_tmp,
    na.rm = TRUE
  )

  measure_sd <- sqrt(extract_svyvar_scalar(measure_var))

  missing_mean <- as.numeric(
    stats::coef(
      survey::svymean(~n_missing_teeth_analysis, des_tmp, na.rm = TRUE)
    )[1]
  )

  missing_var <- survey::svyvar(
    ~n_missing_teeth_analysis,
    des_tmp,
    na.rm = TRUE
  )

  missing_sd <- sqrt(extract_svyvar_scalar(missing_var))

  if (!is.finite(measure_sd) || measure_sd <= 0) {
    warning("Skipping ", measure, ": measure SD is zero or not finite.")
    next
  }

  if (!is.finite(missing_sd) || missing_sd <= 0) {
    warning("Skipping ", measure, ": missing-teeth SD is zero or not finite.")
    next
  }

  df_m <- df_m |>
    dplyr::mutate(
      measure_value = .data[[measure]],
      measure_z = (measure_value - measure_mean) / measure_sd,
      n_missing_teeth_z = (n_missing_teeth_analysis - missing_mean) / missing_sd
    )

  des_m <- survey::svydesign(
    ids = ~sdmvpsu_analysis,
    strata = ~sdmvstra_analysis,
    weights = ~wtmec6yr_analysis,
    nest = TRUE,
    data = df_m
  )

  correlation_results[[measure]] <- tibble::tibble(
    measure = measure,
    measure_label = measure_label,
    n_unweighted = nrow(df_m),
    weighted_mean_measure = measure_mean,
    weighted_sd_measure = measure_sd,
    weighted_mean_missing_teeth = missing_mean,
    weighted_sd_missing_teeth = missing_sd,
    weighted_cor_hba1c = weighted_cor(
      x = df_m$measure_value,
      y = df_m$hba1c_analysis,
      w = df_m$wtmec6yr_analysis
    ),
    weighted_cor_missing_teeth = weighted_cor(
      x = df_m$measure_value,
      y = df_m$n_missing_teeth_analysis,
      w = df_m$wtmec6yr_analysis
    )
  )

  fit_unadjusted <- survey::svyglm(
    hba1c_analysis ~ measure_z,
    design = des_m
  )

  fit_age_adjusted <- survey::svyglm(
    hba1c_analysis ~ measure_z + age_c + age_c_sq,
    design = des_m
  )

  fit_age_missing_adjusted <- survey::svyglm(
    hba1c_analysis ~ measure_z + age_c + age_c_sq + n_missing_teeth_z,
    design = des_m
  )

  coefficient_results[[paste0(measure, "_unadjusted")]] <-
    extract_measure_coefficient(
      fit = fit_unadjusted,
      measure = measure,
      measure_label = measure_label,
      model_label = "Unadjusted"
    )

  coefficient_results[[paste0(measure, "_age_adjusted")]] <-
    extract_measure_coefficient(
      fit = fit_age_adjusted,
      measure = measure,
      measure_label = measure_label,
      model_label = "Age-adjusted quadratic"
    )

  coefficient_results[[paste0(measure, "_age_missing_adjusted")]] <-
    extract_measure_coefficient(
      fit = fit_age_missing_adjusted,
      measure = measure,
      measure_label = measure_label,
      model_label = "Age + missing-teeth adjusted"
    )

  # Prediction curve from the 5th to 95th percentile of the observed measure.
  q_vals <- stats::quantile(
    df_m$measure_value,
    probs = c(0.05, 0.95),
    na.rm = TRUE,
    names = FALSE,
    type = 8
  )

  if (is.finite(q_vals[1]) && is.finite(q_vals[2]) && q_vals[2] > q_vals[1]) {
    measure_values <- seq(q_vals[1], q_vals[2], length.out = 50)

    prediction_results[[paste0(measure, "_age_adjusted")]] <- make_prediction_curve(
      fit = fit_age_adjusted,
      model_label = "Age-adjusted quadratic",
      measure_values = measure_values,
      measure_mean = measure_mean,
      measure_sd = measure_sd,
      measure = measure,
      measure_label = measure_label,
      prediction_age = 60
    )

    prediction_results[[paste0(measure, "_age_missing_adjusted")]] <- make_prediction_curve(
      fit = fit_age_missing_adjusted,
      model_label = "Age + missing-teeth adjusted",
      measure_values = measure_values,
      measure_mean = measure_mean,
      measure_sd = measure_sd,
      measure = measure,
      measure_label = measure_label,
      prediction_age = 60,
      n_missing_teeth_z_value = 0
    )
  }
}

correlation_table <- dplyr::bind_rows(correlation_results)
coefficient_table <- dplyr::bind_rows(coefficient_results) |>
  dplyr::left_join(
    measure_lookup,
    by = c("measure", "measure_label")
  )
prediction_table <- dplyr::bind_rows(prediction_results) |>
  dplyr::left_join(
    measure_lookup,
    by = c("measure", "measure_label")
  )

readr::write_csv(
  correlation_table,
  "outputs/tables/mce_hba1c_weighted_correlations.csv"
)

readr::write_csv(
  coefficient_table,
  "outputs/tables/mce_hba1c_model_coefficients.csv"
)

readr::write_csv(
  prediction_table,
  "outputs/tables/mce_hba1c_prediction_curves_age60.csv"
)

# ------------------------------------------------------------
# 6. HbA1c by MCE quintile
# ------------------------------------------------------------

mce_quintile_data <- dat |>
  dplyr::filter(
    is.finite(mce),
    is.finite(hba1c_analysis),
    is.finite(n_missing_teeth_analysis),
    is.finite(wtmec6yr_analysis),
    !is.na(sdmvpsu_analysis),
    !is.na(sdmvstra_analysis),
    wtmec6yr_analysis > 0
  ) |>
  dplyr::mutate(
    mce_quintile = dplyr::ntile(mce, 5),
    mce_quintile = factor(
      mce_quintile,
      levels = 1:5,
      labels = c("Q1 lowest", "Q2", "Q3", "Q4", "Q5 highest")
    )
  )

des_quintile <- survey::svydesign(
  ids = ~sdmvpsu_analysis,
  strata = ~sdmvstra_analysis,
  weights = ~wtmec6yr_analysis,
  nest = TRUE,
  data = mce_quintile_data
)

hba1c_by_mce_quintile_raw <- survey::svyby(
  ~hba1c_analysis,
  ~mce_quintile,
  design = des_quintile,
  FUN = survey::svymean,
  na.rm = TRUE,
  vartype = c("se", "ci")
) |>
  as.data.frame() |>
  tibble::as_tibble()

# Robust column names across survey package versions.
if ("hba1c_analysis" %in% names(hba1c_by_mce_quintile_raw)) {
  names(hba1c_by_mce_quintile_raw)[names(hba1c_by_mce_quintile_raw) == "hba1c_analysis"] <- "mean_hba1c"
}

se_col <- grep("^se", names(hba1c_by_mce_quintile_raw), value = TRUE)[1]
if (!is.na(se_col)) {
  names(hba1c_by_mce_quintile_raw)[names(hba1c_by_mce_quintile_raw) == se_col] <- "se"
}

ci_l_col <- grep("^ci_l", names(hba1c_by_mce_quintile_raw), value = TRUE)[1]
ci_u_col <- grep("^ci_u", names(hba1c_by_mce_quintile_raw), value = TRUE)[1]

if (!is.na(ci_l_col)) {
  names(hba1c_by_mce_quintile_raw)[names(hba1c_by_mce_quintile_raw) == ci_l_col] <- "ci_low"
}
if (!is.na(ci_u_col)) {
  names(hba1c_by_mce_quintile_raw)[names(hba1c_by_mce_quintile_raw) == ci_u_col] <- "ci_high"
}

hba1c_by_mce_quintile <- hba1c_by_mce_quintile_raw |>
  dplyr::mutate(
    ci_low = if ("ci_low" %in% names(hba1c_by_mce_quintile_raw)) {
      ci_low
    } else {
      mean_hba1c - 1.96 * se
    },
    ci_high = if ("ci_high" %in% names(hba1c_by_mce_quintile_raw)) {
      ci_high
    } else {
      mean_hba1c + 1.96 * se
    }
  ) |>
  dplyr::left_join(
    mce_quintile_data |>
      dplyr::count(mce_quintile, name = "n_unweighted"),
    by = "mce_quintile"
  ) |>
  dplyr::left_join(
    mce_quintile_data |>
      dplyr::group_by(mce_quintile) |>
      dplyr::summarise(
        mce_min = min(mce, na.rm = TRUE),
        mce_median = stats::median(mce, na.rm = TRUE),
        mce_max = max(mce, na.rm = TRUE),
        mean_missing_teeth = mean(n_missing_teeth_analysis, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "mce_quintile"
  )

readr::write_csv(
  hba1c_by_mce_quintile,
  "outputs/tables/mce_hba1c_by_mce_quintile.csv"
)

# ------------------------------------------------------------
# 7. Figures
# ------------------------------------------------------------

coef_plot_data <- coefficient_table |>
  dplyr::filter(model %in% c("Age-adjusted quadratic", "Age + missing-teeth adjusted")) |>
  dplyr::mutate(
    model = factor(
      model,
      levels = c("Age-adjusted quadratic", "Age + missing-teeth adjusted")
    ),
    measure_label = factor(
      measure_label,
      levels = rev(measure_lookup$measure_label)
    )
  )

coef_plot <- ggplot2::ggplot(
  coef_plot_data,
  ggplot2::aes(x = estimate, y = measure_label, shape = model)
) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.4) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = conf.low, xmax = conf.high),
    height = 0.20,
    linewidth = 0.5,
    position = ggplot2::position_dodge(width = 0.45)
  ) +
  ggplot2::geom_point(
    size = 2,
    position = ggplot2::position_dodge(width = 0.45)
  ) +
  ggplot2::labs(
    title = "HbA1c alignment with and without missing-tooth adjustment",
    subtitle = "Expected difference in HbA1c for a one-SD higher periodontal measure",
    x = "Difference in HbA1c percentage points",
    y = NULL,
    shape = "Model"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  filename = "outputs/figures/mce_hba1c_coefficients.png",
  plot = coef_plot,
  width = 9,
  height = 6,
  dpi = 300
)

prediction_plot_data <- prediction_table |>
  dplyr::filter(model %in% c("Age-adjusted quadratic", "Age + missing-teeth adjusted")) |>
  dplyr::mutate(
    model = factor(
      model,
      levels = c("Age-adjusted quadratic", "Age + missing-teeth adjusted")
    )
  )

prediction_plot <- ggplot2::ggplot(
  prediction_plot_data,
  ggplot2::aes(x = measure_value, y = predicted_hba1c, linetype = model)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = ci_low, ymax = ci_high, group = model),
    alpha = 0.10
  ) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::facet_wrap(
    ~measure_label,
    scales = "free_x"
  ) +
  ggplot2::labs(
    title = "Model-predicted HbA1c by periodontal measure",
    subtitle = "Predictions shown at age 60 years; missing-teeth-adjusted model shown at mean missing teeth",
    x = "Periodontal measure value",
    y = "Predicted HbA1c",
    linetype = "Model"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  filename = "outputs/figures/mce_hba1c_prediction_curves_age60.png",
  plot = prediction_plot,
  width = 11,
  height = 7,
  dpi = 300
)

quintile_plot <- ggplot2::ggplot(
  hba1c_by_mce_quintile,
  ggplot2::aes(x = mce_quintile, y = mean_hba1c, group = 1)
) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = ci_low, ymax = ci_high),
    width = 0.15,
    linewidth = 0.5
  ) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2) +
  ggplot2::labs(
    title = "Survey-weighted mean HbA1c by MCE quintile",
    subtitle = "MCE quintiles are formed from the analytic sample",
    x = "MCE quintile",
    y = "Mean HbA1c"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
  )

ggplot2::ggsave(
  filename = "outputs/figures/mce_hba1c_by_mce_quintile.png",
  plot = quintile_plot,
  width = 9,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 8. Console summary
# ------------------------------------------------------------

message("HbA1c alignment analysis complete.")
message("Rows with non-missing HbA1c and MCE: ", nrow(mce_quintile_data))
message("Saved tables to outputs/tables/")
message("Saved figures to outputs/figures/")

message("\nCoefficients for one-SD higher measure:")
print(
  coef_plot_data |>
    dplyr::select(measure_label, model, estimate, conf.low, conf.high, p.value),
  n = Inf
)
