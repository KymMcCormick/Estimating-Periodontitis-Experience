# scripts/05_tooth_count_calibration.R
# ============================================================
# Tooth-count dependence and ordinal calibration of MCE
# ============================================================
#
# Purpose:
#   Evaluate how MCE and comparison periodontal measures behave as the
#   retained dentition contracts, and assess whether MCE categories are
#   ordinally calibrated against mean CAL within tooth-count groups.
#
# Inputs:
#   data/processed/analytic_dataset_mce.rds
#
# Outputs:
#   Tables:
#     outputs/tables/mce_tooth_count_distribution.csv
#     outputs/tables/mce_by_tooth_count_group.csv
#     outputs/tables/mce_by_exact_tooth_count.csv
#     outputs/tables/mce_tooth_count_age_adjusted_models.csv
#     outputs/tables/mce_tooth_count_age_adjusted_predictions.csv
#     outputs/tables/mce_ordinal_calibration_mean_cal.csv
#     outputs/tables/mce_ordinal_calibration_by_tooth_group.csv
#     outputs/tables/mce_ordinal_calibration_interaction_model.csv
#     outputs/tables/mce_ordinal_calibration_interaction_test.csv
#
#   Figures:
#     outputs/figures/mce_tooth_count_distribution.png
#     outputs/figures/mce_by_tooth_count_group.png
#     outputs/figures/mce_comparator_by_tooth_count_group.png
#     outputs/figures/mce_age_adjusted_by_tooth_count.png
#     outputs/figures/mce_mean_cal_by_mce_category.png
#     outputs/figures/mce_ordinal_calibration_by_tooth_group.png
#
# Notes:
#   - MCE is a retained-dentition measure. These analyses are intended to
#     show how its interpretation depends on observable periodontal support.
#   - Age-adjusted tooth-count models are descriptive measurement-behaviour
#     models, not causal models of tooth retention or periodontal disease.
# ============================================================

source("scripts/00_setup.R")

# ------------------------------------------------------------
# 1. Load analytic dataset
# ------------------------------------------------------------

input_file <- "data/processed/analytic_dataset_mce.rds"

if (!file.exists(input_file)) {
  stop(
    "Missing input file: ", input_file, "\n",
    "Run scripts/02_build_mce_measure.R before this script."
  )
}

dat <- readRDS(input_file)

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Standardise required variables
# ------------------------------------------------------------

if (!"mce" %in% names(dat)) {
  stop("Expected variable 'mce' is missing from analytic_dataset_mce.rds.")
}

if (!"mean_CAL" %in% names(dat)) {
  stop("Expected variable 'mean_CAL' is missing from analytic_dataset_mce.rds.")
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

# Survey weights
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

# Retained tooth count
if ("n_present_teeth" %in% names(dat)) {
  dat$n_present_teeth_analysis <- as.numeric(dat$n_present_teeth)
} else if ("approx_teeth_observed" %in% names(dat)) {
  dat$n_present_teeth_analysis <- as.numeric(dat$approx_teeth_observed)
} else if ("n_sites_observed" %in% names(dat)) {
  dat$n_present_teeth_analysis <- as.numeric(dat$n_sites_observed) / 6
} else {
  stop("Could not find retained tooth count. Expected n_present_teeth, approx_teeth_observed, or n_sites_observed.")
}

# Missing tooth count. For a 28-tooth modelled dentition, this is deterministic
# from retained tooth count, but we keep it explicitly for interpretation.
if ("n_missing_teeth" %in% names(dat)) {
  dat$n_missing_teeth_analysis <- as.numeric(dat$n_missing_teeth)
} else {
  dat$n_missing_teeth_analysis <- 28 - dat$n_present_teeth_analysis
}

# Age terms centred at 30 years.
dat <- dat |>
  dplyr::mutate(
    age_c = age_analysis - 30,
    age_c_sq = age_c^2,
    n_present_teeth_analysis = pmin(pmax(n_present_teeth_analysis, 0), 28),
    n_missing_teeth_analysis = pmin(pmax(n_missing_teeth_analysis, 0), 28),

    # Coarse tooth-count groups used in the original ordinal calibration.
    tooth_group = dplyr::case_when(
      n_present_teeth_analysis > 25 ~ ">25 teeth",
      n_present_teeth_analysis >= 15 & n_present_teeth_analysis <= 25 ~ "15–25 teeth",
      n_present_teeth_analysis < 15 ~ "<15 teeth",
      TRUE ~ NA_character_
    ),
    tooth_group = factor(
      tooth_group,
      levels = c(">25 teeth", "15–25 teeth", "<15 teeth")
    ),

    # Finer tooth-count groups used for descriptive compression plots.
    tooth_count_group = dplyr::case_when(
      n_present_teeth_analysis >= 1 & n_present_teeth_analysis <= 9 ~ "1–9",
      n_present_teeth_analysis >= 10 & n_present_teeth_analysis <= 19 ~ "10–19",
      n_present_teeth_analysis >= 20 & n_present_teeth_analysis <= 24 ~ "20–24",
      n_present_teeth_analysis >= 25 & n_present_teeth_analysis <= 27 ~ "25–27",
      n_present_teeth_analysis == 28 ~ "28",
      TRUE ~ NA_character_
    ),
    tooth_count_group = factor(
      tooth_count_group,
      levels = c("1–9", "10–19", "20–24", "25–27", "28")
    ),

    # Ordered MCE categories used for ordinal calibration.
    mce_category = dplyr::case_when(
      is.na(mce) ~ NA_character_,
      mce == 0 ~ "0",
      mce > 0 & mce <= 25 ~ ">0–25",
      mce > 25 & mce <= 50 ~ ">25–50",
      mce > 50 & mce <= 75 ~ ">50–75",
      mce > 75 & mce <= 100 ~ ">75–100",
      TRUE ~ NA_character_
    ),
    mce_category = factor(
      mce_category,
      levels = c("0", ">0–25", ">25–50", ">50–75", ">75–100"),
      ordered = TRUE
    )
  )

# ------------------------------------------------------------
# 3. Analysis dataset and survey design
# ------------------------------------------------------------

analysis_dat <- dat |>
  dplyr::filter(
    is.finite(age_analysis),
    is.finite(wtmec6yr_analysis),
    wtmec6yr_analysis > 0,
    !is.na(sdmvpsu_analysis),
    !is.na(sdmvstra_analysis),
    is.finite(n_present_teeth_analysis),
    is.finite(mce),
    is.finite(mean_CAL)
  )

des <- survey::svydesign(
  ids = ~sdmvpsu_analysis,
  strata = ~sdmvstra_analysis,
  weights = ~wtmec6yr_analysis,
  nest = TRUE,
  data = analysis_dat
)

# ------------------------------------------------------------
# 4. Helpers
# ------------------------------------------------------------

normalise_svyby_mean <- function(x, var, group_vars) {
  out <- as.data.frame(x) |>
    tibble::as_tibble()

  if (var %in% names(out)) {
    names(out)[names(out) == var] <- "mean"
  }

  se_col <- grep("^se", names(out), value = TRUE)[1]
  if (!is.na(se_col)) {
    names(out)[names(out) == se_col] <- "se"
  }

  ci_l_col <- grep("^ci_l", names(out), value = TRUE)[1]
  ci_u_col <- grep("^ci_u", names(out), value = TRUE)[1]

  if (!is.na(ci_l_col)) {
    names(out)[names(out) == ci_l_col] <- "ci_low"
  }
  if (!is.na(ci_u_col)) {
    names(out)[names(out) == ci_u_col] <- "ci_high"
  }

  out |>
    dplyr::mutate(
      variable = var,
      ci_low = if ("ci_low" %in% names(out)) {
        ci_low
      } else {
        mean - 1.96 * se
      },
      ci_high = if ("ci_high" %in% names(out)) {
        ci_high
      } else {
        mean + 1.96 * se
      }
    ) |>
    dplyr::select(dplyr::all_of(group_vars), variable, mean, se, ci_low, ci_high, dplyr::everything())
}

svy_mean_by <- function(var, group, design) {
  f <- stats::as.formula(paste0("~", var))
  by <- stats::as.formula(paste0("~", group))

  raw <- survey::svyby(
    f,
    by,
    design,
    FUN = survey::svymean,
    na.rm = TRUE,
    vartype = c("se", "ci"),
    keep.names = FALSE
  )

  normalise_svyby_mean(raw, var = var, group_vars = group)
}

svy_mean_by_two_groups <- function(var, group1, group2, design) {
  f <- stats::as.formula(paste0("~", var))
  by <- stats::as.formula(paste0("~", group1, "+", group2))

  raw <- survey::svyby(
    f,
    by,
    design,
    FUN = survey::svymean,
    na.rm = TRUE,
    vartype = c("se", "ci"),
    keep.names = FALSE
  )

  normalise_svyby_mean(raw, var = var, group_vars = c(group1, group2))
}

extract_tooth_model_row <- function(fit, measure, measure_label) {
  beta <- stats::coef(fit)
  vc <- stats::vcov(fit)

  term <- "n_present_teeth_analysis"

  est <- unname(beta[term])
  se <- sqrt(unname(vc[term, term]))
  t_value <- est / se

  tibble::tibble(
    measure = measure,
    measure_label = measure_label,
    term = term,
    estimate = est,
    std.error = se,
    statistic = t_value,
    p.value = 2 * stats::pt(abs(t_value), df = fit$df.residual, lower.tail = FALSE),
    conf.low = est - stats::qt(0.975, df = fit$df.residual) * se,
    conf.high = est + stats::qt(0.975, df = fit$df.residual) * se
  )
}

make_tooth_prediction_curve <- function(fit, measure, measure_label, at_age = 60) {
  nd <- tibble::tibble(
    n_present_teeth_analysis = seq(1, 28, by = 1),
    age_c = at_age - 30,
    age_c_sq = (at_age - 30)^2
  )

  pr <- stats::predict(fit, newdata = nd, se.fit = TRUE)
  fit_values <- as.numeric(pr)
  se_values <- sqrt(as.numeric(attr(pr, "var")))

  tibble::tibble(
    measure = measure,
    measure_label = measure_label,
    prediction_age = at_age,
    n_present_teeth = nd$n_present_teeth_analysis,
    predicted = fit_values,
    se = se_values,
    ci_low = fit_values - 1.96 * se_values,
    ci_high = fit_values + 1.96 * se_values,

    # Plot-only truncation for bounded percentage measures. The model-implied
    # values are retained above in predicted/ci_low/ci_high.
    predicted_plot = dplyr::if_else(measure == "mean_CAL", predicted, pmax(predicted, 0)),
    ci_low_plot = dplyr::if_else(measure == "mean_CAL", ci_low, pmax(ci_low, 0)),
    ci_high_plot = dplyr::if_else(measure == "mean_CAL", ci_high, pmax(ci_high, 0))
  )
}

# ------------------------------------------------------------
# 5. Tooth-count distribution
# ------------------------------------------------------------

tooth_count_distribution <- analysis_dat |>
  dplyr::count(n_present_teeth_analysis, name = "n_unweighted") |>
  dplyr::left_join(
    analysis_dat |>
      dplyr::group_by(n_present_teeth_analysis) |>
      dplyr::summarise(
        weighted_total = sum(wtmec6yr_analysis, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "n_present_teeth_analysis"
  ) |>
  dplyr::mutate(
    weighted_percent = 100 * weighted_total / sum(weighted_total, na.rm = TRUE)
  ) |>
  dplyr::rename(n_present_teeth = n_present_teeth_analysis)

readr::write_csv(
  tooth_count_distribution,
  "outputs/tables/mce_tooth_count_distribution.csv"
)

# ------------------------------------------------------------
# 6. MCE and comparator measures by tooth-count group
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

missing_measures <- measure_lookup$measure[!measure_lookup$measure %in% names(analysis_dat)]
if (length(missing_measures) > 0) {
  stop(
    "The following measure variable(s) are missing from analytic_dataset_mce.rds:\n",
    paste0("  - ", missing_measures, collapse = "\n")
  )
}

by_tooth_group <- purrr::map_dfr(measure_lookup$measure, function(v) {
  svy_mean_by(v, "tooth_count_group", des)
}) |>
  dplyr::left_join(measure_lookup, by = c("variable" = "measure")) |>
  dplyr::rename(measure = variable)

by_tooth_group_counts <- analysis_dat |>
  dplyr::filter(!is.na(tooth_count_group)) |>
  dplyr::count(tooth_count_group, name = "n_unweighted") |>
  dplyr::left_join(
    analysis_dat |>
      dplyr::filter(!is.na(tooth_count_group)) |>
      dplyr::group_by(tooth_count_group) |>
      dplyr::summarise(
        weighted_total = sum(wtmec6yr_analysis, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "tooth_count_group"
  ) |>
  dplyr::mutate(
    weighted_percent = 100 * weighted_total / sum(weighted_total, na.rm = TRUE)
  )

by_tooth_group <- by_tooth_group |>
  dplyr::left_join(by_tooth_group_counts, by = "tooth_count_group")

readr::write_csv(
  by_tooth_group,
  "outputs/tables/mce_by_tooth_count_group.csv"
)

by_exact_tooth <- purrr::map_dfr(measure_lookup$measure, function(v) {
  svy_mean_by(v, "n_present_teeth_analysis", des)
}) |>
  dplyr::left_join(measure_lookup, by = c("variable" = "measure")) |>
  dplyr::rename(
    measure = variable,
    n_present_teeth = n_present_teeth_analysis
  ) |>
  dplyr::left_join(tooth_count_distribution, by = "n_present_teeth")

readr::write_csv(
  by_exact_tooth,
  "outputs/tables/mce_by_exact_tooth_count.csv"
)

# ------------------------------------------------------------
# 7. Age-adjusted retained-tooth-count models
# ------------------------------------------------------------

tooth_model_results <- list()
tooth_prediction_results <- list()

for (i in seq_len(nrow(measure_lookup))) {
  measure <- measure_lookup$measure[i]
  measure_label <- measure_lookup$measure_label[i]

  df_m <- analysis_dat |>
    dplyr::filter(is.finite(.data[[measure]]))

  des_m <- survey::svydesign(
    ids = ~sdmvpsu_analysis,
    strata = ~sdmvstra_analysis,
    weights = ~wtmec6yr_analysis,
    nest = TRUE,
    data = df_m
  )

  f <- stats::as.formula(
    paste0(measure, " ~ n_present_teeth_analysis + age_c + age_c_sq")
  )

  fit <- survey::svyglm(f, design = des_m)

  tooth_model_results[[measure]] <- extract_tooth_model_row(
    fit = fit,
    measure = measure,
    measure_label = measure_label
  )

  tooth_prediction_results[[measure]] <- make_tooth_prediction_curve(
    fit = fit,
    measure = measure,
    measure_label = measure_label,
    at_age = 60
  )
}

tooth_model_table <- dplyr::bind_rows(tooth_model_results) |>
  dplyr::left_join(measure_lookup, by = c("measure", "measure_label"))

tooth_prediction_table <- dplyr::bind_rows(tooth_prediction_results) |>
  dplyr::left_join(measure_lookup, by = c("measure", "measure_label"))

readr::write_csv(
  tooth_model_table,
  "outputs/tables/mce_tooth_count_age_adjusted_models.csv"
)

readr::write_csv(
  tooth_prediction_table,
  "outputs/tables/mce_tooth_count_age_adjusted_predictions.csv"
)

# ------------------------------------------------------------
# 8. Ordinal calibration of MCE against mean CAL
# ------------------------------------------------------------

calibration_dat <- analysis_dat |>
  dplyr::filter(
    !is.na(mce_category),
    !is.na(tooth_group),
    is.finite(mean_CAL)
  ) |>
  droplevels()

calibration_des <- survey::svydesign(
  ids = ~sdmvpsu_analysis,
  strata = ~sdmvstra_analysis,
  weights = ~wtmec6yr_analysis,
  nest = TRUE,
  data = calibration_dat
)

calibration_overall <- svy_mean_by("mean_CAL", "mce_category", calibration_des) |>
  dplyr::rename(mean_CAL = mean, mean_CAL_se = se, mean_CAL_ci_low = ci_low, mean_CAL_ci_high = ci_high) |>
  dplyr::left_join(
    calibration_dat |>
      dplyr::count(mce_category, name = "n_unweighted"),
    by = "mce_category"
  ) |>
  dplyr::left_join(
    calibration_dat |>
      dplyr::group_by(mce_category) |>
      dplyr::summarise(
        weighted_total = sum(wtmec6yr_analysis, na.rm = TRUE),
        mce_min = min(mce, na.rm = TRUE),
        mce_median = stats::median(mce, na.rm = TRUE),
        mce_max = max(mce, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "mce_category"
  )

readr::write_csv(
  calibration_overall,
  "outputs/tables/mce_ordinal_calibration_mean_cal.csv"
)

calibration_by_tooth_group <- svy_mean_by_two_groups(
  "mean_CAL",
  "mce_category",
  "tooth_group",
  calibration_des
) |>
  dplyr::rename(mean_CAL = mean, mean_CAL_se = se, mean_CAL_ci_low = ci_low, mean_CAL_ci_high = ci_high) |>
  dplyr::left_join(
    calibration_dat |>
      dplyr::count(mce_category, tooth_group, name = "n_unweighted"),
    by = c("mce_category", "tooth_group")
  ) |>
  dplyr::left_join(
    calibration_dat |>
      dplyr::group_by(mce_category, tooth_group) |>
      dplyr::summarise(
        weighted_total = sum(wtmec6yr_analysis, na.rm = TRUE),
        mce_min = min(mce, na.rm = TRUE),
        mce_median = stats::median(mce, na.rm = TRUE),
        mce_max = max(mce, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("mce_category", "tooth_group")
  )

readr::write_csv(
  calibration_by_tooth_group,
  "outputs/tables/mce_ordinal_calibration_by_tooth_group.csv"
)

# Interaction model: is calibration of MCE category against mean CAL invariant
# across tooth-count groups?
calibration_model_table <- tibble::tibble()
calibration_interaction_test <- tibble::tibble()

try({
  fit_main <- survey::svyglm(
    mean_CAL ~ mce_category + tooth_group,
    design = calibration_des
  )

  fit_interaction <- survey::svyglm(
    mean_CAL ~ mce_category * tooth_group,
    design = calibration_des
  )

  calibration_model_table <- broom::tidy(fit_interaction, conf.int = TRUE)

  lrt <- survey::regTermTest(
    fit_interaction,
    ~mce_category:tooth_group
  )

  calibration_interaction_test <- tibble::tibble(
    test = "MCE category x tooth-count group interaction",
    F_statistic = as.numeric(lrt$Ftest),
    df_num = as.numeric(lrt$df),
    df_denom = as.numeric(lrt$ddf),
    p.value = as.numeric(lrt$p)
  )
}, silent = TRUE)

readr::write_csv(
  calibration_model_table,
  "outputs/tables/mce_ordinal_calibration_interaction_model.csv"
)

readr::write_csv(
  calibration_interaction_test,
  "outputs/tables/mce_ordinal_calibration_interaction_test.csv"
)

# ------------------------------------------------------------
# 9. Figures
# ------------------------------------------------------------

p_tooth_distribution <- ggplot2::ggplot(
  tooth_count_distribution,
  ggplot2::aes(x = n_present_teeth, y = n_unweighted)
) +
  ggplot2::geom_col() +
  ggplot2::scale_x_continuous(breaks = seq(0, 28, by = 4)) +
  ggplot2::labs(
    title = "Distribution of retained tooth count",
    subtitle = "Modelled permanent teeth, third molars excluded",
    x = "Number of retained teeth",
    y = "Number of participants"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  filename = "outputs/figures/mce_tooth_count_distribution.png",
  plot = p_tooth_distribution,
  width = 9,
  height = 6,
  dpi = 300
)

p_mce_tooth_group <- by_tooth_group |>
  dplyr::filter(measure == "mce") |>
  ggplot2::ggplot(
    ggplot2::aes(x = tooth_count_group, y = mean, group = 1)
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = ci_low, ymax = ci_high),
    width = 0.15,
    linewidth = 0.5
  ) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2) +
  ggplot2::labs(
    title = "MCE by retained tooth count",
    subtitle = "Survey-weighted mean with 95% confidence interval",
    x = "Number of retained teeth",
    y = "MCE"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  filename = "outputs/figures/mce_by_tooth_count_group.png",
  plot = p_mce_tooth_group,
  width = 9,
  height = 6,
  dpi = 300
)

p_comparator_tooth_group <- by_tooth_group |>
  dplyr::filter(plot_set == "extent_scale") |>
  dplyr::mutate(
    measure_label = factor(
      measure_label,
      levels = c(
        "Extent CAL >=3 mm",
        "Extent CAL >=4 mm",
        "MCE",
        "Extent CAL >=5 mm",
        "Extent CAL >=6 mm"
      )
    )
  ) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = tooth_count_group,
      y = mean,
      group = measure_label,
      linetype = measure_label
    )
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = ci_low, ymax = ci_high),
    width = 0.15,
    linewidth = 0.35,
    alpha = 0.75
  ) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 1.8) +
  ggplot2::labs(
    title = "MCE and extent measures by retained tooth count",
    subtitle = "Survey-weighted means with 95% confidence intervals",
    x = "Number of retained teeth",
    y = "Measure value",
    linetype = "Measure"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  filename = "outputs/figures/mce_comparator_by_tooth_count_group.png",
  plot = p_comparator_tooth_group,
  width = 11,
  height = 7,
  dpi = 300
)

p_age_adjusted <- tooth_prediction_table |>
  dplyr::filter(plot_set == "extent_scale") |>
  dplyr::mutate(
    measure_label = factor(
      measure_label,
      levels = c(
        "Extent CAL >=3 mm",
        "Extent CAL >=4 mm",
        "MCE",
        "Extent CAL >=5 mm",
        "Extent CAL >=6 mm"
      )
    )
  ) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = n_present_teeth,
      y = predicted_plot,
      group = measure_label,
      linetype = measure_label
    )
  ) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = ci_low_plot, ymax = ci_high_plot),
    alpha = 0.12,
    linetype = 0
  ) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::scale_x_continuous(breaks = seq(4, 28, by = 4)) +
  ggplot2::labs(
    title = "Age-adjusted association with retained tooth count",
    subtitle = "Survey-weighted quadratic age-adjusted models; predictions shown at age 60 years",
    x = "Number of retained teeth",
    y = "Predicted measure value",
    linetype = "Measure"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  filename = "outputs/figures/mce_age_adjusted_by_tooth_count.png",
  plot = p_age_adjusted,
  width = 11,
  height = 7,
  dpi = 300
)

p_calibration_overall <- calibration_overall |>
  ggplot2::ggplot(
    ggplot2::aes(x = mce_category, y = mean_CAL, group = 1)
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = mean_CAL_ci_low, ymax = mean_CAL_ci_high),
    width = 0.15,
    linewidth = 0.5
  ) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2) +
  ggplot2::labs(
    title = "Ordinal calibration of MCE against mean CAL",
    subtitle = "Survey-weighted mean CAL by ordered MCE category",
    x = "MCE category",
    y = "Mean CAL, mm"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  filename = "outputs/figures/mce_mean_cal_by_mce_category.png",
  plot = p_calibration_overall,
  width = 9,
  height = 6,
  dpi = 300
)

p_calibration_tooth_group <- calibration_by_tooth_group |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = mce_category,
      y = mean_CAL,
      group = tooth_group,
      linetype = tooth_group
    )
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = mean_CAL_ci_low, ymax = mean_CAL_ci_high),
    width = 0.15,
    linewidth = 0.45,
    alpha = 0.85
  ) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2) +
  ggplot2::labs(
    title = "Ordinal calibration of MCE by retained tooth count",
    subtitle = "Survey-weighted mean CAL within ordered MCE categories and tooth-count groups",
    x = "MCE category",
    y = "Mean CAL, mm",
    linetype = "Observed tooth-count group"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  filename = "outputs/figures/mce_ordinal_calibration_by_tooth_group.png",
  plot = p_calibration_tooth_group,
  width = 10,
  height = 6.5,
  dpi = 300
)

# ------------------------------------------------------------
# 10. Console summary
# ------------------------------------------------------------

message("Tooth-count dependence and ordinal calibration analysis complete.")
message("Rows analysed: ", nrow(analysis_dat))
message("Saved tables to outputs/tables/")
message("Saved figures to outputs/figures/")

message("\nAge-adjusted tooth-count coefficients:")
print(
  tooth_model_table |>
    dplyr::select(measure_label, estimate, conf.low, conf.high, p.value),
  n = Inf
)

message("\nOrdinal calibration interaction test:")
print(calibration_interaction_test, n = Inf)
