# scripts/08_mce_variant_sensitivity.R
# ============================================================
# Sensitivity analysis: alternative MCE threshold ranges
# ============================================================
#
# Purpose:
#   Test whether the primary MCE definition based on CAL thresholds 3-6 mm
#   is robust to including lower CAL thresholds such as 1 mm and 2 mm.
#
#   This script constructs alternative threshold-averaged extent measures:
#     - mce_1to6: average extent CAL >=1, >=2, >=3, >=4, >=5, >=6 mm
#     - mce_2to6: average extent CAL >=2, >=3, >=4, >=5, >=6 mm
#     - mce_3to6: average extent CAL >=3, >=4, >=5, >=6 mm (primary MCE)
#     - mce_4to6: average extent CAL >=4, >=5, >=6 mm
#
# Inputs:
#   data/processed/analytic_dataset_mce.rds
#   data/raw/OHXPER_F.xpt
#   data/raw/OHXPER_G.xpt
#   data/raw/OHXPER_H.xpt
#
# Outputs:
#   data/processed/analytic_dataset_mce_variants.rds
#   data/processed/analytic_dataset_mce_variants.csv
#
#   outputs/tables/mce_variant_overall_summary.csv
#   outputs/tables/mce_variant_by_age_group.csv
#   outputs/tables/mce_variant_correlations.csv
#   outputs/tables/mce_variant_pca_threshold_1to6_loadings.csv
#   outputs/tables/mce_variant_pca_threshold_1to6_variance.csv
#   outputs/tables/mce_variant_hba1c_coefficients.csv
#   outputs/tables/mce_variant_age_slopes.csv
#
#   outputs/figures/mce_variant_distributions.png
#   outputs/figures/mce_variant_age_trajectories.png
#   outputs/figures/mce_variant_hba1c_coefficients.png
#
# Notes:
#   - Lower thresholds may increase sensitivity but may also capture very mild
#     attachment loss, gingival recession, and measurement noise. Treat this
#     as a measurement sensitivity analysis, not automatic improvement.
# ============================================================

source("scripts/00_setup.R")

options(survey.lonely.psu = "adjust")

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Load analytic MCE dataset
# ------------------------------------------------------------

input_file <- "data/processed/analytic_dataset_mce.rds"

if (!file.exists(input_file)) {
  stop(
    "Missing input file: ", input_file, "\n",
    "Run scripts/02_build_mce_measure.R before this script."
  )
}

analytic <- readRDS(input_file)

# ------------------------------------------------------------
# 2. Helpers
# ------------------------------------------------------------

peek_header <- function(path, n = 60) {
  rawToChar(readBin(path, what = "raw", n = n))
}

stop_if_html <- function(path) {
  hdr <- peek_header(path, n = 60)
  if (grepl("^\\s*<!DOCTYPE html>|^\\s*<html", hdr, ignore.case = TRUE)) {
    stop("File is HTML, not XPT: ", path)
  }
  invisible(TRUE)
}

read_raw_xpt <- function(filename) {
  path <- file.path("data/raw", filename)
  if (!file.exists(path)) stop("Missing file: ", path)
  stop_if_html(path)
  haven::read_xpt(path)
}

cycle_from_filename <- function(filename) {
  out <- sub("^.*_([A-Z])\\.XPT$", "\\1", toupper(basename(filename)))
  if (identical(out, toupper(basename(filename)))) {
    stop("Could not extract NHANES cycle from filename: ", filename)
  }
  out
}

weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (sum(ok) == 0) return(NA_real_)
  sum(w[ok] * x[ok]) / sum(w[ok])
}

weighted_cor <- function(x, y, w) {
  ok <- is.finite(x) & is.finite(y) & is.finite(w) & w > 0
  x <- x[ok]
  y <- y[ok]
  w <- w[ok]
  if (length(x) < 3) return(NA_real_)

  mx <- sum(w * x) / sum(w)
  my <- sum(w * y) / sum(w)
  cov_xy <- sum(w * (x - mx) * (y - my)) / sum(w)
  var_x <- sum(w * (x - mx)^2) / sum(w)
  var_y <- sum(w * (y - my)^2) / sum(w)

  if (var_x <= 0 || var_y <= 0) return(NA_real_)
  cov_xy / sqrt(var_x * var_y)
}

extract_svyvar_scalar <- function(x) {
  x_num <- as.numeric(x)
  if (length(x_num) < 1 || !is.finite(x_num[1])) return(NA_real_)
  x_num[1]
}

normalise_cycle <- function(x) {
  x <- toupper(as.character(x))
  dplyr::case_when(
    grepl("_F\\.", x) | grepl("F$", x) ~ "F",
    grepl("_G\\.", x) | grepl("G$", x) ~ "G",
    grepl("_H\\.", x) | grepl("H$", x) ~ "H",
    TRUE ~ x
  )
}

calc_extent_thresholds <- function(df, thresholds = 1:6) {
  cal_cols <- names(df)[grepl("^OHX\\d{2}LA[A-Z]$", names(df))]

  if (length(cal_cols) == 0) {
    stop("No site-level CAL columns found. Expected names like OHX02LAD.")
  }

  tooth_number <- as.integer(sub("^OHX(\\d{2}).*$", "\\1", cal_cols))
  cal_cols <- cal_cols[!tooth_number %in% c(1, 16, 17, 32)]

  cal_mat <- as.matrix(df[, cal_cols, drop = FALSE])
  storage.mode(cal_mat) <- "numeric"

  # NHANES missing/invalid CAL codes and implausible values.
  cal_mat[cal_mat %in% c(99, 999)] <- NA_real_
  cal_mat[cal_mat < 0 | cal_mat > 30] <- NA_real_

  n_sites_observed <- rowSums(!is.na(cal_mat))
  mean_CAL_site <- rowMeans(cal_mat, na.rm = TRUE)
  mean_CAL_site[n_sites_observed == 0] <- NA_real_

  out <- tibble::tibble(
    SEQN = df$SEQN,
    mean_CAL_site_recalc = mean_CAL_site,
    n_sites_observed_recalc = n_sites_observed
  )

  for (thr in thresholds) {
    nm <- paste0("extent_ge_", thr, "mm")
    val <- rowSums(cal_mat >= thr, na.rm = TRUE) / n_sites_observed * 100
    val[n_sites_observed == 0] <- NA_real_
    out[[nm]] <- val
  }

  out
}

safe_row_mean <- function(df, cols) {
  out <- rowMeans(as.matrix(df[, cols, drop = FALSE]), na.rm = TRUE)
  out[is.nan(out)] <- NA_real_
  out
}

get_analysis_vars <- function(dat) {
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

  # HbA1c
  if ("hba1c" %in% names(dat)) {
    dat$hba1c_analysis <- as.numeric(dat$hba1c)
  } else if ("LBXGH" %in% names(dat)) {
    dat$hba1c_analysis <- as.numeric(dat$LBXGH)
  }

  # Survey design variables
  if ("WTMEC6YR" %in% names(dat)) {
    dat$wtmec6yr_analysis <- as.numeric(dat$WTMEC6YR)
  } else if ("wtmec6yr" %in% names(dat)) {
    dat$wtmec6yr_analysis <- as.numeric(dat$wtmec6yr)
  } else if ("WTMEC2YR" %in% names(dat)) {
    dat$wtmec6yr_analysis <- as.numeric(dat$WTMEC2YR) / 3
  } else if ("wtmec2yr" %in% names(dat)) {
    dat$wtmec6yr_analysis <- as.numeric(dat$wtmec2yr) / 3
  } else {
    stop("Could not find NHANES MEC weight.")
  }

  if ("SDMVPSU" %in% names(dat)) {
    dat$sdmvpsu_analysis <- dat$SDMVPSU
  } else if ("sdmvpsu" %in% names(dat)) {
    dat$sdmvpsu_analysis <- dat$sdmvpsu
  } else {
    stop("Could not find PSU variable.")
  }

  if ("SDMVSTRA" %in% names(dat)) {
    dat$sdmvstra_analysis <- dat$SDMVSTRA
  } else if ("sdmvstra" %in% names(dat)) {
    dat$sdmvstra_analysis <- dat$sdmvstra
  } else {
    stop("Could not find strata variable.")
  }

  if (!"n_missing_teeth" %in% names(dat)) {
    stop("Could not find n_missing_teeth. Re-run scripts/02_build_mce_measure.R.")
  }

  dat |>
    dplyr::mutate(
      n_missing_teeth_analysis = as.numeric(n_missing_teeth),
      age_c = age_analysis - 30,
      age_c_sq = age_c^2,
      AgeGroup_Table2 = if ("AgeGroup_Table2" %in% names(dat)) {
        AgeGroup_Table2
      } else {
        cut(
          age_analysis,
          breaks = c(30, 40, 50, 60, 70, 80, Inf),
          labels = c("30-39", "40-49", "50-59", "60-69", "70-79", "80+"),
          right = FALSE
        )
      }
    )
}

svy_sd <- function(formula, design) {
  sqrt(extract_svyvar_scalar(survey::svyvar(formula, design, na.rm = TRUE)))
}

extract_measure_coef <- function(fit, measure, measure_label, model_label) {
  beta <- stats::coef(fit)
  V <- stats::vcov(fit)
  est <- unname(beta["measure_z"])
  se <- sqrt(unname(V["measure_z", "measure_z"]))

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

# ------------------------------------------------------------
# 3. Recalculate CAL extents for thresholds 1-6 mm
# ------------------------------------------------------------

ohx_files <- c("OHXPER_F.xpt", "OHXPER_G.xpt", "OHXPER_H.xpt")

extent_1to6 <- purrr::map_dfr(ohx_files, function(f) {
  read_raw_xpt(f) |>
    calc_extent_thresholds(thresholds = 1:6) |>
    dplyr::mutate(source_cycle = cycle_from_filename(f))
}) |>
  dplyr::mutate(source_cycle = normalise_cycle(source_cycle))

analytic <- analytic |>
  dplyr::mutate(
    source_cycle = if ("source_cycle" %in% names(analytic)) {
      normalise_cycle(source_cycle)
    } else {
      NA_character_
    }
  )

join_vars <- if ("source_cycle" %in% names(analytic) && !all(is.na(analytic$source_cycle))) {
  c("SEQN", "source_cycle")
} else {
  "SEQN"
}

# Drop old extent variables before joining recalculated 1-6 variables, so the
# sensitivity dataset uses a single internally consistent set of extents.
extent_names <- paste0("extent_ge_", 1:6, "mm")
analytic_variants <- analytic |>
  dplyr::select(-dplyr::any_of(c(extent_names, "mean_CAL_site_recalc", "n_sites_observed_recalc"))) |>
  dplyr::left_join(extent_1to6, by = join_vars) |>
  dplyr::mutate(
    mce_1to6 = safe_row_mean(dplyr::pick(dplyr::all_of(extent_names)), extent_names),
    mce_2to6 = safe_row_mean(
      dplyr::pick(dplyr::all_of(paste0("extent_ge_", 2:6, "mm"))),
      paste0("extent_ge_", 2:6, "mm")
    ),
    mce_3to6 = safe_row_mean(
      dplyr::pick(dplyr::all_of(paste0("extent_ge_", 3:6, "mm"))),
      paste0("extent_ge_", 3:6, "mm")
    ),
    mce_4to6 = safe_row_mean(
      dplyr::pick(dplyr::all_of(paste0("extent_ge_", 4:6, "mm"))),
      paste0("extent_ge_", 4:6, "mm")
    ),
    # Keep primary name aligned with previous scripts.
    mce = mce_3to6,
    mean_CAL = dplyr::coalesce(
      if ("mean_CAL" %in% names(dplyr::pick(dplyr::everything()))) mean_CAL else NA_real_,
      mean_CAL_site_recalc
    ),
    n_sites_observed = dplyr::coalesce(
      if ("n_sites_observed" %in% names(dplyr::pick(dplyr::everything()))) n_sites_observed else NA_real_,
      n_sites_observed_recalc
    )
  ) |>
  get_analysis_vars()

# Save the extended analytic dataset.
saveRDS(
  analytic_variants,
  "data/processed/analytic_dataset_mce_variants.rds"
)

readr::write_csv(
  analytic_variants,
  "data/processed/analytic_dataset_mce_variants.csv"
)

# ------------------------------------------------------------
# 4. Survey design
# ------------------------------------------------------------

survey_data <- analytic_variants |>
  dplyr::filter(
    is.finite(wtmec6yr_analysis),
    wtmec6yr_analysis > 0,
    !is.na(sdmvpsu_analysis),
    !is.na(sdmvstra_analysis)
  )

des <- survey::svydesign(
  ids = ~sdmvpsu_analysis,
  strata = ~sdmvstra_analysis,
  weights = ~wtmec6yr_analysis,
  nest = TRUE,
  data = survey_data
)

variant_lookup <- tibble::tribble(
  ~measure,   ~measure_label,
  "mce_1to6", "MCE 1-6 mm",
  "mce_2to6", "MCE 2-6 mm",
  "mce_3to6", "MCE 3-6 mm (primary)",
  "mce_4to6", "MCE 4-6 mm"
)

# ------------------------------------------------------------
# 5. Overall and age-group summaries
# ------------------------------------------------------------

summary_list <- purrr::map_dfr(seq_len(nrow(variant_lookup)), function(i) {
  measure <- variant_lookup$measure[i]
  label <- variant_lookup$measure_label[i]
  f <- stats::as.formula(paste0("~", measure))

  tibble::tibble(
    measure = measure,
    measure_label = label,
    n_unweighted = sum(is.finite(survey_data[[measure]])),
    weighted_mean = as.numeric(stats::coef(survey::svymean(f, des, na.rm = TRUE))[1]),
    weighted_sd = svy_sd(f, des),
    unweighted_median = stats::median(survey_data[[measure]], na.rm = TRUE),
    unweighted_q25 = stats::quantile(survey_data[[measure]], 0.25, na.rm = TRUE, names = FALSE),
    unweighted_q75 = stats::quantile(survey_data[[measure]], 0.75, na.rm = TRUE, names = FALSE)
  )
})

readr::write_csv(
  summary_list,
  "outputs/tables/mce_variant_overall_summary.csv"
)

age_group_summary <- purrr::map_dfr(seq_len(nrow(variant_lookup)), function(i) {
  measure <- variant_lookup$measure[i]
  label <- variant_lookup$measure_label[i]
  f <- stats::as.formula(paste0("~", measure))

  out <- survey::svyby(
    f,
    ~AgeGroup_Table2,
    design = des,
    FUN = survey::svymean,
    na.rm = TRUE,
    vartype = c("se", "ci")
  ) |>
    as.data.frame() |>
    tibble::as_tibble()

  if (measure %in% names(out)) names(out)[names(out) == measure] <- "estimate"
  se_col <- grep("^se", names(out), value = TRUE)[1]
  if (!is.na(se_col)) names(out)[names(out) == se_col] <- "se"
  ci_l_col <- grep("^ci_l", names(out), value = TRUE)[1]
  ci_u_col <- grep("^ci_u", names(out), value = TRUE)[1]
  if (!is.na(ci_l_col)) names(out)[names(out) == ci_l_col] <- "ci_low"
  if (!is.na(ci_u_col)) names(out)[names(out) == ci_u_col] <- "ci_high"

  out |>
    dplyr::mutate(
      measure = measure,
      measure_label = label,
      ci_low = if ("ci_low" %in% names(out)) ci_low else estimate - 1.96 * se,
      ci_high = if ("ci_high" %in% names(out)) ci_high else estimate + 1.96 * se
    )
})

readr::write_csv(
  age_group_summary,
  "outputs/tables/mce_variant_by_age_group.csv"
)

# ------------------------------------------------------------
# 6. Correlations among variants and with comparators
# ------------------------------------------------------------

cor_vars <- c(
  variant_lookup$measure,
  "extent_ge_3mm", "extent_ge_4mm", "extent_ge_5mm", "extent_ge_6mm",
  "mean_CAL", "n_missing_teeth"
)
cor_vars <- unique(cor_vars[cor_vars %in% names(survey_data)])

cor_table <- tidyr::crossing(var1 = cor_vars, var2 = cor_vars) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    weighted_cor = weighted_cor(
      x = survey_data[[var1]],
      y = survey_data[[var2]],
      w = survey_data$wtmec6yr_analysis
    )
  ) |>
  dplyr::ungroup()

readr::write_csv(
  cor_table,
  "outputs/tables/mce_variant_correlations.csv"
)

# ------------------------------------------------------------
# 7. PCA of thresholds 1-6 mm
# ------------------------------------------------------------

threshold_vars <- paste0("extent_ge_", 1:6, "mm")

pca_df <- survey_data |>
  dplyr::filter(
    dplyr::if_all(dplyr::all_of(threshold_vars), is.finite),
    is.finite(wtmec6yr_analysis),
    wtmec6yr_analysis > 0
  )

weighted_cor_matrix <- matrix(
  NA_real_,
  nrow = length(threshold_vars),
  ncol = length(threshold_vars),
  dimnames = list(threshold_vars, threshold_vars)
)

for (a in threshold_vars) {
  for (b in threshold_vars) {
    weighted_cor_matrix[a, b] <- weighted_cor(
      pca_df[[a]],
      pca_df[[b]],
      pca_df$wtmec6yr_analysis
    )
  }
}

pca <- eigen(weighted_cor_matrix, symmetric = TRUE)

pca_variance <- tibble::tibble(
  component = paste0("PC", seq_along(pca$values)),
  eigenvalue = pca$values,
  proportion_variance = pca$values / sum(pca$values),
  cumulative_variance = cumsum(pca$values / sum(pca$values))
)

pca_loadings <- as.data.frame(pca$vectors) |>
  tibble::as_tibble() |>
  dplyr::mutate(threshold = threshold_vars) |>
  dplyr::rename_with(
    .cols = dplyr::starts_with("V"),
    .fn = ~paste0("PC", seq_along(.x))
  ) |>
  dplyr::select(threshold, dplyr::everything())

# Orient PC1 so loadings are positive.
if (mean(pca_loadings$PC1, na.rm = TRUE) < 0) {
  pca_loadings$PC1 <- -pca_loadings$PC1
}

readr::write_csv(
  pca_variance,
  "outputs/tables/mce_variant_pca_threshold_1to6_variance.csv"
)

readr::write_csv(
  pca_loadings,
  "outputs/tables/mce_variant_pca_threshold_1to6_loadings.csv"
)

# ------------------------------------------------------------
# 8. HbA1c alignment for MCE variants
# ------------------------------------------------------------

coef_results <- list()

for (i in seq_len(nrow(variant_lookup))) {
  measure <- variant_lookup$measure[i]
  label <- variant_lookup$measure_label[i]

  df_m <- survey_data |>
    dplyr::filter(
      is.finite(.data[[measure]]),
      is.finite(hba1c_analysis),
      is.finite(age_analysis),
      is.finite(n_missing_teeth_analysis),
      is.finite(wtmec6yr_analysis),
      wtmec6yr_analysis > 0,
      !is.na(sdmvpsu_analysis),
      !is.na(sdmvstra_analysis)
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
  measure_mean <- as.numeric(stats::coef(survey::svymean(measure_formula, des_tmp, na.rm = TRUE))[1])
  measure_sd <- svy_sd(measure_formula, des_tmp)

  missing_mean <- as.numeric(stats::coef(survey::svymean(~n_missing_teeth_analysis, des_tmp, na.rm = TRUE))[1])
  missing_sd <- svy_sd(~n_missing_teeth_analysis, des_tmp)

  df_m <- df_m |>
    dplyr::mutate(
      measure_z = (.data[[measure]] - measure_mean) / measure_sd,
      n_missing_teeth_z = (n_missing_teeth_analysis - missing_mean) / missing_sd
    )

  des_m <- survey::svydesign(
    ids = ~sdmvpsu_analysis,
    strata = ~sdmvstra_analysis,
    weights = ~wtmec6yr_analysis,
    nest = TRUE,
    data = df_m
  )

  fit_unadj <- survey::svyglm(
    hba1c_analysis ~ measure_z,
    design = des_m
  )

  fit_age <- survey::svyglm(
    hba1c_analysis ~ measure_z + age_c + age_c_sq,
    design = des_m
  )

  fit_age_missing <- survey::svyglm(
    hba1c_analysis ~ measure_z + age_c + age_c_sq + n_missing_teeth_z,
    design = des_m
  )

  coef_results[[paste0(measure, "_unadjusted")]] <- extract_measure_coef(
    fit_unadj, measure, label, "Unadjusted"
  )

  coef_results[[paste0(measure, "_age")]] <- extract_measure_coef(
    fit_age, measure, label, "Age-adjusted quadratic"
  )

  coef_results[[paste0(measure, "_age_missing")]] <- extract_measure_coef(
    fit_age_missing, measure, label, "Age + missing-teeth adjusted"
  )
}

coef_table <- dplyr::bind_rows(coef_results)

readr::write_csv(
  coef_table,
  "outputs/tables/mce_variant_hba1c_coefficients.csv"
)

# ------------------------------------------------------------
# 9. Age trajectory slopes for variants
# ------------------------------------------------------------

age_slopes <- list()
selected_ages <- c(40, 50, 60, 70, 80)

for (i in seq_len(nrow(variant_lookup))) {
  measure <- variant_lookup$measure[i]
  label <- variant_lookup$measure_label[i]

  df_m <- survey_data |>
    dplyr::filter(
      is.finite(.data[[measure]]),
      is.finite(age_c),
      is.finite(age_c_sq)
    )

  des_m <- survey::svydesign(
    ids = ~sdmvpsu_analysis,
    strata = ~sdmvstra_analysis,
    weights = ~wtmec6yr_analysis,
    nest = TRUE,
    data = df_m
  )

  fit <- survey::svyglm(
    stats::as.formula(paste0(measure, " ~ age_c + age_c_sq")),
    design = des_m
  )

  b <- stats::coef(fit)

  age_slopes[[measure]] <- tibble::tibble(
    measure = measure,
    measure_label = label,
    age = selected_ages,
    age_c = age - 30,
    annual_slope = unname(b["age_c"] + 2 * b["age_c_sq"] * age_c)
  )
}

age_slope_table <- dplyr::bind_rows(age_slopes)

readr::write_csv(
  age_slope_table,
  "outputs/tables/mce_variant_age_slopes.csv"
)

# ------------------------------------------------------------
# 10. Figures
# ------------------------------------------------------------

variant_long <- survey_data |>
  dplyr::select(SEQN, dplyr::all_of(variant_lookup$measure)) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(variant_lookup$measure),
    names_to = "measure",
    values_to = "value"
  ) |>
  dplyr::left_join(variant_lookup, by = "measure")

dist_plot <- ggplot2::ggplot(
  variant_long,
  ggplot2::aes(x = value)
) +
  ggplot2::geom_histogram(bins = 40, boundary = 0) +
  ggplot2::facet_wrap(~measure_label, scales = "free_y") +
  ggplot2::labs(
    title = "Distributions of alternative MCE threshold ranges",
    x = "Threshold-averaged extent score",
    y = "Unweighted count"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

ggplot2::ggsave(
  filename = "outputs/figures/mce_variant_distributions.png",
  plot = dist_plot,
  width = 11,
  height = 7,
  dpi = 300
)

age_plot <- ggplot2::ggplot(
  age_group_summary,
  ggplot2::aes(x = AgeGroup_Table2, y = estimate, group = measure_label)
) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = ci_low, ymax = ci_high),
    width = 0.15,
    linewidth = 0.4
  ) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2) +
  ggplot2::facet_wrap(~measure_label) +
  ggplot2::labs(
    title = "Survey-weighted MCE variants by age group",
    x = "Age group",
    y = "Mean threshold-averaged extent score"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
  )

ggplot2::ggsave(
  filename = "outputs/figures/mce_variant_age_trajectories.png",
  plot = age_plot,
  width = 11,
  height = 7,
  dpi = 300
)

coef_plot_data <- coef_table |>
  dplyr::filter(model %in% c("Age-adjusted quadratic", "Age + missing-teeth adjusted")) |>
  dplyr::mutate(
    measure_label = factor(
      measure_label,
      levels = rev(variant_lookup$measure_label)
    )
  )

coef_plot <- ggplot2::ggplot(
  coef_plot_data,
  ggplot2::aes(x = estimate, y = measure_label)
) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.4) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(xmin = conf.low, xmax = conf.high),
    height = 0.18,
    linewidth = 0.5
  ) +
  ggplot2::geom_point(size = 2) +
  ggplot2::facet_wrap(~model) +
  ggplot2::labs(
    title = "HbA1c alignment for alternative MCE threshold ranges",
    subtitle = "Expected difference in HbA1c for a one-SD higher score",
    x = "Difference in HbA1c percentage points",
    y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

ggplot2::ggsave(
  filename = "outputs/figures/mce_variant_hba1c_coefficients.png",
  plot = coef_plot,
  width = 11,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------
# 11. Console summary
# ------------------------------------------------------------

message("MCE variant sensitivity analysis complete.")
message("Saved extended dataset: data/processed/analytic_dataset_mce_variants.rds")
message("Saved tables to outputs/tables/")
message("Saved figures to outputs/figures/")

message("\nOverall weighted means:")
print(summary_list, n = Inf)

message("\nPCA variance for thresholds 1-6 mm:")
print(pca_variance, n = Inf)

message("\nHbA1c coefficients:")
print(coef_table, n = Inf)
