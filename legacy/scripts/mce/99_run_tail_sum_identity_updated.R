# scripts/99_run_tail_sum_identity.R
# ============================================================
# Verify the discrete tail-sum identity for CAL extent measures
# ============================================================
#
# For a non-negative integer-valued CAL measurement X:
#
#   E(X) = sum_{k >= 1} P(X >= k)
#
# Applied to the observed CAL sites of participant i:
#
#   mean_CAL_i = sum_k extent_i(CAL >= k) / 100
#
# where extent is expressed as a percentage of observed sites.
#
# The truncated 1-6 mm version is also verified:
#
#   mean[min(CAL, 6)]
#     = sum_{k=1}^6 extent(CAL >= k) / 100
#     = 6 * MCE_1to6 / 100,
#
# because MCE_1to6 is the arithmetic mean of the six percentage extents.
#
# Inputs:
#   data/derived/nhanes_perio_site_core_2009_2014.rds
#
# If that file is absent, the script runs:
#   scripts/01_prepare_nhanes_site_core.R
#
# Outputs:
#   outputs/tables/tail_sum_identity_summary.csv
#   outputs/tables/tail_sum_identity_weighted_check.csv
#   outputs/tables/tail_sum_identity_participant_preview.csv
#   outputs/tables/tail_sum_identity_qq_quantiles.csv
#   outputs/figures/tail_sum_qq_mean_cal_vs_mce_1to6.png
#   data/derived/nhanes_tail_sum_identity_check.rds
# ============================================================

source("scripts/00_setup.R")

options(survey.lonely.psu = "adjust")

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("data/derived", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Load the unified observed-site CAL dataset
# ------------------------------------------------------------

core_path <- "data/derived/nhanes_perio_site_core_2009_2014.rds"

if (!file.exists(core_path)) {
  message("Site-level core not found; running scripts/01_prepare_nhanes_site_core.R")
  source("scripts/01_prepare_nhanes_site_core.R")
}

if (!file.exists(core_path)) {
  stop("The site-level core was not created: ", core_path)
}

dat <- readRDS(core_path)

# ------------------------------------------------------------
# 2. Identify and clean eligible CAL site variables
# ------------------------------------------------------------

cal_vars <- grep(
  "^OHX\\d{2}LA(D|M|S|P|L|A)$",
  names(dat),
  value = TRUE
)

if (length(cal_vars) == 0) {
  stop("No site-level CAL variables were found in the site-level core.")
}

tooth_number <- as.integer(sub("^OHX(\\d{2}).*$", "\\1", cal_vars))
cal_vars <- cal_vars[!tooth_number %in% c(1, 16, 17, 32)]

cal_mat <- as.matrix(dat[, cal_vars, drop = FALSE])
storage.mode(cal_mat) <- "numeric"

cal_mat[cal_mat %in% c(99, 999)] <- NA_real_
cal_mat[cal_mat < 0 | cal_mat > 30] <- NA_real_

n_sites_observed <- rowSums(!is.na(cal_mat))

if (any(n_sites_observed == 0)) {
  stop("At least one retained participant has no observed CAL sites.")
}

# The discrete identity requires integer-valued CAL.
noninteger <- is.finite(cal_mat) & abs(cal_mat - round(cal_mat)) > 1e-10

if (any(noninteger)) {
  stop(
    "The CAL data contain non-integer values. The discrete tail-sum identity ",
    "cannot be applied exactly without defining an appropriate discretisation."
  )
}

max_cal <- max(cal_mat, na.rm = TRUE)
max_threshold <- max(6L, as.integer(max_cal))
thresholds <- seq_len(max_threshold)

message("Observed CAL range: 0 to ", max_cal, " mm")
message("Verifying thresholds 1 to ", max_threshold, " mm")

# ------------------------------------------------------------
# 3. Direct means and threshold extents
# ------------------------------------------------------------

mean_cal_direct <- rowMeans(cal_mat, na.rm = TRUE)

extent_pct_mat <- vapply(
  thresholds,
  function(k) {
    100 * rowSums(cal_mat >= k, na.rm = TRUE) / n_sites_observed
  },
  numeric(nrow(cal_mat))
)

if (is.null(dim(extent_pct_mat))) {
  extent_pct_mat <- matrix(extent_pct_mat, ncol = 1)
}

colnames(extent_pct_mat) <- paste0("extent_ge_", thresholds, "mm")

# Full discrete tail sum: exact identity for observed integer-valued CAL.
mean_cal_tail_sum <- rowSums(extent_pct_mat) / 100

# Truncated identity through 6 mm.
extent_1to6_names <- paste0("extent_ge_", 1:6, "mm")
extent_1to6_mat <- extent_pct_mat[, extent_1to6_names, drop = FALSE]

tail_sum_1to6 <- rowSums(extent_1to6_mat) / 100
mce_1to6_pct <- rowMeans(extent_1to6_mat)
mean_cal_capped_6 <- rowMeans(pmin(cal_mat, 6), na.rm = TRUE)
mce_1to6_as_capped_mean <- 6 * mce_1to6_pct / 100

participant_check <- tibble::tibble(
  SEQN = dat$SEQN,
  source_cycle = dat$source_cycle,
  n_sites_observed = n_sites_observed,
  mean_cal_direct = mean_cal_direct,
  mean_cal_tail_sum = mean_cal_tail_sum,
  full_identity_difference = mean_cal_direct - mean_cal_tail_sum,
  mean_cal_capped_6 = mean_cal_capped_6,
  tail_sum_1to6 = tail_sum_1to6,
  mce_1to6_pct = mce_1to6_pct,
  mce_1to6_as_capped_mean = mce_1to6_as_capped_mean,
  truncated_identity_difference = mean_cal_capped_6 - tail_sum_1to6,
  mce_identity_difference = mean_cal_capped_6 - mce_1to6_as_capped_mean
)

# ------------------------------------------------------------
# 4. Participant-level verification
# ------------------------------------------------------------

tolerance <- 1e-10

identity_summary <- tibble::tibble(
  check = c(
    "Mean CAL = sum of all threshold extents / 100",
    "Mean min(CAL, 6) = sum of extents 1-6 / 100",
    "Mean min(CAL, 6) = 6 x MCE_1to6 / 100"
  ),
  n_participants = nrow(participant_check),
  maximum_absolute_difference = c(
    max(abs(participant_check$full_identity_difference), na.rm = TRUE),
    max(abs(participant_check$truncated_identity_difference), na.rm = TRUE),
    max(abs(participant_check$mce_identity_difference), na.rm = TRUE)
  ),
  number_exceeding_tolerance = c(
    sum(abs(participant_check$full_identity_difference) > tolerance, na.rm = TRUE),
    sum(abs(participant_check$truncated_identity_difference) > tolerance, na.rm = TRUE),
    sum(abs(participant_check$mce_identity_difference) > tolerance, na.rm = TRUE)
  ),
  tolerance = tolerance
) |>
  dplyr::mutate(identity_verified = number_exceeding_tolerance == 0)

readr::write_csv(
  identity_summary,
  "outputs/tables/tail_sum_identity_summary.csv"
)

readr::write_csv(
  participant_check |>
    dplyr::slice_head(n = 100),
  "outputs/tables/tail_sum_identity_participant_preview.csv"
)

saveRDS(
  participant_check,
  "data/derived/nhanes_tail_sum_identity_check.rds"
)

# ------------------------------------------------------------
# 5. Survey-weighted verification
# ------------------------------------------------------------

get_first_existing <- function(data, candidates, label) {
  hit <- candidates[candidates %in% names(data)]
  if (length(hit) == 0) {
    stop("Could not find ", label, ". Tried: ", paste(candidates, collapse = ", "))
  }
  data[[hit[1]]]
}

identity_data <- dplyr::bind_cols(
  participant_check,
  tibble::as_tibble(extent_pct_mat)
) |>
  dplyr::mutate(
    survey_weight = as.numeric(get_first_existing(dat, c("WTMEC6YR", "wtmec6yr"), "6-year MEC weight")),
    survey_psu = get_first_existing(dat, c("SDMVPSU", "sdmvpsu"), "survey PSU"),
    survey_strata = get_first_existing(dat, c("SDMVSTRA", "sdmvstra"), "survey strata")
  ) |>
  dplyr::filter(
    is.finite(survey_weight),
    survey_weight > 0,
    !is.na(survey_psu),
    !is.na(survey_strata)
  )

nhanes_design <- survey::svydesign(
  ids = ~survey_psu,
  strata = ~survey_strata,
  weights = ~survey_weight,
  nest = TRUE,
  data = identity_data
)

core_means <- survey::svymean(
  ~mean_cal_direct + mean_cal_tail_sum + mean_cal_capped_6 +
    tail_sum_1to6 + mce_1to6_pct + mce_1to6_as_capped_mean,
  design = nhanes_design,
  na.rm = TRUE
)

core_coef <- stats::coef(core_means)
core_se <- sqrt(diag(stats::vcov(core_means)))

extent_formula <- stats::as.formula(
  paste("~", paste(colnames(extent_pct_mat), collapse = " + "))
)

weighted_extent_means <- survey::svymean(
  extent_formula,
  design = nhanes_design,
  na.rm = TRUE
)

weighted_extent_coef <- stats::coef(weighted_extent_means)
weighted_sum_all_extents <- sum(weighted_extent_coef) / 100
weighted_sum_1to6_extents <- sum(weighted_extent_coef[extent_1to6_names]) / 100
weighted_mce_1to6_relation <- 6 * core_coef["mce_1to6_pct"] / 100

weighted_check <- tibble::tibble(
  quantity = c(
    "Survey-weighted mean CAL, calculated directly",
    paste0("Sum of survey-weighted mean extents, thresholds 1-", max_threshold, ", divided by 100"),
    "Difference: direct mean CAL minus full tail sum",
    "Survey-weighted mean min(CAL, 6), calculated directly",
    "Sum of survey-weighted mean extents, thresholds 1-6, divided by 100",
    "Six times survey-weighted mean MCE_1to6 divided by 100",
    "Difference: capped mean CAL minus 1-6 tail sum",
    "Difference: capped mean CAL minus MCE expression"
  ),
  estimate = c(
    core_coef["mean_cal_direct"],
    weighted_sum_all_extents,
    core_coef["mean_cal_direct"] - weighted_sum_all_extents,
    core_coef["mean_cal_capped_6"],
    weighted_sum_1to6_extents,
    weighted_mce_1to6_relation,
    core_coef["mean_cal_capped_6"] - weighted_sum_1to6_extents,
    core_coef["mean_cal_capped_6"] - weighted_mce_1to6_relation
  ),
  standard_error = c(
    core_se["mean_cal_direct"],
    core_se["mean_cal_tail_sum"],
    NA_real_,
    core_se["mean_cal_capped_6"],
    core_se["tail_sum_1to6"],
    6 * core_se["mce_1to6_pct"] / 100,
    NA_real_,
    NA_real_
  )
)

readr::write_csv(
  weighted_check,
  "outputs/tables/tail_sum_identity_weighted_check.csv"
)

# ------------------------------------------------------------
# 6. Survey-weighted Q-Q comparison of mean CAL and MCE_1to6
# ------------------------------------------------------------
#
# MCE_1to6 is expressed as a percentage, whereas mean CAL is expressed in mm.
# Under the exact truncated identity:
#
#   mean[min(CAL, 6)] = 0.06 * MCE_1to6.
#
# The reference line in the Q-Q plot is therefore y = 0.06x. The capped mean
# CAL quantiles should follow this line exactly. Unrestricted mean CAL may rise
# above it in the upper tail because CAL values above 6 mm remain in the mean.

weighted_empirical_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]
  w <- w[ok]

  if (length(x) == 0) {
    stop("No valid observations were available for weighted quantiles.")
  }

  ord <- order(x)
  x <- x[ord]
  w <- w[ord]

  cumulative_weight <- cumsum(w) / sum(w)

  vapply(
    probs,
    function(p) {
      if (p <= 0) {
        return(x[1])
      }

      if (p >= 1) {
        return(x[length(x)])
      }

      x[which(cumulative_weight >= p)[1]]
    },
    numeric(1)
  )
}

qq_probs <- seq(0.01, 0.99, by = 0.01)

qq_quantiles <- tibble::tibble(
  probability = qq_probs,
  mce_1to6_pct_quantile = weighted_empirical_quantile(
    identity_data$mce_1to6_pct,
    identity_data$survey_weight,
    qq_probs
  ),
  mean_cal_direct_quantile = weighted_empirical_quantile(
    identity_data$mean_cal_direct,
    identity_data$survey_weight,
    qq_probs
  ),
  mean_cal_capped_6_quantile = weighted_empirical_quantile(
    identity_data$mean_cal_capped_6,
    identity_data$survey_weight,
    qq_probs
  )
) |>
  dplyr::mutate(
    mce_1to6_mm_equivalent_quantile = 0.06 * mce_1to6_pct_quantile,
    unrestricted_minus_mce_equivalent =
      mean_cal_direct_quantile - mce_1to6_mm_equivalent_quantile,
    capped_minus_mce_equivalent =
      mean_cal_capped_6_quantile - mce_1to6_mm_equivalent_quantile
  )

qq_table_path <- "outputs/tables/tail_sum_identity_qq_quantiles.csv"

readr::write_csv(
  qq_quantiles,
  qq_table_path
)

qq_plot_data <- qq_quantiles |>
  dplyr::select(
    probability,
    mce_1to6_pct_quantile,
    mean_cal_direct_quantile,
    mean_cal_capped_6_quantile
  ) |>
  tidyr::pivot_longer(
    cols = c(mean_cal_direct_quantile, mean_cal_capped_6_quantile),
    names_to = "cal_measure",
    values_to = "cal_quantile_mm"
  ) |>
  dplyr::mutate(
    cal_measure = dplyr::recode(
      cal_measure,
      mean_cal_direct_quantile = "Mean CAL",
      mean_cal_capped_6_quantile = "Mean CAL capped at 6 mm"
    ),
    cal_measure = factor(
      cal_measure,
      levels = c("Mean CAL capped at 6 mm", "Mean CAL")
    )
  )

qq_plot <- ggplot2::ggplot(
  qq_plot_data,
  ggplot2::aes(
    x = mce_1to6_pct_quantile,
    y = cal_quantile_mm,
    linetype = cal_measure,
    shape = cal_measure
  )
) +
  ggplot2::geom_abline(
    intercept = 0,
    slope = 0.06,
    linewidth = 0.6
  ) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 1.8, alpha = 0.85) +
  ggplot2::scale_x_continuous(
    name = expression(paste("Survey-weighted MCE"[1-6], " quantile (%)")),
    sec.axis = ggplot2::sec_axis(
      trans = ~ . * 0.06,
      name = "Equivalent capped mean CAL quantile (mm)"
    )
  ) +
  ggplot2::labs(
    title = expression(paste("Q-Q comparison of mean CAL and MCE"[1-6])),
    subtitle = "Survey-weighted quantiles; the reference line y = 0.06x is the exact capped-mean identity",
    y = "Survey-weighted mean CAL quantile (mm)",
    linetype = NULL,
    shape = NULL,
    caption = "The capped mean follows the identity line; unrestricted mean CAL can diverge where CAL exceeds 6 mm."
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "bottom"
  )

qq_figure_path <- "outputs/figures/tail_sum_qq_mean_cal_vs_mce_1to6.png"

print(qq_plot)

ggplot2::ggsave(
  filename = qq_figure_path,
  plot = qq_plot,
  width = 9,
  height = 6.5,
  dpi = 300
)

if (!file.exists(qq_table_path) || file.info(qq_table_path)$size <= 0) {
  stop("The Q-Q quantile table was not written successfully: ", qq_table_path)
}

if (!file.exists(qq_figure_path) || file.info(qq_figure_path)$size <= 0) {
  stop("The Q-Q figure was not written successfully: ", qq_figure_path)
}

message("Saved: ", qq_table_path)
message("Saved: ", qq_figure_path)

# ------------------------------------------------------------
# 7. Hard verification and console report
# ------------------------------------------------------------

stopifnot(
  all(identity_summary$identity_verified),
  abs(core_coef["mean_cal_direct"] - weighted_sum_all_extents) <= tolerance,
  abs(core_coef["mean_cal_capped_6"] - weighted_sum_1to6_extents) <= tolerance,
  abs(core_coef["mean_cal_capped_6"] - weighted_mce_1to6_relation) <= tolerance
)

message("")
message("Discrete tail-sum identity verified.")
message("")
message("Participant-level checks:")
print(identity_summary, n = Inf)
message("")
message("Survey-weighted checks:")
print(weighted_check, n = Inf)
message("")
message("Interpretation:")
message("  mean CAL = sum of all integer-threshold extents / 100")
message("  mean min(CAL, 6) = sum of extents 1-6 / 100")
message("                   = 6 x MCE_1to6 / 100")
message("")
message("Completed successfully.")
