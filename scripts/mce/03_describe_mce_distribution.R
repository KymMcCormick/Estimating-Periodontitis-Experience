# scripts/03_describe_mce_distribution.R
# Describe the distribution of MCE and comparison periodontal measures.
#
# Inputs:
#   data/processed/analytic_dataset_mce.rds
#
# Outputs:
#   outputs/tables/mce_distribution_overall.csv
#   outputs/tables/mce_distribution_by_age_group.csv
#   outputs/tables/mce_measure_correlations.csv
#   outputs/figures/mce_distribution_histogram.png
#   outputs/figures/mce_by_age_group_boxplot.png
#   outputs/figures/mce_comparator_density.png

source("scripts/00_setup.R")

# ------------------------------------------------------------
# 1. Load analytic dataset
# ------------------------------------------------------------

mce_path <- "data/processed/analytic_dataset_mce.rds"

if (!file.exists(mce_path)) {
  stop(
    "Missing analytic dataset: ", mce_path, "\n",
    "Run scripts/02_build_mce_measure.R before this script."
  )
}

analytic_mce <- readRDS(mce_path)

required_vars <- c(
  "SEQN", "Age", "mce", "mean_CAL",
  "extent_ge_3mm", "extent_ge_4mm", "extent_ge_5mm", "extent_ge_6mm",
  "WTMEC6YR", "SDMVPSU", "SDMVSTRA", "AgeGroup_Table2"
)

missing_vars <- setdiff(required_vars, names(analytic_mce))

if (length(missing_vars) > 0) {
  stop(
    "The analytic MCE dataset is missing required variable(s):\n",
    paste0("  - ", missing_vars, collapse = "\n")
  )
}

# Keep rows eligible for survey summaries.
analysis_data <- analytic_mce |>
  dplyr::filter(
    !is.na(mce),
    !is.na(WTMEC6YR),
    !is.na(SDMVPSU),
    !is.na(SDMVSTRA)
  )

message("Rows in MCE analytic dataset: ", nrow(analytic_mce))
message("Rows available for survey-weighted MCE summaries: ", nrow(analysis_data))

if (nrow(analysis_data) == 0) {
  stop("No complete rows available for survey-weighted MCE summaries.")
}

# ------------------------------------------------------------
# 2. Survey design
# ------------------------------------------------------------

nhanes_design <- survey::svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC6YR,
  nest = TRUE,
  data = analysis_data
)

# ------------------------------------------------------------
# 3. Helpers
# ------------------------------------------------------------

svy_mean_one <- function(var, design) {
  f <- stats::as.formula(paste0("~", var))
  est <- survey::svymean(f, design = design, na.rm = TRUE)

  tibble::tibble(
    measure = var,
    mean = as.numeric(stats::coef(est)[1]),
    se = as.numeric(SE(est)[1]),
    ci_low = mean - 1.96 * se,
    ci_high = mean + 1.96 * se
  )
}

unweighted_quantiles <- function(data, var) {
  x <- data[[var]]
  qs <- stats::quantile(
    x,
    probs = c(0, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 1),
    na.rm = TRUE,
    names = FALSE
  )

  tibble::tibble(
    measure = var,
    min = qs[1],
    p05 = qs[2],
    p10 = qs[3],
    p25 = qs[4],
    median = qs[5],
    p75 = qs[6],
    p90 = qs[7],
    p95 = qs[8],
    max = qs[9],
    n_non_missing = sum(!is.na(x))
  )
}

# ------------------------------------------------------------
# 4. Overall distribution
# ------------------------------------------------------------

comparison_measures <- c(
  "mce",
  "mean_CAL",
  "extent_ge_3mm",
  "extent_ge_4mm",
  "extent_ge_5mm",
  "extent_ge_6mm"
)

weighted_means <- purrr::map_dfr(
  comparison_measures,
  svy_mean_one,
  design = nhanes_design
)

quantile_summary <- purrr::map_dfr(
  comparison_measures,
  unweighted_quantiles,
  data = analysis_data
)

distribution_overall <- weighted_means |>
  dplyr::left_join(quantile_summary, by = "measure") |>
  dplyr::mutate(
    measure_label = dplyr::recode(
      measure,
      mce = "MCE: averaged extent CAL >=3-6 mm",
      mean_CAL = "Mean CAL",
      extent_ge_3mm = "Extent CAL >=3 mm",
      extent_ge_4mm = "Extent CAL >=4 mm",
      extent_ge_5mm = "Extent CAL >=5 mm",
      extent_ge_6mm = "Extent CAL >=6 mm"
    )
  ) |>
  dplyr::select(measure, measure_label, dplyr::everything())

readr::write_csv(
  distribution_overall,
  "outputs/tables/mce_distribution_overall.csv"
)

message("Saved: outputs/tables/mce_distribution_overall.csv")

# ------------------------------------------------------------
# 5. Distribution by age group
# ------------------------------------------------------------

age_group_design <- subset(nhanes_design, !is.na(AgeGroup_Table2))

age_group_mean_raw <- survey::svyby(
  ~mce + mean_CAL + extent_ge_3mm + extent_ge_4mm + extent_ge_5mm + extent_ge_6mm,
  ~AgeGroup_Table2,
  design = age_group_design,
  FUN = survey::svymean,
  na.rm = TRUE,
  vartype = c("se", "ci")
)

age_group_distribution <- age_group_mean_raw |>
  as.data.frame() |>
  tibble::as_tibble()

readr::write_csv(
  age_group_distribution,
  "outputs/tables/mce_distribution_by_age_group.csv"
)

message("Saved: outputs/tables/mce_distribution_by_age_group.csv")

# ------------------------------------------------------------
# 6. Correlations among measures
# ------------------------------------------------------------
# These are unweighted descriptive correlations. They are intended only to
# describe how the retained-dentition summaries co-vary in the analytic sample.

correlation_data <- analysis_data |>
  dplyr::select(dplyr::all_of(comparison_measures))

correlation_matrix <- stats::cor(
  correlation_data,
  use = "pairwise.complete.obs",
  method = "pearson"
)

correlation_table <- correlation_matrix |>
  as.data.frame() |>
  tibble::rownames_to_column("measure_1") |>
  tidyr::pivot_longer(
    cols = -measure_1,
    names_to = "measure_2",
    values_to = "pearson_r"
  )

readr::write_csv(
  correlation_table,
  "outputs/tables/mce_measure_correlations.csv"
)

message("Saved: outputs/tables/mce_measure_correlations.csv")

# ------------------------------------------------------------
# 7. Figures
# ------------------------------------------------------------

mce_histogram <- ggplot2::ggplot(
  analysis_data,
  ggplot2::aes(x = mce)
) +
  ggplot2::geom_histogram(bins = 40, boundary = 0, closed = "left") +
  ggplot2::scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20)
  ) +
  ggplot2::labs(
    title = "Distribution of MCE",
    subtitle = "MCE is the averaged extent of observed sites with CAL >=3, >=4, >=5, and >=6 mm",
    x = "MCE",
    y = "Number of participants"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  filename = "outputs/figures/mce_distribution_histogram.png",
  plot = mce_histogram,
  width = 8,
  height = 5,
  dpi = 300
)

mce_age_boxplot <- ggplot2::ggplot(
  analysis_data |>
    dplyr::filter(!is.na(AgeGroup_Table2)),
  ggplot2::aes(x = AgeGroup_Table2, y = mce)
) +
  ggplot2::geom_boxplot(outlier.alpha = 0.25) +
  ggplot2::scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20)
  ) +
  ggplot2::labs(
    title = "MCE by age group",
    x = "Age group",
    y = "MCE"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  filename = "outputs/figures/mce_by_age_group_boxplot.png",
  plot = mce_age_boxplot,
  width = 8,
  height = 5,
  dpi = 300
)

comparator_long <- analysis_data |>
  dplyr::select(dplyr::all_of(comparison_measures)) |>
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "measure",
    values_to = "value"
  ) |>
  dplyr::filter(!is.na(value)) |>
  dplyr::mutate(
    measure = dplyr::recode(
      measure,
      mce = "MCE",
      mean_CAL = "Mean CAL",
      extent_ge_3mm = "Extent CAL >=3 mm",
      extent_ge_4mm = "Extent CAL >=4 mm",
      extent_ge_5mm = "Extent CAL >=5 mm",
      extent_ge_6mm = "Extent CAL >=6 mm"
    )
  )

mce_comparator_density <- ggplot2::ggplot(
  comparator_long,
  ggplot2::aes(x = value)
) +
  ggplot2::geom_density(na.rm = TRUE) +
  ggplot2::facet_wrap(~measure, scales = "free") +
  ggplot2::labs(
    title = "Distribution of MCE and comparison periodontal measures",
    x = "Measure value",
    y = "Density"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  filename = "outputs/figures/mce_comparator_density.png",
  plot = mce_comparator_density,
  width = 10,
  height = 7,
  dpi = 300
)

message("Saved: outputs/figures/mce_distribution_histogram.png")
message("Saved: outputs/figures/mce_by_age_group_boxplot.png")
message("Saved: outputs/figures/mce_comparator_density.png")

# ------------------------------------------------------------
# 8. Console summary
# ------------------------------------------------------------

message("\nOverall survey-weighted means:")
print(distribution_overall |>
        dplyr::select(measure_label, mean, se, ci_low, ci_high, median, p25, p75),
      n = Inf)

message("\nCompleted MCE distribution analysis.")
