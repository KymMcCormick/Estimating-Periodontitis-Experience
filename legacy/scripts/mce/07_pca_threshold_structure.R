# scripts/07_pca_threshold_structure.R
# ============================================================
# PCA of ordered CAL extent thresholds for MCE project
# ============================================================
#
# Purpose:
#   Assess whether the ordered extent thresholds (CAL >=3, >=4, >=5, >=6 mm)
#   behave as a single underlying observable-burden dimension, and compare that
#   empirical dimension with MCE and extent CAL >=4 mm.
#
# Inputs:
#   data/processed/analytic_dataset_mce.rds
#
# Outputs:
#   outputs/tables/mce_pca_threshold_loadings.csv
#   outputs/tables/mce_pca_threshold_variance.csv
#   outputs/tables/mce_pca_score_correlations.csv
#   outputs/tables/mce_pca_hba1c_model_coefficients.csv
#   outputs/tables/mce_pca_incremental_hba1c_models.csv
#   outputs/tables/mce_pca_scores_preview.csv
#
#   data/processed/analytic_dataset_mce_with_pca.rds
#   data/processed/analytic_dataset_mce_with_pca.csv
#
#   outputs/figures/mce_pca_threshold_loadings.png
#   outputs/figures/mce_pca_score_correlations.png
#   outputs/figures/mce_pca_vs_mce.png
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

# ------------------------------------------------------------
# 2. Standardise variable names robustly
# ------------------------------------------------------------

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
} else {
  dat$hba1c_analysis <- NA_real_
  warning("Could not find HbA1c variable. HbA1c models will be skipped if all values are missing.")
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
  warning("Could not find n_missing_teeth or n_present_teeth. Missing-teeth-adjusted models will be skipped.")
  dat$n_missing_teeth_analysis <- NA_real_
}

dat <- dat |>
  dplyr::mutate(
    age_c = age_analysis - 30,
    age_c_sq = age_c^2
  )

# ------------------------------------------------------------
# 3. PCA variables
# ------------------------------------------------------------

extent_vars <- c(
  "extent_ge_3mm",
  "extent_ge_4mm",
  "extent_ge_5mm",
  "extent_ge_6mm"
)

extent_labels <- tibble::tribble(
  ~variable,         ~variable_label,
  "extent_ge_3mm",  "Extent CAL >=3 mm",
  "extent_ge_4mm",  "Extent CAL >=4 mm",
  "extent_ge_5mm",  "Extent CAL >=5 mm",
  "extent_ge_6mm",  "Extent CAL >=6 mm"
)

missing_extent_vars <- extent_vars[!extent_vars %in% names(dat)]
if (length(missing_extent_vars) > 0) {
  stop(
    "The following extent variable(s) are missing from analytic_dataset_mce.rds:\n",
    paste0("  - ", missing_extent_vars, collapse = "\n")
  )
}

comparison_vars <- c("mce", "mean_CAL", extent_vars)
comparison_vars <- comparison_vars[comparison_vars %in% names(dat)]

# ------------------------------------------------------------
# 4. Helpers
# ------------------------------------------------------------

weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (sum(ok) == 0) return(NA_real_)
  sum(w[ok] * x[ok]) / sum(w[ok])
}

weighted_var <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (sum(ok) < 2) return(NA_real_)
  x <- x[ok]
  w <- w[ok]
  mu <- sum(w * x) / sum(w)
  sum(w * (x - mu)^2) / sum(w)
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

weighted_cor_matrix <- function(df, vars, w) {
  mat <- matrix(NA_real_, nrow = length(vars), ncol = length(vars))
  dimnames(mat) <- list(vars, vars)

  for (i in seq_along(vars)) {
    for (j in seq_along(vars)) {
      mat[i, j] <- weighted_cor(df[[vars[i]]], df[[vars[j]]], w)
    }
  }

  mat
}

extract_term <- function(fit, term, model_label) {
  beta <- stats::coef(fit)
  vc <- stats::vcov(fit)

  if (!term %in% names(beta)) {
    stop("Model does not contain term: ", term)
  }

  est <- unname(beta[term])
  se <- sqrt(unname(vc[term, term]))

  tibble::tibble(
    model = model_label,
    term = term,
    estimate = est,
    std.error = se,
    statistic = est / se,
    p.value = 2 * stats::pt(abs(est / se), df = fit$df.residual, lower.tail = FALSE),
    conf.low = est - stats::qt(0.975, df = fit$df.residual) * se,
    conf.high = est + stats::qt(0.975, df = fit$df.residual) * se
  )
}

# ------------------------------------------------------------
# 5. Weighted PCA on extent thresholds
# ------------------------------------------------------------

pca_dat <- dat |>
  dplyr::filter(
    dplyr::if_all(dplyr::all_of(extent_vars), is.finite),
    is.finite(wtmec6yr_analysis),
    wtmec6yr_analysis > 0
  )

if (nrow(pca_dat) < 100) {
  stop("Fewer than 100 complete cases for threshold PCA.")
}

w <- pca_dat$wtmec6yr_analysis

extent_means <- vapply(extent_vars, function(v) weighted_mean(pca_dat[[v]], w), numeric(1))
extent_sds <- sqrt(vapply(extent_vars, function(v) weighted_var(pca_dat[[v]], w), numeric(1)))

if (any(!is.finite(extent_sds) | extent_sds <= 0)) {
  stop("At least one extent variable has zero or non-finite weighted SD.")
}

X_scaled <- sweep(as.matrix(pca_dat[, extent_vars]), 2, extent_means, FUN = "-")
X_scaled <- sweep(X_scaled, 2, extent_sds, FUN = "/")

R_w <- weighted_cor_matrix(pca_dat, extent_vars, w)

eig <- eigen(R_w, symmetric = TRUE)

loadings <- eig$vectors
rownames(loadings) <- extent_vars
colnames(loadings) <- paste0("PC", seq_len(ncol(loadings)))

# Orient PC1 so that all threshold loadings are positive on average.
if (sum(loadings[, "PC1"]) < 0) {
  loadings[, "PC1"] <- -loadings[, "PC1"]
}

scores <- X_scaled %*% loadings
colnames(scores) <- paste0("pca_extent_pc", seq_len(ncol(scores)), "_score")

# Standardise PC1 to weighted mean 0, weighted SD 1 for regression comparison.
pc1_raw <- as.numeric(scores[, "pca_extent_pc1_score"])
pc1_mean <- weighted_mean(pc1_raw, w)
pc1_sd <- sqrt(weighted_var(pc1_raw, w))
pca_dat$pca_extent_pc1_z <- (pc1_raw - pc1_mean) / pc1_sd

# Attach PC1 scores back to full analytic dataset.
score_df <- tibble::tibble(
  row_id_pca = as.integer(rownames(pca_dat)),
  pca_extent_pc1_z = pca_dat$pca_extent_pc1_z
)

dat_with_id <- dat |>
  dplyr::mutate(row_id_pca = dplyr::row_number())

dat_with_pca <- dat_with_id |>
  dplyr::left_join(score_df, by = "row_id_pca") |>
  dplyr::select(-row_id_pca)

pca_variance <- tibble::tibble(
  component = paste0("PC", seq_along(eig$values)),
  eigenvalue = eig$values,
  proportion_variance = eig$values / sum(eig$values),
  cumulative_variance = cumsum(eig$values / sum(eig$values))
)

pca_loadings <- as.data.frame(loadings) |>
  tibble::rownames_to_column("variable") |>
  tibble::as_tibble() |>
  dplyr::left_join(extent_labels, by = "variable") |>
  dplyr::relocate(variable_label, .after = variable)

# ------------------------------------------------------------
# 6. Correlations with MCE, extent >=4 mm, and comparators
# ------------------------------------------------------------

correlation_rows <- list()

for (v in comparison_vars) {
  correlation_rows[[v]] <- tibble::tibble(
    variable = v,
    weighted_correlation_with_pc1 = weighted_cor(
      dat_with_pca[[v]],
      dat_with_pca$pca_extent_pc1_z,
      dat_with_pca$wtmec6yr_analysis
    )
  )
}

pca_correlations <- dplyr::bind_rows(correlation_rows) |>
  dplyr::mutate(
    variable_label = dplyr::case_when(
      variable == "mce" ~ "MCE",
      variable == "mean_CAL" ~ "Mean CAL",
      variable == "extent_ge_3mm" ~ "Extent CAL >=3 mm",
      variable == "extent_ge_4mm" ~ "Extent CAL >=4 mm",
      variable == "extent_ge_5mm" ~ "Extent CAL >=5 mm",
      variable == "extent_ge_6mm" ~ "Extent CAL >=6 mm",
      TRUE ~ variable
    )
  ) |>
  dplyr::arrange(dplyr::desc(abs(weighted_correlation_with_pc1)))

# ------------------------------------------------------------
# 7. HbA1c alignment for PCA score and incremental comparisons
# ------------------------------------------------------------

hba1c_model_results <- tibble::tibble()
incremental_results <- tibble::tibble()

hba1c_dat <- dat_with_pca |>
  dplyr::filter(
    is.finite(pca_extent_pc1_z),
    is.finite(hba1c_analysis),
    is.finite(age_analysis),
    is.finite(wtmec6yr_analysis),
    !is.na(sdmvpsu_analysis),
    !is.na(sdmvstra_analysis),
    wtmec6yr_analysis > 0
  )

if (nrow(hba1c_dat) >= 100) {

  # Standardise MCE, extent >=4 mm, and missing teeth on the same complete-case sample.
  vars_to_z <- c("mce", "extent_ge_4mm", "n_missing_teeth_analysis")
  vars_to_z <- vars_to_z[vars_to_z %in% names(hba1c_dat)]

  for (v in vars_to_z) {
    mu_v <- weighted_mean(hba1c_dat[[v]], hba1c_dat$wtmec6yr_analysis)
    sd_v <- sqrt(weighted_var(hba1c_dat[[v]], hba1c_dat$wtmec6yr_analysis))
    hba1c_dat[[paste0(v, "_z")]] <- (hba1c_dat[[v]] - mu_v) / sd_v
  }

  des_hba1c <- survey::svydesign(
    ids = ~sdmvpsu_analysis,
    strata = ~sdmvstra_analysis,
    weights = ~wtmec6yr_analysis,
    nest = TRUE,
    data = hba1c_dat
  )

  fit_pc1_unadj <- survey::svyglm(
    hba1c_analysis ~ pca_extent_pc1_z,
    design = des_hba1c
  )

  fit_pc1_age <- survey::svyglm(
    hba1c_analysis ~ pca_extent_pc1_z + age_c + age_c_sq,
    design = des_hba1c
  )

  hba1c_model_results <- dplyr::bind_rows(
    extract_term(fit_pc1_unadj, "pca_extent_pc1_z", "PC1 unadjusted"),
    extract_term(fit_pc1_age, "pca_extent_pc1_z", "PC1 age-adjusted quadratic")
  )

  if ("n_missing_teeth_analysis_z" %in% names(hba1c_dat) &&
      all(is.finite(hba1c_dat$n_missing_teeth_analysis_z))) {

    fit_pc1_age_missing <- survey::svyglm(
      hba1c_analysis ~ pca_extent_pc1_z + age_c + age_c_sq + n_missing_teeth_analysis_z,
      design = des_hba1c
    )

    hba1c_model_results <- dplyr::bind_rows(
      hba1c_model_results,
      extract_term(fit_pc1_age_missing, "pca_extent_pc1_z", "PC1 age + missing-teeth adjusted")
    )
  }

  # Incremental test: does PC1 add anything beyond extent >=4 mm?
  # This asks whether the empirical common threshold dimension contributes
  # HbA1c information after the usual 4 mm extent threshold is already in the model.
  if (all(c("extent_ge_4mm_z", "n_missing_teeth_analysis_z") %in% names(hba1c_dat))) {

    fit_ext4_base <- survey::svyglm(
      hba1c_analysis ~ extent_ge_4mm_z + age_c + age_c_sq + n_missing_teeth_analysis_z,
      design = des_hba1c
    )

    fit_ext4_plus_pc1 <- survey::svyglm(
      hba1c_analysis ~ extent_ge_4mm_z + pca_extent_pc1_z + age_c + age_c_sq + n_missing_teeth_analysis_z,
      design = des_hba1c
    )

    incremental_results <- dplyr::bind_rows(
      extract_term(fit_ext4_base, "extent_ge_4mm_z", "Base: extent >=4 mm + age + missing teeth"),
      extract_term(fit_ext4_plus_pc1, "extent_ge_4mm_z", "Extent >=4 mm plus PC1"),
      extract_term(fit_ext4_plus_pc1, "pca_extent_pc1_z", "PC1 beyond extent >=4 mm")
    )
  }

} else {
  warning("Skipping HbA1c PCA models: fewer than 100 complete cases.")
}

# ------------------------------------------------------------
# 8. Save tables and updated dataset
# ------------------------------------------------------------

readr::write_csv(
  pca_loadings,
  "outputs/tables/mce_pca_threshold_loadings.csv"
)

readr::write_csv(
  pca_variance,
  "outputs/tables/mce_pca_threshold_variance.csv"
)

readr::write_csv(
  pca_correlations,
  "outputs/tables/mce_pca_score_correlations.csv"
)

readr::write_csv(
  hba1c_model_results,
  "outputs/tables/mce_pca_hba1c_model_coefficients.csv"
)

readr::write_csv(
  incremental_results,
  "outputs/tables/mce_pca_incremental_hba1c_models.csv"
)

readr::write_csv(
  dat_with_pca |>
    dplyr::select(
      dplyr::any_of(c(
        "SEQN", "source_cycle", "age_analysis", "mce", "mean_CAL",
        extent_vars, "n_missing_teeth", "pca_extent_pc1_z"
      ))
    ) |>
    utils::head(100),
  "outputs/tables/mce_pca_scores_preview.csv"
)

saveRDS(
  dat_with_pca,
  "data/processed/analytic_dataset_mce_with_pca.rds"
)

readr::write_csv(
  dat_with_pca,
  "data/processed/analytic_dataset_mce_with_pca.csv"
)

# ------------------------------------------------------------
# 9. Figures
# ------------------------------------------------------------

loading_plot_data <- pca_loadings |>
  dplyr::select(variable, variable_label, PC1) |>
  dplyr::mutate(
    variable_label = factor(variable_label, levels = extent_labels$variable_label)
  )

loading_plot <- ggplot2::ggplot(
  loading_plot_data,
  ggplot2::aes(x = variable_label, y = PC1)
) +
  ggplot2::geom_col(width = 0.65) +
  ggplot2::labs(
    title = "PC1 loadings for ordered CAL extent thresholds",
    subtitle = "Weighted PCA using the survey-weighted correlation matrix",
    x = NULL,
    y = "PC1 loading"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
  )

ggplot2::ggsave(
  filename = "outputs/figures/mce_pca_threshold_loadings.png",
  plot = loading_plot,
  width = 8,
  height = 6,
  dpi = 300
)

cor_plot_data <- pca_correlations |>
  dplyr::mutate(
    variable_label = factor(variable_label, levels = rev(variable_label))
  )

cor_plot <- ggplot2::ggplot(
  cor_plot_data,
  ggplot2::aes(x = weighted_correlation_with_pc1, y = variable_label)
) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.4) +
  ggplot2::geom_col(width = 0.65) +
  ggplot2::labs(
    title = "Correlation with PCA threshold score",
    subtitle = "Survey-weighted correlations with PC1 from CAL extent thresholds",
    x = "Weighted correlation with PC1",
    y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  filename = "outputs/figures/mce_pca_score_correlations.png",
  plot = cor_plot,
  width = 8,
  height = 6,
  dpi = 300
)

scatter_dat <- dat_with_pca |>
  dplyr::filter(is.finite(mce), is.finite(pca_extent_pc1_z))

pca_vs_mce_plot <- ggplot2::ggplot(
  scatter_dat,
  ggplot2::aes(x = mce, y = pca_extent_pc1_z)
) +
  ggplot2::geom_point(alpha = 0.20, size = 0.7) +
  ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = TRUE, linewidth = 0.8) +
  ggplot2::labs(
    title = "MCE compared with PCA threshold score",
    subtitle = "PC1 is the empirical common dimension among extent CAL >=3, >=4, >=5, and >=6 mm",
    x = "MCE",
    y = "PCA threshold score, PC1 z"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  filename = "outputs/figures/mce_pca_vs_mce.png",
  plot = pca_vs_mce_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 10. Console summary
# ------------------------------------------------------------

message("PCA threshold structure analysis complete.")
message("Complete cases for PCA: ", nrow(pca_dat))
message("PC1 variance explained: ", scales::percent(pca_variance$proportion_variance[1], accuracy = 0.1))
message("Saved tables to outputs/tables/")
message("Saved figures to outputs/figures/")
message("Saved dataset with PC1 score to data/processed/")

message("\nPC1 loadings:")
print(
  pca_loadings |>
    dplyr::select(variable_label, PC1),
  n = Inf
)

message("\nCorrelations with PC1:")
print(pca_correlations, n = Inf)

if (nrow(hba1c_model_results) > 0) {
  message("\nHbA1c coefficients for PC1:")
  print(hba1c_model_results, n = Inf)
}

if (nrow(incremental_results) > 0) {
  message("\nIncremental HbA1c model terms:")
  print(incremental_results, n = Inf)
}
