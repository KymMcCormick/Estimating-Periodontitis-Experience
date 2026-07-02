# functions/calc_severity_extent.R
# Conventional observed-site periodontal CAL measures.

calc_observed_site_measures <- function(df, exclude_third_molars = TRUE) {
  if (!"SEQN" %in% names(df)) {
    stop("Input data must contain SEQN.")
  }

  cal_cols <- names(df)[grepl("^OHX\\d{2}LA[A-Z]$", names(df))]

  if (length(cal_cols) == 0) {
    stop("No site-level CAL columns found. Expected names like OHX02LAD.")
  }

  tooth_number <- as.integer(sub("^OHX(\\d{2}).*$", "\\1", cal_cols))

  if (exclude_third_molars) {
    cal_cols <- cal_cols[!tooth_number %in% c(1, 16, 17, 32)]
  }

  cal_mat <- as.matrix(df[, cal_cols, drop = FALSE])
  storage.mode(cal_mat) <- "numeric"

  cal_mat[cal_mat %in% c(99, 999)] <- NA_real_
  cal_mat[cal_mat < 0 | cal_mat > 30] <- NA_real_

  n_sites_observed <- rowSums(!is.na(cal_mat))

  mean_cal_site <- rowMeans(cal_mat, na.rm = TRUE)
  mean_cal_site[n_sites_observed == 0] <- NA_real_

  extent_at <- function(threshold) {
    out <- rowSums(cal_mat >= threshold, na.rm = TRUE) / n_sites_observed * 100
    out[n_sites_observed == 0] <- NA_real_
    out
  }

  extent_ge_3mm <- extent_at(3)
  extent_ge_4mm <- extent_at(4)
  extent_ge_5mm <- extent_at(5)
  extent_ge_6mm <- extent_at(6)

  mean_extent_3to6mm <- rowMeans(
    cbind(extent_ge_3mm, extent_ge_4mm, extent_ge_5mm, extent_ge_6mm),
    na.rm = TRUE
  )
  mean_extent_3to6mm[is.nan(mean_extent_3to6mm)] <- NA_real_

  tibble::tibble(
    SEQN = df$SEQN,
    mean_CAL_site = mean_cal_site,
    extent_ge_3mm = extent_ge_3mm,
    extent_ge_4mm = extent_ge_4mm,
    extent_ge_5mm = extent_ge_5mm,
    extent_ge_6mm = extent_ge_6mm,
    mean_extent_3to6mm = mean_extent_3to6mm,
    n_sites_observed = n_sites_observed,

    # Compatibility names used by earlier scripts.
    pct_sites_ge_3mm = extent_ge_3mm,
    pct_sites_ge_4mm = extent_ge_4mm,
    pct_sites_ge_5mm = extent_ge_5mm,
    pct_sites_ge_6mm = extent_ge_6mm
  )
}

# Backwards-compatible wrapper name used by the uploaded cumulative script.
calc_severity_extent_site <- calc_observed_site_measures
