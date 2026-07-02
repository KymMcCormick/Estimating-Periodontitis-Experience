# R/make_max_cal_by_tooth.R
# Construct tooth-level maximum CAL from NHANES site-level periodontal CAL columns.

make_max_cal_by_tooth <- function(df, exclude_third_molars = TRUE) {
  if (!"SEQN" %in% names(df)) {
    stop("Input data must contain SEQN.")
  }

  cal_cols <- names(df)[grepl("^OHX\\d{2}LA[A-Z]$", names(df))]

  if (length(cal_cols) == 0) {
    stop("No site-level CAL columns found. Expected names like OHX02LAD.")
  }

  tooth_numbers <- sort(unique(as.integer(sub("^OHX(\\d{2}).*$", "\\1", cal_cols))))

  if (exclude_third_molars) {
    tooth_numbers <- setdiff(tooth_numbers, c(1, 16, 17, 32))
  }

  out <- tibble::tibble(SEQN = df$SEQN)

  for (tooth in tooth_numbers) {
    tooth_prefix <- sprintf("OHX%02dLA", tooth)
    tooth_cols <- cal_cols[startsWith(cal_cols, tooth_prefix)]

    cal_mat <- as.matrix(df[, tooth_cols, drop = FALSE])
    storage.mode(cal_mat) <- "numeric"

    # NHANES periodontal missing / invalid codes should not be treated as CAL.
    cal_mat[cal_mat %in% c(99, 999)] <- NA_real_
    cal_mat[cal_mat < 0 | cal_mat > 30] <- NA_real_

    max_cal <- apply(cal_mat, 1, function(x) {
      if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
    })

    out[[sprintf("CAL_%02d", tooth)]] <- max_cal
  }

  out
}
