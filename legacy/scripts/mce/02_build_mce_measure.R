# scripts/02_build_mce_measure.R
# Construct the Multithreshold CAL Extent (MCE) analytic dataset.
#
# Purpose:
#   Create a retained-dentition periodontal burden dataset for the MCE project.
#   MCE is defined as the average of ordered single-threshold extent measures:
#
#     MCE_site = mean(extent of observed CAL sites >= 3, 4, 5, and 6 mm)
#
#   This script is intentionally parallel to scripts/02_build_epe_iteration01.R,
#   but removes the EPE-specific missing-tooth attribution model. Missing teeth
#   are counted for observability/tooth-count analyses, but they are not assigned
#   expected CAL burden in MCE.
#
# Inputs:
#   data/raw/OHXPER_F.xpt
#   data/raw/OHXPER_G.xpt
#   data/raw/OHXPER_H.xpt
#   data/processed/nhanes_participant_core.rds      preferred, if available
#     OR
#   data/processed/nhanes_demographics_hba1c.rds    fallback from earlier project
#
# Outputs:
#   data/processed/analytic_dataset_mce.rds
#   data/processed/analytic_dataset_mce.csv
#   outputs/tables/mce_measure_summary.csv
#   outputs/tables/mce_age_group_counts.csv

source("scripts/00_setup.R")

# ------------------------------------------------------------
# 0. Directories and settings
# ------------------------------------------------------------

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

options(survey.lonely.psu = "adjust")


# ------------------------------------------------------------
# 1. Helpers
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

normalise_cycle <- function(x) {
  x_chr <- toupper(as.character(x))

  dplyr::case_when(
    grepl("_F\\.", x_chr) | x_chr == "F" | grepl("F$", x_chr) ~ "F",
    grepl("_G\\.", x_chr) | x_chr == "G" | grepl("G$", x_chr) ~ "G",
    grepl("_H\\.", x_chr) | x_chr == "H" | grepl("H$", x_chr) ~ "H",
    TRUE ~ NA_character_
  )
}

na_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

row_max_na <- function(mat) {
  apply(mat, 1, function(x) {
    if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  })
}

first_existing <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) {
    stop(
      "None of the expected input files exists:\n",
      paste0("  - ", paths, collapse = "\n")
    )
  }
  existing[1]
}

standardise_participant_core <- function(df) {
  if (!"SEQN" %in% names(df)) stop("Participant core is missing SEQN.")
  if (!"source_cycle" %in% names(df)) stop("Participant core is missing source_cycle.")

  df <- df |>
    dplyr::mutate(source_cycle = normalise_cycle(source_cycle))

  if (anyNA(df$source_cycle)) {
    stop("Could not normalise source_cycle for all participant records.")
  }

  if (!"age" %in% names(df)) {
    if (!"RIDAGEYR" %in% names(df)) stop("Participant core needs age or RIDAGEYR.")
    df$age <- as.numeric(df$RIDAGEYR)
  }

  if (!"hba1c" %in% names(df)) {
    if ("LBXGH" %in% names(df)) {
      df$hba1c <- as.numeric(df$LBXGH)
    } else {
      df$hba1c <- NA_real_
    }
  }

  if (!"wtmec2yr" %in% names(df)) {
    if (!"WTMEC2YR" %in% names(df)) stop("Participant core needs wtmec2yr or WTMEC2YR.")
    df$wtmec2yr <- as.numeric(df$WTMEC2YR)
  }

  if (!"wtmec6yr" %in% names(df)) {
    df$wtmec6yr <- df$wtmec2yr / 3
  }

  if (!"sdmvpsu" %in% names(df)) {
    if (!"SDMVPSU" %in% names(df)) stop("Participant core needs sdmvpsu or SDMVPSU.")
    df$sdmvpsu <- df$SDMVPSU
  }

  if (!"sdmvstra" %in% names(df)) {
    if (!"SDMVSTRA" %in% names(df)) stop("Participant core needs sdmvstra or SDMVSTRA.")
    df$sdmvstra <- df$SDMVSTRA
  }

  # Retain both lower-case standard names and the original NHANES-style names
  # expected by some downstream scripts from the EPE project.
  if (!"RIDAGEYR" %in% names(df)) df$RIDAGEYR <- df$age
  if (!"LBXGH" %in% names(df)) df$LBXGH <- df$hba1c
  if (!"WTMEC2YR" %in% names(df)) df$WTMEC2YR <- df$wtmec2yr
  if (!"WTMEC6YR" %in% names(df)) df$WTMEC6YR <- df$wtmec6yr
  if (!"SDMVPSU" %in% names(df)) df$SDMVPSU <- df$sdmvpsu
  if (!"SDMVSTRA" %in% names(df)) df$SDMVSTRA <- df$sdmvstra

  df |>
    dplyr::filter(!is.na(age), age >= 30) |>
    dplyr::arrange(source_cycle, SEQN)
}


# ------------------------------------------------------------
# 2. MCE measure construction function
# ------------------------------------------------------------

calc_mce_measures <- function(
    df,
    thresholds = c(3, 4, 5, 6),
    exclude_third_molars = TRUE
) {
  if (!"SEQN" %in% names(df)) stop("Periodontal file is missing SEQN.")

  cal_cols <- names(df)[grepl("^OHX\\d{2}LA[A-Z]$", names(df))]

  if (length(cal_cols) == 0) {
    stop("No site-level CAL columns found. Expected names like OHX02LAD.")
  }

  tooth_number <- as.integer(sub("^OHX(\\d{2}).*$", "\\1", cal_cols))

  if (exclude_third_molars) {
    keep <- !tooth_number %in% c(1, 16, 17, 32)
    cal_cols <- cal_cols[keep]
    tooth_number <- tooth_number[keep]
  }

  eligible_teeth <- sort(unique(tooth_number))
  n_teeth_modelled <- length(eligible_teeth)

  cal_mat <- as.matrix(df[, cal_cols, drop = FALSE])
  storage.mode(cal_mat) <- "numeric"

  # NHANES periodontal missing / invalid CAL codes.
  cal_mat[cal_mat %in% c(99, 999)] <- NA_real_
  cal_mat[cal_mat < 0 | cal_mat > 30] <- NA_real_

  n_sites_observed <- rowSums(!is.na(cal_mat))

  mean_CAL_site <- rowMeans(cal_mat, na.rm = TRUE)
  mean_CAL_site[n_sites_observed == 0] <- NA_real_

  site_extent_at <- function(threshold) {
    out <- rowSums(cal_mat >= threshold, na.rm = TRUE) /
      n_sites_observed * 100

    out[n_sites_observed == 0] <- NA_real_
    out
  }

  site_extent_mat <- sapply(thresholds, site_extent_at)
  colnames(site_extent_mat) <- paste0("extent_ge_", thresholds, "mm")

  mce_site <- rowMeans(site_extent_mat, na.rm = TRUE)
  mce_site[is.nan(mce_site)] <- NA_real_

  # Tooth-level companion measures are included for retained-tooth-count and
  # sensitivity analyses. Each tooth contributes its maximum observed CAL.
  tooth_max_mat <- sapply(eligible_teeth, function(tooth) {
    cols_this_tooth <- which(tooth_number == tooth)
    row_max_na(cal_mat[, cols_this_tooth, drop = FALSE])
  })

  if (is.null(dim(tooth_max_mat))) {
    tooth_max_mat <- matrix(tooth_max_mat, ncol = 1)
  }

  colnames(tooth_max_mat) <- paste0("CAL_", sprintf("%02d", eligible_teeth))

  n_present_teeth <- rowSums(!is.na(tooth_max_mat))
  n_missing_teeth <- n_teeth_modelled - n_present_teeth

  mean_CAL_tooth <- rowMeans(tooth_max_mat, na.rm = TRUE)
  mean_CAL_tooth[n_present_teeth == 0] <- NA_real_

  max_CAL_tooth <- row_max_na(tooth_max_mat)

  tooth_extent_at <- function(threshold) {
    out <- rowSums(tooth_max_mat >= threshold, na.rm = TRUE) /
      n_present_teeth * 100

    out[n_present_teeth == 0] <- NA_real_
    out
  }

  tooth_extent_mat <- sapply(thresholds, tooth_extent_at)
  colnames(tooth_extent_mat) <- paste0("tooth_extent_ge_", thresholds, "mm")

  mce_tooth <- rowMeans(tooth_extent_mat, na.rm = TRUE)
  mce_tooth[is.nan(mce_tooth)] <- NA_real_

  tibble::tibble(
    SEQN = df$SEQN,

    # Site-level observed periodontal measures.
    mean_CAL_site = mean_CAL_site,
    n_sites_observed = n_sites_observed,
    !!!as.data.frame(site_extent_mat),
    mean_extent_3to6mm = mce_site,
    mce_site = mce_site,

    # Tooth-level observed periodontal companion measures.
    n_teeth_modelled = n_teeth_modelled,
    n_present_teeth = n_present_teeth,
    n_missing_teeth = n_missing_teeth,
    mean_CAL_tooth = mean_CAL_tooth,
    max_CAL_tooth = max_CAL_tooth,
    !!!as.data.frame(tooth_extent_mat),
    mce_tooth = mce_tooth
  ) |>
    dplyr::mutate(
      # Primary shorthand for downstream analyses. The primary MCE is the
      # site-level averaged extent, matching the conventional observed-site
      # extent measures already used in the EPE comparison dataset.
      mce = mce_site,
      mean_CAL = mean_CAL_site,
      measure_model = "Observed multithreshold extent"
    )
}


# ------------------------------------------------------------
# 3. Load participant-level core data
# ------------------------------------------------------------

participant_core_path <- first_existing(c(
  "data/processed/nhanes_participant_core.rds",
  "data/processed/nhanes_demographics_hba1c.rds"
))

message("Reading participant core: ", participant_core_path)

participant_core <- readRDS(participant_core_path) |>
  standardise_participant_core()


# ------------------------------------------------------------
# 4. Load periodontal examination files and calculate MCE
# ------------------------------------------------------------

ohx_files <- c("OHXPER_F.xpt", "OHXPER_G.xpt", "OHXPER_H.xpt")

ohx_list <- purrr::set_names(ohx_files) |>
  purrr::map(read_raw_xpt)

mce_measures <- purrr::imap_dfr(ohx_list, function(df, nm) {
  calc_mce_measures(df) |>
    dplyr::mutate(source_cycle = cycle_from_filename(nm))
}) |>
  dplyr::mutate(source_cycle = normalise_cycle(source_cycle)) |>
  dplyr::arrange(source_cycle, SEQN) |>
  dplyr::filter(n_sites_observed > 0, n_present_teeth > 0)

if (anyNA(mce_measures$source_cycle)) {
  stop("Could not normalise source_cycle for all periodontal records.")
}

message("Observed periodontal site summary:")
print(summary(mce_measures$n_sites_observed))

message("Retained tooth-count summary:")
print(summary(mce_measures$n_present_teeth))


# ------------------------------------------------------------
# 5. Merge participant core and MCE measures
# ------------------------------------------------------------

analytic_mce <- mce_measures |>
  dplyr::inner_join(
    participant_core,
    by = c("SEQN", "source_cycle")
  ) |>
  dplyr::mutate(
    Age = age,
    DiabetesStatus = dplyr::case_when(
      is.na(hba1c) ~ NA_character_,
      hba1c < 5.7 ~ "Normal",
      hba1c >= 5.7 & hba1c < 6.5 ~ "Prediabetes",
      hba1c >= 6.5 ~ "Diabetes"
    ),
    DiabetesStatus = factor(
      DiabetesStatus,
      levels = c("Normal", "Prediabetes", "Diabetes")
    ),
    diabetes_binary = dplyr::case_when(
      is.na(DiabetesStatus) ~ NA_real_,
      DiabetesStatus == "Diabetes" ~ 1,
      DiabetesStatus %in% c("Normal", "Prediabetes") ~ 0
    ),
    AgeGroup_Table2 = cut(
      Age,
      breaks = c(30, 40, 50, 60, 70, 80, Inf),
      labels = c("30\u201339", "40\u201349", "50\u201359", "60\u201369", "70\u201379", "80+"),
      right = FALSE
    ),
    AgeGroup_5yr = cut(
      Age,
      breaks = c(seq(30, 80, by = 5), Inf),
      labels = c(
        "30\u201334", "35\u201339", "40\u201344", "45\u201349", "50\u201354",
        "55\u201359", "60\u201364", "65\u201369", "70\u201374", "75\u201379", "80+"
      ),
      right = FALSE
    ),
    retained_tooth_count = n_present_teeth,
    missing_tooth_count = n_missing_teeth,
    observed_site_count = n_sites_observed
  ) |>
  dplyr::select(
    SEQN,
    source_cycle,
    Age,
    age,
    AgeGroup_Table2,
    AgeGroup_5yr,
    hba1c,
    LBXGH,
    DiabetesStatus,
    diabetes_binary,
    WTMEC2YR,
    WTMEC6YR,
    SDMVPSU,
    SDMVSTRA,
    wtmec2yr,
    wtmec6yr,
    sdmvpsu,
    sdmvstra,
    mce,
    mce_site,
    mce_tooth,
    mean_extent_3to6mm,
    mean_CAL,
    mean_CAL_site,
    mean_CAL_tooth,
    max_CAL_tooth,
    extent_ge_3mm,
    extent_ge_4mm,
    extent_ge_5mm,
    extent_ge_6mm,
    tooth_extent_ge_3mm,
    tooth_extent_ge_4mm,
    tooth_extent_ge_5mm,
    tooth_extent_ge_6mm,
    n_teeth_modelled,
    n_present_teeth,
    n_missing_teeth,
    retained_tooth_count,
    missing_tooth_count,
    n_sites_observed,
    observed_site_count,
    measure_model,
    dplyr::everything()
  )


# ------------------------------------------------------------
# 6. Basic checks and summaries
# ------------------------------------------------------------

message("Rows in analytic_mce: ", nrow(analytic_mce))
message("Rows with non-missing MCE: ", sum(!is.na(analytic_mce$mce)))
message("Rows with non-missing HbA1c: ", sum(!is.na(analytic_mce$hba1c)))
message("Rows with non-missing diabetes status: ", sum(!is.na(analytic_mce$DiabetesStatus)))

measure_summary <- analytic_mce |>
  dplyr::summarise(
    n_persons = dplyr::n(),
    n_nonmissing_mce = sum(!is.na(mce)),
    n_nonmissing_mean_cal = sum(!is.na(mean_CAL)),
    n_nonmissing_hba1c = sum(!is.na(hba1c)),
    mean_mce = mean(mce, na.rm = TRUE),
    sd_mce = sd(mce, na.rm = TRUE),
    median_mce = stats::median(mce, na.rm = TRUE),
    q25_mce = stats::quantile(mce, probs = 0.25, na.rm = TRUE),
    q75_mce = stats::quantile(mce, probs = 0.75, na.rm = TRUE),
    mean_mean_CAL = mean(mean_CAL, na.rm = TRUE),
    mean_retained_teeth = mean(n_present_teeth, na.rm = TRUE),
    mean_observed_sites = mean(n_sites_observed, na.rm = TRUE)
  )

age_group_counts <- analytic_mce |>
  dplyr::count(AgeGroup_Table2, DiabetesStatus, name = "n_unweighted") |>
  dplyr::arrange(AgeGroup_Table2, DiabetesStatus)

readr::write_csv(
  measure_summary,
  "outputs/tables/mce_measure_summary.csv"
)

readr::write_csv(
  age_group_counts,
  "outputs/tables/mce_age_group_counts.csv"
)


# ------------------------------------------------------------
# 7. Save analytic dataset
# ------------------------------------------------------------

readr::write_csv(
  analytic_mce,
  "data/processed/analytic_dataset_mce.csv"
)

saveRDS(
  analytic_mce,
  "data/processed/analytic_dataset_mce.rds"
)

message("Saved: data/processed/analytic_dataset_mce.csv")
message("Saved: data/processed/analytic_dataset_mce.rds")
message("Saved: outputs/tables/mce_measure_summary.csv")
message("Saved: outputs/tables/mce_age_group_counts.csv")
