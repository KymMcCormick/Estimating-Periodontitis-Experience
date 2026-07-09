# scripts/01_prepare_nhanes_site_core.R
# ============================================================
# Prepare unified NHANES 2009-2014 periodontal site-level core
# ============================================================
#
# Purpose:
#   Build one reusable person-level, wide site-level dataset for MCE,
#   EPE, and related retained-dentition periodontal measurement analyses.
#
# Inputs expected in data/raw/:
#   DEMO_F.XPT, DEMO_G.XPT, DEMO_H.XPT
#   OHXPER_F.XPT, OHXPER_G.XPT, OHXPER_H.XPT
#   GHB_F.XPT, GHB_G.XPT, GHB_H.XPT       optional but recommended
#   SMQ_F.XPT, SMQ_G.XPT, SMQ_H.XPT       optional
#
# Outputs:
#   data/derived/nhanes_perio_site_core_2009_2014.rds
#   data/derived/nhanes_perio_site_core_2009_2014.csv
#
#   Compatibility aliases:
#   data/derived/nhanes_perio_site_2009_2014.rds
#   data/derived/nhanes_perio_site_2009_2014.csv
#
#   outputs/tables/nhanes_site_core_sample_flow.csv
#   outputs/tables/nhanes_site_core_variable_check.csv
#
# Notes:
#   - Adults aged >=30 years are retained.
#   - Third molars are excluded when selecting periodontal site variables.
#   - Edentulous participants are excluded using CAL observability.
#   - CAL and PD values coded 99 or 999, or outside plausible bounds, are set to NA.
#   - WTMEC6YR is created as WTMEC2YR / 3 for pooled NHANES 2009-2014 analyses.
# ============================================================

source("scripts/00_setup.R")

# ------------------------------------------------------------
# 1. Helpers
# ------------------------------------------------------------

find_raw_file <- function(filename) {
  candidates <- unique(c(
    file.path("data/raw", filename),
    file.path("data/raw", toupper(filename)),
    file.path("data/raw", tolower(filename))
  ))

  hit <- candidates[file.exists(candidates)]

  if (length(hit) == 0) {
    stop(
      "Missing raw NHANES file: ", filename, "\n",
      "Looked for:\n  ", paste(candidates, collapse = "\n  ")
    )
  }

  hit[1]
}

peek_header <- function(path, n = 80) {
  rawToChar(readBin(path, what = "raw", n = n))
}

stop_if_html <- function(path) {
  hdr <- peek_header(path)

  if (grepl("^\\s*<!DOCTYPE html>|^\\s*<html", hdr, ignore.case = TRUE)) {
    stop(
      "File appears to be HTML rather than an XPT file: ", path, "\n",
      "This can happen when an NHANES download link saved an error page."
    )
  }

  invisible(TRUE)
}

read_nhanes_xpt <- function(filename) {
  path <- find_raw_file(filename)
  stop_if_html(path)
  haven::read_xpt(path)
}

read_cycle_files <- function(files_by_cycle, required = TRUE) {
  available <- vapply(
    files_by_cycle,
    function(f) {
      any(file.exists(file.path("data/raw", c(f, toupper(f), tolower(f)))))
    },
    logical(1)
  )

  if (required && !all(available)) {
    missing_files <- files_by_cycle[!available]
    stop(
      "Missing required NHANES file(s):\n  ",
      paste(missing_files, collapse = "\n  ")
    )
  }

  if (!required && !all(available)) {
    missing_files <- files_by_cycle[!available]
    warning(
      "Skipping unavailable optional file(s):\n  ",
      paste(missing_files, collapse = "\n  "),
      call. = FALSE
    )
  }

  files_to_read <- files_by_cycle[available]

  if (length(files_to_read) == 0) {
    return(tibble::tibble())
  }

  purrr::imap_dfr(
    files_to_read,
    function(filename, cycle) {
      read_nhanes_xpt(filename) |>
        dplyr::mutate(source_cycle = toupper(cycle))
    }
  )
}

standardise_cycle <- function(x) {
  x <- toupper(as.character(x))
  dplyr::case_when(
    grepl("F$|_F\\.", x) ~ "F",
    grepl("G$|_G\\.", x) ~ "G",
    grepl("H$|_H\\.", x) ~ "H",
    TRUE ~ x
  )
}

clean_site_values <- function(x) {
  x <- as.numeric(x)
  x[x %in% c(99, 999)] <- NA_real_
  x[x < 0 | x > 30] <- NA_real_
  x
}

site_tooth_number <- function(site_vars) {
  as.integer(sub("^OHX(\\d{2}).*$", "\\1", site_vars))
}

count_present_teeth_from_sites <- function(df, site_vars) {
  tooth_numbers <- sort(unique(site_tooth_number(site_vars)))

  present_mat <- vapply(
    tooth_numbers,
    function(tooth) {
      cols <- site_vars[site_tooth_number(site_vars) == tooth]
      rowSums(!is.na(df[, cols, drop = FALSE])) > 0
    },
    logical(nrow(df))
  )

  rowSums(present_mat, na.rm = TRUE)
}

na_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

# ------------------------------------------------------------
# 2. File lists
# ------------------------------------------------------------

cycles <- c("F", "G", "H")

demo_files <- c(
  F = "DEMO_F.XPT",
  G = "DEMO_G.XPT",
  H = "DEMO_H.XPT"
)

perio_files <- c(
  F = "OHXPER_F.XPT",
  G = "OHXPER_G.XPT",
  H = "OHXPER_H.XPT"
)

ghb_files <- c(
  F = "GHB_F.XPT",
  G = "GHB_G.XPT",
  H = "GHB_H.XPT"
)

smq_files <- c(
  F = "SMQ_F.XPT",
  G = "SMQ_G.XPT",
  H = "SMQ_H.XPT"
)

eligible_teeth <- setdiff(1:32, c(1, 16, 17, 32))

# ------------------------------------------------------------
# 3. Demographics and core survey variables
# ------------------------------------------------------------

message("Reading DEMO files...")

demo_raw <- read_cycle_files(demo_files, required = TRUE)

demo <- demo_raw |>
  dplyr::select(
    dplyr::any_of(c(
      "SEQN", "source_cycle", "RIDAGEYR", "RIAGENDR",
      "RIDRETH1", "RIDRETH3", "DMDEDUC2", "INDFMPIR",
      "SDMVPSU", "SDMVSTRA", "WTMEC2YR"
    ))
  ) |>
  dplyr::mutate(
    source_cycle = standardise_cycle(source_cycle),
    age = as.numeric(RIDAGEYR),
    wtmec2yr = as.numeric(WTMEC2YR),
    wtmec6yr = wtmec2yr / 3,
    sdmvpsu = SDMVPSU,
    sdmvstra = SDMVSTRA,
    sex = dplyr::case_when(
      RIAGENDR == 1 ~ "Male",
      RIAGENDR == 2 ~ "Female",
      TRUE ~ NA_character_
    ),
    sex = factor(sex, levels = c("Male", "Female")),
    education = dplyr::case_when(
      DMDEDUC2 %in% c(1, 2) ~ "<High school",
      DMDEDUC2 == 3 ~ "High school/GED",
      DMDEDUC2 %in% c(4, 5) ~ ">High school",
      TRUE ~ NA_character_
    ),
    education = factor(
      education,
      levels = c("<High school", "High school/GED", ">High school")
    ),
    poverty_income_ratio = as.numeric(INDFMPIR)
  ) |>
  dplyr::filter(!is.na(age), age >= 30)

# Prefer RIDRETH3 where available, but keep RIDRETH1 as well.
if ("RIDRETH3" %in% names(demo)) {
  demo <- demo |>
    dplyr::mutate(
      race_ethnicity_6 = dplyr::case_when(
        RIDRETH3 == 1 ~ "Mexican American",
        RIDRETH3 == 2 ~ "Other Hispanic",
        RIDRETH3 == 3 ~ "Non-Hispanic White",
        RIDRETH3 == 4 ~ "Non-Hispanic Black",
        RIDRETH3 == 6 ~ "Non-Hispanic Asian",
        RIDRETH3 == 7 ~ "Other/multiracial",
        TRUE ~ NA_character_
      )
    )
}

if ("RIDRETH1" %in% names(demo)) {
  demo <- demo |>
    dplyr::mutate(
      race_ethnicity_5 = dplyr::case_when(
        RIDRETH1 == 1 ~ "Mexican American",
        RIDRETH1 == 2 ~ "Other Hispanic",
        RIDRETH1 == 3 ~ "Non-Hispanic White",
        RIDRETH1 == 4 ~ "Non-Hispanic Black",
        RIDRETH1 == 5 ~ "Other/multiracial",
        TRUE ~ NA_character_
      )
    )
}

# ------------------------------------------------------------
# 4. HbA1c
# ------------------------------------------------------------

message("Reading GHB files...")

ghb_raw <- read_cycle_files(ghb_files, required = FALSE)

if (nrow(ghb_raw) > 0) {
  ghb <- ghb_raw |>
    dplyr::select(dplyr::any_of(c("SEQN", "source_cycle", "LBXGH"))) |>
    dplyr::mutate(
      source_cycle = standardise_cycle(source_cycle),
      hba1c = as.numeric(LBXGH)
    ) |>
    dplyr::group_by(SEQN, source_cycle) |>
    dplyr::summarise(
      hba1c = na_mean(hba1c),
      .groups = "drop"
    )
} else {
  ghb <- tibble::tibble(SEQN = numeric(), source_cycle = character(), hba1c = numeric())
}

# ------------------------------------------------------------
# 5. Smoking
# ------------------------------------------------------------

message("Reading SMQ files if available...")

smq_raw <- read_cycle_files(smq_files, required = FALSE)

if (nrow(smq_raw) > 0) {
  smq <- smq_raw |>
    dplyr::select(dplyr::any_of(c("SEQN", "source_cycle", "SMQ020", "SMQ040"))) |>
    dplyr::mutate(
      source_cycle = standardise_cycle(source_cycle),
      ever_smoker = dplyr::case_when(
        SMQ020 == 1 ~ "Ever smoker",
        SMQ020 == 2 ~ "Never smoker",
        TRUE ~ NA_character_
      ),
      ever_smoker = factor(ever_smoker, levels = c("Never smoker", "Ever smoker")),
      smoking_status = dplyr::case_when(
        SMQ020 == 2 ~ "Never smoker",
        SMQ020 == 1 & SMQ040 %in% c(1, 2) ~ "Current smoker",
        SMQ020 == 1 & SMQ040 == 3 ~ "Former smoker",
        TRUE ~ NA_character_
      ),
      smoking_status = factor(
        smoking_status,
        levels = c("Never smoker", "Former smoker", "Current smoker")
      )
    ) |>
    dplyr::select(SEQN, source_cycle, SMQ020, SMQ040, ever_smoker, smoking_status)
} else {
  smq <- tibble::tibble(
    SEQN = numeric(),
    source_cycle = character(),
    SMQ020 = numeric(),
    SMQ040 = numeric(),
    ever_smoker = factor(levels = c("Never smoker", "Ever smoker")),
    smoking_status = factor(levels = c("Never smoker", "Former smoker", "Current smoker"))
  )
}

# ------------------------------------------------------------
# 6. Periodontal site-level CAL and PD
# ------------------------------------------------------------

message("Reading OHXPER files...")

ohx_raw <- read_cycle_files(perio_files, required = TRUE) |>
  dplyr::mutate(source_cycle = standardise_cycle(source_cycle))

cal_site_vars_all <- grep(
  "^OHX\\d{2}LA(D|M|S|P|L|A)$",
  names(ohx_raw),
  value = TRUE
)

pd_site_vars_all <- grep(
  "^OHX\\d{2}PC(D|M|S|P|L|A)$",
  names(ohx_raw),
  value = TRUE
)

if (length(cal_site_vars_all) == 0) {
  stop("No CAL site variables detected. Expected names like OHX02LAD.")
}

if (length(pd_site_vars_all) == 0) {
  stop("No PD site variables detected. Expected names like OHX02PCD.")
}

cal_site_vars <- cal_site_vars_all[
  site_tooth_number(cal_site_vars_all) %in% eligible_teeth
]

pd_site_vars <- pd_site_vars_all[
  site_tooth_number(pd_site_vars_all) %in% eligible_teeth
]

if (length(cal_site_vars) == 0) {
  stop("No eligible non-third-molar CAL site variables remained after exclusions.")
}

if (length(pd_site_vars) == 0) {
  stop("No eligible non-third-molar PD site variables remained after exclusions.")
}

perio_site <- ohx_raw |>
  dplyr::select(
    SEQN,
    source_cycle,
    dplyr::all_of(c(cal_site_vars, pd_site_vars))
  ) |>
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(c(cal_site_vars, pd_site_vars)),
      clean_site_values
    )
  ) |>
  dplyr::mutate(
    n_cal_sites_observed = rowSums(!is.na(dplyr::across(dplyr::all_of(cal_site_vars)))),
    n_pd_sites_observed = rowSums(!is.na(dplyr::across(dplyr::all_of(pd_site_vars)))),
    n_possible_cal_sites = length(cal_site_vars),
    n_possible_pd_sites = length(pd_site_vars),
    pct_cal_sites_observed = 100 * n_cal_sites_observed / n_possible_cal_sites,
    pct_pd_sites_observed = 100 * n_pd_sites_observed / n_possible_pd_sites,
    n_present_teeth = count_present_teeth_from_sites(
      dplyr::pick(dplyr::all_of(cal_site_vars)),
      cal_site_vars
    ),
    n_modelled_teeth = length(eligible_teeth),
    n_missing_teeth = n_modelled_teeth - n_present_teeth,
    edentulous_by_cal = n_cal_sites_observed == 0
  )

# ------------------------------------------------------------
# 7. Merge and save
# ------------------------------------------------------------

sample_flow <- tibble::tibble(
  step = c(
    "DEMO adults aged >=30 years",
    "OHXPER records read",
    "After DEMO-OHXPER inner join",
    "After excluding no observed CAL sites"
  ),
  n_rows = c(
    nrow(demo),
    nrow(perio_site),
    nrow(demo |> dplyr::inner_join(perio_site, by = c("SEQN", "source_cycle"))),
    NA_integer_
  ),
  n_unique_seqn = c(
    dplyr::n_distinct(demo$SEQN),
    dplyr::n_distinct(perio_site$SEQN),
    dplyr::n_distinct((demo |> dplyr::inner_join(perio_site, by = c("SEQN", "source_cycle")))$SEQN),
    NA_integer_
  )
)

nhanes_perio_site_core <- demo |>
  dplyr::inner_join(perio_site, by = c("SEQN", "source_cycle")) |>
  dplyr::filter(!edentulous_by_cal) |>
  dplyr::left_join(ghb, by = c("SEQN", "source_cycle")) |>
  dplyr::left_join(smq, by = c("SEQN", "source_cycle")) |>
  dplyr::mutate(
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
    age_group_5yr = cut(
      age,
      breaks = c(30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, Inf),
      right = FALSE,
      labels = c(
        "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
        "60-64", "65-69", "70-74", "75-79", "80+"
      )
    ),
    age_group_10yr = cut(
      age,
      breaks = c(30, 40, 50, 60, 70, 80, Inf),
      right = FALSE,
      labels = c("30-39", "40-49", "50-59", "60-69", "70-79", "80+")
    )
  ) |>
  dplyr::arrange(source_cycle, SEQN)

sample_flow$n_rows[sample_flow$step == "After excluding no observed CAL sites"] <- nrow(nhanes_perio_site_core)
sample_flow$n_unique_seqn[sample_flow$step == "After excluding no observed CAL sites"] <-
  dplyr::n_distinct(nhanes_perio_site_core$SEQN)

variable_check <- tibble::tibble(
  item = c(
    "CAL site variables retained",
    "PD site variables retained",
    "Modelled teeth",
    "Third molars excluded",
    "Rows in analytic site core",
    "Unique participants",
    "Rows with non-missing HbA1c",
    "Rows with ever-smoker variable",
    "Rows with education variable"
  ),
  value = c(
    length(cal_site_vars),
    length(pd_site_vars),
    length(eligible_teeth),
    paste(c(1, 16, 17, 32), collapse = ", "),
    nrow(nhanes_perio_site_core),
    dplyr::n_distinct(nhanes_perio_site_core$SEQN),
    sum(!is.na(nhanes_perio_site_core$hba1c)),
    sum(!is.na(nhanes_perio_site_core$ever_smoker)),
    sum(!is.na(nhanes_perio_site_core$education))
  )
)

readr::write_csv(
  sample_flow,
  "outputs/tables/nhanes_site_core_sample_flow.csv"
)

readr::write_csv(
  variable_check,
  "outputs/tables/nhanes_site_core_variable_check.csv"
)

saveRDS(
  nhanes_perio_site_core,
  "data/derived/nhanes_perio_site_core_2009_2014.rds"
)

readr::write_csv(
  nhanes_perio_site_core,
  "data/derived/nhanes_perio_site_core_2009_2014.csv"
)

# Backward-compatible aliases for earlier exploratory scripts.
saveRDS(
  nhanes_perio_site_core,
  "data/derived/nhanes_perio_site_2009_2014.rds"
)

readr::write_csv(
  nhanes_perio_site_core,
  "data/derived/nhanes_perio_site_2009_2014.csv"
)

message("Saved: data/derived/nhanes_perio_site_core_2009_2014.rds")
message("Saved: data/derived/nhanes_perio_site_core_2009_2014.csv")
message("Saved compatibility aliases: data/derived/nhanes_perio_site_2009_2014.*")
message("Rows: ", nrow(nhanes_perio_site_core))
message("Unique participants: ", dplyr::n_distinct(nhanes_perio_site_core$SEQN))
message("CAL site variables retained: ", length(cal_site_vars))
message("PD site variables retained: ", length(pd_site_vars))
