# ============================================================
# 00_setup.R
# Shared settings, constants, and helper functions for the EPE study.
# Site-definition correction: EPE uses four interproximal sites per tooth;
# conventional mean/extent/observability measures use all six sites.
# ============================================================
# This file contains no analysis and writes no reporting outputs.
# Every numbered script sources this file so that variable definitions,
# model specifications, formatting, and empirical attribution schedules
# have a single source of truth.

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

required_packages <- c(
  "dplyr", "tidyr", "purrr", "readr", "stringr", "tibble",
  "haven", "survey", "ggplot2", "patchwork", "scales"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running the analysis: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

options(survey.lonely.psu = "adjust")
set.seed(20260730)


# ------------------------------------------------------------
# 2. Project paths
# ------------------------------------------------------------

DIR_RAW <- file.path("data", "raw")
DIR_DERIVED <- file.path("data", "derived")
DIR_REPORT <- file.path("outputs", "reporting")

PARTICIPANT_CORE_PATH <- file.path(DIR_DERIVED, "nhanes_participant_core.rds")
ANALYTIC_DATA_PATH <- file.path(DIR_DERIVED, "analytic_dataset_epe.rds")
CALIBRATION_DATA_PATH <- file.path(DIR_DERIVED, "calibration_primary.rds")

invisible(lapply(
  c(DIR_RAW, DIR_DERIVED, DIR_REPORT),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))


# ------------------------------------------------------------
# 3. Study constants
# ------------------------------------------------------------

NHANES_CYCLES <- c("F", "G", "H")
N_ELIGIBLE_TEETH <- 28L
N_SITES_PER_TOOTH_ALL <- 6L
N_SITES_PER_TOOTH_EPE <- 4L
N_POSSIBLE_SITES <- N_ELIGIBLE_TEETH * N_SITES_PER_TOOTH_ALL
N_EPE_INTERPROXIMAL_SITES <- N_ELIGIBLE_TEETH * N_SITES_PER_TOOTH_EPE

# Site definitions are deliberately different for the proposed and
# conventional measures:
#   - EPE retained-tooth values use the maximum CAL across four
#     interproximal sites per tooth (112 possible sites across 28 teeth).
#   - Mean CAL, extent CAL >=4 mm, and periodontal observability use all
#     six NHANES sites per tooth (168 possible sites across 28 teeth).
#
# NHANES suffixes used for the EPE retained-tooth component:
#   LAS = mesiofacial; LAD = distofacial;
#   LAA = mesiolingual; LAP = distolingual.
# Mid-facial (LAM) and mid-lingual (LAL) sites are excluded from EPE.
EPE_INTERPROXIMAL_CAL_SUFFIXES <- c("LAS", "LAD", "LAA", "LAP")

OPEN_AGE_UPPER <- 80L
PRIMARY_MISSING_CAL_FLOOR <- 10

MISSING_CAL_FLOORS <- c(6, 8, 10, 12, 14)

AGE_GROUP_LEVELS <- c("30–39", "40–49", "50–59", "60–69", "70–79", "80+")
GLYCAEMIC_STATUS_LEVELS <- c("Normal", "Prediabetes", "Diabetes")
GLYCAEMIC_STATUS_LABELS <- c(
  Normal = "Normal glycaemia",
  Prediabetes = "Prediabetes",
  Diabetes = "Diabetes"
)

# The values below reproduce main-manuscript Table 1 and are the sole
# empirical inputs used by the primary and leave-one-source-out models.
# Keeping them here prevents the previous divergence between scripts that
# used different Schedule A sample sizes.
ATTRIBUTION_SCHEDULES <- tibble::tribble(
  ~schedule, ~author, ~country, ~age_low, ~age_high, ~n, ~p,
  "Schedule A", "Chrysanthakopoulos (2011)", "Greece", 18, 24, 180, .022,
  "Schedule A", "Chrysanthakopoulos (2011)", "Greece", 25, 34, 446, .100,
  "Schedule A", "Chrysanthakopoulos (2011)", "Greece", 35, 44, 576, .222,
  "Schedule A", "Chrysanthakopoulos (2011)", "Greece", 45, 54, 454, .370,
  "Schedule A", "Chrysanthakopoulos (2011)", "Greece", 55, 64, 529, .500,
  "Schedule A", "Chrysanthakopoulos (2011)", "Greece", 65, NA, 277, .650,
  "Schedule B", "Murray et al. (1996)", "Canada", 1, 12, 50, .000,
  "Schedule B", "Murray et al. (1996)", "Canada", 13, 19, 175, .000,
  "Schedule B", "Murray et al. (1996)", "Canada", 20, 39, 535, .193,
  "Schedule B", "Murray et al. (1996)", "Canada", 40, 59, 488, .606,
  "Schedule B", "Murray et al. (1996)", "Canada", 60, NA, 462, .465,
  "Schedule C", "Al-Shammari et al. (2006)", "Kuwait", 12, 20, 223, .009,
  "Schedule C", "Al-Shammari et al. (2006)", "Kuwait", 21, 30, 506, .061,
  "Schedule C", "Al-Shammari et al. (2006)", "Kuwait", 31, 40, 633, .207,
  "Schedule C", "Al-Shammari et al. (2006)", "Kuwait", 41, 50, 495, .451,
  "Schedule C", "Al-Shammari et al. (2006)", "Kuwait", 51, 60, 592, .666,
  "Schedule C", "Al-Shammari et al. (2006)", "Kuwait", 61, NA, 334, .775
) |>
  dplyr::mutate(
    age_high_plot = dplyr::if_else(is.na(age_high), OPEN_AGE_UPPER, age_high),
    age_mid = (age_low + age_high_plot) / 2,
    age_band = dplyr::if_else(
      is.na(age_high),
      paste0(age_low, "+"),
      paste0(age_low, "–", age_high)
    ),
    schedule = factor(schedule, levels = c("Schedule A", "Schedule B", "Schedule C"))
  )


# ------------------------------------------------------------
# 4. General validation and I/O helpers
# ------------------------------------------------------------

assert_has_columns <- function(data, required, object_name = deparse(substitute(data))) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(
      object_name, " is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

find_raw_file <- function(filename) {
  exact_path <- file.path(DIR_RAW, filename)
  if (file.exists(exact_path)) return(exact_path)

  candidates <- list.files(DIR_RAW, full.names = TRUE)
  matched <- candidates[tolower(basename(candidates)) == tolower(filename)]

  if (length(matched) == 1) return(matched)
  if (length(matched) > 1) {
    stop("Multiple case-insensitive matches found for ", filename, call. = FALSE)
  }

  stop("Missing raw NHANES file: ", exact_path, call. = FALSE)
}

peek_header <- function(path, n = 80L) {
  rawToChar(readBin(path, what = "raw", n = n))
}

stop_if_html <- function(path) {
  header <- peek_header(path)
  if (grepl("^\\s*<!DOCTYPE html>|^\\s*<html", header, ignore.case = TRUE)) {
    stop("Expected an XPT file but found HTML: ", path, call. = FALSE)
  }
  invisible(TRUE)
}

read_nhanes_xpt <- function(filename) {
  path <- find_raw_file(filename)
  stop_if_html(path)
  haven::read_xpt(path)
}

cycle_from_filename <- function(filename) {
  match <- stringr::str_match(toupper(basename(filename)), "_([FGH])\\.XPT$")[, 2]
  if (is.na(match)) {
    stop("Could not determine NHANES cycle from filename: ", filename, call. = FALSE)
  }
  match
}

save_report_csv <- function(data, filename) {
  path <- file.path(DIR_REPORT, filename)
  readr::write_csv(data, path, na = "")
  message("Saved: ", path)
  invisible(path)
}

save_report_plot <- function(plot, filename, width, height, dpi = 300) {
  path <- file.path(DIR_REPORT, filename)
  ggplot2::ggsave(path, plot = plot, width = width, height = height, dpi = dpi)
  message("Saved: ", path)
  invisible(path)
}


# ------------------------------------------------------------
# 5. Periodontal data helpers
# ------------------------------------------------------------

find_all_cal_site_columns <- function(data, exclude_third_molars = TRUE) {
  # Conventional mean CAL, extent CAL >=4 mm, and observability use all
  # six NHANES periodontal sites. CAL variables have names OHX##LA<site>.
  columns <- names(data)[grepl("^OHX[0-9]{2}LA[A-Z]$", toupper(names(data)))]

  if (length(columns) == 0) {
    stop(
      "No site-level CAL columns were found (expected names matching OHX##LA<site>).",
      call. = FALSE
    )
  }

  tooth_number <- as.integer(sub("^OHX([0-9]{2}).*$", "\\1", toupper(columns)))

  if (exclude_third_molars) {
    columns <- columns[!tooth_number %in% c(1, 16, 17, 32)]
  }

  if (exclude_third_molars && length(columns) != N_POSSIBLE_SITES) {
    stop(
      "Expected ", N_POSSIBLE_SITES,
      " CAL columns for all six sites on 28 eligible teeth, but found ",
      length(columns), ". Check the OHXPER variable names and site selection.",
      call. = FALSE
    )
  }

  columns
}

find_epe_interproximal_cal_columns <- function(
    data,
    exclude_third_molars = TRUE,
    suffixes = EPE_INTERPROXIMAL_CAL_SUFFIXES
) {
  # EPE uses the maximum CAL across four interproximal sites per retained
  # tooth. Mid-facial and mid-lingual sites are intentionally excluded so
  # that an isolated non-interproximal value cannot determine the complete
  # tooth-level maximum used by EPE.
  suffixes <- unique(toupper(suffixes))
  if (length(suffixes) != N_SITES_PER_TOOTH_EPE) {
    stop(
      "The EPE site definition must contain exactly ",
      N_SITES_PER_TOOTH_EPE, " unique suffixes.",
      call. = FALSE
    )
  }

  suffix_pattern <- paste0("(", paste(suffixes, collapse = "|"), ")")
  pattern <- paste0("^OHX[0-9]{2}", suffix_pattern, "$")
  columns <- names(data)[grepl(pattern, toupper(names(data)))]

  if (length(columns) == 0) {
    stop(
      "No EPE interproximal CAL columns were found. Expected suffixes: ",
      paste(suffixes, collapse = ", "), ".",
      call. = FALSE
    )
  }

  tooth_number <- as.integer(sub("^OHX([0-9]{2}).*$", "\\1", toupper(columns)))

  if (exclude_third_molars) {
    columns <- columns[!tooth_number %in% c(1, 16, 17, 32)]
  }

  if (exclude_third_molars && length(columns) != N_EPE_INTERPROXIMAL_SITES) {
    stop(
      "Expected ", N_EPE_INTERPROXIMAL_SITES,
      " EPE CAL columns (four interproximal sites on 28 eligible teeth), but found ",
      length(columns), ". Check the OHXPER variable names and EPE site definition.",
      call. = FALSE
    )
  }

  columns
}

clean_cal_matrix <- function(data, columns) {
  matrix <- as.matrix(data[, columns, drop = FALSE])
  storage.mode(matrix) <- "numeric"
  matrix[matrix %in% c(99, 999)] <- NA_real_
  matrix[matrix < 0 | matrix > 30] <- NA_real_
  matrix
}

make_epe_max_cal_by_tooth <- function(data) {
  assert_has_columns(data, "SEQN")
  columns <- find_epe_interproximal_cal_columns(
    data,
    exclude_third_molars = TRUE
  )
  tooth_number <- as.integer(sub("^OHX([0-9]{2}).*$", "\\1", toupper(columns)))
  matrix <- clean_cal_matrix(data, columns)

  out <- tibble::tibble(SEQN = data$SEQN)
  eligible_teeth <- setdiff(1:32, c(1, 16, 17, 32))

  for (tooth in eligible_teeth) {
    selected <- which(tooth_number == tooth)
    name <- sprintf("CAL_%02d", tooth)

    if (length(selected) != N_SITES_PER_TOOTH_EPE) {
      stop(
        "Tooth ", sprintf("%02d", tooth), " has ", length(selected),
        " selected EPE sites; expected ",
        N_SITES_PER_TOOTH_EPE, ".",
        call. = FALSE
      )
    }

    values <- apply(matrix[, selected, drop = FALSE], 1, function(x) {
      if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
    })
    out[[name]] <- as.numeric(values)
  }

  out
}

calc_observed_site_measures <- function(data) {
  assert_has_columns(data, "SEQN")
  columns <- find_all_cal_site_columns(data, exclude_third_molars = TRUE)
  matrix <- clean_cal_matrix(data, columns)

  n_sites_observed <- rowSums(!is.na(matrix))
  mean_cal <- rowMeans(matrix, na.rm = TRUE)
  mean_cal[n_sites_observed == 0] <- NA_real_

  extent_at <- function(threshold) {
    value <- rowSums(matrix >= threshold, na.rm = TRUE) / n_sites_observed * 100
    value[n_sites_observed == 0] <- NA_real_
    value
  }

  tibble::tibble(
    SEQN = data$SEQN,
    mean_CAL = mean_cal,
    extent_ge_4mm = extent_at(4),
    n_sites_observed = n_sites_observed
  )
}


# ------------------------------------------------------------
# 6. Survey and formatting helpers
# ------------------------------------------------------------

make_survey_design <- function(data) {
  assert_has_columns(data, c("SDMVPSU", "SDMVSTRA", "WTMEC6YR"))
  survey::svydesign(
    ids = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~WTMEC6YR,
    nest = TRUE,
    data = data
  )
}

weighted_mean_safe <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  sum(x[keep] * w[keep]) / sum(w[keep])
}

weighted_sd_safe <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  mean_value <- weighted_mean_safe(x[keep], w[keep])
  variance <- sum(w[keep] * (x[keep] - mean_value)^2) / sum(w[keep])
  if (!is.finite(variance) || variance <= 0) return(NA_real_)
  sqrt(variance)
}

z_weighted <- function(x, w) {
  mean_value <- weighted_mean_safe(x, w)
  sd_value <- weighted_sd_safe(x, w)
  if (!is.finite(mean_value) || !is.finite(sd_value) || sd_value <= 0) {
    return(rep(NA_real_, length(x)))
  }
  (x - mean_value) / sd_value
}

survey_mean_ci <- function(design, variable) {
  estimate <- survey::svymean(
    stats::as.formula(paste0("~", variable)),
    design = design,
    na.rm = TRUE
  )
  interval <- stats::confint(estimate)

  tibble::tibble(
    estimate = as.numeric(stats::coef(estimate)[1]),
    conf_low = as.numeric(interval[1, 1]),
    conf_high = as.numeric(interval[1, 2])
  )
}

format_number <- function(x, digits = 2) {
  formatC(x, format = "f", digits = digits)
}

format_ci <- function(estimate, conf_low, conf_high, digits = 2) {
  paste0(
    format_number(estimate, digits),
    " (", format_number(conf_low, digits), ", ",
    format_number(conf_high, digits), ")"
  )
}

format_percent_value <- function(x, digits = 1) {
  paste0(format_number(100 * x, digits), "%")
}

survey_rank_correlation <- function(data, x, y = "LBXGH") {
  assert_has_columns(data, c(x, y, "SDMVPSU", "SDMVSTRA", "WTMEC6YR"))

  rank_data <- data |>
    dplyr::mutate(
      .rank_x = rank(.data[[x]], ties.method = "average"),
      .rank_y = rank(.data[[y]], ties.method = "average")
    )

  design <- make_survey_design(rank_data)
  covariance <- survey::svyvar(~.rank_x + .rank_y, design = design, na.rm = TRUE)
  as.numeric(stats::cov2cor(as.matrix(covariance))[1, 2])
}

tidy_svy_term <- function(fit, term) {
  coefficient_table <- as.data.frame(summary(fit)$coefficients)
  coefficient_table$term <- rownames(coefficient_table)
  rownames(coefficient_table) <- NULL
  names(coefficient_table)[1:4] <- c("estimate", "std_error", "statistic", "p_value")

  intervals <- as.data.frame(stats::confint(fit))
  intervals$term <- rownames(intervals)
  rownames(intervals) <- NULL
  names(intervals)[1:2] <- c("conf_low", "conf_high")

  coefficient_table |>
    dplyr::left_join(intervals, by = "term") |>
    dplyr::filter(.data$term == !!term) |>
    dplyr::select(term, estimate, conf_low, conf_high, std_error, p_value)
}


# ------------------------------------------------------------
# 7. Attribution-model helpers
# ------------------------------------------------------------

fit_incident_attribution <- function(excluded_schedule = NULL) {
  data <- ATTRIBUTION_SCHEDULES

  if (!is.null(excluded_schedule)) {
    data <- data |>
      dplyr::filter(as.character(schedule) != excluded_schedule) |>
      dplyr::mutate(schedule = droplevels(schedule))
  }

  fit <- stats::glm(
    p ~ splines::ns(age_mid, df = 3) + schedule,
    data = data,
    weights = n,
    family = stats::quasibinomial(link = "logit")
  )

  prediction_grid <- tidyr::crossing(
    age_mid = 18:OPEN_AGE_UPPER,
    schedule = levels(data$schedule)
  ) |>
    dplyr::mutate(schedule = factor(schedule, levels = levels(data$schedule)))

  prediction_grid$q_source <- as.numeric(stats::predict(
    fit,
    newdata = prediction_grid,
    type = "response"
  ))

  curve <- prediction_grid |>
    dplyr::group_by(age = age_mid) |>
    dplyr::summarise(
      q_incident = mean(q_source),
      q_min_source = min(q_source),
      q_max_source = max(q_source),
      .groups = "drop"
    )

  list(
    excluded_schedule = excluded_schedule,
    data = data,
    fit = fit,
    prediction_grid = prediction_grid,
    curve = curve
  )
}

weighted_pava <- function(y, w = rep(1, length(y))) {
  stopifnot(length(y) == length(w))

  blocks <- tibble::tibble(
    start = seq_along(y),
    end = seq_along(y),
    value = y,
    weight = w
  )

  index <- 1L
  while (index < nrow(blocks)) {
    if (blocks$value[index] > blocks$value[index + 1L]) {
      pooled_weight <- blocks$weight[index] + blocks$weight[index + 1L]
      pooled_value <- (
        blocks$value[index] * blocks$weight[index] +
          blocks$value[index + 1L] * blocks$weight[index + 1L]
      ) / pooled_weight

      blocks$value[index] <- pooled_value
      blocks$weight[index] <- pooled_weight
      blocks$end[index] <- blocks$end[index + 1L]
      blocks <- blocks[-(index + 1L), , drop = FALSE]
      if (index > 1L) index <- index - 1L
    } else {
      index <- index + 1L
    }
  }

  fitted <- numeric(length(y))
  for (row in seq_len(nrow(blocks))) {
    fitted[blocks$start[row]:blocks$end[row]] <- blocks$value[row]
  }
  fitted
}

make_missing_accumulation <- function(data) {
  assert_has_columns(
    data,
    c("Age", "n_missing_teeth", "WTMEC6YR", "SDMVPSU", "SDMVSTRA")
  )

  analysis_data <- data |>
    dplyr::transmute(
      Age = as.numeric(Age),
      age = as.integer(Age),
      n_missing_teeth = as.numeric(n_missing_teeth),
      WTMEC6YR = as.numeric(WTMEC6YR),
      SDMVPSU,
      SDMVSTRA
    ) |>
    dplyr::filter(
      age >= 30,
      age <= OPEN_AGE_UPPER,
      is.finite(n_missing_teeth),
      is.finite(WTMEC6YR),
      WTMEC6YR > 0,
      !is.na(SDMVPSU),
      !is.na(SDMVSTRA)
    )

  design <- make_survey_design(analysis_data)
  by_age <- survey::svyby(
    ~n_missing_teeth,
    ~age,
    design = design,
    FUN = survey::svymean,
    na.rm = TRUE,
    vartype = "se"
  ) |>
    as.data.frame() |>
    tibble::as_tibble()

  if ("n_missing_teeth" %in% names(by_age)) {
    names(by_age)[names(by_age) == "n_missing_teeth"] <- "mean_missing"
  }

  se_column <- grep("^se", names(by_age), value = TRUE)[1]
  if (is.na(se_column)) stop("Could not identify the survey standard-error column.")
  names(by_age)[names(by_age) == se_column] <- "se_missing"

  by_age <- by_age |>
    dplyr::left_join(
      analysis_data |> dplyr::count(age, name = "n_persons"),
      by = "age"
    )

  fit <- stats::lm(
    mean_missing ~ splines::ns(age, df = 4),
    data = by_age,
    weights = n_persons
  )

  curve <- tibble::tibble(age = 30:OPEN_AGE_UPPER) |>
    dplyr::mutate(
      mean_missing_smooth = as.numeric(stats::predict(fit, newdata = tibble::tibble(age = age))),
      mean_missing_smooth = pmin(pmax(mean_missing_smooth, 0), N_ELIGIBLE_TEETH)
    ) |>
    dplyr::left_join(by_age |> dplyr::select(age, n_persons), by = "age") |>
    dplyr::mutate(
      n_persons = dplyr::coalesce(as.integer(n_persons), 1L),
      mean_missing_ordered = weighted_pava(mean_missing_smooth, n_persons)
    ) |>
    dplyr::arrange(age) |>
    dplyr::mutate(
      delta_missing = mean_missing_ordered - dplyr::lag(mean_missing_ordered),
      delta_missing = dplyr::if_else(
        is.na(delta_missing),
        mean_missing_ordered,
        delta_missing
      ),
      delta_missing = pmax(delta_missing, 0)
    )

  list(points = by_age, curve = curve, fit = fit)
}

make_cumulative_attribution <- function(incident_curve, missing_curve) {
  q_curve <- incident_curve |>
    dplyr::filter(age >= 30, age <= OPEN_AGE_UPPER) |>
    dplyr::select(age, q_incident)

  missing_curve |>
    dplyr::left_join(q_curve, by = "age") |>
    dplyr::mutate(
      periodontal_missing_increment = q_incident * delta_missing,
      cumulative_missing = cumsum(delta_missing),
      cumulative_periodontal_missing = cumsum(periodontal_missing_increment),
      pi_raw = cumulative_periodontal_missing / cumulative_missing,
      pi_raw = dplyr::if_else(
        is.nan(pi_raw) | is.infinite(pi_raw),
        0,
        pi_raw
      ),
      pi_cumulative = cummax(pi_raw)
    )
}

lookup_curve_by_age <- function(age, curve, value) {
  age_clipped <- pmin(pmax(age, min(curve$age)), max(curve$age))
  stats::approx(
    x = curve$age,
    y = curve[[value]],
    xout = age_clipped,
    rule = 2
  )$y
}

calculate_epe <- function(tooth_cal_data, age, cumulative_curve, floor = 10) {
  # tooth_cal_data must contain the retained-tooth maxima produced by
  # make_epe_max_cal_by_tooth(), i.e., maxima across the four interproximal CAL sites only. 
  # Missing tooth positions are assigned a fixed CAL severity,
  # weighted by the age-specific cumulative attribution function.
  cal_columns <- names(tooth_cal_data)[grepl("^CAL_\\d{2}$", names(tooth_cal_data))]
  if (length(cal_columns) != N_ELIGIBLE_TEETH) {
    stop(
      "Expected ", N_ELIGIBLE_TEETH, " tooth-level CAL columns, found ",
      length(cal_columns), "."
    )
  }

  matrix <- as.matrix(tooth_cal_data[, cal_columns, drop = FALSE])
  storage.mode(matrix) <- "numeric"

  observed_retained_burden <- rowSums(matrix, na.rm = TRUE)
  n_present_teeth <- rowSums(!is.na(matrix))
  n_missing_teeth <- N_ELIGIBLE_TEETH - n_present_teeth

  max_observed_cal <- apply(matrix, 1, function(x) {
    if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  })

  assigned_missing_cal = as.numeric(floor)

  pi_cumulative <- lookup_curve_by_age(age, cumulative_curve, "pi_cumulative")
  expected_missing_burden <- n_missing_teeth * pi_cumulative * assigned_missing_cal
  epe_total <- observed_retained_burden + expected_missing_burden
  epe <- epe_total / N_ELIGIBLE_TEETH

  tibble::tibble(
    observed_retained_burden,
    n_present_teeth,
    n_missing_teeth,
    max_observed_cal,
    assigned_missing_cal,
    pi_cumulative,
    expected_missing_burden,
    epe_total,
    epe
  )
}

recalculate_epe_from_components <- function(data, cumulative_curve, floor = 10) {
  # observed_retained_burden and max_observed_cal originate from the
  # four-interproximal-site retained-tooth maxima in the primary dataset.
  data |>
    dplyr::mutate(
      pi_model = lookup_curve_by_age(Age, cumulative_curve, "pi_cumulative"),
      assigned_missing_cal_model = as.numeric(floor),
      assigned_missing_cal_model = dplyr::if_else(
        is.finite(assigned_missing_cal_model),
        assigned_missing_cal_model,
        as.numeric(floor)
      ),
      expected_missing_burden_model =
        n_missing_teeth * pi_model * assigned_missing_cal_model,
      observed_retained_component_model =
        observed_retained_burden / N_ELIGIBLE_TEETH,
      expected_missing_component_model =
        expected_missing_burden_model / N_ELIGIBLE_TEETH,
      epe_model = observed_retained_component_model + expected_missing_component_model
    )
}


# ------------------------------------------------------------
# 8. Model-based late-life contrasts
# ------------------------------------------------------------

fit_late_life_contrasts <- function(data, outcome, status, z_critical = 1.96) {
  group_data <- data |>
    dplyr::filter(
      DiabetesStatus == status,
      is.finite(Age),
      is.finite(.data[[outcome]]),
      is.finite(WTMEC6YR),
      WTMEC6YR > 0,
      !is.na(SDMVPSU),
      !is.na(SDMVSTRA)
    )

  if (nrow(group_data) == 0) {
    stop("No usable observations for ", outcome, " in ", status, ".")
  }

  design <- make_survey_design(group_data)
  fit <- survey::svyglm(
    stats::as.formula(paste0(outcome, " ~ Age + I(Age^2)")),
    design = design
  )

  prediction_ages <- tibble::tibble(Age = c(60, 70, 80))
  matrix <- stats::model.matrix(stats::delete.response(stats::terms(fit)), prediction_ages)
  matrix <- matrix[, names(stats::coef(fit)), drop = FALSE]

  contrasts <- list(
    change_60_70 = matrix[2, ] - matrix[1, ],
    change_70_80 = matrix[3, ] - matrix[2, ],
    late_life_curvature = matrix[3, ] - 2 * matrix[2, ] + matrix[1, ]
  )

  purrr::imap_dfr(contrasts, function(vector, quantity) {
    estimate <- drop(vector %*% stats::coef(fit))
    variance <- drop(vector %*% stats::vcov(fit) %*% vector)
    standard_error <- sqrt(max(variance, 0))

    tibble::tibble(
      DiabetesStatus = status,
      outcome = outcome,
      quantity = quantity,
      estimate = as.numeric(estimate),
      standard_error = as.numeric(standard_error),
      conf_low = as.numeric(estimate - z_critical * standard_error),
      conf_high = as.numeric(estimate + z_critical * standard_error),
      n_unweighted = nrow(group_data)
    )
  })
}


# ------------------------------------------------------------
# 9. Common analysis-data preparations
# ------------------------------------------------------------

add_reporting_variables <- function(data) {
  data |>
    dplyr::mutate(
      WTMEC6YR = WTMEC2YR / 3,
      pct_sites_observed = n_sites_observed / N_POSSIBLE_SITES * 100,
      DiabetesStatus = dplyr::case_when(
        is.na(LBXGH) ~ NA_character_,
        LBXGH < 5.7 ~ "Normal",
        LBXGH >= 5.7 & LBXGH < 6.5 ~ "Prediabetes",
        LBXGH >= 6.5 ~ "Diabetes",
        TRUE ~ NA_character_
      ),
      DiabetesStatus = factor(DiabetesStatus, levels = GLYCAEMIC_STATUS_LEVELS),
      AgeGroup = cut(
        Age,
        breaks = c(30, 40, 50, 60, 70, 80, Inf),
        labels = AGE_GROUP_LEVELS,
        right = FALSE
      )
    )
}

prepare_hba1c_complete_case <- function(data) {
  data |>
    add_reporting_variables() |>
    dplyr::filter(
      is.finite(LBXGH),
      is.finite(Age),
      is.finite(epe),
      is.finite(mean_CAL),
      is.finite(extent_ge_4mm),
      is.finite(pct_sites_observed),
      is.finite(WTMEC6YR),
      WTMEC6YR > 0,
      !is.na(SDMVPSU),
      !is.na(SDMVSTRA)
    ) |>
    dplyr::mutate(
      age_c = Age - 55,
      z_epe = z_weighted(epe, WTMEC6YR),
      z_mean_CAL = z_weighted(mean_CAL, WTMEC6YR),
      z_extent_ge_4mm = z_weighted(extent_ge_4mm, WTMEC6YR),
      z_sites_observed = z_weighted(pct_sites_observed, WTMEC6YR)
    )
}

study_theme <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_line(linewidth = 0.25)
    )
}
