# scripts/02_build_epe_iteration01.R
# Construct periodontal measures using a cumulative missing-tooth attribution model.
#
# This script is intentionally separate from 02_construct_measures.R.
# It leaves the original three-schedule ELB pipeline unchanged and creates a new
# cumulative EPE analytic dataset for the major revision.
#
# Conceptual model:
#   q(a)  = P(periodontal reason | extraction at age a)
#   pi(A) = cumulative periodontal-attribution weight for a tooth already
#           missing at current age A
#
# The extraction-reason papers estimate q(a). NHANES observes teeth already
# missing at examination, so EPE uses pi(A), not q(A).

source("scripts/00_setup.R")
source("R/make_max_cal_by_tooth.R")
source("R/calc_observed_site_measures.R")

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

# Needed for NHANES design-based standard errors in the tooth-loss curve.
if (!requireNamespace("survey", quietly = TRUE)) {
  stop(
    "The 'survey' package is required for this script. ",
    "Install it with install.packages('survey') and rerun."
  )
}

options(survey.lonely.psu = "adjust")


# ============================================================
# 0. Helpers
# ============================================================

peek_header <- function(path, n = 60) {
  rawToChar(readBin(path, "raw", n = n))
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
  sub("^.*_([A-Z])\\.xpt$", "\\1", toupper(filename))
}

weighted_pava <- function(y, w = rep(1, length(y))) {
  # Weighted pooled adjacent violators algorithm.
  # Returns the closest non-decreasing sequence to y under weights w.

  stopifnot(length(y) == length(w))

  blocks <- tibble::tibble(
    start = seq_along(y),
    end = seq_along(y),
    value = y,
    weight = w
  )

  i <- 1

  while (i < nrow(blocks)) {
    if (blocks$value[i] > blocks$value[i + 1]) {
      pooled_weight <- blocks$weight[i] + blocks$weight[i + 1]
      pooled_value <- (
        blocks$value[i] * blocks$weight[i] +
          blocks$value[i + 1] * blocks$weight[i + 1]
      ) / pooled_weight

      blocks$value[i] <- pooled_value
      blocks$weight[i] <- pooled_weight
      blocks$end[i] <- blocks$end[i + 1]
      blocks <- blocks[-(i + 1), ]

      if (i > 1) i <- i - 1
    } else {
      i <- i + 1
    }
  }

  fitted <- numeric(length(y))
  for (b in seq_len(nrow(blocks))) {
    fitted[blocks$start[b]:blocks$end[b]] <- blocks$value[b]
  }

  fitted
}


calc_competing_site_measures <- function(df) {
  # Conventional comparison measures from observed periodontal CAL sites.
  #
  # Outputs:
  #   mean_CAL_site       = mean observed site-level CAL
  #   extent_ge_3mm       = % observed sites with CAL >= 3 mm
  #   extent_ge_4mm       = % observed sites with CAL >= 4 mm
  #   extent_ge_5mm       = % observed sites with CAL >= 5 mm
  #   extent_ge_6mm       = % observed sites with CAL >= 6 mm
  #   mean_extent_3to6mm  = average of the four extent measures
  #
  # Important:
  #   NHANES missing / invalid CAL codes such as 99 must be set to NA
  #   before calculating means or extent measures.
  
  cal_cols <- names(df)[grepl("^OHX\\d{2}LA[A-Z]$", names(df))]
  
  if (length(cal_cols) == 0) {
    stop("No site-level CAL columns found. Expected names like OHX02LAD.")
  }
  
  tooth_number <- as.integer(sub("^OHX(\\d{2}).*$", "\\1", cal_cols))
  
  # Exclude third molars
  cal_cols <- cal_cols[!tooth_number %in% c(1, 16, 17, 32)]
  
  cal_mat <- as.matrix(df[, cal_cols, drop = FALSE])
  storage.mode(cal_mat) <- "numeric"
  
  # Critical cleaning step
  cal_mat[cal_mat %in% c(99, 999)] <- NA_real_
  cal_mat[cal_mat < 0 | cal_mat > 30] <- NA_real_
  
  n_sites_observed <- rowSums(!is.na(cal_mat))
  
  mean_CAL_site <- rowMeans(cal_mat, na.rm = TRUE)
  mean_CAL_site[n_sites_observed == 0] <- NA_real_
  
  extent_at <- function(threshold) {
    out <- rowSums(cal_mat >= threshold, na.rm = TRUE) /
      n_sites_observed * 100
    
    out[n_sites_observed == 0] <- NA_real_
    out
  }
  
  extent_ge_3mm <- extent_at(3)
  extent_ge_4mm <- extent_at(4)
  extent_ge_5mm <- extent_at(5)
  extent_ge_6mm <- extent_at(6)
  
  mean_extent_3to6mm <- rowMeans(
    cbind(
      extent_ge_3mm,
      extent_ge_4mm,
      extent_ge_5mm,
      extent_ge_6mm
    ),
    na.rm = TRUE
  )
  
  mean_extent_3to6mm[
    is.nan(mean_extent_3to6mm)
  ] <- NA_real_
  
  tibble::tibble(
    SEQN = df$SEQN,
    mean_CAL_site = mean_CAL_site,
    extent_ge_3mm = extent_ge_3mm,
    extent_ge_4mm = extent_ge_4mm,
    extent_ge_5mm = extent_ge_5mm,
    extent_ge_6mm = extent_ge_6mm,
    mean_extent_3to6mm = mean_extent_3to6mm,
    n_sites_observed = n_sites_observed
  )
}


# ============================================================
# 1. Load participant-level demographics and HbA1c
# ============================================================

demo_ghb <- readRDS("data/processed/nhanes_demographics_hba1c.rds")


# ============================================================
# 2. Load periodontal examination files and construct tooth-level CAL
# ============================================================

ohx_files <- c("OHXPER_F.xpt", "OHXPER_G.xpt", "OHXPER_H.xpt")

ohx_list <- purrr::set_names(ohx_files) |>
  purrr::map(read_raw_xpt)

CAL_wide <- purrr::imap_dfr(ohx_list, function(df, nm) {
  df_out <- make_max_cal_by_tooth(df)
  df_out |>
    dplyr::mutate(source_cycle = cycle_from_filename(nm))
}) |>
  dplyr::arrange(source_cycle, SEQN)

# Drop participants with all tooth-level CAL missing.
cal_only_tmp <- CAL_wide |>
  dplyr::select(dplyr::starts_with("CAL_"))

CAL_wide <- CAL_wide |>
  dplyr::filter(rowSums(!is.na(cal_only_tmp)) > 0)


# ============================================================
# 3. Align demographics and periodontal data
# ============================================================

stopifnot("source_cycle" %in% names(demo_ghb))

CAL_wide <- CAL_wide |>
  dplyr::mutate(
    source_cycle = toupper(source_cycle),
    source_cycle = sub("^.*_([FGH])\\..*$", "\\1", source_cycle)
  )

demo_ghb <- demo_ghb |>
  dplyr::mutate(
    source_cycle = toupper(source_cycle),
    source_cycle = dplyr::case_when(
      grepl("_F\\.", source_cycle) | grepl("F$", source_cycle) ~ "F",
      grepl("_G\\.", source_cycle) | grepl("G$", source_cycle) ~ "G",
      grepl("_H\\.", source_cycle) | grepl("H$", source_cycle) ~ "H",
      TRUE ~ NA_character_
    )
  )

stopifnot(!anyNA(demo_ghb$source_cycle))

demo_ghb_sub <- demo_ghb |>
  dplyr::semi_join(CAL_wide |> dplyr::select(SEQN), by = "SEQN") |>
  dplyr::arrange(SEQN)

CAL_wide <- CAL_wide |>
  dplyr::semi_join(demo_ghb_sub |> dplyr::select(SEQN), by = "SEQN") |>
  dplyr::arrange(SEQN)

stopifnot(identical(CAL_wide$SEQN, demo_ghb_sub$SEQN))


# ============================================================
# 4. Construct conventional comparison measures
# ============================================================
# These are the observed-site measures used as competitors/comparators
# for Expected Periodontitis Experience (EPE):
#   - mean observed site-level CAL
#   - extent of observed sites with CAL >= 3, 4, 5, and 6 mm
#   - mean extent across thresholds 3-6 mm

competing_site_measures <- purrr::imap_dfr(ohx_list, function(df, nm) {
  calc_competing_site_measures(df) |>
    dplyr::mutate(source_cycle = cycle_from_filename(nm))
}) |>
  dplyr::mutate(source_cycle = toupper(source_cycle)) |>
  dplyr::arrange(source_cycle, SEQN)

competing_site_measures <- competing_site_measures |>
  dplyr::semi_join(
    CAL_wide |> dplyr::select(SEQN) |> dplyr::distinct(),
    by = "SEQN"
  )

message("Observed periodontal site summary:")
print(summary(competing_site_measures$n_sites_observed))


# ============================================================
# 5. Incident periodontal-attribution model q(a)
# ============================================================
# The primary calibration model uses three empirical extraction-attribution
# schedules. Hiltunen et al. (2023) is not included here because it describes
# late-life incident extractions in older adults and is better used to justify
# the distinction between incident extraction reasons and accumulated missing
# teeth.

open_age_upper <- 80

prob_bands <- tibble::tribble(
  ~schedule, ~author, ~country, ~age_low, ~age_high, ~n, ~p,

  # Schedule A: Chrysanthakopoulos 2011, Greece
  "Schedule A", "Chrysanthakopoulos (2011)", "Greece", 18, 24,  91, .022,
  "Schedule A", "Chrysanthakopoulos (2011)", "Greece", 25, 34, 200, .100,
  "Schedule A", "Chrysanthakopoulos (2011)", "Greece", 35, 44, 305, .222,
  "Schedule A", "Chrysanthakopoulos (2011)", "Greece", 45, 54, 198, .370,
  "Schedule A", "Chrysanthakopoulos (2011)", "Greece", 55, 64, 130, .500,
  "Schedule A", "Chrysanthakopoulos (2011)", "Greece", 65, NA,  94, .650,

  # Schedule B: Murray et al. 1996, Canada
  "Schedule B", "Murray et al. (1996)", "Canada",  1, 12,  50, .000,
  "Schedule B", "Murray et al. (1996)", "Canada", 13, 19, 175, .000,
  "Schedule B", "Murray et al. (1996)", "Canada", 20, 39, 535, .193,
  "Schedule B", "Murray et al. (1996)", "Canada", 40, 59, 488, .606,
  "Schedule B", "Murray et al. (1996)", "Canada", 60, NA, 462, .465,

  # Schedule C: Al-Shammari et al. 2006, Kuwait
  "Schedule C", "Al-Shammari et al. (2006)", "Kuwait", 12, 20, 223, .009,
  "Schedule C", "Al-Shammari et al. (2006)", "Kuwait", 21, 30, 506, .061,
  "Schedule C", "Al-Shammari et al. (2006)", "Kuwait", 31, 40, 633, .207,
  "Schedule C", "Al-Shammari et al. (2006)", "Kuwait", 41, 50, 495, .451,
  "Schedule C", "Al-Shammari et al. (2006)", "Kuwait", 51, 60, 592, .666,
  "Schedule C", "Al-Shammari et al. (2006)", "Kuwait", 61, NA, 334, .775
) |>
  dplyr::mutate(
    age_high_plot = dplyr::if_else(is.na(age_high), open_age_upper, age_high),
    age_mid = (age_low + age_high_plot) / 2,
    age_band = dplyr::if_else(
      is.na(age_high),
      paste0(age_low, "+"),
      paste0(age_low, "\u2013", age_high)
    ),
    schedule = factor(schedule),
    author = factor(author),
    approx_perio_n = n * p
  )

fit_incident <- glm(
  p ~ splines::ns(age_mid, df = 3) + schedule,
  data = prob_bands,
  weights = n,
  family = quasibinomial(link = "logit")
)

message("Incident attribution model:")
print(summary(fit_incident))

prediction_ages <- 18:open_age_upper

pred_grid <- tidyr::crossing(
  age_mid = prediction_ages,
  schedule = levels(prob_bands$schedule)
) |>
  dplyr::mutate(
    schedule = factor(schedule, levels = levels(prob_bands$schedule))
  )

pred_grid$q_source <- predict(
  fit_incident,
  newdata = pred_grid,
  type = "response"
)

incident_attribution_lookup <- pred_grid |>
  dplyr::group_by(age = age_mid) |>
  dplyr::summarise(
    q_incident = mean(q_source),
    q_min_source = min(q_source),
    q_max_source = max(q_source),
    .groups = "drop"
  )

readr::write_csv(
  incident_attribution_lookup,
  "data/processed/incident_periodontal_attribution_lookup.csv"
)

# Optional diagnostic figure for the incident model.
incident_model_plot <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = prob_bands,
    ggplot2::aes(
      x = age_low,
      xend = age_high_plot,
      y = p,
      yend = p,
      colour = schedule
    ),
    linewidth = 1,
    alpha = 0.60
  ) +
  ggplot2::geom_point(
    data = prob_bands,
    ggplot2::aes(
      x = age_mid,
      y = p,
      colour = schedule,
      shape = schedule,
      size = n
    ),
    alpha = 0.85
  ) +
  ggplot2::geom_line(
    data = pred_grid,
    ggplot2::aes(
      x = age_mid,
      y = q_source,
      colour = schedule,
      group = schedule
    ),
    linewidth = 0.7,
    alpha = 0.35
  ) +
  ggplot2::geom_ribbon(
    data = incident_attribution_lookup,
    ggplot2::aes(
      x = age,
      ymin = q_min_source,
      ymax = q_max_source
    ),
    alpha = 0.12
  ) +
  ggplot2::geom_line(
    data = incident_attribution_lookup,
    ggplot2::aes(
      x = age,
      y = q_incident,
      linetype = "Source-adjusted incident model"
    ),
    colour = "black",
    linewidth = 1.3
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 0.85),
    breaks = seq(0, 0.8, 0.2)
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(10, open_age_upper, 10),
    limits = c(10, open_age_upper)
  ) +
  ggplot2::scale_size_continuous(range = c(2, 6)) +
  ggplot2::labs(
    title = "Incident periodontal-attribution model",
    subtitle = "Empirical age-band probabilities with source-adjusted fitted curve",
    x = "Age",
    y = "Probability extraction was attributed to periodontitis",
    colour = "Empirical schedule",
    shape = "Empirical schedule",
    size = "Age-band n",
    linetype = "Model"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.box = "vertical",
    plot.title = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  filename = "outputs/figures/incident_periodontal_attribution_model.png",
  plot = incident_model_plot,
  width = 10,
  height = 7,
  dpi = 300
)


# ============================================================
# 6. Cumulative missing-tooth attribution model pi(A)
# ============================================================
# Convert q(a) into pi(A) by weighting q(a) by the age-pattern of
# accumulated missing teeth in the analytic sample.

CAL_only_df <- CAL_wide |>
  dplyr::select(dplyr::starts_with("CAL_"))

age_vec <- demo_ghb_sub$RIDAGEYR
n_teeth_modelled <- ncol(CAL_only_df)

missing_tooth_data <- demo_ghb_sub |>
  dplyr::select(SEQN, RIDAGEYR, WTMEC2YR, SDMVPSU, SDMVSTRA) |>
  dplyr::bind_cols(
    tibble::tibble(
      n_present_teeth = rowSums(!is.na(CAL_only_df)),
      n_missing_teeth = rowSums(is.na(CAL_only_df))
    )
  ) |>
  dplyr::mutate(
    age = as.integer(RIDAGEYR),
    n_missing_teeth = pmin(pmax(n_missing_teeth, 0), n_teeth_modelled),

    # NHANES 2009-2014 combines three 2-year cycles.
    # Dividing by 3 gives the pooled 6-year MEC weight.
    # This constant rescaling does not change weighted means, but it is
    # the correct pooled-cycle weight for design-based SEs.
    WTMEC6YR = WTMEC2YR / 3
  )

min_cumulative_age <- max(
  min(missing_tooth_data$age, na.rm = TRUE),
  min(incident_attribution_lookup$age, na.rm = TRUE)
)

max_cumulative_age <- min(
  max(missing_tooth_data$age, na.rm = TRUE),
  max(incident_attribution_lookup$age, na.rm = TRUE),
  80
)

q_curve <- tibble::tibble(
  age = min_cumulative_age:max_cumulative_age
) |>
  dplyr::mutate(
    q_incident = approx(
      x = incident_attribution_lookup$age,
      y = incident_attribution_lookup$q_incident,
      xout = age,
      rule = 2
    )$y
  )

# Estimate age-specific mean missing teeth and SEs using the NHANES
# complex survey design. These observed points are used for plotting and
# to support the smoothed accumulation curve used in the cumulative model.
missing_tooth_design_data <- missing_tooth_data |>
  dplyr::filter(
    age >= min_cumulative_age,
    age <= max_cumulative_age,
    !is.na(n_missing_teeth),
    !is.na(WTMEC6YR),
    !is.na(SDMVPSU),
    !is.na(SDMVSTRA)
  )

nhanes_missing_design <- survey::svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC6YR,
  nest = TRUE,
  data = missing_tooth_design_data
)

age_missing_svy <- survey::svyby(
  ~n_missing_teeth,
  ~age,
  design = nhanes_missing_design,
  FUN = survey::svymean,
  na.rm = TRUE,
  vartype = c("se", "ci")
)

age_missing_raw <- age_missing_svy |>
  as.data.frame() |>
  tibble::as_tibble()

# Make the column names robust across survey package versions.
if ("n_missing_teeth" %in% names(age_missing_raw)) {
  names(age_missing_raw)[names(age_missing_raw) == "n_missing_teeth"] <- "mean_missing"
}

se_col <- grep("^se", names(age_missing_raw), value = TRUE)[1]
if (is.na(se_col)) {
  stop("Could not find the standard error column returned by survey::svyby().")
}
names(age_missing_raw)[names(age_missing_raw) == se_col] <- "se_missing"

ci_l_col <- grep("^ci_l", names(age_missing_raw), value = TRUE)[1]
ci_u_col <- grep("^ci_u", names(age_missing_raw), value = TRUE)[1]

if (!is.na(ci_l_col)) {
  names(age_missing_raw)[names(age_missing_raw) == ci_l_col] <- "ci_low"
}
if (!is.na(ci_u_col)) {
  names(age_missing_raw)[names(age_missing_raw) == ci_u_col] <- "ci_high"
}

age_missing_raw <- age_missing_raw |>
  dplyr::mutate(
    ci_low = if ("ci_low" %in% names(age_missing_raw)) {
      pmax(ci_low, 0)
    } else {
      pmax(mean_missing - 1.96 * se_missing, 0)
    },
    ci_high = if ("ci_high" %in% names(age_missing_raw)) {
      pmin(ci_high, n_teeth_modelled)
    } else {
      pmin(mean_missing + 1.96 * se_missing, n_teeth_modelled)
    }
  ) |>
  dplyr::left_join(
    missing_tooth_design_data |>
      dplyr::count(age, name = "n_persons"),
    by = "age"
  ) |>
  dplyr::left_join(
    missing_tooth_design_data |>
      dplyr::group_by(age) |>
      dplyr::summarise(
        weight_sum = sum(WTMEC6YR, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "age"
  )

fit_missing <- lm(
  mean_missing ~ splines::ns(age, df = 4),
  data = age_missing_raw,
  weights = n_persons
)

missing_curve <- tibble::tibble(
  age = min_cumulative_age:max_cumulative_age
) |>
  dplyr::mutate(
    mean_missing_smooth = predict(
      fit_missing,
      newdata = tibble::tibble(age = age)
    ),
    mean_missing_smooth = pmin(
      pmax(mean_missing_smooth, 0),
      n_teeth_modelled
    )
  ) |>
  dplyr::left_join(
    age_missing_raw |> dplyr::select(age, n_persons),
    by = "age"
  ) |>
  dplyr::mutate(
    n_persons = dplyr::if_else(is.na(n_persons), 1L, n_persons),
    mean_missing_ordered = weighted_pava(
      y = mean_missing_smooth,
      w = n_persons
    )
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

cumulative_attribution_lookup <- missing_curve |>
  dplyr::left_join(q_curve, by = "age") |>
  dplyr::mutate(
    periodontal_missing_increment = q_incident * delta_missing
  ) |>
  dplyr::arrange(age) |>
  dplyr::mutate(
    cumulative_missing = cumsum(delta_missing),
    cumulative_periodontal_missing = cumsum(periodontal_missing_increment),
    pi_raw = cumulative_periodontal_missing / cumulative_missing,
    pi_raw = dplyr::if_else(
      is.nan(pi_raw) | is.infinite(pi_raw),
      0,
      pi_raw
    ),
    # Final cumulative attribution weight, order-restricted.
    pi_cumulative = cummax(pi_raw)
  )

cumulative_summary <- cumulative_attribution_lookup |>
  dplyr::filter(age %in% c(30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80)) |>
  dplyr::transmute(
    age,
    q_incident = scales::percent(q_incident, accuracy = 0.1),
    mean_missing = round(mean_missing_ordered, 2),
    delta_missing = round(delta_missing, 3),
    pi_raw = scales::percent(pi_raw, accuracy = 0.1),
    pi_cumulative = scales::percent(pi_cumulative, accuracy = 0.1)
  )

message("Cumulative attribution summary:")
print(cumulative_summary, n = Inf)

readr::write_csv(
  cumulative_attribution_lookup,
  "data/processed/cumulative_periodontal_attribution_lookup.csv"
)

predict_pi_cumulative_periodontal <- function(age) {
  age_clipped <- pmin(
    pmax(age, min(cumulative_attribution_lookup$age)),
    max(cumulative_attribution_lookup$age)
  )

  approx(
    x = cumulative_attribution_lookup$age,
    y = cumulative_attribution_lookup$pi_cumulative,
    xout = age_clipped,
    rule = 2
  )$y
}

# Cumulative attribution figure.
cumulative_probability_plot <- cumulative_attribution_lookup |>
  dplyr::select(age, q_incident, pi_cumulative) |>
  tidyr::pivot_longer(
    cols = c(q_incident, pi_cumulative),
    names_to = "curve",
    values_to = "probability"
  ) |>
  dplyr::mutate(
    curve = dplyr::recode(
      curve,
      q_incident = "Incident extraction attribution q(a)",
      pi_cumulative = "Cumulative missing-tooth attribution pi(A)"
    )
  ) |>
  ggplot2::ggplot(ggplot2::aes(x = age, y = probability, linetype = curve)) +
  ggplot2::geom_line(linewidth = 1.2) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 0.85),
    breaks = seq(0, 0.8, 0.2)
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(30, max_cumulative_age, 10)
  ) +
  ggplot2::labs(
    title = "Incident and cumulative periodontal-attribution probabilities",
    subtitle = "Cumulative attribution weights incident probabilities by accumulated missing teeth",
    x = "Age",
    y = "Probability / attribution weight",
    linetype = "Curve"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "bottom",
    plot.title = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  filename = "outputs/figures/cumulative_periodontal_attribution_probability.png",
  plot = cumulative_probability_plot,
  width = 9,
  height = 6,
  dpi = 300
)

missing_teeth_plot <- ggplot2::ggplot() +
  ggplot2::geom_errorbar(
    data = age_missing_raw,
    ggplot2::aes(
      x = age,
      ymin = pmax(mean_missing - se_missing, 0),
      ymax = mean_missing + se_missing
    ),
    width = 0.30,
    linewidth = 0.35,
    alpha = 0.55
  ) +
  ggplot2::geom_point(
    data = age_missing_raw,
    ggplot2::aes(x = age, y = mean_missing),
    alpha = 0.70,
    size = 1.8
  ) +
  ggplot2::geom_line(
    data = missing_curve,
    ggplot2::aes(x = age, y = mean_missing_ordered),
    linewidth = 1.2
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(30, max_cumulative_age, 10)
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 15, 5)
  ) +
  ggplot2::coord_cartesian(
    ylim = c(0, 15)
  ) +
  ggplot2::labs(
    title = "Estimated accumulation of missing teeth by age",
    subtitle = "Points show survey-weighted age-specific means with SE; line shows smoothed non-decreasing curve",
    x = "Age",
    y = "Mean number of missing teeth"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  filename = "outputs/figures/missing_teeth_accumulation_curve.png",
  plot = missing_teeth_plot,
  width = 9,
  height = 6,
  dpi = 300
)


# ============================================================
# 7. Expected Periodontitis Experience calculation
# ============================================================
# Retained teeth contribute observed tooth-level maximum CAL.
# Missing teeth contribute their expected periodontal burden:
#   pi(A) * assigned_missing_CAL
#
# This is expectation-based. No individual missing tooth is classified as
# definitively periodontal.
#
# EPE_total is the direct replacement for the previous ELB total-burden
# measure. EPE_mean_tooth standardises EPE_total by the number of modelled
# teeth, making it easier to compare with mean CAL.

calc_EPE_cumulative <- function(
    CAL_df,
    age_vec,
    pi_fun,
    advanced_CAL = 10
) {
  CAL_mat <- as.matrix(CAL_df)

  observed_retained_burden <- rowSums(CAL_mat, na.rm = TRUE)
  n_present_teeth <- rowSums(!is.na(CAL_mat))
  n_missing_teeth <- rowSums(is.na(CAL_mat))
  n_teeth_modelled <- ncol(CAL_mat)

  max_observed_CAL <- apply(
    CAL_mat,
    1,
    function(x) {
      if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
    }
  )

  assigned_missing_CAL <- pmax(
    advanced_CAL,
    max_observed_CAL,
    na.rm = TRUE
  )

  assigned_missing_CAL[is.na(assigned_missing_CAL)] <- advanced_CAL

  pi_cumulative <- pi_fun(age_vec)

  expected_missing_burden <-
    n_missing_teeth * pi_cumulative * assigned_missing_CAL

  EPE_total <-
    observed_retained_burden + expected_missing_burden

  EPE_mean_tooth <-
    EPE_total / n_teeth_modelled

  tibble::tibble(
    observed_retained_burden = observed_retained_burden,
    n_present_teeth = n_present_teeth,
    n_missing_teeth = n_missing_teeth,
    max_observed_CAL = max_observed_CAL,
    assigned_missing_CAL = assigned_missing_CAL,
    pi_cumulative = pi_cumulative,
    expected_missing_burden = expected_missing_burden,
    EPE_total = EPE_total,
    EPE_mean_tooth = EPE_mean_tooth,

    # Primary shorthand used downstream. Change to EPE_total if you want
    # the unstandardised total-burden scale as the default EPE variable.
    EPE = EPE_mean_tooth
  )
}

epe_cumulative <- calc_EPE_cumulative(
  CAL_df = CAL_only_df,
  age_vec = age_vec,
  pi_fun = predict_pi_cumulative_periodontal,
  advanced_CAL = 10
) |>
  dplyr::mutate(measure_model = "Cumulative attribution model")

id_df <- CAL_wide |>
  dplyr::select(SEQN, source_cycle)

epe_cumulative <- dplyr::bind_cols(id_df, epe_cumulative)


# ============================================================
# 8. Build analytic dataset
# ============================================================

analytic_cumulative <- epe_cumulative |>
  dplyr::left_join(
    demo_ghb_sub |>
      dplyr::select(SEQN, RIDAGEYR, LBXGH, WTMEC2YR, SDMVPSU, SDMVSTRA),
    by = "SEQN"
  ) |>
  dplyr::mutate(
    Age = RIDAGEYR,
    DiabetesStatus = dplyr::case_when(
      is.na(LBXGH) ~ NA_character_,
      LBXGH < 5.7 ~ "Normal",
      LBXGH >= 5.7 & LBXGH < 6.5 ~ "Prediabetes",
      LBXGH >= 6.5 ~ "Diabetes"
    ),
    DiabetesStatus = factor(
      DiabetesStatus,
      levels = c("Normal", "Prediabetes", "Diabetes")
    ),
    AgeGroup_Table2 = cut(
      Age,
      breaks = c(30, 40, 50, 60, 70, 80, Inf),
      labels = c("30\u201339", "40\u201349", "50\u201359", "60\u201369", "70\u201379", "80+"),
      right = FALSE
    )
  ) 

analytic_cumulative <- analytic_cumulative |>
  dplyr::left_join(
    competing_site_measures |>
      dplyr::select(
        SEQN,
        mean_CAL_site,
        extent_ge_3mm,
        extent_ge_4mm,
        extent_ge_5mm,
        extent_ge_6mm,
        mean_extent_3to6mm,
        n_sites_observed
      ),
    by = "SEQN"
  ) |>
  dplyr::mutate(
    mean_CAL = mean_CAL_site
  ) |>
  # Standardised lower-case names for derived EPE variables.
  # NHANES source variables are deliberately left unchanged.
  dplyr::rename(
    max_observed_cal = max_observed_CAL,
    assigned_missing_cal = assigned_missing_CAL,
    epe_total = EPE_total,
    epe_mean_tooth = EPE_mean_tooth,
    epe = EPE
  )

message("Rows in analytic_cumulative before site-measure join: ", nrow(analytic_cumulative))
message("Rows with non-missing HbA1c: ", sum(!is.na(analytic_cumulative$LBXGH)))
message("Rows with non-missing diabetes status: ", sum(!is.na(analytic_cumulative$DiabetesStatus)))

readr::write_csv(
  analytic_cumulative,
  "data/processed/analytic_dataset_epe.csv"
)

saveRDS(
  analytic_cumulative,
  "data/processed/analytic_dataset_epe.rds"
)

message("Saved: data/processed/analytic_dataset_epe.csv")
message("Saved: data/processed/analytic_dataset_epe.rds")
message("Saved figures to outputs/figures/")

