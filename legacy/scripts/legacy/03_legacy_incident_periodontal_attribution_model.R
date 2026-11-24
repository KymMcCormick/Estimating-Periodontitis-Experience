# ============================================================
# Incident periodontal-attribution model
# ============================================================
#
# Purpose:
#   Estimate q(a), the age-specific probability that an extraction
#   occurring at age a was attributed to periodontitis.
#
# Important:
#   This is an incident extraction-reason model.
#   It is NOT the final probability assigned to teeth already missing
#   at NHANES examination.
#
# Later step:
#   q(a) will be integrated over the age-pattern of accumulated
#   missing teeth to estimate pi(A), the cumulative missing-tooth
#   attribution weight.
#
# Conceptual distinction:
#   q(a)  = probability current extraction was periodontal
#   pi(A) = cumulative periodontal contribution to teeth already missing
# ============================================================


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

library(tidyverse)
library(splines)
library(scales)


# ------------------------------------------------------------
# 1. Enter primary empirical calibration schedules
# ------------------------------------------------------------

open_age_upper <- 80

prob_bands <- tribble(
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
  mutate(
    age_high_plot = if_else(is.na(age_high), open_age_upper, age_high),
    age_mid = (age_low + age_high_plot) / 2,
    age_band = if_else(
      is.na(age_high),
      paste0(age_low, "+"),
      paste0(age_low, "\u2013", age_high)
    ),
    schedule = factor(schedule),
    author = factor(author),
    approx_perio_n = n * p
  )


print(prob_bands)


# ------------------------------------------------------------
# 2. Fit source-adjusted incident attribution model
# ------------------------------------------------------------
# q(a) = P(periodontal reason | extraction at age a)
#
# The model is weighted by age-band n.
# quasibinomial is used because the inputs are grouped proportions
# from heterogeneous empirical sources.

fit_incident <- glm(
  p ~ ns(age_mid, df = 3) + schedule,
  data = prob_bands,
  weights = n,
  family = quasibinomial(link = "logit")
)

summary(fit_incident)


# ------------------------------------------------------------
# 3. Predict q(a) across the adult age range
# ------------------------------------------------------------

prediction_ages <- 18:open_age_upper

pred_grid <- crossing(
  age_mid = prediction_ages,
  schedule = levels(prob_bands$schedule)
) |>
  mutate(
    schedule = factor(schedule, levels = levels(prob_bands$schedule))
  )

pred_grid$q_source <- predict(
  fit_incident,
  newdata = pred_grid,
  type = "response"
)

incident_curve <- pred_grid |>
  group_by(age = age_mid) |>
  summarise(
    q_incident = mean(q_source),
    q_min_source = min(q_source),
    q_max_source = max(q_source),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 4. Inspect fitted incident probabilities
# ------------------------------------------------------------

incident_summary <- incident_curve |>
  filter(age %in% c(18, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80)) |>
  transmute(
    age,
    q_incident = percent(q_incident, accuracy = 0.1),
    q_min_source = percent(q_min_source, accuracy = 0.1),
    q_max_source = percent(q_max_source, accuracy = 0.1)
  )

print(incident_summary, n = Inf)


# ------------------------------------------------------------
# 5. Plot empirical schedules and incident attribution model
# ------------------------------------------------------------

incident_model_plot <- ggplot() +
  
  # Empirical age-band probabilities
  geom_segment(
    data = prob_bands,
    aes(
      x = age_low,
      xend = age_high_plot,
      y = p,
      yend = p,
      colour = schedule
    ),
    linewidth = 1,
    alpha = 0.60
  ) +
  
  # Empirical midpoint probabilities, sized by n
  geom_point(
    data = prob_bands,
    aes(
      x = age_mid,
      y = p,
      colour = schedule,
      shape = schedule,
      size = n
    ),
    alpha = 0.85
  ) +
  
  # Source-specific fitted curves
  geom_line(
    data = pred_grid,
    aes(
      x = age_mid,
      y = q_source,
      colour = schedule,
      group = schedule
    ),
    linewidth = 0.7,
    alpha = 0.35
  ) +
  
  # Range across schedule-specific predictions
  geom_ribbon(
    data = incident_curve,
    aes(
      x = age,
      ymin = q_min_source,
      ymax = q_max_source
    ),
    alpha = 0.12
  ) +
  
  # Pooled source-adjusted incident model
  geom_line(
    data = incident_curve,
    aes(
      x = age,
      y = q_incident,
      linetype = "Source-adjusted incident model"
    ),
    colour = "black",
    linewidth = 1.3
  ) +
  
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 0.85),
    breaks = seq(0, 0.8, 0.2)
  ) +
  
  scale_x_continuous(
    breaks = seq(10, open_age_upper, 10),
    limits = c(10, open_age_upper)
  ) +
  
  scale_size_continuous(
    range = c(2, 6)
  ) +
  
  labs(
    title = "Incident periodontal-attribution model",
    subtitle = "Empirical age-band probabilities with source-adjusted fitted curve",
    x = "Age",
    y = "Probability extraction was attributed to periodontitis",
    colour = "Empirical schedule",
    shape = "Empirical schedule",
    size = "Age-band n",
    linetype = "Model"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    plot.title = element_text(face = "bold"),
    plot.margin = margin(10, 20, 10, 10)
  )

print(incident_model_plot)


# Optional save
# ggsave(
#   filename = "incident_periodontal_attribution_model.png",
#   plot = incident_model_plot,
#   width = 10,
#   height = 7,
#   dpi = 300
# )


# ------------------------------------------------------------
# 6. Leave-one-source-out sensitivity check
# ------------------------------------------------------------
# This checks whether the incident curve is strongly driven by
# any single empirical schedule.

fit_leave_one_out <- function(excluded_schedule = NULL) {
  
  dat <- prob_bands
  
  if (!is.null(excluded_schedule)) {
    dat <- dat |>
      filter(schedule != excluded_schedule) |>
      mutate(schedule = droplevels(schedule))
  }
  
  fit <- glm(
    p ~ ns(age_mid, df = 3) + schedule,
    data = dat,
    weights = n,
    family = quasibinomial(link = "logit")
  )
  
  pg <- crossing(
    age_mid = prediction_ages,
    schedule = levels(dat$schedule)
  ) |>
    mutate(
      schedule = factor(schedule, levels = levels(dat$schedule))
    )
  
  pg$q_source <- predict(fit, newdata = pg, type = "response")
  
  pg |>
    group_by(age = age_mid) |>
    summarise(
      q_incident = mean(q_source),
      .groups = "drop"
    ) |>
    mutate(
      model = if_else(
        is.null(excluded_schedule),
        "All schedules",
        paste0("Without ", excluded_schedule)
      )
    )
}


sensitivity_curves <- bind_rows(
  fit_leave_one_out(NULL),
  fit_leave_one_out("Schedule A"),
  fit_leave_one_out("Schedule B"),
  fit_leave_one_out("Schedule C")
)


sensitivity_summary <- sensitivity_curves |>
  filter(age %in% c(30, 40, 50, 60, 70, 80)) |>
  transmute(
    model,
    age,
    q_incident = percent(q_incident, accuracy = 0.1)
  ) |>
  pivot_wider(
    names_from = age,
    values_from = q_incident,
    names_prefix = "Age "
  )

print(sensitivity_summary, n = Inf)


sensitivity_plot <- ggplot(sensitivity_curves) +
  geom_line(
    aes(
      x = age,
      y = q_incident,
      linetype = model
    ),
    colour = "black",
    linewidth = 1
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 0.85),
    breaks = seq(0, 0.8, 0.2)
  ) +
  scale_x_continuous(
    breaks = seq(10, open_age_upper, 10),
    limits = c(10, open_age_upper)
  ) +
  labs(
    title = "Leave-one-source-out sensitivity of incident attribution model",
    x = "Age",
    y = "Incident periodontal-attribution probability",
    linetype = "Model"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

print(sensitivity_plot)


# ------------------------------------------------------------
# 7. Prediction function for q(a)
# ------------------------------------------------------------
# Returns the incident periodontal-attribution probability.
# This is the input for the later cumulative model.

predict_q_incident_periodontal <- function(age) {
  
  age_clipped <- pmin(
    pmax(age, min(incident_curve$age)),
    max(incident_curve$age)
  )
  
  approx(
    x = incident_curve$age,
    y = incident_curve$q_incident,
    xout = age_clipped,
    rule = 2
  )$y
}


# Quick test
tibble(
  age = c(30, 40, 50, 60, 70, 80),
  q_incident_periodontal = predict_q_incident_periodontal(age)
) |>
  mutate(
    q_incident_periodontal =
      percent(q_incident_periodontal, accuracy = 0.1)
  ) |>
  print()


# ------------------------------------------------------------
# 8. Save lookup table for cumulative model
# ------------------------------------------------------------

incident_attribution_lookup <- incident_curve |>
  transmute(
    age,
    q_incident,
    q_min_source,
    q_max_source
  )

# Optional save
# write_csv(
#   incident_attribution_lookup,
#   "incident_periodontal_attribution_lookup.csv"
# )