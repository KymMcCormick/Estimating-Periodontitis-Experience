# code/functions/make_max_cal_by_tooth.R
# Construct tooth-level maximum clinical attachment loss (CAL) from NHANES OHXPER files.
#
# - NHANES periodontal exam includes multiple sites per tooth.
# - CAL is recorded at sites with variable names beginning with OHX<tooth><site>.
# - This function:
#   1) selects CAL variables for retained teeth (excludes 01,16,17,32),
#   2) recodes 99 and negative values to NA (NHANES missing codes),
#   3) computes per-tooth max CAL across sites,
#   4) returns a wide dataset with columns CAL_02 ... CAL_31 plus SEQN.

make_max_cal_by_tooth <- function(
    ohxper_df,
    patterns_CAL = c("LAS", "LAD", "LAP", "LAA"),
    teeth_keep = sprintf("%02d", setdiff(1:32, c(1, 16, 17, 32)))
) {
  
  df <- dplyr::as_tibble(ohxper_df)
  
  cal_cols <- names(df)[
    grepl(paste0(patterns_CAL, collapse = "|"), names(df)) &
      grepl("^OHX\\d{2}", names(df))
  ]
  
  df_cal <- df %>%
    dplyr::select(SEQN, dplyr::all_of(cal_cols)) %>%
    dplyr::mutate(
      dplyr::across(
        -SEQN,
        ~ {
          x <- as.numeric(.x)
          x[x == 99 | x < 0] <- NA_real_
          x
        }
      )
    )
  
  out <- df_cal %>% dplyr::select(SEQN)
  
  for (tooth in teeth_keep) {
    tooth_cols <- names(df_cal)[
      grepl(paste0("^OHX", tooth), names(df_cal)) &
        grepl(paste0(patterns_CAL, collapse = "|"), names(df_cal))
    ]
    
    colname <- paste0("CAL_", tooth)
    
    if (length(tooth_cols) == 0) {
      out[[colname]] <- NA_real_
    } else {
      out[[colname]] <- do.call(pmax, c(df_cal[tooth_cols], na.rm = TRUE))
      out[[colname]][is.infinite(out[[colname]])] <- NA_real_
    }
  }
  
  out
}

