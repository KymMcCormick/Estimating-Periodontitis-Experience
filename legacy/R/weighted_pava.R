# R/weighted_pava.R
# Weighted pooled adjacent violators algorithm.

weighted_pava <- function(y, w = rep(1, length(y))) {
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
