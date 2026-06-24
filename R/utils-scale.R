get_minmax_params <- function(X) {
  X <- as.matrix(X)

  list(
    min = apply(X, 2, min, na.rm = TRUE),
    max = apply(X, 2, max, na.rm = TRUE)
  )
}

apply_minmax_params <- function(X, params) {
  X <- as.matrix(X)

  for (j in seq_len(ncol(X))) {
    den <- params$max[j] - params$min[j]

    if (den == 0) {
      X[, j] <- 0
    } else {
      X[, j] <- (X[, j] - params$min[j]) / den
    }
  }

  X
}
