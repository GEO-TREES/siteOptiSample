#' Calculate euclidean distance between two matrices
#'
#' @param x numeric matrix of matrix 1
#' @param y numeric matrix of matrix 2
#' @param w optional, numeric vector of weights (one per column). Larger values
#'     increase a variable’s contribution; if NULL, all variables are equally
#'     weighted.
#'
#' @noRd
#'

euclideanDist <- function(x, y, w = NULL) {
  # Check matrices have the same number of columns
  if (ncol(x) != ncol(y)) {
    stop("'x' and 'y' must have the same number of columns")
  }

  # Check weights equal the number of columns as matrices
  if (!is.null(w) && ncol(x) != length(w)) {
    stop("The length of 'w' must equal the number of columns in 'x'")
  }

  # If no weights supplied, use equal weights
  if (is.null(w)) {
    w <- rep(1, ncol(x))
  }

  if (any(w < 0)) {
    stop("'w' must be non-negative")
  }

  # Weights multiply squared differences, so scale columns by sqrt(w)
  x <- sweep(as.matrix(x), 2, sqrt(w), "*")
  y <- sweep(as.matrix(y), 2, sqrt(w), "*")

  return(distMat(x, y))
}

#' Calculate unweighted euclidean distance between rows of two matrices
#'
#' @param x numeric matrix
#' @param y numeric matrix with the same number of columns as `x`
#'
#' @return matrix with `nrow(x)` rows and `nrow(y)` columns
#'
#' @noRd
#'
distMat <- function(x, y) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  d2 <- outer(rowSums(x^2), rowSums(y^2), "+") - 2 * tcrossprod(x, y)
  sqrt(pmax(d2, 0))
}

#' Calculate euclidean distance from each row of a matrix to its nearest row
#'     in a second matrix
#'
#' @param x numeric matrix
#' @param y numeric matrix with the same number of columns as `x`
#' @param chunk_size maximum number of distances held in memory at once
#'
#' @return vector of length `nrow(x)`
#'
#' @noRd
#'
nnDist <- function(x, y, chunk_size = 1e7) {
  out <- rep(Inf, nrow(x))
  step <- max(1, floor(chunk_size / nrow(x)))
  for (s in seq(1, nrow(y), by = step)) {
    idx <- s:min(nrow(y), s + step - 1)
    out <- pmin(out, apply(distMat(x, y[idx, , drop = FALSE]), 1, min))
  }
  out
}
