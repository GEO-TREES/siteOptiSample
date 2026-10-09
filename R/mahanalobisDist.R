#' Calculate mahalanobis distance between two matrices
#'
#' @param x numeric matrix of matrix 1
#' @param y numeric matrix of matrix 2
#' @param w optional, numeric vector of weights (one per column). Larger values
#'     increase a variable’s contribution; if NULL, all variables are equally
#'     weighted.
#'
#' @noRd
#'

mahalanobisDist <- function(x, y, w = NULL) {
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

  x <- as.matrix(x)
  y <- as.matrix(y)
  S_inv <- solve(stats::cov(rbind(x, y)))

  # Mahalanobis distance is euclidean distance after transforming by the
  # Cholesky factor of the inverse covariance: S_inv = t(L) %*% L
  L <- chol(S_inv)
  x_t <- sweep(x, 2, w, "*") %*% t(L)
  y_t <- sweep(y, 2, w, "*") %*% t(L)

  return(distMat(x_t, y_t))
}
