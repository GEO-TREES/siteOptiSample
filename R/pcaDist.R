#' Get nearest neighbour distances between two PCAs
#'
#' @param x PCA scores
#' @param y PCA scores
#' @param w optional, numeric vector of weights (one per column). Larger values
#'     increase a variable’s contribution; if NULL, all variables are equally
#'     weighted.
#' @param n_pca number of PCA axes used in analysis
#' @param k number of nearest neighbours in y to return for each row in x
#' @param method distance method, either "euclidean" or "mahalanobis"
#'
#' @return if `k = 1` a vector of distances to the nearest neighbour in `y` for
#'     each row in `x`. If `k > 1` a matrix with k columns with ordered nearest
#'     neighbour distances
#'
#' @export
#'
pcaDist <- function(x, y, w = NULL, n_pca = 3, k = 1, method = "euclidean") {

  # Check input
  if (n_pca > ncol(x)) {
    stop("n_pca must not be greater than the number of principal components in x")
  }

  if (ncol(x) != ncol(y)) {
    stop("The number of principal components in x must be equal to the number in y")
  }

  if (k > nrow(y)) {
    stop("k must not be greater than the number of rows in y")
  }

  # Calculate nearest neighbor distances
  # For each landscape pixel, find distance to nearest plot
  if (method == "mahalanobis") {
    dists_mat <- mahalanobisDist(x[,1:n_pca, drop = FALSE], y[,1:n_pca, drop = FALSE], w)
  } else if (method == "euclidean") {
    dists_mat <- euclideanDist(x[,1:n_pca, drop = FALSE], y[,1:n_pca, drop = FALSE], w)
  } else {
    stop("method must be either 'euclidean' or 'mahalanobis'")
  }

  if (k == 1) {
    out <- apply(dists_mat, 1, min)
  } else {
    out <- t(apply(dists_mat, 1, function(i) {
      sort(i)[seq_len(k)]
    }))
  }

  # Return
  return(out)
}
