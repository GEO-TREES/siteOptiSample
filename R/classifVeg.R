#' Classify pixels into structural types
#'
#' @param r_pca PCA scores of pixels in structural space, as returned by
#'     `PCALandscape()$r_pca`
#' @param k integer, number of groups (clusters)
#' @param method either "kmeans" or "gmm" (Gaussian Mixture Model)
#' @param n_pca number of PCA axes used in analysis
#' @param nstart integer, number of random sets for K-means (`method` = "kmeans"). Ignored if `method` = "gmm"
#' @param spatial_weight weight applied to spatial coordinates to encourage 
#'     contiguous patches. 0 ignores coordinates, >0 increases spatial grouping. Kmeans only
#' @param coords matrix or dataframe of X, Y coordinates corresponding to the rows of r_pca. 
#'     Required if spat_weight > 0.
#' 
#' @return vector containing cluster assignments for each pixel
#' 
#' @import mclust
#' @importFrom stats kmeans
#'
#' @export
#' 
classifVeg <- function(r_pca, k, method = c("kmeans", "gmm"), n_pca = 3, nstart = 25, spatial_weight = 0, coords = NULL) {
  # Match the method argument to allow shorthand inputs
  method <- match.arg(method)

  # Extract PCA scores from chosen PCs
  r_df <- as.data.frame(r_pca)[,paste0("PC", 1:n_pca), drop = FALSE]

  # Incorporate spatial coordinates if spat_weight > 0
  if (spatial_weight > 0) {
    if (method == "gmm") {
      warning("method is 'gmm', spatial_weight will be ignored")
    } else {
      if (is.null(coords)) {
        stop("'coords' must be provided when spatial_weight > 0.")
      }

      # Ensure coords is a matrix
      coords <- as.matrix(coords)

      if (nrow(coords) != nrow(r_df)) {
        stop("Number of rows in 'coords' does not match the number of rows in 'r_pca'.")
      }

      # Calculate average variance of structural PCA axes
      pca_mean_var <- mean(apply(r_df, 2, stats::var))

      # Scale coordinates to unit variance, then apply the spatial weight
      xy_scaled <- scale(coords) * sqrt(pca_mean_var) * spatial_weight
      
      # Bind coordinates to the PCA scores
      r_df <- cbind(r_df, xy_scaled)
    }
  }

  if (method == "kmeans") {
    # Run K-means clustering
    km <- stats::kmeans(r_df, centers = k, nstart = nstart)
    out <- as.factor(unname(km$cluster))
    
  } else if (method == "gmm") {
    gmm_mod <- mclust::Mclust(r_df, G = k, verbose = FALSE)
    
    # Extract hard classifications (highest probability cluster)
    out <- as.factor(unname(gmm_mod$classification))
  }
  
  # Return
  return(out)
}


