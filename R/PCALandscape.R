#' Run a PCA on landscape structural metrics, and optionally place plots in
#' that space
#' 
#' @param r either a SpatRaster or dataframe with structural metrics
#' @param p optional dataframe containing structural metrics for each plot, e.g.
#'     as returned by `extractPlotMetrics()`
#' @param ... additional arguments passed to `prcomp()`  
#'
#' @return list containing: `r_pca`: PCA object from pixel values, and
#'     optionally `p_pca`: plot values in PCA space
#' 
#' @export
#' 
PCALandscape <- function(r, p = NULL, ...) {

  # Convert raster to data frame for landscape
  if (inherits(r, "SpatRaster")) {
    r_df <- as.data.frame(r, xy = FALSE, na.rm = TRUE)
  } else {
    r_df <- as.data.frame(r)
  }
  
  # Select only metric columns 
  metric_cols <- names(r_df)
  
  # PCA to reduce dimensionality
  r_pca <- prcomp(r_df[,metric_cols], ...)
  
  # Project plots into PCA space
  if (!is.null(p)) { 
    p_pca <- scale(p[,metric_cols], 
      center = r_pca$center, scale = r_pca$scale) %*% r_pca$rotation
  } else {
    p_pca <- NULL
  }

  # Return
  return(list(
    "r_pca" = r_pca,
    "p_pca" = p_pca))
}

