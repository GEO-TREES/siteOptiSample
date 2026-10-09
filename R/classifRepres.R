#' Classify pixels by whether they are well-represented by existing (and
#'     proposed) plots
#'
#' @details
#' Representativeness is measured relative to the centroid of the plots in
#'     PCA space. The cutoff is the `ci` quantile of the distances of the plots
#'     from their own centroid. Pixels further than this cutoff from the
#'     centroid are classified as poorly represented. If `p_new` is supplied,
#'     the centroid and cutoff are recalculated using both existing and
#'     proposed plots. 
#'
#' @return dataframe with one row per pixel, containing `dist`: euclidean
#'     distance from the centroid of existing plots, and `group`: the
#'     classification of each pixel
#'
#' @param r_pca PCA scores of pixels in structural space, as returned by
#'     `PCALandscape()$r_pca`
#' @param p PCA scores of candidate plots in structural space, as returned
#'     by `PCALandscape()$p_pca`
#' @param p_new optional PCA scores of candidate plots in structural space, as
#'     returned by `PCALandscape()$p_pca`
#' @param n_pca number of PCA axes used in analysis
#' @param ci quantile threshold for distances, e.g. 0.95 = 95th percentile of
#'     distances among plots 
#' 
#' @importFrom stats sd
#'
#' @export
#' 
classifRepres <- function(r_pca, p, p_new = NULL, 
  n_pca = 3, ci = 0.95) {

  # Extract PCA scores from chosen PCs
  r_df <- as.data.frame(r_pca)[,paste0("PC", 1:n_pca)]
  p <- p[,paste0("PC", 1:n_pca)]

  if (!is.null(p_new)) { 
    p_new <- p_new[,paste0("PC", 1:n_pca)]
  }

  # Compute centroid (mean vector) of the plots
  p_cent <- colMeans(p)

  if (!is.null(p_new)) { 
    p_all <- rbind(p, p_new)
    p_all_cent <- colMeans(p_all)
  }

  # Compute distances
  # Compute Euclidean distances for plots and pixels using existing and suggested plots
  p_dist <- sqrt(rowSums((p - matrix(p_cent, nrow(p), ncol(p), byrow = TRUE))^2))
  r_dist <- sqrt(rowSums((r_df - matrix(p_cent, nrow(r_df), ncol(r_df), byrow = TRUE))^2))

  if (!is.null(p_new)) { 
    p_all_dist <- sqrt(rowSums((p_all - matrix(p_all_cent, nrow(p_all), ncol(p_all), byrow = TRUE))^2))
    r_all_dist <- sqrt(rowSums((r_df - matrix(p_all_cent, nrow(r_df), ncol(r_df), byrow = TRUE))^2))
  }

  # Choose cutoff — e.g. 95th percentile of plot distances
  p_cutoff <- stats::quantile(p_dist, ci)

  if (!is.null(p_new)) { 
    p_all_cutoff <- stats::quantile(p_all_dist, ci)
  }

  # Flag pixels that are too far from plots
  px_out <- r_dist > p_cutoff

  if (!is.null(p_new)) { 
    px_all_out <- r_all_dist > p_all_cutoff
  }

  # Create output dataframe
  r_df$px_out <- px_out
  r_df$dist <- r_dist

  r_df$group <- ifelse(r_df$px_out, 
    "Poorly represented", "Well-represented by existing plots")

  if (!is.null(p_new)) {
    r_df$px_all_out <- px_all_out
    r_df$group[r_df$px_out & !r_df$px_all_out] <- 
      "Well-represented by existing and proposed plots"
  }
  r_df$group <- factor(r_df$group, 
    levels = c(
      "Poorly represented", 
      "Well-represented by existing plots", 
      "Well-represented by existing and proposed plots"))

  out <- r_df[,c("dist", "group")]

  # Return
  return(out)
}

