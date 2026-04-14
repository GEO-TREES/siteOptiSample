#' Select candidate plots using K-Means clustering 
#' 
#' @param r_pca `SpatRaster` of PCA scores or raw structural metrics. 
#' @param old_ind optional, pixel IDs in `r_pca` specifying existing plot locations. 
#' @param new_ind optional, pixel IDs in `r_pca` specifying candidate plot locations. 
#'     If not supplied, the full set of valid pixels in `r_pca` is used.
#' @param n_plots maximum number of new plots to add.
#' @param p_new_dim dimensions of new plots in the same coordinate system as `r`. 
#'      A vector of two values. Should be perfectly divisible by the resolution 
#'      of `r_pca`. 
#' @param het_q optional, numeric value between 0 and 1. Defines a threshold 
#'     to filter out internally heterogeneous candidate locations, preventing plots 
#'     being placed on sharp structural transitions. Represents maximum 
#'     quantile for total focal standard deviation of metrics within 
#'     candidate footprint (e.g., 0.8 excludes top 20% most heterogeneous 
#'     areas). 
#' 
#' @details 
#' The K-means algorithm aims to capture the full structural diversity of the
#'      landscape through stratified sampling. It first partitions available
#'      candidate pixels into structural clusters, equal to the number
#'      of proposed plots, to identify "ideal" representative centroids in the
#'      feature space. For each cluster, it iteratively evaluates the available
#'      candidate pixels and selects the location whose structural attributes are
#'      closest to the cluster's centroid. As a result, rather than pushing plots
#'      toward structural extremes, the proposed plots are distributed
#'      representatively across the dominant structural conditions of the landscape.
#' 
#' @return list of `sf` polygons for proposed new plots. 
#' 
#' @import terra
#' @import sf
#' @importFrom stats kmeans
#' 
#' @export
#' 
kmeansSelect <- function(r_pca, old_ind, new_ind, n_plots, p_new_dim, het_q = 0.8) { 

  p_list <- list()

  # Extract data 
  v_pca <- terra::values(r_pca)
  valid_idx <- which(complete.cases(v_pca))
  
  # Build focal window 
  res_x <- terra::res(r_pca)[1]
  res_y <- terra::res(r_pca)[2]

  p_x <- round(p_new_dim[1] / res_x)
  p_y <- round(p_new_dim[2] / res_y)

  n_x <- ceiling(p_x / 2)
  n_y <- ceiling(p_y / 2)
  full_cols <- (2 * n_x) + 1 
  full_rows <- (2 * n_y) + 1

  pad_left_x  <- floor((full_cols - p_x) / 2)
  pad_right_x <- ceiling((full_cols - p_x) / 2)
  pad_left_y  <- floor((full_rows - p_y) / 2)
  pad_right_y <- ceiling((full_rows - p_y) / 2)

  w_matrix <- matrix(NA, nrow = full_rows, ncol = full_cols)
  w_matrix[(pad_left_y + 1):(full_rows - pad_right_y), 
           (pad_left_x + 1):(full_cols - pad_right_x)] <- 1

  # Pre-compute cell offsets for the selected window
  pos <- which(w_matrix == 1, arr.ind = TRUE)
  row_offsets <- pos[,1] - ceiling(nrow(w_matrix) / 2)
  col_offsets <- pos[,2] - ceiling(ncol(w_matrix) / 2)
  
  # Calculate exactly how many pixels constitute a complete, valid plot
  target_sum <- sum(w_matrix == 1, na.rm = TRUE)

  # Save mask boundaries outside the loop
  base_candidates <- new_ind

  # Heterogeneity filter 
  if (!is.null(het_q)) {
    if (target_sum == 1) {
      message("Plot size is a single pixel. Internal heterogeneity is 0. Skipping heterogeneity filter...")
      het_safe_centers <- seq_len(terra::ncell(r_pca))
    } else {
      r_sd <- terra::focal(r_pca, w = w_matrix, fun = "sd", na.rm = TRUE)
      r_het <- sum(r_sd, na.rm = TRUE)
      
      het_vals <- terra::values(r_het)[base_candidates]
      het_cutoff <- quantile(het_vals, probs = het_q, na.rm = TRUE)
      het_safe_centers <- which(terra::values(r_het) <= het_cutoff)
    }
  } else {
    het_safe_centers <- seq_len(terra::ncell(r_pca))
  }

  # Find initial pool of pixels to cluster
  r_avail <- r_pca[[1]]
  terra::values(r_avail) <- 0
  initial_avail <- setdiff(base_candidates, old_ind)
  terra::values(r_avail)[initial_avail] <- 1
  
  r_safe <- terra::focal(r_avail, w = w_matrix, fun = "sum", na.rm = TRUE)
  initial_safe <- which(terra::values(r_safe) == target_sum)
  initial_candidates <- intersect(initial_avail, initial_safe)
  initial_candidates <- intersect(initial_candidates, het_safe_centers)

  if (length(initial_candidates) < n_plots) {
    message("No possible locations for plot ", i, "/", n_plots, ". Stopping ...")
    break
  }

  # Run K-means on the valid candidate space
  km <- stats::kmeans(v_pca[initial_candidates, , drop = FALSE], centers = n_plots, iter.max = 100, nstart = 10)
  centroids <- km$centers

  # Selection loop
  for (i in seq_len(n_plots)) { 
    
    # Ensure entire footprint of new plots falls within unoccupied pixels
    r_avail <- r_pca[[1]]
    terra::values(r_avail) <- 0
    avail_idx <- setdiff(base_candidates, old_ind)
    terra::values(r_avail)[avail_idx] <- 1
    
    r_safe <- terra::focal(r_avail, w = w_matrix, fun = "sum", na.rm = TRUE)
    safe_centers <- which(terra::values(r_safe) == target_sum)
    
    # Restrict candidate centers to those with safe footprints
    current_candidates <- intersect(avail_idx, safe_centers)
    current_candidates <- intersect(current_candidates, het_safe_centers)

    if (length(current_candidates) == 0) { 
      message("No possible locations remaining for plot ", i, ". Stopping ...")
      break
    }

    # Extract PCA values for the currently valid candidates
    cand_pca <- v_pca[current_candidates, , drop = FALSE]
    target_centroid <- centroids[i, ]

    # Calculate distance between candidates and the target cluster centroid
    # Using fast vector recycling matrix math
    dists <- sqrt(colSums((t(cand_pca) - target_centroid)^2))

    # Pick the pixel closest to the centroid
    best_idx <- which.min(dists)
    sel_center <- current_candidates[best_idx]

    # Generate polygon coordinates
    rc <- terra::rowColFromCell(r_pca, sel_center)
    sel_rows <- rc[1] + row_offsets
    sel_cols <- rc[2] + col_offsets

    sel_id <- terra::cellFromRowCol(r_pca, sel_rows, sel_cols)
    sel_id <- sel_id[!is.na(sel_id)]

    # Generate polygon
    r_sub <- r_pca[[1]]
    terra::values(r_sub) <- NA
    r_sub[sel_id] <- 1
    p_sub <- terra::as.polygons(r_sub, dissolve = TRUE)
    
    p_list[[i]] <- sf::st_geometry(sf::st_as_sf(p_sub))
    
    # Update old_ind so the next plot can't overlap it
    old_ind <- c(old_ind, sel_id)
  }

  return(p_list)
}
