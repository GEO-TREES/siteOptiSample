#' Iteratively add candidate plots using the meanmin algorithm
#' 
#' @param r_pca `SpatRaster` of PCA scores or raw structural metrics. 
#' @param p_pca optional, PCA values or raw structural metrics of existing 
#'     plots in structural space
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
#' The mean-min algorithm aims to minimize the mean distance of all candidate 
#'      locations (pixels or rows) to their nearest plot in structural/PCA space. 
#'      It iteratively places new plots in locations that most reduce the average 
#'      location-to-nearest-plot distance. For each candidate, it computes how 
#'      adding a plot there would decrease the overall mean distance, choosing 
#'      the candidate that yields the largest decrease. As a result, selected plots 
#'      tend to concentrate in the most common/representative areas of structural space.
#'
#' @return list of `sf` polygons for proposed new plots. 
#' 
#' @importFrom proxy dist
#' @import terra
#' @import sf
#' 
#' @export
#'
meanminSelect <- function(r_pca, p_pca, old_ind, new_ind, n_plots, p_new_dim, het_q) { 

  # Initialise empty list to store new plot polygons
  p_list <- list()

  # Extract data 
  v_pca <- terra::values(r_pca)
  valid_idx <- which(complete.cases(v_pca))
  v_pca_valid <- v_pca[valid_idx, , drop = FALSE]
  t_v_pca_valid <- t(v_pca_valid) 

  # Create template raster for focal operations
  r_cand <- r_pca[[1]] 
  terra::values(r_cand) <- NA_real_
  
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

  # Calculate how many pixels constitute a valid plot
  target_sum <- sum(w_matrix == 1, na.rm = TRUE)

  # Save mask boundaries outside the loop
  base_candidates <- new_ind

  # Heterogeneity filter 
  if (!is.null(het_q)) {
    if (target_sum == 1) {
      message("Plot size is a single pixel. Internal heterogeneity is 0. Skipping heterogeneity filter...")
      het_safe_centers <- seq_len(terra::ncell(r_pca))
    } else {
      # Calculate focal standard deviation for each PCA layer
      r_sd <- terra::focal(r_pca, w = w_matrix, fun = "sd", na.rm = TRUE)
      
      # Sum standard deviations across all PCA layers to get a Total Heterogeneity Index
      r_het <- sum(r_sd, na.rm = TRUE)
      
      # Find the threshold value based ONLY on pixels inside the user's mask
      het_vals <- terra::values(r_het)[base_candidates]
      het_cutoff <- quantile(het_vals, probs = het_q, na.rm = TRUE)
      
      # Identify all center pixels whose footprint falls below the variance threshold
      het_safe_centers <- which(terra::values(r_het) <= het_cutoff)
    }
  } else {
    het_safe_centers <- seq_len(terra::ncell(r_pca)) 
  }

  # Selection loop
  for (i in seq_len(n_plots)) { 
    message(i, "/", n_plots)
    
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
      message("No possible locations for plot ", i, "/", n_plots, ". Stopping ...")
      break
    }

    # Compile current set of plots in structural space
    valid_old_ind <- intersect(old_ind, valid_idx)
    current_p <- p_pca
    if (length(valid_old_ind) > 0) {
      current_p <- rbind(current_p, v_pca[valid_old_ind, , drop = FALSE])
    }

    # Find minimum distance to any existing plot for all pixels
    if (!is.null(current_p) && nrow(current_p) > 0) {
      min_dists <- apply(as.matrix(proxy::dist(v_pca_valid, current_p)), 1, min)
    } else {
      min_dists <- rep(Inf, length(valid_idx))
    }

    # Calculate reduction of mean distance 
    mean_after <- sapply(current_candidates, function(j) {
      cand_vec <- v_pca[j, ]
      cand_dist <- sqrt(colSums((t_v_pca_valid - cand_vec)^2))
      new_min <- pmin(min_dists, cand_dist)
      mean(new_min, na.rm = TRUE) 
    })
    
    # Populate template raster and run focal smoothing
    r_cand_vals <- rep(NA_real_, terra::ncell(r_cand))
    r_cand_vals[current_candidates] <- mean_after
    terra::values(r_cand) <- r_cand_vals
    
    # Run focal operation
    r_rs <- terra::focal(r_cand, w = w_matrix, fun = "mean", na.rm = TRUE)

    # Mask invalid centres
    rs_vals <- as.numeric(terra::values(r_rs))
    mask_idx <- setdiff(seq_along(rs_vals), current_candidates)
    rs_vals[mask_idx] <- Inf

    if (all(!is.finite(rs_vals))) {
      message("No possible locations for plot ", i, "/", n_plots, ". Stopping ...")
    break
    }

    # Extract central coordinate of selected plot (first occurrence if tied)
    min_val <- min(rs_vals, na.rm = TRUE)
    which_min_val <- which(rs_vals == min_val)[1] 
    
    rc <- terra::rowColFromCell(r_rs, which_min_val)
    sel_rows <- rc[1] + row_offsets
    sel_cols <- rc[2] + col_offsets

    sel_id <- terra::cellFromRowCol(r_rs, sel_rows, sel_cols)
    sel_id <- sel_id[!is.na(sel_id)]

    # Generate polygon of selected plot directly using raster cells
    r_sub <- r_pca[[1]]
    terra::values(r_sub) <- NA
    r_sub[sel_id] <- 1
    p_sub <- terra::as.polygons(r_sub, dissolve = TRUE)
    
    # Save geometry and update old indices
    p_list[[i]] <- sf::st_geometry(sf::st_as_sf(p_sub))
    old_ind <- c(old_ind, sel_id)
  }

  # Return
  return(p_list)
}
