#' Iteratively add candidate plots using the minimax algorithm
#' 
#' @param r_pca `SpatRaster`, `data.frame`, or `matrix` of PCA scores or raw 
#'     structural metrics. If spatial, values represent pixels across a landscape. 
#'     If non-spatial, each row represents a candidate plot 
#' @param p_pca optional, PCA values or raw structural metrics of existing 
#'     plots in structural space
#' @param old_ind optional, pixel IDs (if spatial) or row indices (if non-spatial) 
#'      in `r_pca` specifying existing plot locations. 
#' @param new_ind optiona, pixel IDs (if spatial) or row indices (if non-spatial) 
#'      in `r_pca` specifying candidate plot locations. If not supplied, the full 
#'      set of valid landscape pixels or dataframe rows in `r_pca` is used.
#' @param n_plots maximum number of new plots to add.
#' @param p_new_dim dimensions of new plots in the same coordinate system as `r`. 
#'      A vector of two values. Should be perfectly divisible by the resolution 
#'      of `r_pca`. **Note:** Ignored if `r_pca` is a `data.frame` or `matrix`.
#'
#' @details 
#' The minimax algorithm aims to minimise the maximum distance between
#'     proposed plots and existing plots, iteratively placing new plots in
#'     pixels with structural attributes most dissimilar to existing plots. For
#'     each candidate pixel, it computes the distance to the nearest existing plot,
#'     then chooses the pixel with the highest distance value. As a result, plots
#'     occupy structural extremes.
#'
#' @return If `r_pca` is a `SpatRaster`, returns a list of `sf` polygons for the 
#'      proposed new plots. If `r_pca` is a `data.frame` or `matrix`, returns a 
#'      numeric vector of row indices corresponding to the selected plots.
#' 
#' @importFrom proxy dist
#' @import terra
#' @import sf
#' 
#' @export
#' 
minimaxSelect <- function(r_pca, p_pca, old_ind, new_ind, n_plots, p_new_dim) { 

  # Initialise empty list to store new plot polygons
  p_list <- list()

  if (!inherits(r_pca, "SpatRaster")) {
    r_mat <- as.matrix(r_pca)
    
    if (!is.null(p_pca)) {
      p_mat_ext <- as.matrix(p_pca)
    } else {
      p_mat_ext <- NULL
    } 

    # Selection Loop
    for (i in seq_len(n_plots)) { 
      message(i, "/", n_plots)
      
      # Exclude already selected plots from candidates
      new_ind <- setdiff(new_ind, old_ind)
      
      if (length(new_ind) == 0) { 
        message("No possible locations for plot ", i, "/", n_plots, ". Stopping ...")
        break
      }
      
      # Compile the current set of all "existing" plots
      current_p_mat <- p_mat_ext
      if (length(old_ind) > 0) {
        current_p_mat <- rbind(current_p_mat, r_mat[old_ind, , drop = FALSE])
      }
      
      if (is.null(current_p_mat) || nrow(current_p_mat) == 0) {
        # If there are absolutely no existing plots, pick the first available candidate
        sel_id <- new_ind[1]
      } else {
        # Calculate minimum distance from all remaining candidates to existing plots
        cand_mat <- r_mat[new_ind, , drop = FALSE]
        min_dists <- apply(as.matrix(proxy::dist(cand_mat, current_p_mat)), 1, min)
        
        # Pick the candidate that has the MAXIMUM minimum distance
        best_idx <- which.max(min_dists)
        sel_id <- new_ind[best_idx]
      }
      
      # Store result and update existing internal plots for the next iteration
      p_list <- c(p_list, sel_id)
      old_ind <- c(old_ind, sel_id)
    }
    
    return(p_list)
  } else {

    # Safely store external existing plots for the spatial loop
    p_mat_ext <- NULL
    if (!is.null(p_pca)) {
      p_mat_ext <- as.matrix(p_pca)
    }

    # For each plot 
    for (i in seq_len(n_plots)) { 
      message(i, "/", n_plots)
      # Get index values in PCA for candidate plots
      # Excluding pixels occupied by existing plots 
      new_ind <- setdiff(new_ind, old_ind)

      # Stop if no potential locations available
      if (length(new_ind) == 0) { 
        message("No possible locations for plot ", i, "/", n_plots, ". Stopping ...")
        break
      }

      # Compile the current set of all "existing" plots
      current_p_mat <- p_mat_ext
      if (length(old_ind) > 0) {
        current_p_mat <- rbind(current_p_mat, terra::values(r_pca)[old_ind, , drop = FALSE])
      }

      # Initialise minimum distances between each pixel and existing plots
      r_cand <- r_pca[[1]]
      if (!is.null(current_p_mat) && nrow(current_p_mat) > 0) {
        terra::values(r_cand) <- apply(as.matrix(
          proxy::dist(terra::values(r_pca), current_p_mat)), 1, min)
      } else {
        terra::values(r_cand)[complete.cases(terra::values(r_cand))] <- Inf
      }

      # Create a moving window to compute mean distance of other subplots
      p_width_x <- p_new_dim[1]
      p_height_y <- p_new_dim[2]
      res_x <- terra::res(r_pca)[1]
      res_y <- terra::res(r_pca)[2]

      # Number of cells needed in half-window
      n_x <- ceiling((p_width_x / res_x) / 2)
      n_y <- ceiling((p_height_y / res_y) / 2)

      # Full window dimensions
      full_cols <- (2 * n_x) + 1 
      full_rows <- (2 * n_y) + 1

      # Size of the new plot in raster cells
      p_x <- round(p_new_dim[1] / res_x)
      p_y <- round(p_new_dim[2] / res_y)
      P <- c(p_x, p_y)

      # Compute symmetric padding (handles odd differences)
      pad_left_x  <- floor((full_cols - P[1]) / 2)
      pad_right_x <- ceiling((full_cols - P[1]) / 2)
      pad_left_y  <- floor((full_rows - P[2]) / 2)
      pad_right_y <- ceiling((full_rows - P[2]) / 2)

      # Create centered window
      w_matrix <- matrix(NA, nrow = full_rows, ncol = full_cols)

      # Placement indices
      start_row <- pad_left_y + 1
      end_row   <- full_rows - pad_right_y
      start_col <- pad_left_x + 1
      end_col   <- full_cols - pad_right_x

      # Insert centered block
      w_matrix[start_row:end_row, start_col:end_col] <- 1

      # Run focal operation
      r_rs <- terra::focal(r_cand, w = w_matrix, fun = "mean")

      # Find pixel with largest value (Minimax behavior)
      max_val <- max(terra::values(r_rs, na.rm = TRUE))
      which_max_val <- which(terra::values(r_rs) == max_val)

      # Extract coordinates of selected plot
      rc <- terra::rowColFromCell(r_rs, which_max_val)
      row0 <- rc[1]
      col0 <- rc[2]

      # Dimensions of window
      h <- nrow(w_matrix)
      w <- ncol(w_matrix)
      pos <- which(w_matrix == 1, arr.ind = TRUE)

      # The cell in the window that corresponds to the focal raster cell:
      # It is ALWAYS the middle cell of the full window matrix
      row_center <- ceiling(h / 2)
      col_center <- ceiling(w / 2)

      # Compute row/col offsets FROM THE CENTER
      row_offsets <- pos[,1] - row_center
      col_offsets <- pos[,2] - col_center

      # Apply offsets to the focal cell
      sel_rows <- row0 + row_offsets
      sel_cols <- col0 + col_offsets

      # Convert to raster cell numbers
      sel_id <- terra::cellFromRowCol(r_rs, sel_rows, sel_cols)
      sel_id <- sel_id[!is.na(sel_id)]

      # Extract polygon of selected plot
      r_sub <- r_pca[[1]]
      terra::values(r_sub) <- NA
      r_sub[sel_id] <- 1
      p_sub <- terra::as.polygons(r_sub, dissolve = TRUE)
      p_sf <- sf::st_geometry(sf::st_as_sf(p_sub))

      # Add polygon of selected plot to list
      p_list[[i]] <- p_sf

      # Update cell IDs of selected plots
      old_ind <- c(old_ind, sel_id)
    }

    # Return
    return(p_list)
  }
}
