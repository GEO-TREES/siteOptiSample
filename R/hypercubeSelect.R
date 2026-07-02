#' Select candidate plots using Multi-dimensional Quantiles (Latin Hypercube)
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
#' The Conditioned Latin Hypercube Sampling (cLHS) algorithm aims to capture the 
#'      full multi-dimensional structural gradient of the landscape evenly.
#'      It divides the cumulative distribution of each structural metric (or principal
#'      component) into equal-probability quantiles, corresponding to the
#'      number of proposed plots. It then randomly combines these marginal
#'      quantiles to generate a set of theoretical target coordinates in the
#'      multi-dimensional feature space. For each theoretical target, it
#'      iteratively evaluates the physically available candidate pixels and
#'      selects the location whose structural attributes minimize the Euclidean
#'      distance to the target. As a result, the proposed plots are forced to
#'      spread evenly across the entire range of every structural gradient,
#'      ensuring that both average conditions and rare structural combinations
#'      are sampled representatively.
#' 
#' @return list of `sf` polygons for proposed new plots. 
#' 
#' @import terra
#' @import sf
#' @importFrom stats quantile
#' 
#' @export
#'
hypercubeSelect <- function(r_pca, old_ind, new_ind, n_plots, p_new_dim, het_q = 0.8) { 

  p_list <- list()
  res <- terra::res(r_pca)
  
  # 1. Precise Window Generation
  p_n_cells_x <- round(p_new_dim[1] / res[1])
  p_n_cells_y <- round(p_new_dim[2] / res[2])
  
  # Ensure focal window is odd-sized for terra::focal
  focal_rows <- if (p_n_cells_y %% 2 == 0) p_n_cells_y + 1 else p_n_cells_y
  focal_cols <- if (p_n_cells_x %% 2 == 0) p_n_cells_x + 1 else p_n_cells_x
  
  # Create weight matrix to mask padding pixels
  weights <- matrix(0, nrow = focal_rows, ncol = focal_cols)
  y_start <- floor((focal_rows - p_n_cells_y) / 2) + 1
  x_start <- floor((focal_cols - p_n_cells_x) / 2) + 1
  weights[y_start:(y_start + p_n_cells_y - 1), x_start:(x_start + p_n_cells_x - 1)] <- 1
  
  # Offsets for polygon generation
  row_offsets <- rep(seq(0, p_n_cells_y - 1), times = p_n_cells_x) - floor(p_n_cells_y / 2)
  col_offsets <- rep(seq(0, p_n_cells_x - 1), each = p_n_cells_y) - floor(p_n_cells_x / 2)

  # 2. Candidate Pool Filtering
  v_pca <- terra::values(r_pca)
  base_candidates <- new_ind
  
  if (!is.null(het_q)) {
    r_sd <- terra::focal(r_pca, w = weights, fun = "sd", na.rm = TRUE)
    r_het <- sum(r_sd, na.rm = TRUE)
    het_cutoff <- quantile(terra::values(r_het)[base_candidates], probs = het_q, na.rm = TRUE)
    het_safe_centers <- which(terra::values(r_het) <= het_cutoff)
  } else {
    het_safe_centers <- seq_len(terra::ncell(r_pca))
  }

  # 3. LHS Target Generation
  n_dim <- ncol(v_pca)
  target_matrix <- matrix(NA, nrow = n_plots, ncol = n_dim)
  quant_probs <- (seq_len(n_plots) - 0.5) / n_plots 
  
  # Use full valid space for quantile calculation
  for (d in seq_len(n_dim)) {
    target_vals <- quantile(v_pca[base_candidates, d], probs = quant_probs, na.rm = TRUE)
    target_matrix[, d] <- sample(target_vals)
  }

  # 4. Selection Loop
  for (i in seq_len(n_plots)) { 
    r_avail <- r_pca[[1]]
    terra::values(r_avail) <- 0
    terra::values(r_avail)[setdiff(base_candidates, old_ind)] <- 1
    
    r_safe <- terra::focal(r_avail, w = weights, fun = "sum", na.rm = TRUE)
    current_candidates <- intersect(which(terra::values(r_safe) == sum(weights)), base_candidates)
    current_candidates <- intersect(current_candidates, het_safe_centers)
    current_candidates <- setdiff(current_candidates, old_ind)

    if (length(current_candidates) == 0) break

    dists <- sqrt(colSums((t(v_pca[current_candidates, , drop = FALSE]) - target_matrix[i, ])^2))
    sel_center <- current_candidates[which.min(dists)]

    # Generate polygon
    rc <- terra::rowColFromCell(r_pca, sel_center)
    sel_rows <- rc[1] + row_offsets
    sel_cols <- rc[2] + col_offsets
    
    valid_idx <- sel_rows >= 1 & sel_rows <= terra::nrow(r_pca) & 
                 sel_cols >= 1 & sel_cols <= terra::ncol(r_pca)
    sel_id <- terra::cellFromRowCol(r_pca, sel_rows[valid_idx], sel_cols[valid_idx])
    
    r_sub <- r_pca[[1]]
    terra::values(r_sub) <- NA
    r_sub[sel_id] <- 1
    p_sub <- terra::as.polygons(r_sub, dissolve = TRUE)
    
    p_list[[i]] <- sf::st_geometry(sf::st_as_sf(p_sub))
    names(p_list)[[i]] <- paste(sel_id, collapse = ":")
    old_ind <- c(old_ind, sel_id)
  }

  return(p_list)
}
