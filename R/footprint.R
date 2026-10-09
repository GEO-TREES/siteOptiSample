#' Build a focal window describing a plot footprint
#'
#' @param r `SpatRaster`
#' @param p_new_dim dimensions of new plots in the same coordinate system as `r`.
#'     A vector of two values.
#'
#' @return list containing: `w`: focal weights matrix, with 1 for cells in the
#'     footprint and NA elsewhere, `row_offsets` and `col_offsets`: row and
#'     column offsets of footprint cells relative to the focal centre, and
#'     `n_cells`: number of cells in the footprint
#'
#' @noRd
#'
plotFootprint <- function(r, p_new_dim) {
  res_x <- terra::res(r)[1]
  res_y <- terra::res(r)[2]

  p_x <- round(p_new_dim[1] / res_x)
  p_y <- round(p_new_dim[2] / res_y)

  n_x <- ceiling(p_x / 2)
  n_y <- ceiling(p_y / 2)
  raw_cols <- (2 * n_x) + 1
  raw_rows <- (2 * n_y) + 1
  max_cols <- max(1, (2 * terra::ncol(r)) - 1)
  max_rows <- max(1, (2 * terra::nrow(r)) - 1)
  full_cols <- min(raw_cols, max_cols)
  full_rows <- min(raw_rows, max_rows)

  pad_left_x  <- floor((full_cols - p_x) / 2)
  pad_right_x <- ceiling((full_cols - p_x) / 2)
  pad_left_y  <- floor((full_rows - p_y) / 2)
  pad_right_y <- ceiling((full_rows - p_y) / 2)

  w <- matrix(NA, nrow = full_rows, ncol = full_cols)
  w[(pad_left_y + 1):(full_rows - pad_right_y),
    (pad_left_x + 1):(full_cols - pad_right_x)] <- 1

  # Cell offsets of the footprint relative to the focal centre
  pos <- which(w == 1, arr.ind = TRUE)

  list(
    w = w,
    row_offsets = pos[,1] - ceiling(nrow(w) / 2),
    col_offsets = pos[,2] - ceiling(ncol(w) / 2),
    n_cells = nrow(pos))
}

#' Get cell IDs of a footprint centred on a cell
#'
#' @param r `SpatRaster`
#' @param fp footprint, as returned by `plotFootprint()`
#' @param center cell ID of the focal centre
#'
#' @noRd
#'
footprintCells <- function(r, fp, center) {
  rc <- terra::rowColFromCell(r, center)
  ids <- terra::cellFromRowCol(r, rc[1] + fp$row_offsets, rc[2] + fp$col_offsets)
  ids[!is.na(ids)]
}

#' Calculate mean values within the footprint centred on each cell
#'
#' @param r_pca `SpatRaster` of PCA scores or raw structural metrics
#' @param fp footprint, as returned by `plotFootprint()`
#'
#' @return matrix with one row per cell in `r_pca` and one column per layer.
#'     NA where any cell in the footprint is NA.
#'
#' @noRd
#'
footprintMeans <- function(r_pca, fp) {
  if (fp$n_cells == 1) {
    return(terra::values(r_pca))
  }
  terra::values(terra::focal(r_pca, w = fp$w, fun = "mean", na.rm = FALSE))
}

#' Identify cells whose footprint is not internally heterogeneous
#'
#' @param r_pca `SpatRaster` of PCA scores or raw structural metrics
#' @param fp footprint, as returned by `plotFootprint()`
#' @param base_candidates cell IDs of candidate plot locations
#' @param het_q optional, maximum quantile of total focal standard deviation
#'
#' @return vector of cell IDs
#'
#' @noRd
#'
hetSafeCenters <- function(r_pca, fp, base_candidates, het_q) {
  if (is.null(het_q)) {
    return(seq_len(terra::ncell(r_pca)))
  }

  if (fp$n_cells == 1) {
    message("Plot size is a single pixel. Internal heterogeneity is 0. Skipping heterogeneity filter...")
    return(seq_len(terra::ncell(r_pca)))
  }

  # Sum focal standard deviations across layers to get a Total Heterogeneity Index
  # Only complete footprints are used, so partial windows at edges don't
  # affect the cutoff
  r_sd <- terra::focal(r_pca, w = fp$w, fun = "sd", na.rm = FALSE)
  het <- terra::values(sum(r_sd))[,1]

  # Find the threshold value based ONLY on pixels inside the user's mask
  het_cutoff <- stats::quantile(het[base_candidates], probs = het_q, na.rm = TRUE)

  which(het <= het_cutoff)
}

#' Build a focal window marking cells within a minimum distance of a cell
#'
#' @param r `SpatRaster` 
#' @param min_dist minimum distance between the edges of plots, in the units 
#'     of `r`
#'
#' @return focal weights matrix, with 1 for cells whose edge is closer than
#'     `min_dist` to the edge of the central cell, and NA elsewhere
#'
#' @noRd
#' 
bufferWindow <- function(r, min_dist) {
  res <- terra::res(r)
  n_x <- ceiling(min_dist / res[1])
  n_y <- ceiling(min_dist / res[2])

  # Gaps between cell edges along each axis
  gap_x <- pmax(abs(-n_x:n_x) - 1, 0) * res[1]
  gap_y <- pmax(abs(-n_y:n_y) - 1, 0) * res[2]
  gap <- sqrt(outer(gap_y^2, gap_x^2, "+"))

  ifelse(gap < min_dist, 1, NA)
}

#' Identify cells where a plot footprint fits entirely within unoccupied
#'     candidate cells
#'
#' @param r `SpatRaster` 
#' @param fp footprint, as returned by `plotFootprint()`
#' @param base_candidates cell IDs of candidate plot locations
#' @param old_ind cell IDs of occupied locations
#' @param het_safe_centers cell IDs passing heterogeneity filter, as returned
#'     by `hetSafeCenters()`
#' @param buffer_w optional, focal window of cells too close to occupied
#'     cells, as returned by `bufferWindow()`
#'
#' @return vector of cell IDs
#'
#' @noRd
#' 
availCenters <- function(r, fp, base_candidates, old_ind, het_safe_centers, 
  buffer_w = NULL) {

  blocked <- old_ind

  # Block cells too close to occupied cells
  if (!is.null(buffer_w) && length(old_ind) > 0) {
    r_occ <- r[[1]]
    occ_vals <- rep(0, terra::ncell(r_occ))
    occ_vals[old_ind] <- 1
    terra::values(r_occ) <- occ_vals
    r_near <- terra::focal(r_occ, w = buffer_w, fun = "max", na.rm = TRUE)
    blocked <- union(blocked, which(terra::values(r_near)[,1] > 0))
  }

  avail_idx <- setdiff(base_candidates, blocked)

  r_avail <- r[[1]]
  avail_vals <- rep(0, terra::ncell(r_avail))
  avail_vals[avail_idx] <- 1
  terra::values(r_avail) <- avail_vals

  r_safe <- terra::focal(r_avail, w = fp$w, fun = "sum", na.rm = TRUE)
  safe_centers <- which(terra::values(r_safe)[,1] == fp$n_cells)

  out <- intersect(avail_idx, safe_centers)
  intersect(out, het_safe_centers)
}

#' Create a polygon from a set of raster cells
#'
#' @param r `SpatRaster`
#' @param ids cell IDs
#'
#' @return `sfc` polygon
#'
#' @noRd
#'
footprintPolygon <- function(r, ids) {
  r_sub <- r[[1]]
  terra::values(r_sub) <- NA
  r_sub[ids] <- 1
  p_sub <- terra::as.polygons(r_sub, dissolve = TRUE)
  sf::st_geometry(sf::st_as_sf(p_sub))
}

#' Prepare common inputs for selection algorithms
#'
#' @param r_pca `SpatRaster` of PCA scores or raw structural metrics
#' @param p_pca optional, structural values of existing plots
#' @param new_ind optional, cell IDs of candidate plot locations
#' @param p_new_dim dimensions of new plots
#' @param het_q optional, heterogeneity quantile threshold
#' @param min_dist optional, minimum distance between plot edges
#' @param r_cost optional, `SpatRaster` of access costs
#'
#' @return list containing: 
#'     * `p_pca`: matrix of existing plots with complete values, or NULL
#'     * `base_candidates`: candidate cell IDs
#'     * `fp`: footprint, as returned by `plotFootprint()`
#'     * `v_fp`: matrix of footprint means for each cell
#'     * `pop`: matrix of footprint means of all complete footprints in the 
#'       landscape, representing the population of possible plots
#'     * `het_safe`: cell IDs passing heterogeneity and cost filters
#'     * `buffer_w`: focal window for minimum distance, or NULL
#'     * `cost`: vector of mean cost within the footprint of each cell, or NULL
#'
#' @noRd
#' 
selectContext <- function(r_pca, p_pca, new_ind, p_new_dim, het_q = NULL, 
  min_dist = NULL, r_cost = NULL) {

  if (is.null(new_ind)) {
    new_ind <- which(stats::complete.cases(terra::values(r_pca)))
  }

  if (!is.null(p_pca)) {
    p_pca <- as.matrix(p_pca)
    p_complete <- stats::complete.cases(p_pca)
    if (!all(p_complete)) {
      warning(sum(!p_complete), " existing plot(s) have missing structural values and will be ignored")
    }
    p_pca <- p_pca[p_complete, , drop = FALSE]
    if (nrow(p_pca) == 0) {
      p_pca <- NULL
    }
  }

  fp <- plotFootprint(r_pca, p_new_dim)
  v_fp <- footprintMeans(r_pca, fp)
  pop <- v_fp[stats::complete.cases(v_fp), , drop = FALSE]
  het_safe <- hetSafeCenters(r_pca, fp, new_ind, het_q)

  buffer_w <- NULL
  if (!is.null(min_dist) && min_dist > 0) {
    buffer_w <- bufferWindow(r_pca, min_dist)
  }

  cost <- NULL
  if (!is.null(r_cost)) {
    if (!terra::compareGeom(r_pca, r_cost, stopOnError = FALSE)) {
      stop("`r_cost` must have the same extent and resolution as `r_pca`")
    }
    cost <- footprintMeans(r_cost[[1]], fp)[,1]
    # Candidates with unknown cost are excluded
    het_safe <- intersect(het_safe, which(!is.na(cost)))
  }

  list(p_pca = p_pca, base_candidates = new_ind, fp = fp, v_fp = v_fp, 
    pop = pop, het_safe = het_safe, buffer_w = buffer_w, cost = cost)
}

#' Find available candidate centres given occupied cells
#'
#' @param r_pca `SpatRaster` 
#' @param ctx context, as returned by `selectContext()`
#' @param old_ind cell IDs of occupied locations
#'
#' @noRd
#' 
ctxCandidates <- function(r_pca, ctx, old_ind) {
  availCenters(r_pca, ctx$fp, ctx$base_candidates, old_ind, ctx$het_safe, 
    ctx$buffer_w)
}

#' Choose a candidate, optionally trading off score against access cost
#'
#' @param score vector of candidate scores, where lower is better
#' @param acceptable logical vector, candidates whose score is close enough
#'     to the best to be considered when `cost` is supplied
#' @param cost optional, vector of candidate costs
#'
#' @return index of chosen candidate. Without `cost`, the candidate with the
#'     best score. With `cost`, the cheapest acceptable candidate, with ties
#'     broken by score.
#'
#' @noRd
#' 
pickCandidate <- function(score, acceptable, cost = NULL) {
  if (is.null(cost)) {
    return(which.min(score))
  }
  idx <- which(acceptable)
  idx[order(cost[idx], score[idx])[1]]
}

#' Calculate mean structural metrics within every possible plot footprint
#'
#' @param r `SpatRaster` of structural metrics or PCA scores
#' @param p_new_dim dimensions of plots in the same coordinate system as `r`.
#'     Either a single value for square plots, or a vector of two values for
#'     rectangular plots.
#'
#' @details
#' Plots are compared to the landscape at the scale of a plot, rather than
#'     the scale of a pixel. The landscape is represented by the mean values
#'     within every possible plot footprint (a moving window), which can then
#'     be compared to plots using e.g. `pcaDist()`. Footprints containing any
#'     missing values are excluded.
#'
#' @return matrix with one row per complete footprint and one column per
#'     layer in `r`
#'
#' @export
#' 
footprintMetrics <- function(r, p_new_dim) {
  if (!inherits(r, "SpatRaster")) {
    stop("`r` must be a SpatRaster")
  }
  fp <- plotFootprint(r, rep_len(p_new_dim, 2))
  v_fp <- footprintMeans(r, fp)
  colnames(v_fp) <- names(r)
  v_fp[stats::complete.cases(v_fp), , drop = FALSE]
}

#' Check cost tolerance is valid
#'
#' @param cost_tol cost tolerance
#'
#' @noRd
#' 
checkCostTol <- function(cost_tol) {
  if (!is.numeric(cost_tol) || length(cost_tol) != 1 || cost_tol < 0 || cost_tol >= 1) {
    stop("`cost_tol` must be a single numeric value >= 0 and < 1")
  }
}
