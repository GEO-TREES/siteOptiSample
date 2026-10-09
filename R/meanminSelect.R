#' Iteratively add candidate plots using the meanmin algorithm
#' 
#' @param r_pca `SpatRaster` of PCA scores or raw structural metrics. 
#' @param p_pca optional, PCA values or raw structural metrics of existing 
#'     plots in structural space
#' @param old_ind optional, pixel IDs in `r_pca` specifying existing plot locations. 
#' @param new_ind optional, pixel IDs in `r_pca` specifying candidate plot locations. 
#'     If not supplied, the full set of valid pixels in `r_pca` is used.
#' @param n_plots maximum number of new plots to add.
#' @param p_new_dim dimensions of new plots in the same coordinate system as `r_pca`. 
#'      A vector of two values. Should be perfectly divisible by the resolution 
#'      of `r_pca`. 
#' @param het_q optional, numeric value between 0 and 1. Defines a threshold 
#'     to filter out internally heterogeneous candidate locations, preventing plots 
#'     being placed on sharp structural transitions. Represents maximum 
#'     quantile for total focal standard deviation of metrics within 
#'     candidate footprint (e.g., 0.8 excludes top 20% most heterogeneous 
#'     areas). 
#' @param min_dist optional, minimum distance between the edges of a new plot
#'     and any other plot, new or existing, in the units of `r_pca`.
#' @param r_cost optional, `SpatRaster` with the same geometry as `r_pca`
#'     describing the cost of accessing each pixel, e.g. travel time. Candidate
#'     locations where any pixel has a missing cost are excluded. See
#'     `cost_tol`.
#' @param cost_tol numeric value between 0 and 1, used if `r_cost` is
#'     supplied. Each plot is placed at the cheapest candidate location whose
#'     score is within a tolerance of the best candidate's score. Larger
#'     values give more weight to cost. See Details for how the tolerance is
#'     defined.
#' @param refine logical, if TRUE, the greedy selection is followed by a
#'     refinement step which swaps selected plots for other candidate
#'     locations while doing so reduces the mean distance. See Details.
#'
#' @details
#' The mean-min algorithm aims to minimize the mean distance between every
#'      possible plot location in the landscape and its nearest plot, in
#'      structural/PCA space. Representativeness is assessed at the scale of a
#'      plot: the landscape is represented by the mean structural values
#'      within every possible plot footprint (see `footprintMetrics()`), and
#'      each plot is represented by its mean structural values: existing plots
#'      by `p_pca`, and new plots by the mean of the pixels within their
#'      footprint. The algorithm iteratively places new plots in locations
#'      that most reduce the mean distance. As a result, selected plots tend
#'      to concentrate in the most common/representative areas of structural
#'      space.
#'
#' Greedy selection is not guaranteed to find the set of plots with the lowest
#'      possible mean distance, because early plots are placed without
#'      knowledge of later plots. The first plot, for example, is placed at the
#'      most central location in structural space, which is rarely part of the
#'      best set of several plots. If `refine = TRUE`, after greedy selection
#'      each new plot is in turn removed and replaced with the candidate
#'      location which gives the lowest mean distance, given all other plots.
#'      This is repeated until no swap reduces the mean distance. Existing
#'      plots are never moved. The result is never worse than greedy
#'      selection, but takes longer to compute. Refined plots are returned in
#'      the order they were originally placed by greedy selection. If `r_cost`
#'      is supplied, a plot is only swapped to a location with an equal or
#'      lower cost.
#'
#' If `r_cost` is supplied, a candidate's score is the reduction in mean
#'      distance it would give, and candidates are accepted if they achieve at
#'      least `1 - cost_tol` times the best possible reduction, e.g. 0.1
#'      accepts candidates that achieve at least 90% of the best reduction.
#'
#' @return list of `sf` polygons for proposed new plots. 
#' 
#' @import terra
#' @import sf
#' 
#' @export
#'
meanminSelect <- function(r_pca, p_pca = NULL, old_ind = NULL, new_ind = NULL, 
  n_plots, p_new_dim, het_q = NULL, min_dist = NULL, r_cost = NULL, 
  cost_tol = 0.1, refine = FALSE) { 

  ctx <- selectContext(r_pca, p_pca, new_ind, p_new_dim, het_q, min_dist, r_cost)
  checkCostTol(cost_tol)
  v_fp <- ctx$v_fp
  pop <- ctx$pop

  # Find minimum distance from all possible plots to existing plots
  if (!is.null(ctx$p_pca)) {
    base_dists <- nnDist(pop, ctx$p_pca)
  } else {
    base_dists <- rep(Inf, nrow(pop))
  }

  # Mean distance to nearest plot after adding each candidate, processed in
  # chunks to limit memory use
  chunk_size <- max(1, floor(1e7 / nrow(pop)))
  meanAfter <- function(candidates, min_dists) {
    out <- numeric(length(candidates))
    for (s in seq(1, length(candidates), by = chunk_size)) {
      idx <- s:min(length(candidates), s + chunk_size - 1)
      cand_dist <- distMat(pop, v_fp[candidates[idx], , drop = FALSE])
      out[idx] <- colMeans(pmin(cand_dist, min_dists))
    }
    out
  }

  # Distance from all possible plots to a candidate plot
  candDist <- function(center) {
    distMat(pop, v_fp[center, , drop = FALSE])[,1]
  }

  sel_centers <- integer(0)
  sel_ids <- list()
  min_dists <- base_dists

  # Greedy selection loop
  for (i in seq_len(n_plots)) { 
    message(i, "/", n_plots)
    
    # Ensure entire footprint of new plots falls within unoccupied pixels
    current_candidates <- ctxCandidates(r_pca, ctx, c(old_ind, unlist(sel_ids)))

    if (length(current_candidates) == 0) { 
      message("No possible locations for plot ", i, "/", n_plots, ". Stopping ...")
      break
    }

    mean_after <- meanAfter(current_candidates, min_dists)

    # Candidates close enough to the best to be traded off against cost
    current_mean <- mean(min_dists)
    if (is.finite(current_mean)) {
      gain <- current_mean - mean_after
      if (max(gain) > 0) {
        acceptable <- gain >= (1 - cost_tol) * max(gain)
      } else {
        acceptable <- mean_after <= min(mean_after)
      }
    } else {
      acceptable <- mean_after <= min(mean_after) * (1 + cost_tol)
    }

    # Select candidate (first occurrence if tied)
    best <- pickCandidate(mean_after, acceptable, ctx$cost[current_candidates])
    sel_centers[i] <- current_candidates[best]
    sel_ids[[i]] <- footprintCells(r_pca, ctx$fp, sel_centers[i])

    # Update minimum distances with new plot
    min_dists <- pmin(min_dists, candDist(sel_centers[i]))
  }

  # Refinement by swapping
  if (refine && length(sel_centers) > 0) {
    sel_dists <- sapply(sel_centers, candDist)
    sel_dists <- matrix(sel_dists, nrow = nrow(pop))
    current_mean <- mean(pmin(base_dists, apply(sel_dists, 1, min)))
    n_pass <- 0
    max_pass <- 100

    repeat {
      n_pass <- n_pass + 1
      improved <- FALSE

      for (i in seq_along(sel_centers)) {
        # Minimum distances to all plots except plot i
        others <- sel_dists[, -i, drop = FALSE]
        min_others <- if (ncol(others) > 0) {
          pmin(base_dists, apply(others, 1, min))
        } else {
          base_dists
        }

        # Candidates which don't overlap existing plots or other new plots
        swap_candidates <- ctxCandidates(r_pca, ctx, 
          c(old_ind, unlist(sel_ids[-i])))

        # Only swap to locations which are no more expensive
        if (!is.null(ctx$cost)) {
          swap_candidates <- swap_candidates[
            ctx$cost[swap_candidates] <= ctx$cost[sel_centers[i]]]
        }

        if (length(swap_candidates) == 0) {
          next
        }

        mean_after <- meanAfter(swap_candidates, min_others)
        gain <- current_mean - mean_after
        best_gain <- max(gain)

        if (best_gain > 1e-12 * max(1, current_mean)) {
          acceptable <- gain >= (1 - cost_tol) * best_gain
          best <- pickCandidate(mean_after, acceptable, 
            ctx$cost[swap_candidates])
          sel_centers[i] <- swap_candidates[best]
          sel_ids[[i]] <- footprintCells(r_pca, ctx$fp, sel_centers[i])
          sel_dists[, i] <- candDist(sel_centers[i])
          current_mean <- mean_after[best]
          improved <- TRUE
        }
      }

      if (!improved) {
        break
      }

      if (n_pass >= max_pass) {
        message("Refinement stopped after ", max_pass, " passes without converging")
        break
      }
    }
  }

  # Generate plot polygons
  p_list <- lapply(sel_ids, function(ids) footprintPolygon(r_pca, ids))
  names(p_list) <- vapply(sel_ids, paste, character(1), collapse = ":")

  # Return
  return(p_list)
}
