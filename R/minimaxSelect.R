#' Iteratively add candidate plots using the minimax algorithm
#' 
#' @inheritParams meanminSelect
#' @param cost_tol numeric value between 0 and 1, used if `r_cost` is
#'     supplied. Each plot is placed at the cheapest candidate location whose
#'     distance to its nearest plot is at least `1 - cost_tol` times the
#'     largest distance among candidates.
#'
#' @details 
#' The minimax algorithm aims to minimise the maximum distance between
#'     proposed plots and existing plots, iteratively placing new plots in
#'     locations with structural attributes most dissimilar to existing plots.
#'     Each plot is represented by its mean structural values: existing plots
#'     by `p_pca`, and new plots by the mean of the pixels within their
#'     footprint. For each candidate footprint, it computes the distance to the
#'     nearest plot, then chooses the candidate with the highest distance
#'     value. If there are no existing plots, the first plot is placed at the
#'     candidate most dissimilar to the landscape mean. As a result, plots
#'     occupy structural extremes.
#'
#' @return list of `sf` polygons for proposed new plots. 
#' 
#' @import terra
#' @import sf
#' 
#' @export
#' 
minimaxSelect <- function(r_pca, p_pca = NULL, old_ind = NULL, new_ind = NULL, 
  n_plots, p_new_dim, het_q = NULL, min_dist = NULL, r_cost = NULL, 
  cost_tol = 0.1) { 

  p_list <- list()

  ctx <- selectContext(r_pca, p_pca, new_ind, p_new_dim, het_q, min_dist, r_cost)
  checkCostTol(cost_tol)
  v_fp <- ctx$v_fp

  # Current set of plots in structural space
  current_p <- ctx$p_pca

  # Selection loop
  for (i in seq_len(n_plots)) { 
    message(i, "/", n_plots)
    
    # Ensure entire footprint of new plots falls within unoccupied pixels
    current_candidates <- ctxCandidates(r_pca, ctx, old_ind)

    if (length(current_candidates) == 0) {  
      message("No possible locations for plot ", i, "/", n_plots, ". Stopping ...")
      break
    }

    cand_fp <- v_fp[current_candidates, , drop = FALSE]

    # Calculate distance from each candidate to nearest plot
    if (!is.null(current_p)) {
      cand_dists <- nnDist(cand_fp, current_p)
    } else {
      # If no existing plots, use distance to landscape mean
      v_mean <- colMeans(ctx$pop)
      cand_dists <- distMat(cand_fp, matrix(v_mean, nrow = 1))[,1]
    }

    # Select candidate with highest distance (first occurrence if tied)
    acceptable <- cand_dists >= (1 - cost_tol) * max(cand_dists)
    best <- pickCandidate(-cand_dists, acceptable, ctx$cost[current_candidates])
    sel_center <- current_candidates[best]
    sel_id <- footprintCells(r_pca, ctx$fp, sel_center)

    # Save geometry and update old indices
    p_list[[i]] <- footprintPolygon(r_pca, sel_id)
    names(p_list)[[i]] <- paste(sel_id, collapse = ":")

    old_ind <- c(old_ind, sel_id)
    current_p <- rbind(current_p, v_fp[sel_center, , drop = FALSE])
  }

  # Return
  return(p_list)
}
