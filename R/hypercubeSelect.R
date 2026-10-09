#' Select candidate plots using Multi-dimensional Quantiles (Latin Hypercube)
#' 
#' @inheritParams meanminSelect
#' 
#' @details 
#' The Latin Hypercube algorithm aims to capture the full multi-dimensional
#'      structural gradient of the landscape evenly. The landscape is
#'      represented by the mean structural values within every possible plot
#'      footprint (see `footprintMetrics()`). The distribution of each
#'      structural metric (or principal component) is divided into
#'      equal-probability strata, one per plot, including existing plots.
#'      Existing plots occupy some of these strata, and new plots are targeted
#'      at the empty strata, so that together the existing and new plots
#'      cover each gradient as evenly as possible. Target values are randomly
#'      combined across metrics to generate theoretical target coordinates in
#'      the multi-dimensional feature space. For each target, it selects the
#'      available candidate location whose mean structural attributes are
#'      closest to the target. As a result, the proposed plots are spread
#'      across the entire range of every structural gradient, ensuring that
#'      both average conditions and rare structural combinations are sampled
#'      representatively.
#'
#' If `r_cost` is supplied, candidates are accepted if their distance to
#'      their target is no more than the distance of the closest candidate
#'      plus `cost_tol` times the mean distance between the target and the
#'      landscape locations closer to it than to any other target or existing
#'      plot.
#' 
#' @return list of `sf` polygons for proposed new plots. 
#' 
#' @import terra
#' @import sf
#' 
#' @export
#'
hypercubeSelect <- function(r_pca, p_pca = NULL, old_ind = NULL, new_ind = NULL, 
  n_plots, p_new_dim, het_q = 0.8, min_dist = NULL, r_cost = NULL, 
  cost_tol = 0.1) { 

  p_list <- list()

  ctx <- selectContext(r_pca, p_pca, new_ind, p_new_dim, het_q, min_dist, r_cost)
  checkCostTol(cost_tol)
  v_fp <- ctx$v_fp
  pop <- ctx$pop

  # Divide each dimension into equal-probability strata, one per plot
  n_exist <- if (is.null(ctx$p_pca)) 0 else nrow(ctx$p_pca)
  n_total <- n_exist + n_plots
  breaks_probs <- seq(0, 1, length.out = n_total + 1)

  target_matrix <- matrix(NA_real_, nrow = n_plots, ncol = ncol(pop))
  for (d in seq_len(ncol(pop))) {
    breaks <- stats::quantile(pop[, d], probs = breaks_probs, names = FALSE)

    # Find strata not already occupied by existing plots
    occupied <- if (n_exist > 0) {
      findInterval(ctx$p_pca[, d], breaks, all.inside = TRUE)
    } else {
      integer(0)
    }
    empty <- setdiff(seq_len(n_total), occupied)
    chosen <- empty[sample.int(length(empty), n_plots)]

    # Target the middle of each chosen stratum, in random order
    target_vals <- stats::quantile(pop[, d], probs = (chosen - 0.5) / n_total, 
      names = FALSE)
    target_matrix[, d] <- target_vals
  }

  # Mean distance from each target to the landscape locations nearest to it
  d_targets <- distMat(pop, rbind(ctx$p_pca, target_matrix))
  assign <- max.col(-d_targets, ties.method = "first") - n_exist
  min_d <- apply(d_targets, 1, min)
  radius <- vapply(seq_len(n_plots), function(j) {
    members <- which(assign == j)
    if (length(members) > 0) mean(min_d[members]) else mean(min_d)
  }, numeric(1))

  # Selection loop
  for (i in seq_len(n_plots)) { 
    current_candidates <- ctxCandidates(r_pca, ctx, old_ind)

    if (length(current_candidates) == 0) {
      message("No possible locations remaining for plot ", i, ". Stopping ...")
      break
    }

    dists <- distMat(v_fp[current_candidates, , drop = FALSE], 
      target_matrix[i, , drop = FALSE])[,1]
    acceptable <- dists <= min(dists) + cost_tol * radius[i]
    best <- pickCandidate(dists, acceptable, ctx$cost[current_candidates])
    sel_center <- current_candidates[best]
    sel_id <- footprintCells(r_pca, ctx$fp, sel_center)

    # Generate polygon
    p_list[[i]] <- footprintPolygon(r_pca, sel_id)
    names(p_list)[[i]] <- paste(sel_id, collapse = ":")
    old_ind <- c(old_ind, sel_id)
  }

  return(p_list)
}
