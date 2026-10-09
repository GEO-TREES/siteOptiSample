#' Select candidate plots using K-Means clustering 
#' 
#' @inheritParams meanminSelect
#' 
#' @details 
#' The K-means algorithm aims to capture the full structural diversity of the
#'      landscape through stratified sampling. The landscape is represented by
#'      the mean structural values within every possible plot footprint (see
#'      `footprintMetrics()`). These are partitioned into clusters, with one
#'      cluster per existing plot and one cluster per proposed plot. The
#'      centroids of the existing plots' clusters are fixed at the existing
#'      plots' structural values, so new clusters form in parts of
#'      structural space not already covered by existing plots. For each new
#'      cluster, starting with the largest, it selects the available
#'      candidate location whose mean structural attributes are closest to
#'      the cluster's centroid. As a result, rather than pushing plots toward
#'      structural extremes, the proposed plots are distributed
#'      representatively across the dominant structural conditions of the
#'      landscape.
#'
#' If `r_cost` is supplied, candidates are accepted if their distance to the
#'      cluster centroid is no more than the distance of the closest candidate
#'      plus `cost_tol` times the mean distance between the centroid and the
#'      members of its cluster.
#' 
#' @return list of `sf` polygons for proposed new plots. 
#' 
#' @import terra
#' @import sf
#' 
#' @export
#' 
kmeansSelect <- function(r_pca, p_pca = NULL, old_ind = NULL, new_ind = NULL, 
  n_plots, p_new_dim, het_q = 0.8, min_dist = NULL, r_cost = NULL, 
  cost_tol = 0.1) { 

  p_list <- list()

  ctx <- selectContext(r_pca, p_pca, new_ind, p_new_dim, het_q, min_dist, r_cost)
  checkCostTol(cost_tol)
  v_fp <- ctx$v_fp

  n_distinct <- nrow(unique(ctx$pop))
  if (n_plots > n_distinct) {
    message("Only ", n_distinct, " distinct plot locations in the landscape. ",
      "Reducing number of plots to ", n_distinct, " ...")
    n_plots <- n_distinct
  }

  if (n_plots == 0) {
    return(p_list)
  }

  # Cluster the landscape, holding existing plots fixed as centroids
  km <- fixedKmeans(ctx$pop, ctx$p_pca, n_plots)
  centroids <- km$centers
  radius <- km$radius

  # Select plots for the largest clusters first
  cluster_order <- order(km$size, decreasing = TRUE)

  # Selection loop
  for (i in seq_along(cluster_order)) { 
    j <- cluster_order[i]
    
    # Ensure entire footprint of new plots falls within unoccupied pixels
    current_candidates <- ctxCandidates(r_pca, ctx, old_ind)

    if (length(current_candidates) == 0) { 
      message("No possible locations remaining for plot ", i, ". Stopping ...")
      break
    }

    # Pick the candidate closest to the target cluster centroid
    dists <- distMat(v_fp[current_candidates, , drop = FALSE], 
      centroids[j, , drop = FALSE])[,1]
    acceptable <- dists <= min(dists) + cost_tol * radius[j]
    best <- pickCandidate(dists, acceptable, ctx$cost[current_candidates])
    sel_center <- current_candidates[best]
    sel_id <- footprintCells(r_pca, ctx$fp, sel_center)

    # Generate polygon
    p_list[[i]] <- footprintPolygon(r_pca, sel_id)
    names(p_list)[[i]] <- paste(sel_id, collapse = ":")
    
    # Update old_ind so the next plot can't overlap it
    old_ind <- c(old_ind, sel_id)
  }

  return(p_list)
}

#' K-means clustering with some centroids held fixed
#'
#' @param x numeric matrix of observations
#' @param fixed optional, numeric matrix of fixed centroids
#' @param k number of free centroids
#' @param nstart number of random starts
#' @param iter_max maximum number of iterations per start
#'
#' @return list containing: `centers`: matrix of free centroids, `size`:
#'     number of observations assigned to each free centroid, and `radius`:
#'     mean distance between each free centroid and its observations
#'
#' @noRd
#' 
fixedKmeans <- function(x, fixed = NULL, k, nstart = 10, iter_max = 100) {
  m <- if (is.null(fixed)) 0 else nrow(fixed)
  best <- NULL

  for (s in seq_len(nstart)) {
    # k-means++ seeding, treating fixed centroids as already chosen
    centers <- matrix(NA_real_, nrow = k, ncol = ncol(x))
    if (m > 0) {
      d2 <- nnDist(x, fixed)^2
    } else {
      d2 <- rep(1, nrow(x))
    }
    for (j in seq_len(k)) {
      prob <- if (sum(d2) > 0) d2 else rep(1, nrow(x))
      centers[j, ] <- x[sample.int(nrow(x), 1, prob = prob), ]
      new_d2 <- distMat(x, centers[j, , drop = FALSE])[,1]^2
      d2 <- if (j == 1 && m == 0) new_d2 else pmin(d2, new_d2)
    }

    # Lloyd's algorithm, only updating free centroids
    for (it in seq_len(iter_max)) {
      d <- distMat(x, rbind(fixed, centers))
      assign <- max.col(-d, ties.method = "first") - m
      new_centers <- centers
      for (j in seq_len(k)) {
        members <- which(assign == j)
        if (length(members) > 0) {
          new_centers[j, ] <- colMeans(x[members, , drop = FALSE])
        } else {
          # Reseed empty cluster at the observation furthest from any centroid
          new_centers[j, ] <- x[which.max(apply(d, 1, min)), ]
        }
      }
      converged <- max(abs(new_centers - centers)) < 1e-10
      centers <- new_centers
      if (converged) {
        break
      }
    }

    d <- distMat(x, rbind(fixed, centers))
    min_d <- apply(d, 1, min)
    wss <- sum(min_d^2)
    if (is.null(best) || wss < best$wss) {
      assign <- max.col(-d, ties.method = "first") - m
      best <- list(wss = wss, centers = centers, assign = assign, min_d = min_d)
    }
  }

  size <- tabulate(pmax(best$assign, 0), nbins = k)
  radius <- vapply(seq_len(k), function(j) {
    members <- which(best$assign == j)
    if (length(members) > 0) mean(best$min_d[members]) else mean(best$min_d)
  }, numeric(1))

  list(centers = best$centers, size = size, radius = radius)
}
