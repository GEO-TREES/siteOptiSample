#' Recommend locations for plots
#' 
#' @param r `SpatRaster` or dataframe with structural metrics
#' @param p optional, existing plots, in the format given by `p_type`
#' @param p_type type of `p`, one of:
#'     * "polygon" (default): an `sf` or `SpatVector` object containing 
#'       polygons of existing plots
#'     * "point": an `sf` or `SpatVector` object containing points of 
#'       existing plots
#'     * "xy": a dataframe or matrix with columns `x` and `y` giving the 
#'       coordinates of existing plots in the same coordinate system as `r`. 
#'       Requires `r` to be spatial
#'     * "values": a dataframe or matrix with one row per existing plot and 
#'       one column per structural variable, with column names matching the 
#'       variables in `r`. Other columns are ignored. See Details
#'     * "rows": a vector of row indices in `r` giving locations of existing 
#'       plots. Requires `r` to be a dataframe
#' @param n_plots maximum number of new plots to add
#' @param p_new_dim optional, dimensions of new plots in the same coordinate 
#'     system as `r`. Either a single value for square plots, or a vector of
#'     two values for rectangular plots. Should be perfectly divisible by the
#'     resolution of `r`
#' @param r_mask optional, a raster which defines a mask of potential plot
#'     locations, e.g. based on accessibility. Ignored if `r` is not spatial
#' @param pca logical, if TRUE (default) the variables in `r` are run through a
#'     PCA to reduce issues of collinearity among variables. If TRUE, provide
#'     `n_pca`
#' @param n_pca optional, the number of PCA axes used to calculate the distance
#'     between pixels and plots. If not provided and `pca` is TRUE, all PCA
#'     axes are used
#' @param coord optional, if `r` is a dataframe, an optional character vector
#'     containing two column names in `r` specifying the X and Y coordinates of
#'     the grid cell centres. Coordinates must form a regular grid. If NULL and
#'     `r` is a dataframe, workflow is non-spatial.
#' @param het_q optional, numeric value between 0 and 1. Defines a threshold 
#'     to filter out internally heterogeneous candidate locations, preventing
#'     plots being placed on sharp structural transitions. Represents maximum
#'     quantile for total focal standard deviation of metrics within candidate
#'     footprint (e.g., 0.8 excludes top 20% most heterogeneous areas) 
#' @param min_dist optional, minimum distance between the edges of a new plot
#'     and any other plot, new or existing, in the units of `r`. Ignored if
#'     `r` is not spatial
#' @param r_cost optional, the cost of accessing each location, e.g. travel
#'     time. If `r` is a `SpatRaster`, a `SpatRaster` with the same geometry
#'     as `r`. If `r` is a dataframe, a numeric vector with one value per row
#'     of `r`. Locations with missing costs are excluded. See `cost_tol`
#' @param cost_tol numeric value between 0 and 1, used if `r_cost` is
#'     supplied. Each plot is placed at the cheapest candidate location whose
#'     score is within `cost_tol` of the best candidate's score. 0 ignores
#'     cost except to break ties, while larger values give more weight to
#'     cost. See the documentation of each `method` for how the tolerance
#'     is defined
#' @param distance distance metric used to compare locations in structural
#'     space, either "euclidean" (default) or "mahalanobis". See Details
#' @param method a function to optimally place plots, e.g. `meanminSelect()` or
#'     `minimaxSelect()`
#' @param ... additional arguments passed to `method`, e.g. `refine = TRUE`
#'     for `meanminSelect()`
#' 
#' @details
#' Diagnostics are returned for each new plot, in the order in which plots
#'     were selected. Representativeness is measured at the scale of a plot:
#'     the landscape is represented by the mean structural values within every
#'     possible plot footprint (see `footprintMetrics()`), and each plot by its
#'     mean structural values. `mean_dist` and `max_dist` give the mean and
#'     maximum distance in structural space from every possible plot location
#'     to its nearest plot, including existing plots and all new plots up to
#'     and including this one. Plotting `mean_dist` against `plot_order` shows
#'     how representativeness improves as plots are added, which can help to
#'     decide how many plots are needed. The values before adding any new
#'     plots are stored in the `baseline` attribute.
#'
#' Mahalanobis distance accounts for differences in variance and correlations
#'     among structural variables (or PCA axes), using the covariance of all
#'     pixels. It is implemented by transforming pixels and existing plots so
#'     that euclidean distance in the transformed space equals Mahalanobis
#'     distance in the original space, then running `method` as usual. When
#'     `pca = TRUE`, PCA axes are uncorrelated, so this is equivalent to
#'     dividing each axis by its standard deviation, giving each retained axis
#'     equal weight. Minor axes, which may mostly contain noise, then have as
#'     much influence as the first axis, so it is advisable to set `n_pca` to
#'     retain only the most important axes. When `distance = "mahalanobis"`,
#'     `mean_dist` and `max_dist` are Mahalanobis distances, and the
#'     heterogeneity filter (`het_q`) is calculated in the transformed space.
#'     Values of structural variables or PCA axes are reported in their
#'     original units.
#'
#' Existing plots can be supplied as structural values rather than locations
#'     (`p_type = "values"`), e.g. for plots outside the extent of `r`. Values
#'     should be comparable to the mean of `r` within a plot footprint, i.e.
#'     derived from the same data source, in the same units, and averaged over
#'     a plot-sized area. As the locations of these plots are unknown, new
#'     plots are not prevented from overlapping them, and `min_dist` only
#'     applies among new plots.
#'
#' @return if `r` is a `SpatRaster`, an sf dataframe with polygons of proposed
#'     new plots and columns: `plot_order`: order of selection, one column per
#'     structural variable or PCA axis giving the mean value within the plot,
#'     `cost`: mean cost within the plot, if `r_cost` is supplied, and
#'     `mean_dist` and `max_dist`: see Details. If `r` is a dataframe, a vector
#'     of row indices in `r` corresponding to the selected rows, with the same
#'     diagnostics in the `selection` attribute, which also contains a column
#'     `rows` giving the row indices of each plot.
#'
#' @import terra
#' @import sf
#' 
#' @export
#'
plotSelect <- function(r, p = NULL, n_plots, 
  p_type = c("polygon", "point", "xy", "values", "rows"),
  p_new_dim = NULL, r_mask = NULL, pca = TRUE, n_pca = NULL, coord = NULL, het_q = NULL, min_dist = NULL, 
  r_cost = NULL, cost_tol = 0.1, distance = c("euclidean", "mahalanobis"), 
  method = meanminSelect, ...) {

  distance <- match.arg(distance)
  p_type <- match.arg(p_type)
  
  # Input validation 
  is_rast <- inherits(r, "SpatRaster")
  is_spatial <- is_rast || !is.null(coord)
  
  if (!is.null(het_q) && 
      (!is.numeric(het_q) || length(het_q) != 1 || het_q <= 0 || het_q > 1)) {
    stop("`het_q` must be a single numeric value > 0 and <= 1")
  }
  
  if (!is.null(min_dist) && !is_spatial) {
    message("`min_dist` is ignored for non-spatial inputs.")
    min_dist <- NULL
  }

  if (is_rast && !is.null(r_cost) && !inherits(r_cost, "SpatRaster")) {
    stop("`r_cost` must be a SpatRaster when `r` is a SpatRaster")
  }

  # Data coercion 
  if (!is_rast) {
    if (!is.null(r_cost) && (!is.numeric(r_cost) || length(r_cost) != nrow(r))) {
      stop("`r_cost` must be a numeric vector with one value per row of `r`")
    }
    if (!inherits(r, c("data.frame", "matrix"))) {
      stop("`r` must be a 'SpatRaster', 'data.frame', or 'matrix'")
    }
    if (inherits(r, "sf")) {
      r <- sf::st_drop_geometry(r)
    }
    r <- as.data.frame(r)
    
    if (!is.null(r_mask)) {
      message("`r_mask` is ignored for non-SpatRaster inputs.")
      r_mask <- NULL
    }
    
    if (is_spatial) {
      if (!is.character(coord) || length(coord) != 2 || !all(coord %in% colnames(r))) {
        stop("`coord` must contain two valid column names found in `r`.")
      }
      r_coords <- r[, coord, drop = FALSE]
      r <- tryCatch(
        {
          terra::rast(cbind(r_coords, r[, setdiff(colnames(r), coord), drop = FALSE]), type = "xyz")
        },
        error = function(e) {
          stop("Spatial coordinates do not form a regular grid.\nOriginal error: ", e$message)
        }
      )
    } else {
      r_rast <- terra::rast(nrows = nrow(r), ncols = 1, nlyrs = ncol(r), 
        crs = "", ext = c(0, 1, 0, nrow(r)))
      terra::values(r_rast) <- r
      names(r_rast) <- names(r)
      r <- r_rast
    }
    
    # Coordinates are assumed to share the coordinate system of spatial `p`
    if (is_spatial && inherits(p, c("sf", "sfc", "SpatVector")) &&
        terra::crs(r) == "" && terra::crs(p) != "") {
      terra::crs(r) <- terra::crs(p)
    }

    # Convert costs to a raster aligned with `r`
    if (!is.null(r_cost)) {
      cost_vals <- rep(NA_real_, terra::ncell(r))
      if (is_spatial) {
        cost_vals[terra::cellFromXY(r, as.matrix(r_coords))] <- r_cost
      } else {
        cost_vals <- r_cost
      }
      r_cost <- terra::setValues(r[[1]], cost_vals)
    }
  }

  # Existing plots
  p_vals <- NULL
  if (!is.null(p)) {
    crs_r <- if (terra::crs(r) == "") NA else terra::crs(r)
    if (p_type %in% c("polygon", "point")) {
      sf_types <- if (p_type == "polygon") {
        c("POLYGON", "MULTIPOLYGON") 
      } else {
        c("POINT", "MULTIPOINT")
      }
      vect_type <- if (p_type == "polygon") "polygons" else "points"
      if (!is_spatial) {
        stop("`p_type = \"", p_type, "\"` requires `r` to be spatial")
      }
      if (!(isSFType(p, sf_types) || 
          (inherits(p, "SpatVector") && terra::geomtype(p) == vect_type))) {
        stop("`p` must be an sf or SpatVector object containing ", vect_type, 
          " when `p_type = \"", p_type, "\"`")
      }
    } else if (p_type == "xy") {
      if (!is_spatial) {
        stop("`p_type = \"xy\"` requires `r` to be spatial")
      }
      if (!inherits(p, c("data.frame", "matrix")) || 
          !all(c("x", "y") %in% colnames(p))) {
        stop("`p` must be a dataframe or matrix with columns `x` and `y` ",
          "when `p_type = \"xy\"`")
      }
      if (inherits(p, "sf")) {
        p <- sf::st_drop_geometry(p)
      }
      p <- sf::st_as_sf(as.data.frame(p)[, c("x", "y")], coords = c("x", "y"), 
        crs = crs_r)
    } else if (p_type == "values") {
      # Structural values are used directly, with no locations
      if (!inherits(p, c("data.frame", "matrix"))) {
        stop("`p` must be a dataframe or matrix when `p_type = \"values\"`")
      }
      if (inherits(p, "sf")) {
        p <- sf::st_drop_geometry(p)
      }
      p_vals <- as.data.frame(p)
      p <- NULL
    } else if (p_type == "rows") {
      if (is_rast) {
        stop("`p_type = \"rows\"` requires `r` to be a dataframe")
      }
      max_idx <- if (is_spatial) nrow(r_coords) else terra::ncell(r)
      if (!is.numeric(p) || 
          any(is.na(p) | p < 1 | p > max_idx | p != floor(p))) {
        stop("`p` must be valid row indices in `r` when `p_type = \"rows\"`")
      }
      if (is_spatial) {
        p <- sf::st_as_sf(r_coords[p, , drop = FALSE], coords = coord, 
          crs = crs_r) 
      }
    }
  }
  
  # Dimension and mask constraints
  res_r <- terra::res(r)
  if (is.null(p_new_dim)) {
    p_new_dim <- res_r 
  } else {
    p_new_dim <- rep_len(p_new_dim, 2)
  }
  
  if (any(p_new_dim < res_r)) {
    warning("'p_new_dim' is smaller than raster resolution. Setting to raster resolution.")
    p_new_dim <- res_r
  } else if (any(p_new_dim %% res_r != 0)) {
    p_new_dim <- floor(p_new_dim / res_r) * res_r
    warning("Rounding 'p_new_dim' down to evenly divisible values: ", 
      paste(p_new_dim, collapse = ", "))
  }
  
  r_mask <- if (!is.null(r_mask)) terra::mask(r, r_mask) else r
  
  # Feature extraction and PCA
  if (!is.null(p_vals)) {
    missing_vars <- setdiff(names(r), colnames(p_vals))
    if (length(missing_vars) > 0) {
      stop("`p` is missing structural variables found in `r`: ",
        paste(missing_vars, collapse = ", "))
    }
    old_ext <- p_vals[, names(r), drop = FALSE]
  } else if (!is.null(p)) {
    if (is_spatial) {
      old_ext <- extractPlotMetrics(r, p) 
    } else { 
      old_ext <- as.data.frame(r)[p, , drop = FALSE]
    }
  } else {
    old_ext <- NULL
  }
  
  if (pca && terra::nlyr(r) > 1) {
    n_pca <- if (is.null(n_pca)) terra::nlyr(r) else n_pca
    if (n_pca > terra::nlyr(r)) {
      stop("`n_pca` must not be greater than the number of layers in `r`")
    }
    old_pca <- PCALandscape(r, old_ext, center = TRUE, scale. = TRUE)
    
    r_pca <- rep(r[[1]], n_pca)
    v_pca <- matrix(NA, nrow = terra::ncell(r_pca), ncol = n_pca)
    v_pca[stats::complete.cases(terra::values(r)), ] <- old_pca$r_pca$x[, 1:n_pca, drop = FALSE]
    r_pca <- terra::setValues(r_pca, v_pca)
    names(r_pca) <- colnames(old_pca$r_pca$x[, 1:n_pca, drop = FALSE])
    
    if (!is.null(old_ext)) {
      p_pca <- old_pca$p_pca[, 1:n_pca, drop = FALSE] 
    } else {
      p_pca <- NULL
    }
  } else {
    if (pca) message("Only one variable in `r`. PCA will be skipped.")

    # Scale pixels and plots using the same pixel mean and standard deviation
    v_r <- terra::values(r)
    r_center <- colMeans(v_r, na.rm = TRUE)
    r_scale <- apply(v_r, 2, stats::sd, na.rm = TRUE)
    r_pca <- terra::setValues(r, scale(v_r, center = r_center, scale = r_scale))
    names(r_pca) <- names(r)

    if (!is.null(old_ext)) {
      p_pca <- scale(as.matrix(old_ext)[, names(r), drop = FALSE], 
        center = r_center, scale = r_scale)
    } else {
      p_pca <- NULL
    }
  }
  
  # Values reported in diagnostics, in original units
  r_report <- r_pca

  # Transform so that euclidean distance equals Mahalanobis distance
  if (distance == "mahalanobis") {
    w_mat <- whiteningMatrix(terra::values(r_pca))
    r_pca <- terra::setValues(r_pca, terra::values(r_pca) %*% w_mat)
    names(r_pca) <- names(r_report)
    if (!is.null(p_pca)) {
      p_pca <- as.matrix(p_pca) %*% w_mat
    }
  }
  
  # Define search space and execute selection algorithm
  if (!is.null(p)) {
    if (inherits(p, c("sf", "sfc", "SpatVector")) || is_spatial) {
      p_vect <- if (inherits(p, "SpatVector")) p else terra::vect(p)
      old_ind <- which(stats::complete.cases(terra::values(terra::mask(r, p_vect))))
    } else {
      old_ind <- p 
    } 
  } else {
    old_ind <- NULL
  }
  
  new_ind <- which(stats::complete.cases(terra::values(r_mask)))
  
  cand_args <- list(r_pca = r_pca, p_pca = p_pca, 
    old_ind = old_ind, new_ind = new_ind, 
    n_plots = n_plots, p_new_dim = p_new_dim, het_q = het_q, 
    min_dist = min_dist, r_cost = r_cost, cost_tol = cost_tol)
  cand_args <- cand_args[intersect(names(cand_args), names(formals(method)))]
  p_list <- do.call(method, c(cand_args, list(...)))

  # Diagnostics
  sel_cells <- if (length(p_list) > 0) {
    lapply(strsplit(names(p_list), ":"), as.integer)
  } else {
    list()
  }
  diag <- selectionDiagnostics(r_pca, p_pca, sel_cells, p_new_dim, r_cost, 
    r_report)
  
  # Format output
  if (is_rast) {
    if (terra::crs(r) == "") {
      crs_val <- NA 
    } else {
      crs_val <- terra::crs(r)
    }
    if (length(p_list) == 0) {
      out <- sf::st_sf(diag$selection, geometry = sf::st_sfc(crs = crs_val))
    } else {
      out <- sf::st_sf(diag$selection, geometry = do.call(c, p_list), crs = crs_val)
    }
    attr(out, "baseline") <- diag$baseline
    return(out)
  } 
  
  if (is_spatial) {
    r_coords_sf <- sf::st_as_sf(r_coords, coords = coord, crs = terra::crs(r))
    rows <- lapply(p_list, function(poly) {
      if (terra::crs(r) != "") {
        sf::st_crs(poly) <- terra::crs(r)
      }
      sf::st_intersects(poly, r_coords_sf)[[1]]
    })
  } else {
    # Cell IDs equal row indices for non-spatial input
    rows <- sel_cells
  }

  out <- as.integer(unlist(rows))
  selection <- diag$selection
  selection$rows <- vapply(rows, paste, character(1), collapse = ":")
  attr(out, "selection") <- selection
  attr(out, "baseline") <- diag$baseline
  return(out)
}

#' Calculate representativeness diagnostics for selected plots
#'
#' @param r_pca `SpatRaster` of PCA scores or scaled structural metrics
#' @param p_pca optional, structural values of existing plots
#' @param sel_cells list of cell IDs of each new plot, in order of selection
#' @param p_new_dim dimensions of new plots
#' @param r_cost optional, `SpatRaster` of access costs
#' @param r_report optional, `SpatRaster` of values to report for each plot,
#'     if different from `r_pca`
#'
#' @return list containing: `selection`: dataframe with one row per new plot,
#'     and `baseline`: mean and maximum distance before adding new plots
#'
#' @noRd
#' 
selectionDiagnostics <- function(r_pca, p_pca, sel_cells, p_new_dim, r_cost = NULL,
  r_report = r_pca) {
  v_pca <- terra::values(r_pca)
  v_report <- terra::values(r_report)
  pop <- footprintMetrics(r_pca, p_new_dim)

  plot_vals <- matrix(NA_real_, nrow = length(sel_cells), ncol = ncol(v_pca))
  report_vals <- matrix(NA_real_, nrow = length(sel_cells), ncol = ncol(v_report),
    dimnames = list(NULL, names(r_report)))
  for (i in seq_along(sel_cells)) {
    plot_vals[i, ] <- colMeans(v_pca[sel_cells[[i]], , drop = FALSE])
    report_vals[i, ] <- colMeans(v_report[sel_cells[[i]], , drop = FALSE])
  }

  if (!is.null(p_pca)) {
    p_pca <- as.matrix(p_pca)
    p_pca <- p_pca[stats::complete.cases(p_pca), , drop = FALSE]
  }

  if (!is.null(p_pca) && nrow(p_pca) > 0) {
    min_d <- nnDist(pop, p_pca)
    baseline <- c(mean_dist = mean(min_d), max_dist = max(min_d))
  } else {
    min_d <- rep(Inf, nrow(pop))
    baseline <- c(mean_dist = NA_real_, max_dist = NA_real_)
  }

  mean_dist <- max_dist <- numeric(length(sel_cells))
  for (i in seq_along(sel_cells)) {
    min_d <- pmin(min_d, distMat(pop, plot_vals[i, , drop = FALSE])[,1])
    mean_dist[i] <- mean(min_d)
    max_dist[i] <- max(min_d)
  }

  selection <- data.frame(plot_order = seq_along(sel_cells), report_vals)
  if (!is.null(r_cost)) {
    v_cost <- terra::values(r_cost)[,1]
    selection$cost <- vapply(sel_cells, function(ids) mean(v_cost[ids]), numeric(1))
  }
  selection$mean_dist <- mean_dist
  selection$max_dist <- max_dist

  list(selection = selection, baseline = baseline)
}

#' Calculate a matrix which transforms Mahalanobis distance to euclidean
#'     distance
#'
#' @param x numeric matrix of observations, which may contain missing values
#'
#' @return square matrix `W` such that euclidean distances between rows of
#'     `x %*% W` equal Mahalanobis distances between rows of `x`, using the
#'     covariance of the complete rows of `x`
#'
#' @noRd
#' 
whiteningMatrix <- function(x) {
  x <- x[stats::complete.cases(x), , drop = FALSE]
  S_inv <- tryCatch(solve(stats::cov(x)), error = function(e) {
    stop("The covariance matrix of structural variables is singular, so ",
      "Mahalanobis distance cannot be calculated. Try `pca = TRUE` with fewer ",
      "axes (`n_pca`).", call. = FALSE)
  })
  # S_inv = t(L) %*% L, so euclidean distance after x %*% t(L) is Mahalanobis
  t(chol(S_inv))
}
