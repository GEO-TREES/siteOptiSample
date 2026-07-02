#' Recommend locations for plots
#' 
#' @param r `SpatRaster` or dataframe with structural metrics
#' @param p optional, either an `sf` object containing polygons or points of 
#'     existing plots, a dataframe with two columns `x` and `y` describing 
#'     coordinates in the same coordinate system as `r`, or if `r` is not 
#'     spatial a vector of row indices in `r` specifying locations with 
#'     existing plots
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
#' @param method a function to optimally place plots, e.g. `meanminSelect()` or
#'     `minimaxSelect()`
#' 
#' @return sf dataframe with polygons of proposed new plots, or if `r` is a
#'     dataframe, a vector of row indices in `r` corresponding to the selected
#'     rows
#'
#' @import terra
#' @import sf
#' 
#' @export
#'
plotSelect <- function(r, p = NULL, n_plots, p_new_dim = NULL, r_mask = NULL, 
  pca = TRUE, n_pca = NULL, coord = NULL, het_q = NULL, method = meanminSelect) {
  
  # Input validation 
  is_rast <- inherits(r, "SpatRaster")
  is_spatial <- is_rast || !is.null(coord)
  
  if (!is.null(het_q) && 
      (!is.numeric(het_q) || length(het_q) != 1 || het_q <= 0 || het_q > 1)) {
    stop("`het_q` must be a single numeric value > 0 and <= 1")
  }
  
  # Data coercion 
  if (!is_rast) {
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
    
    if (!is.null(p) && !inherits(p, c("sf", "sfc", "SpatVector"))) {
      p_vec <- suppressWarnings(as.numeric(unlist(p)))
      max_idx <- if (is_spatial) nrow(r_coords) else terra::ncell(r)
      
      if (any(is.na(p_vec) | p_vec < 1 | p_vec > max_idx | p_vec != floor(p_vec))) {
        stop("`p` must be valid row indices in `r` when input is not a SpatRaster and `p` is not a spatial object.")
      }
      if (is_spatial) {
        p <- sf::st_as_sf(r_coords[p_vec, , drop = FALSE], coords = coord) 
      } else { 
        p <- p_vec
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
  if (!is.null(p)) {
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
    old_pca <- PCALandscape(r, old_ext, center = TRUE, scale. = TRUE)
    
    r_pca <- rep(r[[1]], n_pca)
    v_pca <- matrix(NA, nrow = terra::ncell(r_pca), ncol = n_pca)
    v_pca[complete.cases(terra::values(r)), ] <- old_pca$r_pca$x[, 1:n_pca, drop = FALSE]
    terra::setValues(r_pca, v_pca)
    names(r_pca) <- colnames(old_pca$r_pca$x[, 1:n_pca, drop = FALSE])
    
    if (!is.null(p)) {
      p_pca <- old_pca$p_pca[, 1:n_pca, drop = FALSE] 
    } else {
      p_pca <- NULL
    }
  } else {
    if (pca) message("Only one variable in `r`. PCA will be skipped.")
    r_pca <- scale(r)
    if (!is.null(p)) {
      if (inherits(old_ext, c("sf", "sfc"))) {
        p_pca <- sf::st_drop_geometry(old_ext)[, names(r_pca), drop = FALSE] 
      } else {
        p_pca <- old_ext[, names(r_pca), drop = FALSE]
      }
    } else {
      p_pca <- NULL
    }
  }
  
  # Define search space and execute selection algorithm
  if (!is.null(p)) {
    if (inherits(p, c("sf", "sfc", "SpatVector")) || is_spatial) {
      old_ind <- which(complete.cases(terra::values(terra::mask(r, terra::vect(p)))))
    } else {
      old_ind <- p 
    } 
  } else {
    old_ind <- NULL
  }
  
  new_ind <- which(complete.cases(terra::values(r_mask)))
  
  cand_args <- list(r_pca = r_pca, p_pca = p_pca, 
    old_ind = old_ind, new_ind = new_ind, 
    n_plots = n_plots, p_new_dim = p_new_dim, het_q = het_q)
  p_list <- do.call(method, cand_args[intersect(names(cand_args), names(formals(method)))])
  
  # Format output
  if (is_rast) {
    if (terra::crs(r) == "") {
      crs_val <- NA 
    } else {
      crs_val <- terra::crs(r)
    }
    return(sf::st_sf(geometry = do.call(c, p_list), crs = crs_val))
  } 
  
  if (is_spatial) {
    r_coords_sf <- sf::st_as_sf(r_coords, coords = coord, crs = terra::crs(r))
    return(unlist(lapply(p_list, function(poly) {
      if (terra::crs(r) != "") {
        sf::st_crs(poly) <- terra::crs(r)
      }
      sf::st_intersects(poly, r_coords_sf)[[1]]
    })))
  } 
  
  return(as.numeric(names(p_list)))
}
