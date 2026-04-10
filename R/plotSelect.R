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
#'     system as `r`. Either a single value for square plots, or a vector of two
#'     values for rectangular plots. Should be perfectly divisible by the
#'     resolution of `r`
#' @param r_mask optional, a raster which defines a mask of potential plot
#'     locations, e.g. based on accessibility. Ignored if `r` is not spatial
#' @param pca logical, if TRUE (default) the variables in `r` are run through a
#'     PCA to reduce issues of collinearity among variables. If TRUE, provide
#'     `n_pca`
#' @param n_pca optional, the number of PCA axes used to calculate the distance
#'     between pixels and plots. If not provided and `pca` is TRUE, all PCA axes
#'     are used
#' @param coord optional, if `r` is a dataframe, a character vector containing 
#'     two column names in `r` specifying the X and Y coordinates of the grid 
#'     cell centres. Coordinates must form a regular grid
#' @param het_q optional, numeric value between 0 and 1. Defines a threshold 
#'     to filter out internally heterogeneous candidate locations, preventing plots 
#'     being placed on sharp structural transitions. Represents maximum 
#'     quantile for total focal standard deviation of metrics within 
#'     candidate footprint (e.g., 0.8 excludes top 20% most heterogeneous 
#'     areas) 
#' @param method a function to optimally place plots, e.g. `meanminSelect()` or
#'     `minimaxSelect()`
#' 
#' @return sf dataframe with polygons of proposed new plots, or if `r` is a dataframe, a vector of row indices in `r` corresponding to the selected rows
#'
#' @import terra
#' @import sf
#' 
#' @export
#'
plotSelect <- function(r, p = NULL, n_plots, p_new_dim = NULL, r_mask = NULL, 
  pca = TRUE, n_pca = NULL, coord = NULL, het_q = NULL, method = meanminSelect) {

  is_rast <- TRUE
  # Check type of r and check inputs
  if (!inherits(r, "SpatRaster")) {
    is_rast <- FALSE

    if (!inherits(r, c("data.frame", "matrix"))) { 
      stop("if `r` is not a 'SpatRaster' it must be a 'data.frame' or 'matrix'")
    }

    if (!is.character(coord) || length(coord) != 2) {
      stop("If `r` is not a 'SpatRaster', `coord` must be a vector of two column names specifying X and Y coordinates of each grid cell.")
     }

    if (!is.null(r_mask)) { 
      message("if `r` is not a 'SpatRaster', `r_mask` is ignored.")
      r_mask <- NULL
    }

    # Coerce to dataframe
    if (inherits(r, "sf")) {
      r <- sf::st_drop_geometry(r)
    }
    r <- as.data.frame(r)

    if (!all(coord %in% colnames(r))) {
      stop("Columns specified in `coord` not found in `r`.")
    }

    # Rearrange dataframe so X and Y are initial columns 
    r_coords <- r[, coord, drop = FALSE]
    r_metrics <- r[, setdiff(colnames(r), coord), drop = FALSE]
    r <- cbind(r_coords, r_metrics)

    r <- tryCatch({
      terra::rast(r, type = "xyz")
    }, error = function(e) {
      stop(
        "Unable to create 'SpatRaster' from `r` dataframe.\n", 
        "Spatial coordinates in `coord` do not form a perfectly regular grid.\n",
        "Original terra error: ", e$message)
    })

    # Convert p (row indices) to spatial points 
    if (!is.null(p)) { 
      p_vec <- suppressWarnings(as.numeric(unlist(p)))
      
      if (any(is.na(p_vec)) || any(p_vec > nrow(r_coords)) || 
        any(p_vec < 1) || !all(p_vec == floor(p_vec))) {
        stop("if `r` is not a 'SpatRaster', `p` must be a valid vector of row indices in `r`.")
      }
      p_coords <- r_coords[p_vec, , drop = FALSE]
      
      p <- sf::st_as_sf(p_coords, coords = coord)
    }

  }

  # Check heterogeneity quantile
  if (!is.null(het_q) && 
    (!is.numeric(het_q) || length(het_q) != 1 || (het_q > 1 | het_q <= 0))) {
    stop("`het_q` must be a single numeric value > 0 and <= 1")
  }
  
  # Assign number of PCA axes
  if (pca && is.null(n_pca)) {
    message("`n_pca` not provided, including all PCA axes")
    n_pca <- terra::nlyr(r)
  }

  # If new plot dimensions not provided, set to raster resolution
  if (is.null(p_new_dim)) {
    p_new_dim <- terra::res(r)
  }

  # If new plot dimensions are smaller than the raster resolution
  if (any(p_new_dim < terra::res(r))) {
    warning("'p_new_dim' is smaller than the resolution of 'r'.\nSetting 'p_new_dim' to the resolution of 'r'\n",
      immediate. = TRUE)
    p_new_dim <- terra::res(r)
  }

  if (any(p_new_dim %% terra::res(r) != 0)) {
    p_new_dim <- floor(p_new_dim / terra::res(r)) * terra::res(r)
    warning("'p_new_dim' is not perfectly divisible by the resolution of 'r'.\nRounding 'p_new_dim' down to a value which is evenly divisible: ", paste(p_new_dim, collapse = ", "))
  }

  # If new plot dimensions are just a single value, multiply it up
  if (length(p_new_dim) == 1) {
    p_new_dim <- rep(p_new_dim, 2)
  }

  # If mask for new plots is provided, mask the raster
  if (!is.null(r_mask)) { 
    r_mask <- terra::mask(r, r_mask)
  } else {
    r_mask <- r
  }

  # Convert p to spatial vector object 
  # Extract metrics for each plot polygon
  if (!is.null(p)) {
    old_ext <- extractPlotMetrics(r, p)
  } else {
    old_ext <- NULL
  }

  # Optional PCA to reduce dimensionality
  if (pca) { 
    old_pca <- PCALandscape(r, old_ext, center = TRUE, scale. = TRUE)
    r_pca <- rep(r[[1]], n_pca)
    pca_scores <- old_pca$r_pca$x[, 1:n_pca, drop = FALSE]
    v_pca <- matrix(NA, nrow = terra::ncell(r_pca), ncol = n_pca)
    v_pca[complete.cases(terra::values(r)), ] <- pca_scores
    terra::setValues(r_pca, v_pca)
    names(r_pca) <- colnames(pca_scores)
    
    if (!is.null(p)) {
      p_pca <- old_pca$p_pca[,1:n_pca, drop = FALSE]
    } else {
      p_pca <- NULL
    }
  } else {
    r_pca <- scale(r)
    if (!is.null(p)) {
      p_pca <- old_ext[,names(r_pca), drop = FALSE]
      if (inherits(p_pca, "sf")) {
        p_pca <- sf::st_drop_geometry(p_pca)
      }
    } else {
      p_pca <- NULL
    }
  }

  # Get index values in PCA for existing plots and potential new plots
  if (!is.null(p)) {
    old_ind <- which(complete.cases(terra::values(terra::mask(r, terra::vect(p)))))
  } else {
    old_ind <- NULL
  }

  # Initialise PCA values for candidate plots 
  new_ind <- which(complete.cases(terra::values(r_mask)))

  # Implement plot selection algorithm
  p_list <- method(
    r_pca = r_pca,
    p_pca = p_pca,
    old_ind = old_ind,
    new_ind = new_ind,
    n_plots = n_plots,
    p_new_dim = p_new_dim,
    het_q = het_q
  )

  # Prepare output sf object
  if (!is_rast) {
    # Intersect the generated polygons back with original coords to get row indices
    r_coords_sf <- sf::st_as_sf(r_coords, coords = coord, remove = FALSE)
    
    # Apply CRS if the fake raster generated one (rare for raw df, but good practice)
    if (terra::crs(r) != "") {
      sf::st_crs(r_coords_sf) <- terra::crs(r)
    }
    
    out <- unlist(lapply(p_list, function(poly) {
      if (terra::crs(r) != "") {
        sf::st_crs(poly) <- terra::crs(r)
      }
      intersecting_indices <- sf::st_intersects(poly, r_coords_sf)
      intersecting_indices[[1]]
    }))

  } else {
    crs_val <- terra::crs(r)
    if (crs_val == "") {
      crs_val <- NA
    }

    out <- sf::st_sf(geometry = do.call(c, p_list), crs = crs_val)
  }

  # Return
  return(out)
}


