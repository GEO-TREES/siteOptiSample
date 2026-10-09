# Reference implementations used to check the package's algorithms. These
# deliberately avoid the package's own footprint and distance code: every
# possible plot is enumerated as a rectangular block of cells by its top-left
# row and column.

# Create a SpatRaster with 1 x 1 cells from a matrix of values (one column per
# layer, rows in terra cell order)
makeRast <- function(vals, nrow, ncol) {
  vals <- as.matrix(vals)
  r <- terra::rast(nrows = nrow, ncols = ncol, nlyrs = ncol(vals),
    xmin = 0, xmax = ncol, ymin = 0, ymax = nrow, crs = "")
  terra::values(r) <- vals
  names(r) <- paste0("PC", seq_len(ncol(vals)))
  r
}

# Enumerate all blocks of p_y rows x p_x columns that fit within the raster and
# only contain cells in `allowed`
allBlocks <- function(r, p_x, p_y, allowed = NULL) {
  if (is.null(allowed)) {
    allowed <- which(stats::complete.cases(terra::values(r)))
  }
  out <- list()
  for (i in seq_len(terra::nrow(r) - p_y + 1)) {
    for (j in seq_len(terra::ncol(r) - p_x + 1)) {
      rows <- rep(i:(i + p_y - 1), times = p_x)
      cols <- rep(j:(j + p_x - 1), each = p_y)
      cells <- (rows - 1) * terra::ncol(r) + cols
      if (all(cells %in% allowed)) {
        out[[length(out) + 1]] <- sort(cells)
      }
    }
  }
  out
}

blockMean <- function(r, cells) {
  colMeans(terra::values(r)[cells, , drop = FALSE])
}

# Naive euclidean distance between two vectors
d2 <- function(a, b) sqrt(sum((a - b)^2))

# Means of every possible block in the landscape, representing the population
# of possible plots
popMeans <- function(r, p_x, p_y) {
  t(sapply(allBlocks(r, p_x, p_y), function(b) blockMean(r, b)))
}

# Mean distance from every possible plot in the population to its nearest
# plot, where plots are represented by their mean values
meanNearest <- function(pop, plot_means) {
  mean(apply(pop, 1, function(px) {
    min(apply(plot_means, 1, function(pm) d2(px, pm)))
  }))
}

# Cell IDs selected by a selection function, in selection order
selCells <- function(p_list) {
  lapply(strsplit(names(p_list), ":"), function(x) sort(as.integer(x)))
}
