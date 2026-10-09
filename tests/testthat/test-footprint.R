dims_list <- list(c(1, 1), c(2, 2), c(3, 3), c(2, 3), c(3, 2), c(4, 1))

test_that("footprints are rectangles of the requested dimensions", {
  r <- makeRast(matrix(runif(100), ncol = 1), 10, 10)
  for (dims in dims_list) {
    fp <- siteOptiSample:::plotFootprint(r, dims)
    expect_equal(fp$n_cells, prod(dims))
    cells <- siteOptiSample:::footprintCells(r, fp, 55)
    rc <- terra::rowColFromCell(r, cells)
    # p_new_dim is c(x, y): columns then rows
    expect_equal(diff(range(rc[,2])) + 1, dims[1])
    expect_equal(diff(range(rc[,1])) + 1, dims[2])
    expect_length(unique(cells), prod(dims))
  }
})

test_that("footprint means equal the mean of footprint cells", {
  set.seed(1)
  r <- makeRast(matrix(rnorm(200), ncol = 2), 10, 10)
  for (dims in dims_list) {
    fp <- siteOptiSample:::plotFootprint(r, dims)
    v_fp <- siteOptiSample:::footprintMeans(r, fp)
    for (center in which(stats::complete.cases(v_fp))) {
      cells <- siteOptiSample:::footprintCells(r, fp, center)
      expect_equal(unname(v_fp[center, ]), unname(blockMean(r, cells)))
    }
  }
})

test_that("available centres are exactly the blocks within allowed, unoccupied cells", {
  set.seed(2)
  vals <- matrix(rnorm(100), ncol = 1)
  vals[c(5, 23, 77)] <- NA
  r <- makeRast(vals, 10, 10)
  allowed <- setdiff(which(!is.na(vals)), c(40:45))
  old_ind <- c(88, 89)
  for (dims in dims_list) {
    fp <- siteOptiSample:::plotFootprint(r, dims)
    centers <- siteOptiSample:::availCenters(r, fp, allowed, old_ind,
      seq_len(terra::ncell(r)))
    pkg_blocks <- lapply(centers, function(x) {
      sort(siteOptiSample:::footprintCells(r, fp, x))
    })
    ref_blocks <- allBlocks(r, dims[1], dims[2], setdiff(allowed, old_ind))
    expect_setequal(sapply(pkg_blocks, paste, collapse = ":"),
      sapply(ref_blocks, paste, collapse = ":"))
  }
})

test_that("footprints never extend beyond the raster edge", {
  r <- makeRast(matrix(runif(9), ncol = 1), 3, 3)
  expect_message(out <- meanminSelect(r, n_plots = 5, p_new_dim = c(2, 2)), 
    "No possible locations")
  expect_length(out, 1)
  expect_length(selCells(out)[[1]], 4)
})

test_that("heterogeneity filter keeps only centres at or below the quantile cutoff", {
  set.seed(3)
  vals <- matrix(rnorm(400), ncol = 1)
  r <- makeRast(vals, 20, 20)
  fp <- siteOptiSample:::plotFootprint(r, c(2, 2))
  safe <- siteOptiSample:::hetSafeCenters(r, fp, seq_len(400), 0.5)
  sds <- sapply(seq_len(400), function(x) {
    cells <- siteOptiSample:::footprintCells(r, fp, x)
    if (length(cells) < 4) NA else stats::sd(vals[cells])
  })
  cutoff <- stats::quantile(sds, 0.5, na.rm = TRUE)
  expect_true(all(sds[safe] <= cutoff, na.rm = TRUE))
  expect_true(all(which(sds <= cutoff) %in% safe))
})
