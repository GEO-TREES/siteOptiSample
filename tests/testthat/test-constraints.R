sortRows <- function(m) {
  m <- unname(as.matrix(m))
  m[do.call(order, as.data.frame(m)), , drop = FALSE]
}

polys <- function(out) do.call(c, unname(out))

test_that("footprintMetrics returns the mean of every complete footprint", {
  set.seed(1)
  vals <- matrix(rnorm(200), ncol = 2)
  vals[c(15, 60), ] <- NA
  r <- makeRast(vals, 10, 10)
  for (dims in list(c(1, 1), c(2, 2), c(3, 2))) {
    expect_equal(sortRows(footprintMetrics(r, dims)),
      sortRows(popMeans(r, dims[1], dims[2])))
  }
})

test_that("min_dist keeps plots apart from each other and existing plots", {
  set.seed(2)
  r <- makeRast(matrix(rnorm(800), ncol = 2), 20, 20)
  old_ind <- c(190, 191, 210, 211)
  old_poly <- siteOptiSample:::footprintPolygon(r, old_ind)
  for (fn in list(meanminSelect, minimaxSelect, kmeansSelect, hypercubeSelect)) {
    for (md in c(1, 2.5)) {
      out <- suppressMessages(fn(r, old_ind = old_ind, n_plots = 6,
        p_new_dim = c(2, 2), het_q = NULL, min_dist = md))
      expect_length(out, 6)
      geoms <- polys(out)
      d_new <- sf::st_distance(geoms)
      expect_true(all(d_new[upper.tri(d_new)] >= md - 1e-9))
      expect_true(all(sf::st_distance(geoms, old_poly) >= md - 1e-9))
    }
  }
})

test_that("min_dist keeps plots apart from existing plots on real data", {
  r <- terra::unwrap(san_lorenzo_rast)
  p <- san_lorenzo_plots
  out <- suppressMessages(plotSelect(r, p, n_plots = 10, p_new_dim = c(100, 100),
    n_pca = 3, min_dist = 100))
  d_new <- sf::st_distance(out)
  expect_true(all(as.numeric(d_new[upper.tri(d_new)]) >= 100 - 1e-6))
  d_old <- sf::st_distance(out, sf::st_transform(p, sf::st_crs(out)))
  expect_true(all(as.numeric(d_old) >= 100 - 1e-6))
})

test_that("meanmin with cost picks the cheapest candidate within tolerance", {
  set.seed(3)
  r <- makeRast(matrix(rnorm(128), ncol = 2), 8, 8)
  r_cost <- makeRast(matrix(rep(1:8, times = 8), ncol = 1), 8, 8)
  p_pca <- matrix(c(1, -1), ncol = 2)
  pop <- popMeans(r, 2, 2)
  blocks <- allBlocks(r, 2, 2)
  bm <- t(sapply(blocks, function(b) blockMean(r, b)))
  costs <- sapply(blocks, function(b) mean(terra::values(r_cost)[b, 1]))
  for (tol in c(0, 0.3)) {
    out <- suppressMessages(meanminSelect(r, p_pca = p_pca, n_plots = 1,
      p_new_dim = c(2, 2), r_cost = r_cost, cost_tol = tol))
    current <- meanNearest(pop, p_pca)
    after <- sapply(seq_along(blocks), function(i) meanNearest(pop, rbind(p_pca, bm[i, ])))
    gain <- current - after
    ok <- which(gain >= (1 - tol) * max(gain))
    ref <- ok[order(costs[ok], after[ok])[1]]
    expect_equal(selCells(out)[[1]], blocks[[ref]])
  }
})

test_that("cost tolerance trades representativeness for lower cost", {
  set.seed(4)
  r <- makeRast(matrix(rnorm(400), ncol = 2), 10, 20)
  r_cost <- makeRast(matrix(rep(1:20, times = 10), ncol = 1), 10, 20)
  meanCost <- function(out) {
    mean(sapply(selCells(out), function(b) mean(terra::values(r_cost)[b, 1])))
  }
  for (fn in list(meanminSelect, minimaxSelect, kmeansSelect, hypercubeSelect)) {
    set.seed(5)
    lo <- suppressMessages(fn(r, n_plots = 5, p_new_dim = c(2, 2), het_q = NULL,
      r_cost = r_cost, cost_tol = 0))
    set.seed(5)
    hi <- suppressMessages(fn(r, n_plots = 5, p_new_dim = c(2, 2), het_q = NULL,
      r_cost = r_cost, cost_tol = 0.5))
    expect_lte(meanCost(hi), meanCost(lo))
  }
})

test_that("candidates with missing cost are excluded", {
  set.seed(6)
  r <- makeRast(matrix(rnorm(200), ncol = 2), 10, 10)
  cost <- rep(1, 100)
  cost[1:50] <- NA
  r_cost <- makeRast(matrix(cost, ncol = 1), 10, 10)
  out <- suppressMessages(meanminSelect(r, n_plots = 4, p_new_dim = c(2, 2),
    r_cost = r_cost))
  expect_false(any(unlist(selCells(out)) %in% 1:50))
})

test_that("refinement never increases the cost of a plot", {
  set.seed(7)
  r <- makeRast(matrix(rnorm(400), ncol = 2), 10, 20)
  r_cost <- makeRast(matrix(rep(1:20, times = 10), ncol = 1), 10, 20)
  cellCost <- function(out) {
    sapply(selCells(out), function(b) mean(terra::values(r_cost)[b, 1]))
  }
  greedy <- suppressMessages(meanminSelect(r, n_plots = 5, p_new_dim = c(2, 2),
    r_cost = r_cost, cost_tol = 0.2))
  refined <- suppressMessages(meanminSelect(r, n_plots = 5, p_new_dim = c(2, 2),
    r_cost = r_cost, cost_tol = 0.2, refine = TRUE))
  expect_true(all(cellCost(refined) <= cellCost(greedy)))
})

test_that("kmeans places new plots in clusters not covered by existing plots", {
  set.seed(8)
  col <- rep(seq_len(12), times = 12)
  group <- ceiling(col / 4)
  vals <- cbind(group * 5, -group * 5) + matrix(rnorm(288, sd = 0.1), ncol = 2)
  r <- makeRast(vals, 12, 12)
  p_pca <- rbind(c(5, -5), c(10, -10))
  out <- kmeansSelect(r, p_pca = p_pca, n_plots = 1, p_new_dim = c(2, 2), het_q = NULL)
  expect_equal(unique(group[selCells(out)[[1]]]), 3)
})

test_that("hypercube targets strata not occupied by existing plots", {
  set.seed(9)
  r <- makeRast(matrix(rep(1:20, times = 10), ncol = 1), 10, 20)
  # Existing plots occupy the lowest and highest of 6 strata
  p_pca <- matrix(c(1, 20), ncol = 1)
  out <- hypercubeSelect(r, p_pca = p_pca, n_plots = 4, p_new_dim = c(1, 1),
    het_q = NULL)
  sel_vals <- sort(sapply(selCells(out), function(x) terra::values(r)[x, 1]))
  targets <- stats::quantile(rep(1:20, times = 10), (2:5 - 0.5) / 6)
  expect_true(all(abs(sel_vals - targets) <= 0.5))
})

test_that("plotSelect diagnostics describe representativeness at plot scale", {
  r <- terra::unwrap(san_lorenzo_rast)
  p <- san_lorenzo_plots
  r_cost <- terra::init(r[[1]], "x")
  out <- suppressMessages(plotSelect(r, p, n_plots = 5, p_new_dim = c(100, 100),
    n_pca = 3, r_cost = r_cost))
  expect_equal(out$plot_order, 1:5)
  expect_true(all(c("PC1", "PC2", "PC3", "cost", "mean_dist", "max_dist") %in% names(out)))

  # Recompute distances independently
  old_ext <- extractPlotMetrics(r, p)
  old_pca <- PCALandscape(r, old_ext, center = TRUE, scale. = TRUE)
  pop <- stats::predict(old_pca$r_pca, footprintMetrics(r, c(100, 100)))[, 1:3]
  P <- old_pca$p_pca[, 1:3]
  new_vals <- sf::st_drop_geometry(out)[, c("PC1", "PC2", "PC3")]
  expect_equal(attr(out, "baseline")[["mean_dist"]],
    mean(siteOptiSample:::nnDist(pop, P)))
  for (i in 1:5) {
    d <- siteOptiSample:::nnDist(pop, rbind(P, as.matrix(new_vals[1:i, ])))
    expect_equal(out$mean_dist[i], mean(d))
    expect_equal(out$max_dist[i], max(d))
  }
  expect_true(all(diff(c(attr(out, "baseline")[["mean_dist"]], out$mean_dist)) <= 0))

  # Plot values and costs match the plot footprints
  cells <- terra::cells(r, terra::vect(out[1, ]))[, "cell"]
  expect_equal(out$cost[1], mean(terra::values(r_cost)[cells, 1]))
})

test_that("cost vectors for data frame input match cost rasters", {
  r <- terra::unwrap(san_lorenzo_rast)
  p <- san_lorenzo_plots
  r_cost <- terra::init(r[[1]], "y")
  r_df <- terra::values(r)
  valid_rows <- stats::complete.cases(r_df)
  r_df_all <- cbind(r_df[valid_rows, ], terra::crds(r))
  cost_vec <- terra::values(r_cost)[valid_rows, 1]

  out_r <- suppressMessages(plotSelect(r, p, n_plots = 5, p_new_dim = c(100, 100),
    n_pca = 3, r_cost = r_cost, cost_tol = 0.3))
  out_df <- suppressMessages(plotSelect(r_df_all, p, n_plots = 5,
    p_new_dim = c(100, 100), n_pca = 3, coord = c("x", "y"), r_cost = cost_vec,
    cost_tol = 0.3))
  cells_r <- terra::cells(r, terra::vect(out_r))[, "cell"]
  expect_setequal(which(valid_rows)[out_df], cells_r)
  expect_equal(attr(out_df, "selection")$mean_dist, out_r$mean_dist)
})
