test_that("extractPlotMetrics returns raster values at points", {
  r <- makeRast(matrix(1:20, ncol = 1), 4, 5)
  pts <- data.frame(x = c(0.5, 4.5), y = c(3.5, 0.5))
  expect_equal(unname(extractPlotMetrics(r, pts)[, 1]), c(1, 20))
  pts_sf <- sf::st_as_sf(pts, coords = c("x", "y"))
  expect_equal(unname(extractPlotMetrics(r, pts_sf)[, 1]), c(1, 20))
})

test_that("extractPlotMetrics averages over polygons", {
  r <- makeRast(matrix(1:20, ncol = 1), 4, 5)
  poly <- sf::st_as_sf(sf::st_sfc(sf::st_polygon(list(
    matrix(c(0, 2, 2, 2, 2, 4, 0, 4, 0, 2), ncol = 2, byrow = TRUE)))))
  # Top-left 2 x 2 cells: 1, 2, 6, 7
  expect_equal(unname(extractPlotMetrics(r, poly)[1, 1]), mean(c(1, 2, 6, 7)))
})

test_that("classifRepres classifies pixels by distance from plot centroid", {
  p <- matrix(c(-1, 0, 1, 0, 0, 1, 0, -1), ncol = 2, byrow = TRUE,
    dimnames = list(NULL, c("PC1", "PC2")))
  r_pca <- matrix(c(0, 0, 0.5, 0, 5, 5), ncol = 2, byrow = TRUE,
    dimnames = list(NULL, c("PC1", "PC2")))
  out <- classifRepres(r_pca, p, n_pca = 2)
  expect_equal(out$dist, c(0, 0.5, sqrt(50)))
  expect_equal(as.character(out$group), c(
    "Well-represented by existing plots",
    "Well-represented by existing plots",
    "Poorly represented"))

  p_new <- matrix(c(5, 5, 6, 6, 4, 4), ncol = 2, byrow = TRUE,
    dimnames = list(NULL, c("PC1", "PC2")))
  out_new <- classifRepres(r_pca, p, p_new = p_new, n_pca = 2)
  expect_equal(as.character(out_new$group)[3],
    "Well-represented by existing and proposed plots")
  expect_false(any(is.na(out_new$group)))
})
