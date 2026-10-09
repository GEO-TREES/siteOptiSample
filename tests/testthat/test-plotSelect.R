r <- terra::unwrap(san_lorenzo_rast)
p <- san_lorenzo_plots

test_that("plotSelect passes PCA scores to the selection algorithm", {
  old_ext <- extractPlotMetrics(r, p)
  old_pca <- PCALandscape(r, old_ext, center = TRUE, scale. = TRUE)
  captured <- NULL
  spy <- function(r_pca, p_pca, n_plots, p_new_dim) {
    captured <<- list(r_pca = r_pca, p_pca = p_pca)
    list()
  }
  plotSelect(r, p, n_plots = 1, p_new_dim = c(100, 100), n_pca = 3, method = spy)
  v <- terra::values(captured$r_pca)
  expect_equal(unname(v[stats::complete.cases(v), ]),
    unname(old_pca$r_pca$x[, 1:3]))
  expect_equal(unname(captured$p_pca), unname(old_pca$p_pca[, 1:3]))
})

test_that("plotSelect scales pixels and plots identically when pca = FALSE", {
  captured <- NULL
  spy <- function(r_pca, p_pca, n_plots, p_new_dim) {
    captured <<- list(r_pca = r_pca, p_pca = p_pca)
    list()
  }
  plotSelect(r, p, n_plots = 1, p_new_dim = c(100, 100), pca = FALSE, method = spy)
  # Plot values extracted from the scaled raster should equal scaled plot values
  ext_scaled <- extractPlotMetrics(captured$r_pca, p)
  expect_equal(unname(captured$p_pca), unname(ext_scaled), tolerance = 1e-8, 
    ignore_attr = TRUE)
})

test_that("plotSelect meanmin result matches brute-force greedy on PCA scores", {
  old_ext <- extractPlotMetrics(r, p)
  old_pca <- PCALandscape(r, old_ext, center = TRUE, scale. = TRUE)
  X <- old_pca$r_pca$x[, 1:3]
  P <- old_pca$p_pca[, 1:3]
  out <- suppressMessages(plotSelect(r, p, n_plots = 3, p_new_dim = c(100, 100),
    n_pca = 3, method = meanminSelect))
  r_pca <- terra::rast(r, nlyrs = 3)
  v <- matrix(NA, terra::ncell(r), 3)
  v[stats::complete.cases(terra::values(r)), ] <- X
  terra::values(r_pca) <- v
  old_ind <- which(stats::complete.cases(terra::values(terra::mask(r, terra::vect(p)))))

  valid <- which(stats::complete.cases(v))
  blocks <- allBlocks(r_pca, 2, 2, setdiff(valid, old_ind))
  bm <- t(sapply(blocks, function(b) blockMean(r_pca, b)))
  # The population of possible plots is every complete block in the landscape
  X <- popMeans(r_pca, 2, 2)
  cur <- siteOptiSample:::nnDist(X, P)
  used <- NULL
  for (i in 1:3) {
    ok <- !sapply(blocks, function(b) any(b %in% used))
    D <- siteOptiSample:::distMat(X, bm[ok, , drop = FALSE])
    best <- which(ok)[which.min(colMeans(pmin(D, cur)))]
    cur <- pmin(cur, siteOptiSample:::distMat(X, bm[best, , drop = FALSE])[,1])
    used <- c(used, blocks[[best]])
    sel <- terra::cells(r, terra::vect(out[i, ]))[, "cell"]
    expect_setequal(sel, blocks[[best]])
  }
})

test_that("new plots do not overlap existing plots or each other", {
  for (fn in list(meanminSelect, minimaxSelect, kmeansSelect, hypercubeSelect)) {
    set.seed(1)
    out <- suppressMessages(plotSelect(r, p, n_plots = 10,
      p_new_dim = c(100, 100), n_pca = 3, method = fn))
    expect_equal(nrow(out), 10)
    p_t <- sf::st_transform(p, sf::st_crs(out))
    overlap_old <- sf::st_intersection(sf::st_union(out), sf::st_union(p_t))
    overlap_old_area <- if (length(overlap_old) == 0) 0 else sum(as.numeric(sf::st_area(overlap_old)))
    expect_equal(overlap_old_area, 0)
    pairs <- sf::st_overlaps(out, sparse = FALSE)
    expect_false(any(pairs))
  }
})

test_that("new plots fall within the mask", {
  r_mask <- r[[1]]
  r_mask[terra::values(r_mask) < 35] <- NA
  out <- suppressMessages(plotSelect(r, NULL, n_plots = 10,
    p_new_dim = c(100, 100), r_mask = r_mask, n_pca = 3))
  cells <- terra::cells(r_mask, terra::vect(out))[, "cell"]
  expect_false(any(is.na(terra::values(r_mask)[cells])))
})

test_that("raster and data frame inputs select the same locations", {
  r_df <- terra::values(r)
  valid_rows <- stats::complete.cases(r_df)
  r_df_all <- cbind(r_df[valid_rows, ], terra::crds(r))
  for (fn in list(meanminSelect, minimaxSelect)) {
    out_r <- suppressMessages(plotSelect(r, p, n_plots = 5,
      p_new_dim = c(100, 100), n_pca = 3, method = fn))
    out_df <- suppressMessages(plotSelect(r_df_all, p, n_plots = 5,
      p_new_dim = c(100, 100), n_pca = 3, coord = c("x", "y"), method = fn))
    cells_r <- terra::cells(r, terra::vect(out_r))[, "cell"]
    cells_df <- which(valid_rows)[out_df]
    expect_setequal(cells_df, cells_r)
  }
})

test_that("non-spatial input returns row indices matching brute-force greedy", {
  set.seed(60)
  df <- data.frame(a = rnorm(30), b = rnorm(30), c = rnorm(30))
  out <- suppressMessages(plotSelect(df, p = c(3, 7), p_type = "rows", n_plots = 3, pca = FALSE))
  X <- scale(df)
  cur <- siteOptiSample:::nnDist(X, X[c(3, 7), ])
  used <- c(3, 7)
  ref <- integer(0)
  for (i in 1:3) {
    cand <- setdiff(1:30, used)
    best <- cand[which.min(sapply(cand, function(j) {
      mean(pmin(cur, siteOptiSample:::distMat(X, X[j, , drop = FALSE])[,1]))
    }))]
    cur <- pmin(cur, siteOptiSample:::distMat(X, X[best, , drop = FALSE])[,1])
    used <- c(used, best)
    ref <- c(ref, best)
  }
  expect_equal(as.vector(out), ref)
})

test_that("Mahalanobis distance with PCA matches selection on rescaled PCA axes", {
  old_pca <- PCALandscape(r, center = TRUE, scale. = TRUE)$r_pca
  k <- 3
  v <- matrix(NA, terra::ncell(r), k)
  v[stats::complete.cases(terra::values(r)), ] <- sweep(old_pca$x[, 1:k], 2, 
    old_pca$sdev[1:k], "/")
  r_white <- terra::setValues(rep(r[[1]], k), v)
  names(r_white) <- paste0("PC", 1:k)

  for (fn in list(meanminSelect, minimaxSelect)) {
    out_m <- suppressMessages(plotSelect(r, p, n_plots = 5, p_new_dim = c(100, 100), 
      n_pca = k, distance = "mahalanobis", method = fn))
    out_w <- suppressMessages(plotSelect(r_white, p, n_plots = 5, 
      p_new_dim = c(100, 100), pca = FALSE, method = fn))
    expect_equal(sf::st_geometry(out_m), sf::st_geometry(out_w))
    expect_equal(out_m$mean_dist, out_w$mean_dist)
    expect_equal(attr(out_m, "baseline"), attr(out_w, "baseline"))

    # PCA values are reported in original units
    out_e <- suppressMessages(plotSelect(r, p, n_plots = 1, p_new_dim = c(100, 100), 
      n_pca = k, method = fn))
    cells <- terra::cells(r, terra::vect(out_m[1, ]))[, "cell"]
    pc_vals <- matrix(NA, terra::ncell(r), k)
    pc_vals[stats::complete.cases(terra::values(r)), ] <- old_pca$x[, 1:k]
    expect_equal(unname(unlist(sf::st_drop_geometry(out_m)[1, paste0("PC", 1:k)])), 
      unname(colMeans(pc_vals[cells, ])))
  }
})

test_that("Mahalanobis distance without PCA accounts for correlated variables", {
  vars <- names(r)[c(1, 2, 6)]
  r3 <- r[[vars]]
  out <- suppressMessages(plotSelect(r3, p, n_plots = 3, p_new_dim = c(100, 100), 
    pca = FALSE, distance = "mahalanobis"))

  # Recompute baseline with stats::mahalanobis on scaled values
  v <- terra::values(r3)
  ctr <- colMeans(v, na.rm = TRUE)
  sds <- apply(v, 2, stats::sd, na.rm = TRUE)
  r_scaled <- terra::setValues(r3, scale(v, ctr, sds))
  S <- stats::cov(stats::na.omit(scale(v, ctr, sds)))
  pop <- footprintMetrics(r_scaled, c(100, 100))
  P <- scale(extractPlotMetrics(r3, p), ctr, sds)
  d <- apply(sapply(seq_len(nrow(P)), function(j) {
    sqrt(stats::mahalanobis(pop, P[j, ], S))
  }), 1, min)
  expect_equal(attr(out, "baseline")[["mean_dist"]], mean(d))
  expect_equal(attr(out, "baseline")[["max_dist"]], max(d))

  # Euclidean distance differs because the variables are correlated
  out_e <- suppressMessages(plotSelect(r3, p, n_plots = 3, p_new_dim = c(100, 100), 
    pca = FALSE))
  expect_false(isTRUE(all.equal(attr(out, "baseline"), attr(out_e, "baseline"))))
})

test_that("invalid distance arguments are rejected", {
  expect_error(plotSelect(r, p, n_plots = 1, distance = "manhattan"))
  r_dup <- c(r[[1]], r[[1]] * 2)
  names(r_dup) <- c("a", "b")
  expect_error(suppressMessages(plotSelect(r_dup, n_plots = 1, pca = FALSE, 
    distance = "mahalanobis")), "singular")
})

test_that("existing plots can be supplied as structural values", {
  old_ext <- as.data.frame(extractPlotMetrics(r, p))
  spy_capture <- function(...) {
    captured <- NULL
    spy <- function(r_pca, p_pca, old_ind, n_plots, p_new_dim) {
      captured <<- list(p_pca = p_pca, old_ind = old_ind)
      list()
    }
    plotSelect(r, n_plots = 1, p_new_dim = c(100, 100), n_pca = 3, 
      method = spy, ...)
    captured
  }
  loc <- spy_capture(p = p)
  val <- spy_capture(p = old_ext[, rev(names(r))], p_type = "values")
  expect_equal(val$p_pca, loc$p_pca)
  expect_null(val$old_ind)
  
  # Non-spatial input, with an extra column which is ignored
  set.seed(61)
  df <- data.frame(a = rnorm(30), b = rnorm(30), c = rnorm(30))
  out_ind <- suppressMessages(plotSelect(df, p = c(3, 7), p_type = "rows",
    n_plots = 3, pca = FALSE))
  out_val <- suppressMessages(plotSelect(df, p = cbind(df[c(3, 7), ], id = 1:2),
    p_type = "values", n_plots = 3, pca = FALSE))
  expect_equal(as.vector(out_val), as.vector(out_ind))
})

test_that("p_type = 'values' rejects invalid values", {
  df <- data.frame(a = rnorm(10), b = rnorm(10))
  expect_error(plotSelect(df, p = data.frame(a = 1), p_type = "values", 
    n_plots = 1, pca = FALSE), "missing structural variables found in `r`: b")
  expect_error(plotSelect(df, p = c(1, 2), p_type = "values", n_plots = 1,
    pca = FALSE), "must be a dataframe or matrix")
})

test_that("point and xy existing plots give the same result", {
  pts <- suppressWarnings(sf::st_centroid(sf::st_geometry(p)))
  xy <- as.data.frame(sf::st_coordinates(pts))
  names(xy) <- c("x", "y")
  out_pt <- suppressMessages(plotSelect(r, pts, p_type = "point", n_plots = 3,
    p_new_dim = c(100, 100), n_pca = 3))
  out_xy <- suppressMessages(plotSelect(r, xy, p_type = "xy", n_plots = 3,
    p_new_dim = c(100, 100), n_pca = 3))
  expect_equal(sf::st_drop_geometry(out_xy), sf::st_drop_geometry(out_pt))
  expect_equal(attr(out_xy, "baseline"), attr(out_pt, "baseline"))

  # Coordinates of dataframe input
  r_df <- terra::values(r)
  valid_rows <- stats::complete.cases(r_df)
  r_df_all <- cbind(r_df[valid_rows, ], terra::crds(r))
  out_df <- suppressMessages(plotSelect(r_df_all, xy, p_type = "xy", 
    n_plots = 3, p_new_dim = c(100, 100), n_pca = 3, coord = c("x", "y")))
  expect_equal(attr(out_df, "baseline"), attr(out_pt, "baseline"))
})

test_that("p must match p_type", {
  df <- data.frame(a = rnorm(10), b = rnorm(10))
  pts <- suppressWarnings(sf::st_centroid(p))
  expect_error(plotSelect(r, pts, n_plots = 1), "containing polygons")
  expect_error(plotSelect(r, p, p_type = "point", n_plots = 1), 
    "containing points")
  expect_error(plotSelect(r, data.frame(a = 1), p_type = "xy", n_plots = 1), 
    "columns `x` and `y`")
  expect_error(plotSelect(r, c(1, 2), p_type = "rows", n_plots = 1), 
    "requires `r` to be a dataframe")
  expect_error(plotSelect(df, p, n_plots = 1, pca = FALSE), 
    "requires `r` to be spatial")
  expect_error(plotSelect(df, c(1, 11), p_type = "rows", n_plots = 1, 
    pca = FALSE), "valid row indices")
})

test_that("SpatVector existing plots match sf", {
  pts <- suppressWarnings(sf::st_centroid(p))
  for (args in list(list(p, "polygon"), list(pts, "point"))) {
    out_sf <- suppressMessages(plotSelect(r, args[[1]], p_type = args[[2]], 
      n_plots = 3, p_new_dim = c(100, 100), n_pca = 3))
    out_v <- suppressMessages(plotSelect(r, terra::vect(args[[1]]), 
      p_type = args[[2]], n_plots = 3, p_new_dim = c(100, 100), n_pca = 3))
    expect_equal(sf::st_drop_geometry(out_v), sf::st_drop_geometry(out_sf))
  }
})

test_that("existing plots with missing values give a warning", {
  vals <- as.data.frame(extractPlotMetrics(r, p))
  expect_warning(suppressMessages(plotSelect(r, rbind(vals, NA), 
    p_type = "values", n_plots = 1, p_new_dim = c(100, 100), n_pca = 3)),
    "1 existing plot\\(s\\) have missing")
  xy <- data.frame(x = 0, y = 0)
  expect_warning(suppressMessages(plotSelect(r, xy, p_type = "xy", 
    n_plots = 1, p_new_dim = c(100, 100), n_pca = 3)), "missing structural values")
})

test_that("a single existing plot location is accepted", {
  out <- suppressMessages(plotSelect(r, p[1, ], n_plots = 1, 
    p_new_dim = c(100, 100), n_pca = 3))
  expect_false(is.na(attr(out, "baseline")[["mean_dist"]]))
})
