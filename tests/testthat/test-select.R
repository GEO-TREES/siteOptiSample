# Brute-force greedy meanmin: at each step, try every available block and keep
# the one giving the lowest mean pixel-to-nearest-plot distance
refMeanmin <- function(r, p_pca, n_plots, dims, old_ind = NULL) {
  sel <- list()
  occupied <- old_ind
  plot_means <- p_pca
  valid <- which(stats::complete.cases(terra::values(r)))
  pop <- popMeans(r, dims[1], dims[2])
  for (i in seq_len(n_plots)) {
    blocks <- allBlocks(r, dims[1], dims[2], setdiff(valid, occupied))
    if (length(blocks) == 0) break
    scores <- sapply(blocks, function(b) {
      meanNearest(pop, rbind(plot_means, blockMean(r, b)))
    })
    best <- blocks[[which.min(scores)]]
    sel[[i]] <- best
    occupied <- c(occupied, best)
    plot_means <- rbind(plot_means, blockMean(r, best))
  }
  sel
}

# Brute-force greedy minimax: at each step, keep the available block whose
# mean is furthest from its nearest plot
refMinimax <- function(r, p_pca, n_plots, dims) {
  sel <- list()
  occupied <- NULL
  plot_means <- p_pca
  valid <- which(stats::complete.cases(terra::values(r)))
  for (i in seq_len(n_plots)) {
    blocks <- allBlocks(r, dims[1], dims[2], setdiff(valid, occupied))
    if (length(blocks) == 0) break
    scores <- sapply(blocks, function(b) {
      bm <- blockMean(r, b)
      min(apply(plot_means, 1, function(pm) d2(bm, pm)))
    })
    best <- blocks[[which.max(scores)]]
    sel[[i]] <- best
    occupied <- c(occupied, best)
    plot_means <- rbind(plot_means, blockMean(r, best))
  }
  sel
}

test_that("meanmin matches brute-force greedy selection", {
  set.seed(10)
  r <- makeRast(matrix(rnorm(128), ncol = 2), 8, 8)
  p_pca <- matrix(c(1, 1, -1, 0.5), ncol = 2, byrow = TRUE)
  for (dims in list(c(1, 1), c(2, 2), c(2, 3))) {
    for (pp in list(NULL, p_pca)) {
      out <- suppressMessages(meanminSelect(r, p_pca = pp, n_plots = 4,
        p_new_dim = dims))
      expect_equal(selCells(out), refMeanmin(r, pp, 4, dims))
    }
  }
})

test_that("meanmin respects occupied cells", {
  set.seed(11)
  r <- makeRast(matrix(rnorm(128), ncol = 2), 8, 8)
  old_ind <- c(1:4, 9:12)
  out <- suppressMessages(meanminSelect(r, old_ind = old_ind, n_plots = 3,
    p_new_dim = c(2, 2)))
  expect_equal(selCells(out), refMeanmin(r, NULL, 3, c(2, 2), old_ind))
  expect_false(any(unlist(selCells(out)) %in% old_ind))
})

test_that("meanmin is close to the exhaustive optimum and beats most placements", {
  # Greedy selection is not guaranteed to be optimal: across random rasters
  # it is typically 7-10% above the optimum for 3 plots
  set.seed(12)
  r <- makeRast(matrix(rnorm(128), ncol = 2), 8, 8)
  blocks <- allBlocks(r, 2, 2)
  bm <- t(sapply(blocks, function(b) blockMean(r, b)))
  # The population of possible plots is every block
  D <- as.matrix(stats::dist(bm))
  objective <- function(idx) mean(apply(D[, idx, drop = FALSE], 1, min))

  # Exhaustive search over non-overlapping sets of 3 blocks
  combos <- utils::combn(length(blocks), 3)
  no_overlap <- apply(combos, 2, function(idx) {
    !anyDuplicated(unlist(blocks[idx]))
  })
  all_obj <- apply(combos[, no_overlap], 2, objective)

  out <- suppressMessages(meanminSelect(r, n_plots = 3, p_new_dim = c(2, 2)))
  greedy_idx <- match(sapply(selCells(out), paste, collapse = ":"), 
    sapply(blocks, paste, collapse = ":"))
  greedy <- objective(greedy_idx)

  expect_gte(greedy, min(all_obj) - 1e-12)
  expect_lt(greedy, min(all_obj) * 1.12)
  expect_gt(mean(all_obj > greedy), 0.85)

  # Refinement typically brings the result to within ~1.5% of the optimum
  out_ref <- suppressMessages(meanminSelect(r, n_plots = 3, p_new_dim = c(2, 2), 
    refine = TRUE))
  refined_idx <- match(sapply(selCells(out_ref), paste, collapse = ":"), 
    sapply(blocks, paste, collapse = ":"))
  refined <- objective(refined_idx)
  expect_lte(refined, greedy)
  expect_gte(refined, min(all_obj) - 1e-12)
  expect_lt(refined, min(all_obj) * 1.02)
})

test_that("refined meanmin is a local optimum: no single swap improves it", {
  set.seed(13)
  r <- makeRast(matrix(rnorm(200), ncol = 2), 10, 10)
  p_pca <- matrix(c(1, 1), ncol = 2)
  old_ind <- c(1, 2, 11, 12)
  out <- suppressMessages(meanminSelect(r, p_pca = p_pca, old_ind = old_ind, 
    n_plots = 4, p_new_dim = c(2, 2), refine = TRUE))
  sel <- selCells(out)
  expect_false(any(unlist(sel) %in% old_ind))
  expect_false(anyDuplicated(unlist(sel)) > 0)

  sel_means <- t(sapply(sel, function(b) blockMean(r, b)))
  pop <- popMeans(r, 2, 2)
  current <- meanNearest(pop, rbind(p_pca, sel_means))
  for (i in seq_along(sel)) {
    free <- setdiff(seq_len(100), c(old_ind, unlist(sel[-i])))
    for (b in allBlocks(r, 2, 2, free)) {
      swapped <- sel_means
      swapped[i, ] <- blockMean(r, b)
      expect_gte(meanNearest(pop, rbind(p_pca, swapped)), current - 1e-9)
    }
  }
})

test_that("plotSelect passes refine through to meanminSelect", {
  set.seed(14)
  df <- data.frame(a = rnorm(40), b = rnorm(40))
  greedy <- suppressMessages(plotSelect(df, n_plots = 3, pca = FALSE))
  refined <- suppressMessages(plotSelect(df, n_plots = 3, pca = FALSE, refine = TRUE))
  X <- scale(df)
  obj <- function(idx) mean(siteOptiSample:::nnDist(X, X[idx, , drop = FALSE]))
  expect_lte(obj(refined), obj(greedy))
})

test_that("minimax matches brute-force greedy selection", {
  set.seed(20)
  r <- makeRast(matrix(rnorm(128), ncol = 2), 8, 8)
  p_pca <- matrix(c(0, 0), ncol = 2)
  for (dims in list(c(1, 1), c(2, 2), c(3, 2))) {
    out <- suppressMessages(minimaxSelect(r, p_pca = p_pca, n_plots = 4,
      p_new_dim = dims))
    expect_equal(selCells(out), refMinimax(r, p_pca, 4, dims))
  }
})

test_that("minimax picks a planted structural outlier first", {
  set.seed(21)
  vals <- matrix(rnorm(200, sd = 0.1), ncol = 2)
  r <- makeRast(vals, 10, 10)
  outlier <- c(34, 35, 44, 45)
  r[outlier] <- matrix(10, 4, 2)
  for (pp in list(NULL, matrix(c(0, 0), ncol = 2))) {
    out <- suppressMessages(minimaxSelect(r, p_pca = pp, n_plots = 2,
      p_new_dim = c(2, 2)))
    expect_equal(selCells(out)[[1]], outlier)
  }
})

test_that("kmeans places one plot in each planted structural cluster", {
  set.seed(30)
  ncol <- 12
  col <- rep(seq_len(ncol), times = 12)
  group <- ceiling(col / 4)
  vals <- cbind(group * 5, -group * 5) + matrix(rnorm(288, sd = 0.1), ncol = 2)
  r <- makeRast(vals, 12, ncol)
  out <- kmeansSelect(r, n_plots = 3, p_new_dim = c(2, 2), het_q = NULL)
  groups <- sapply(selCells(out), function(x) unique(group[x]))
  expect_type(groups, "double")
  expect_setequal(groups, 1:3)
})

test_that("hypercube targets are spread across quantiles of the gradient", {
  set.seed(40)
  r <- makeRast(matrix(rep(1:20, times = 10), ncol = 1), 10, 20)
  out <- hypercubeSelect(r, n_plots = 4, p_new_dim = c(1, 1), het_q = NULL)
  sel_vals <- sort(sapply(selCells(out), function(x) terra::values(r)[x, 1]))
  targets <- stats::quantile(rep(1:20, times = 10), (1:4 - 0.5) / 4)
  expect_true(all(abs(sel_vals - targets) <= 0.5))
})

test_that("all algorithms return non-overlapping plots of the right size", {
  set.seed(50)
  r <- makeRast(matrix(rnorm(400), ncol = 2), 10, 20)
  for (fn in list(meanminSelect, minimaxSelect, kmeansSelect, hypercubeSelect)) {
    out <- suppressMessages(fn(r, n_plots = 6, p_new_dim = c(3, 2), het_q = NULL))
    cells <- selCells(out)
    expect_length(cells, 6)
    expect_true(all(lengths(cells) == 6))
    expect_false(anyDuplicated(unlist(cells)) > 0)
    areas <- sapply(out, function(x) as.numeric(sf::st_area(x)))
    expect_equal(unname(areas), rep(6, 6))
  }
})
