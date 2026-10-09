test_that("euclidean distances match stats::dist", {
  set.seed(1)
  x <- matrix(rnorm(30), 10, 3)
  y <- matrix(rnorm(12), 4, 3)
  ref <- as.matrix(stats::dist(rbind(x, y)))[1:10, 11:14]
  expect_equal(siteOptiSample:::euclideanDist(x, y), unname(ref))
  expect_equal(siteOptiSample:::nnDist(x, y), unname(apply(ref, 1, min)))
  expect_equal(siteOptiSample:::nnDist(x, y, chunk_size = 10),
    unname(apply(ref, 1, min)))
})

test_that("weighted euclidean distances weight squared differences", {
  x <- matrix(c(0, 0), 1)
  y <- matrix(c(1, 2), 1)
  expect_equal(siteOptiSample:::euclideanDist(x, y, w = c(4, 1))[1, 1], sqrt(4 + 4))
})

test_that("mahalanobis distances match stats::mahalanobis", {
  set.seed(2)
  x <- matrix(rnorm(60), 20, 3)
  y <- matrix(rnorm(9), 3, 3)
  S <- stats::cov(rbind(x, y))
  ref <- sapply(1:3, function(j) sqrt(stats::mahalanobis(x, y[j, ], S)))
  expect_equal(siteOptiSample:::mahalanobisDist(x, y), ref)
})

test_that("pcaDist returns ordered nearest neighbour distances", {
  set.seed(3)
  x <- matrix(rnorm(40), 10, 4)
  y <- matrix(rnorm(20), 5, 4)
  ref <- t(apply(as.matrix(stats::dist(rbind(x[,1:2], y[,1:2])))[1:10, 11:15], 1, sort))
  expect_equal(pcaDist(x, y, n_pca = 2, k = 1), unname(ref[,1]))
  expect_equal(pcaDist(x, y, n_pca = 2, k = 3), unname(ref[,1:3]))
  expect_error(pcaDist(x, y, n_pca = 5))
  expect_error(pcaDist(x, y, k = 6))
})
