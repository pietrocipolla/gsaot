test_that("squared Euclidean cost matches dist for continuous and repeated outputs", {
  set.seed(19)
  for (n_columns in c(1L, 3L)) {
    y <- matrix(rnorm(90L * n_columns), ncol = n_columns)
    y[2L, ] <- y[1L, ]
    actual <- squared_euclidean_cost(y)
    expected <- as.matrix(stats::dist(y))^2

    expect_equal(unname(actual), unname(expected), tolerance = 1e-12)
    expect_true(all(actual >= 0))
    expect_true(all(diag(actual) == 0))
  }
})
