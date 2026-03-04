test_that("irrelevance_threshold returns valid gsaot_indices object for solver sinkhorn", {
  dat <- gaussian_fun(1000)
  result <- irrelevance_threshold(dat[["y"]], M = 10, solver = "sinkhorn")

  expect_s3_class(result, "gsaot_indices")
  expect_named(result, c("method", "indices", "bound", "x", "y",
                         "separation_measures", "partitions", "boot",
                         "solver_optns"))
})

test_that("irrelevance_threshold returns valid gsaot_indices object for solver sinkhorn_log", {
  dat <- gaussian_fun(1000)
  result <- irrelevance_threshold(dat[["y"]], M = 10, solver = "sinkhorn_log")

  expect_s3_class(result, "gsaot_indices")
  expect_named(result, c("method", "indices", "bound", "x", "y",
                         "separation_measures", "partitions", "boot",
                         "solver_optns"))
})

test_that("irrelevance_threshold returns valid gsaot_indices object for solver transport", {
  dat <- gaussian_fun(1000)
  result <- irrelevance_threshold(dat[["y"]], M = 10, solver = "transport")

  expect_s3_class(result, "gsaot_indices")
  expect_named(result, c("method", "indices", "bound", "x", "y",
                         "separation_measures", "partitions", "boot",
                         "solver_optns"))
})

test_that("irrelevance_threshold returns valid gsaot_indices object for solver wb", {
  dat <- gaussian_fun(1000)
  result <- irrelevance_threshold(dat[["y"]], M = 10, solver = "wasserstein-bures")

  expect_s3_class(result, "gsaot_indices")
  expect_named(result, c("method", "indices", "bound", "x", "y",
                         "separation_measures", "partitions", "boot",
                         "adv", "diff"))
})

test_that("irrelevance_threshold returns valid gsaot_indices object for solver 1d", {
  dat <- ishi_homma_fun(1000)
  result <- irrelevance_threshold(dat[["y"]], M = 10, solver = "1d")

  expect_s3_class(result, "gsaot_indices")
  expect_named(result, c("method", "indices", "bound", "x", "y",
                         "separation_measures", "partitions", "boot"))
})
