test_that("ot_indices returns valid gsaot_indices object for solver sinkhorn", {
  dat <- gaussian_fun(1000)
  result <- ot_indices(dat[["x"]], dat[["y"]], M = 10, solver = "sinkhorn")

  expect_s3_class(result, "gsaot_indices")
  expect_named(result, c("method", "indices", "bound", "x", "y",
                         "separation_measures", "partitions", "boot",
                         "solver_optns"))
  expect_true(all(c("X1", "X2", "X3") %in% names(result$indices)))
})

test_that("ot_indices returns valid gsaot_indices object for solver sinkhorn_log", {
  dat <- gaussian_fun(1000)
  result <- ot_indices(dat[["x"]], dat[["y"]], M = 10, solver = "sinkhorn_log")

  expect_s3_class(result, "gsaot_indices")
  expect_named(result, c("method", "indices", "bound", "x", "y",
                         "separation_measures", "partitions", "boot",
                         "solver_optns"))
  expect_true(all(c("X1", "X2", "X3") %in% names(result$indices)))
})

test_that("ot_indices returns valid gsaot_indices object for solver transport", {
  dat <- gaussian_fun(1000)
  result <- ot_indices(dat[["x"]], dat[["y"]], M = 10, solver = "transport")

  expect_s3_class(result, "gsaot_indices")
  expect_named(result, c("method", "indices", "bound", "x", "y",
                         "separation_measures", "partitions", "boot",
                         "solver_optns"))
  expect_true(all(c("X1", "X2", "X3") %in% names(result$indices)))
})

test_that("ot_indices returns valid gsaot_indices object for bootstrap", {
  dat <- gaussian_fun(1000)
  result <- ot_indices(dat[["x"]], dat[["y"]], M = 10, solver = "sinkhorn",
                       boot = TRUE, R = 100)

  expect_s3_class(result, "gsaot_indices")
  expect_named(result, c("method", "indices", "bound", "x", "y",
                         "separation_measures", "partitions", "boot",
                         "solver_optns", "indices_ci", "separation_measures_ci",
                         "bound_ci", "R", "type", "conf", "W_boot"))
  expect_true(all(c("X1", "X2", "X3") %in% names(result$indices)))
})
