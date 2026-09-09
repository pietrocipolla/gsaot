test_that("ot_indices returns valid gsaot_indices object for solver sinkhorn", {
  dat <- gaussian_fun(1000)
  result <- ot_indices(dat[["x"]], dat[["y"]], M = 10, solver = "sinkhorn")

  expect_s3_class(result, "gsaot_indices")
  expect_named(result, c("method", "indices", "bound", "x", "y",
                         "separation_measures", "partitions", "is_L22", "boot",
                         "solver_optns"))
  expect_true(all(c("X1", "X2", "X3") %in% names(result$indices)))
})

test_that("ot_indices returns valid gsaot_indices object for solver sinkhorn_stable", {
  dat <- gaussian_fun(1000)
  result <- ot_indices(dat[["x"]], dat[["y"]], M = 10, solver = "sinkhorn_stable")

  expect_s3_class(result, "gsaot_indices")
  expect_named(result, c("method", "indices", "bound", "x", "y",
                         "separation_measures", "partitions", "is_L22", "boot",
                         "solver_optns"))
  expect_true(all(c("X1", "X2", "X3") %in% names(result$indices)))
})

test_that("sinkhorn_stable returns the full entropic dual objective", {
  a <- c(0.4, 0.6)
  b <- c(0.55, 0.45)
  cost_matrix <- matrix(c(0, 1, 1, 0), nrow = 2)

  sinkhorn_result <- gsaot:::sinkhorn(
    a, b, cost_matrix,
    numIterations = 1000,
    epsilon = 0.1,
    maxErr = 1e-12
  )
  stable_result <- gsaot:::sinkhorn_stable(
    a, b, cost_matrix,
    numIterations = 1000,
    epsilon = 0.1,
    maxErr = 1e-12,
    tau = 2
  )

  expect_equal(stable_result$cost, sinkhorn_result$cost, tolerance = 1e-10)
})

test_that("ot_indices returns valid gsaot_indices object for solver transport", {
  dat <- gaussian_fun(1000)
  result <- ot_indices(dat[["x"]], dat[["y"]], M = 10, solver = "transport")

  expect_s3_class(result, "gsaot_indices")
  expect_named(result, c("method", "indices", "bound", "x", "y",
                         "separation_measures", "partitions", "is_L22", "boot",
                         "solver_optns"))
  expect_true(all(c("X1", "X2", "X3") %in% names(result$indices)))
})

test_that("ot_indices returns valid gsaot_indices object for bootstrap", {
  dat <- gaussian_fun(1000)
  result <- ot_indices(dat[["x"]], dat[["y"]], M = 10, solver = "sinkhorn",
                       boot = TRUE, R = 100)

  expect_s3_class(result, "gsaot_indices")
  expect_named(result, c("method", "indices", "bound", "x", "y",
                         "separation_measures", "partitions", "is_L22", "boot",
                         "solver_optns", "indices_ci", "separation_measures_ci",
                         "bound_ci", "R", "type", "conf", "W_boot"))
  expect_true(all(c("X1", "X2", "X3") %in% names(result$indices)))
})
