test_that("ot_indices_div returns valid gsaot_indices object for solver sinkhorn", {
  set.seed(2026)
  dat <- gaussian_fun(1000)
  result <- ot_indices_div(dat[["x"]], dat[["y"]], M = 10, solver = "sinkhorn",
                           solver_optns = list(epsilon = 1))

  expect_s3_class(result, "gsaot_indices")
  expect_named(result, c("method", "indices", "bound", "x", "y",
                         "separation_measures", "partitions", "is_L22", "boot",
                         "solver_optns"))
  expect_true(all(c("X1", "X2", "X3") %in% names(result$indices)))
})

test_that("ot_indices_div returns valid gsaot_indices object for solver sinkhorn_stable", {
  dat <- gaussian_fun(1000)
  result <- ot_indices_div(dat[["x"]], dat[["y"]], M = 10, solver = "sinkhorn_stable",
                           solver_optns = list(epsilon = 1e3))

  expect_s3_class(result, "gsaot_indices")
  expect_named(result, c("method", "indices", "bound", "x", "y",
                         "separation_measures", "partitions", "is_L22", "boot",
                         "solver_optns"))
  expect_true(all(c("X1", "X2", "X3") %in% names(result$indices)))
})

test_that("ot_indices_div returns valid gsaot_indices object for bootstrap", {
  dat <- gaussian_fun(1000)
  result <- ot_indices_div(dat[["x"]], dat[["y"]], M = 10, solver = "sinkhorn",
                           solver_optns = list(epsilon = 1e3),
                           boot = TRUE, R = 100)

  expect_s3_class(result, "gsaot_indices")
  expect_named(result, c("method", "indices", "bound", "x", "y",
                         "separation_measures", "partitions", "is_L22", "boot",
                         "solver_optns", "indices_ci", "separation_measures_ci",
                         "bound_ci", "R", "type", "conf", "W_boot"))
  expect_true(all(c("X1", "X2", "X3") %in% names(result$indices)))
})
