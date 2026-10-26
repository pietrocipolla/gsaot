test_that("irrelevance_threshold returns valid gsaot_indices object for solver sinkhorn", {
  dat <- gaussian_fun(1000)
  result <- irrelevance_threshold(dat[["y"]], M = 10, solver = "sinkhorn")

  expect_type(result, "double")
})

test_that("irrelevance_threshold returns valid gsaot_indices object for solver sinkhorn_stable", {
  dat <- gaussian_fun(1000)
  result <- irrelevance_threshold(dat[["y"]], M = 10, solver = "sinkhorn_stable")

  expect_type(result, "double")
})

test_that("entropic_bound supports solver sinkhorn_stable", {
  set.seed(22)
  y <- cbind(rnorm(50), rnorm(50))

  result <- entropic_bound(
    y,
    M = 5,
    solver = "sinkhorn_stable",
    solver_optns = list(epsilon = 0.1, maxErr = 1e-7)
  )

  expect_type(result, "double")
  expect_true(is.finite(result))
})

test_that("irrelevance_threshold returns valid gsaot_indices object for solver transport", {
  dat <- gaussian_fun(1000)
  result <- irrelevance_threshold(dat[["y"]], M = 10, solver = "transport")

  expect_type(result, "double")
})

test_that("irrelevance_threshold returns valid gsaot_indices object for solver wb", {
  dat <- gaussian_fun(1000)
  result <- irrelevance_threshold(dat[["y"]], M = 10, solver = "wasserstein-bures")

  expect_type(result, "double")
})

test_that("irrelevance_threshold returns valid gsaot_indices object for solver 1d", {
  dat <- ishi_homma_fun(1000)
  result <- irrelevance_threshold(dat[["y"]], M = 10, solver = "1d")

  expect_type(result, "double")
})

test_that("higher_order_terms returns the difference between ot and wb indices", {
  dat <- gaussian_fun(1000)

  ot_result <- ot_indices(dat[["x"]], dat[["y"]], M = 10)
  wb_result <- ot_indices_wb(dat[["x"]], dat[["y"]], M = 10)
  result <- residual_gap(ot_result, wb_result)

  expect_s3_class(result, "gsaot_indices")
  expect_named(result$indices, names(ot_result$indices))
  expect_equal(result$indices, ot_result$indices - wb_result$indices)
})
