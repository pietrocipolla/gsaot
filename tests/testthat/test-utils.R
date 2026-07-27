test_that("irrelevance_threshold returns valid gsaot_indices object for solver sinkhorn", {
  dat <- gaussian_fun(1000)
  result <- irrelevance_threshold(dat[["y"]], M = 10, solver = "sinkhorn")

  expect_type(result, "double")
})

test_that("irrelevance_threshold returns valid gsaot_indices object for solver sinkhorn_log", {
  dat <- gaussian_fun(1000)
  result <- irrelevance_threshold(dat[["y"]], M = 10, solver = "sinkhorn_log")

  expect_type(result, "double")
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
  result <- higher_order_terms(ot_result, wb_result)

  expect_s3_class(result, "gsaot_indices")
  expect_named(result$indices, names(ot_result$indices))
  expect_equal(result$indices, ot_result$indices - wb_result$indices)
})
