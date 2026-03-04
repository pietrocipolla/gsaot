test_that("ot_indices_wb returns proper structure for univariate output", {
  dat <- sobol_fun(1000)
  result <- ot_indices_wb(dat[["x"]], dat[["y"]], M = 10)

  expect_s3_class(result, "gsaot_indices")
  expect_named(result, c("method", "indices", "bound", "x", "y",
                         "separation_measures", "partitions", "boot", "adv",
                         "diff"))
  expect_true(all(paste0("X", 1:8) %in% names(result$indices)))
})

test_that("ot_indices_wb returns proper structure for multivariate output", {
  dat <- gaussian_fun(1000)
  result <- ot_indices_wb(dat[["x"]], dat[["y"]], M = 10)

  expect_s3_class(result, "gsaot_indices")
  expect_named(result, c("method", "indices", "bound", "x", "y",
                         "separation_measures", "partitions", "boot", "adv",
                         "diff"))
  expect_true(all(paste0("X", 1:3) %in% names(result$indices)))
})

test_that("ot_indices_wb handles bootstrapping", {
  dat <- gaussian_fun(1000)
  result <- ot_indices_wb(dat[["x"]], dat[["y"]], M = 10, boot = T, R = 100)

  expect_s3_class(result, "gsaot_indices")
  expect_named(result, c("method", "indices", "bound", "x", "y",
                         "separation_measures", "partitions", "boot",
                         "adv", "diff",
                         "indices_ci", "separation_measures_ci",
                         "bound_ci", "R", "type", "conf", "W_boot"))
  expect_true(all(paste0("X", 1:3) %in% rownames(result$indices)))
})
