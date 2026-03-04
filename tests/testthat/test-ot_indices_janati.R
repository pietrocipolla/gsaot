test_that("ot_indices_janati returns proper structure for univariate output", {
  dat <- sobol_fun(1000)
  result <- ot_indices_janati(dat[["x"]], as.matrix(dat[["y"]]), M = 10, solver_optns = list(epsilon = 1e3))

  expect_s3_class(result, "gsaot_indices")
  expect_named(result, c("method", "indices", "bound", "x", "y",
                         "separation_measures", "partitions", "boot", "sobol"))
  expect_true(all(paste0("X", 1:8) %in% names(result$indices)))
})

test_that("ot_indices_janati returns proper structure for multivariate output", {
  dat <- gaussian_fun(1000)
  result <- ot_indices_janati(dat[["x"]], dat[["y"]], M = 10, solver_optns = list(epsilon = 1e3))

  expect_s3_class(result, "gsaot_indices")
  expect_named(result, c("method", "indices", "bound", "x", "y",
                         "separation_measures", "partitions", "boot", "sobol"))
  expect_true(all(paste0("X", 1:3) %in% names(result$indices)))
})

test_that("ot_indices_janati handles bootstrapping", {
  dat <- gaussian_fun(1000)
  result <- ot_indices_janati(dat[["x"]], dat[["y"]], M = 10, boot = T, R = 100, solver_optns = list(epsilon = 1e3))

  expect_s3_class(result, "gsaot_indices")
  expect_named(result, c("method", "indices", "bound", "x", "y",
                         "separation_measures", "partitions", "boot", "sobol",
                         "indices_ci", "separation_measures_ci",
                         "bound_ci", "R", "type", "conf", "W_boot"))
  expect_true(all(paste0("X", 1:3) %in% rownames(result$indices)))
})
