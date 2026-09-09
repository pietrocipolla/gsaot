test_that("partitions support all documented input column types", {
  factor_input <- factor(rep(c("low", "high"), 6),
                         levels = c("low", "unused", "high"))
  x <- data.frame(
    continuous = as.double(seq_len(12)),
    integer = rep(as.integer(c(20, 10, 30)), 4),
    character = rep(c("beta", "alpha", "gamma"), 4),
    factor = factor_input
  )

  partitions <- build_partition(x, M = 3)

  expect_equal(sort(unique(partitions[, 1])), 1:3)
  expect_equal(sort(unique(partitions[, 2])), 1:3)
  expect_equal(sort(unique(partitions[, 3])), 1:3)
  expect_equal(sort(unique(partitions[, 4])), 1:2)
  expect_false(anyNA(partitions))

  result <- ot_indices_1d(x, y = as.double(seq_len(12)), M = 3)
  expect_s3_class(result, "gsaot_indices")
  expect_named(result$indices, names(x))
  expect_true(all(is.finite(result$indices)))
})

test_that("unused factor levels do not create empty partitions", {
  x <- factor(rep(c("a", "c"), each = 6), levels = c("a", "b", "c"))

  partition <- build_discrete_partition(x)

  expect_equal(sort(unique(partition)), 1:2)
  expect_equal(tabulate(partition), c(6L, 6L))
})

test_that("sinkhorn_stable tau default is independent of option-list form", {
  defaults <- check_solver_optns("sinkhorn_stable", NULL)
  empty_options <- check_solver_optns("sinkhorn_stable", list())
  partial_options <- check_solver_optns("sinkhorn_stable", list(epsilon = 0.1))

  expect_equal(defaults$tau, empty_options$tau)
  expect_equal(defaults$tau, partial_options$tau)
  expect_equal(defaults$tau, 1e4)
})
