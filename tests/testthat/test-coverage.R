test_that("ot_indices validates inputs", {
  x <- data.frame(x1 = 1:10, x2 = 11:20)
  y <- matrix(rnorm(10), ncol = 1)

  expect_error(ot_indices(x, y, M = 20), "number of partitions")
  expect_error(ot_indices(x, y[-1, , drop = FALSE], M = 2), "same")
  expect_error(ot_indices(x, matrix(letters[1:10]), M = 2), "`y` must")
  expect_error(ot_indices(x, y, M = 2, cost = "bad"), "cost")
  expect_error(ot_indices(x, y, M = 2, boot = TRUE, R = NULL), "Bootstrapping")
  expect_error(ot_indices(x, y, M = 2, boot = FALSE, R = 10), "Bootstrapping")
})

test_that("ot_indices_1d validates inputs", {
  x <- data.frame(x1 = 1:10, x2 = 11:20)
  y <- rnorm(10)

  expect_error(ot_indices_1d(x, y, M = 20), "number of partitions")
  expect_error(ot_indices_1d(x, y[-1], M = 2), "same")
  expect_error(ot_indices_1d(x, matrix(y, ncol = 1), M = 2), "`y` must")
  expect_error(ot_indices_1d(x, y, M = 2, p = 0.5), "p")
  expect_error(ot_indices_1d(x, y, M = 2, boot = TRUE, R = NULL), "Bootstrapping")
})

test_that("ot_indices_wb validates inputs", {
  x <- data.frame(x1 = 1:10, x2 = 11:20)
  y <- matrix(rnorm(20), ncol = 2)

  expect_error(ot_indices_wb(x, y, M = 20), "number of partitions")
  expect_error(ot_indices_wb(x, y[-1, , drop = FALSE], M = 2), "same")
  expect_error(ot_indices_wb(x, y, M = 2, boot = TRUE, R = NULL), "Bootstrapping")
})

test_that("indices rank dominant inputs", {
  set.seed(10)
  N <- 200
  x1 <- rnorm(N)
  x2 <- rnorm(N)
  x <- data.frame(x1 = x1, x2 = x2)

  y1 <- 3 * x1 + 0.1 * rnorm(N)
  res_1d <- ot_indices_1d(x, y1, M = 5)
  expect_gt(res_1d$indices["x1"], res_1d$indices["x2"])

  y2 <- cbind(2 * x1 + 0.1 * rnorm(N), 0.5 * x1 + rnorm(N))
  res_md <- ot_indices(x, y2, M = 5, solver = "transport")
  expect_gt(res_md$indices["x1"], res_md$indices["x2"])
})

test_that("ot_indices supports custom cost and scaling", {
  set.seed(11)
  N <- 1000
  x <- data.frame(x1 = rnorm(N), x2 = rnorm(N))
  y <- cbind(2 * x$x1 + rnorm(N), rnorm(N))

  cost_fun <- function(y) as.matrix(stats::dist(y))^2

  res <- ot_indices(
    x, y, M = 5,
    cost = cost_fun,
    scaling = FALSE,
    solver = "sinkhorn",
    solver_optns = list(epsilon = 1, numIterations = 1e5)
  )

  expect_true(all(is.finite(res$indices)))
  expect_true(all(res$indices >= -1e-8 & res$indices <= 1 + 1e-8))
})

test_that("bootstrap CIs are consistent and confint works", {
  set.seed(12)
  N <- 120
  x <- data.frame(x1 = rnorm(N), x2 = rnorm(N))
  y <- cbind(2 * x$x1 + rnorm(N), rnorm(N))

  res <- ot_indices(
    x, y, M = 5,
    solver = "sinkhorn",
    solver_optns = list(epsilon = 0.1, numIterations = 500),
    boot = TRUE, R = 20
  )

  ci <- res$indices_ci
  idx <- res$indices[ci$input]

  expect_true(all(ci$low.ci <= idx + 1e-8))
  expect_true(all(ci$high.ci >= idx - 1e-8))
  expect_true(all(ci$low.ci <= ci$high.ci))

  ci2 <- confint(res, level = 0.9)
  expect_true(all(c("input", "index", "original", "bias", "low.ci", "high.ci") %in% names(ci2)))
})

test_that("S3 methods run and return expected types", {
  set.seed(13)
  N <- 120
  x <- data.frame(x1 = rnorm(N), x2 = rnorm(N))
  y <- cbind(2 * x$x1 + rnorm(N), rnorm(N))

  res <- ot_indices(x, y, M = 5)

  expect_output(print(res), "Method")
  expect_true(is.list(summary(res)))
  expect_s3_class(plot(res), "ggplot")
  expect_true(inherits(plot_separations(res, ranking = 1), "patchwork"))
})

test_that("Wasserstein-Bures decomposition is coherent", {
  set.seed(14)
  N <- 120
  x <- data.frame(x1 = rnorm(N), x2 = rnorm(N))
  y <- cbind(2 * x$x1 + rnorm(N), rnorm(N))

  res <- ot_indices_wb(x, y, M = 5)

  expect_true(all(is.finite(res$indices)))
  expect_true(all(abs(res$indices - (res$adv + res$diff)) < 1e-3))
})