test_that("threshold solvers forward custom costs and logical options correctly", {
  set.seed(82)
  y <- cbind(rnorm(30), rnorm(30))
  cost_matrix <- as.matrix(dist(y))^2 / 10
  custom_cost <- function(y) cost_matrix
  for (solver in c("transport", "sinkhorn", "sinkhorn_stable", "sinkhorn_div")) {
    options <- if (solver == "transport") NULL else
      list(epsilon = 1, numIterations = 1000, maxErr = 1e-6)
    set.seed(13)
    actual <- irrelevance_threshold(
      y, M = 3, cost = custom_cost, discrete_out = FALSE,
      solver = solver, solver_optns = options, scaling = FALSE, R_irr = 2
    )

    set.seed(13)
    expected <- mean(replicate(2, {
      x <- data.frame(rnorm = rnorm(nrow(y)))
      estimator <- if (solver == "sinkhorn_div") ot_indices_div else ot_indices
      estimator(
        x, y, M = 3, cost = custom_cost, discrete_out = FALSE,
        solver = if (solver == "sinkhorn_div") "sinkhorn" else solver,
        solver_optns = options, scaling = FALSE
      )$indices
    }))

    expect_true(is.finite(actual), info = solver)
    expect_equal(actual, expected, info = solver)
  }
})
