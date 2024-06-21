#' Calculate lower bounds for Optimal Transport sensitivity indices
#'
#' @inheritParams ot_indices
#' @param y An array or a matrix containing the output values.
#' @param bound (default `"dummy"`) A string defining the type of lower bound to
#'   compute. Should be `"dummy"` or `"entropic"`. See `details` for more
#'   information.
#' @param dummy_optns (default `NULL`) A list containing the options on the
#'   distribution of the dummy variable. See `details` for more information.
#'
#' @details The function allows the computation of two different lower bounds.
#'   With `bound="dummy"`, the function samples from a distribution defined in
#'   `dummy_optns` (by default a standard normal), independent from the output
#'   `y` and then computes the indices using the algorithm specified in
#'   `solver`. Under the hood, `lower_bound` calls the other available functions
#'   in the package:
#' * [ot_indices_1d()] (for `solver="1d`)
#' * [ot_indices_wb()] (for `solver="wasserstein-bures"`)
#' * [ot_indices()] (for `solver %in% c("sinkhorn", "sinkhorn_log", "wasserstein")`)
#'   The user can choose the distribution of the dummy variable using the
#'   argument `dummy_optns`. `dummy_optns` should be a named list with at least
#'   a term called `"distr"` defining the sampling function. The other terms in
#'   the list are used as arguments to the sampling function. With
#'   `bound="entropic"`, the function computes the lower bound of the entropic
#'   indices. In this case, `solver` should be either `"sinkhorn"` or
#'   `"sinkhorn_log"`
#'
#' @return An object of class `gsaot_indices` for `bound="dummy"` or a scalar
#'   for `bound="entropic"`.
#' @export
#'
#' @examples
#' N <- 1000
#'
#' mx <- c(1, 1, 1)
#' Sigmax <- matrix(data = c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), nrow = 3)
#'
#' x1 <- rnorm(N)
#' x2 <- rnorm(N)
#' x3 <- rnorm(N)
#'
#' x <- cbind(x1, x2, x3)
#' x <- mx + x %*% chol(Sigmax)
#'
#' A <- matrix(data = c(4, -2, 1, 2, 5, -1), nrow = 2, byrow = TRUE)
#' y <- t(A %*% t(x))
#'
#' M <- 25
#'
#' sink_lb <- lower_bound(y, M, bound = "entropic")
#' dummy_lb <- lower_bound(y, M, bound = "dummy")
#'
#' # Custom sampling funtion and network simplex solver
#' dummy_optns <- list(distr = "rgamma", shape = 3)
#' dummy_lb_cust <- lower_bound(y, M, bound = "dummy",
#'                                   dummy_optns = dummy_optns,
#'                                   solver = "wasserstein")
lower_bound <- function(y,
                        M,
                        bound = "dummy",
                        dummy_optns = NULL,
                        cost = "L2",
                        discrete_out = FALSE,
                        solver = "sinkhorn",
                        solver_optns = NULL,
                        scaling = TRUE) {
  # INPUT CHECKS
  # ----------------------------------------------------------------------------
  # Check if the output is a numerical
  if (!is.numeric(y) | (is.matrix(y) & is.vector(y)))
    stop("`y` must be a vector or a matrix of numerical values!")

  # Check if the number of partitions is lower than the number of samples
  if (nrow(as.matrix(y)) <= M)
    stop("The number of partitions should be lower than the number of samples")

  # Check if the bound is in the pool
  match.arg(bound, c("dummy", "entropic"))

  # Sinkhorn is allowed only for sinkhorn or sinkhorn log solver
  if (bound == "entropic" & !(solver %in% c("sinkhorn", "sinkhorn_log")))
    stop("`sinkhorn` bound is allowed only for `sinkhorn` or `sinkhorn_log` solvers")

  # Check if the cost is defined correctly
  cost_type <- identical(cost, "L2")
  if (!cost_type & !is.function(cost))
    stop("`cost` should be \"L2\" or a function")

  # Check the logical values
  if (!is.logical(discrete_out)) stop("`discrete_out` should be logical")
  if (!is.logical(scaling)) stop("`scaling` should be logical")

  # Check if the solver is present in the pool
  match.arg(solver, c("1d", "wasserstein-bures", "sinkhorn", "sinkhorn_log", "wasserstein"))

  # RETURN THE DUMMY INDICES
  # ----------------------------------------------------------------------------
  if (bound == "dummy") {
    dummy_optns <- init_dummy_optns(dummy_optns)

    N <- nrow(as.matrix(y))
    x <- do.call(dummy_optns[["distr"]], c(n = N, within(dummy_optns, rm("distr"))))
    x <- as.data.frame(x)
    colnames(x) <- dummy_optns[["distr"]]

    dummy_ind <- switch(solver,
                        "1d" = ot_indices_1d(x, y, M),
                        "wasserstein-bures" = ot_indices_wb(x, y, M),
                        "sinkhorn" = ot_indices(x,
                                                y,
                                                M,
                                                cost,
                                                discrete_out,
                                                solver,
                                                solver_optns,
                                                scaling),
                        "sinkhorn_log" = ot_indices(x,
                                                    y,
                                                    M,
                                                    cost,
                                                    discrete_out,
                                                    solver,
                                                    solver_optns,
                                                    scaling),
                        "wasserstein" = ot_indices(x,
                                                   y,
                                                   M,
                                                   cost,
                                                   discrete_out,
                                                   solver,
                                                   solver_optns,
                                                   scaling),
                        default = NULL
    )

    return(dummy_ind)
  }

  # RETURN THE SINKHORN LOWER BOUND
  sink_ind <- entropic_lower_bound(y, cost, solver, solver_optns, discrete_out, scaling)

  return(sink_ind)
}

# Check that the bound_optns contain the required components, otherwise initialize them
init_dummy_optns <- function(dummy_optns) {
  if (is.null(dummy_optns)) {
    dummy_optns <- list(distr = "rnorm")
    cat("Using `rnorm` to generate the dummy\n")
  }

  return(dummy_optns)
}

# Compute the sinkhorn lower bound
entropic_lower_bound <- function(y,
                                 cost,
                                 solver,
                                 solver_optns,
                                 discrete_out,
                                 scaling) {
  # Initialize useful values
  N <- nrow(as.matrix(y))
  cost_type <- identical(cost, "L2")

  # BUILD AN HISTOGRAM IF REQUESTED (otherwise assume equal weights)
  # ----------------------------------------------------------------------------
  # Default values for continuous output
  a <- rep(1 / N, N)
  y_unique <- NULL

  if (discrete_out) {
    # y_unique is always needed, in order to avoid mismatches between the
    # pre-computed cost matrix and the resampled outputs
    y_unique <- unique(y)

    # Function to build the histogram (maybe we can put this outside to make it
    # easier to change?)
    a <- apply(y_unique, MARGIN = 1, FUN = function(u) sum(
      apply(y, 1, FUN = function(row) all(row == u)))) / N
  }

  if (identical(a, rep(1 / N, N)) & discrete_out)
    warning("The output is continuous, consider using `discrete_out=FALSE`\n")

  # BUILD COST MATRIX
  # ----------------------------------------------------------------------------
  if (discrete_out) {
    if (cost_type) {
      C <- as.matrix(stats::dist(y_unique, method = "euclidean")) ^ 2
    } else {
      C <- cost(y_unique)
    }
  } else {
    if (cost_type) {
      C <- as.matrix(stats::dist(y, method = "euclidean")) ^ 2
    } else {
      C <- cost(y)
    }
  }

  # If the scaling is requested, scale the cost matrix by the highest value
  scaling_param <- 1
  if (scaling) {
    scaling_param <- max(C)
    C <- C / scaling_param
  }

  # DEFINE THE OT SOLVER
  # ----------------------------------------------------------------------------
  solver_fun <- switch (
    solver,
    "sinkhorn" = sinkhorn,
    "sinkhorn_log" = sinkhorn_log,
    default = NULL
  )

  # Check the consistency of the options provided
  solver_optns <- check_solver_optns(solver, solver_optns)

  # Evaluate the upper bound of the indices
  # ----------------------------------------------------------------------------
  V <- higher_bound(C) * scaling_param

  # COMPUTE THE LOWER BOUND
  # ----------------------------------------------------------------------------
  ret <- do.call(solver_fun,
                 c(list(a = a, b = a, costm = C),
                   solver_optns))

  sink_ind <- ret$cost * scaling_param / V

  return(sink_ind)
}
