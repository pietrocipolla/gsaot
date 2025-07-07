#' Entropic lower bounds for entropic optimal transport sensitivity indices
#'
#' Calculate entropic lower bounds for entropic Optimal Transport sensitivity indices
#'
#' @inheritParams ot_indices
#' @param y An array or a matrix containing the output values.
#' @param solver Solver for the Optimal Transport problem. Currently supported
#'   options are:
#' * `"sinkhorn"` (default), the Sinkhorn's solver \insertCite{cuturi2013sinkhorn}{gsaot}.
#' * `"sinkhorn_log"`, the Sinkhorn's solver in log scale \insertCite{peyre2019computational}{gsaot}.
#'
#' @details The function allows the computation of the entropic lower bounds.
#' `solver` should be either `"sinkhorn"` or `"sinkhorn_log"`.
#'
#' @return A scalar representing the entropic lower bound.
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
#' sink_lb <- entropic_bound(y, M)
#'
entropic_bound <- function(y,
                           M,
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

  # Check if the cost is defined correctly
  cost_type <- identical(cost, "L2")
  if (!cost_type & !is.function(cost))
    stop("`cost` should be \"L2\" or a function")

  # Check the logical values
  if (!is.logical(discrete_out)) stop("`discrete_out` should be logical")
  if (!is.logical(scaling)) stop("`scaling` should be logical")

  # Check if the solver is present in the pool
  match.arg(solver, c("sinkhorn", "sinkhorn_log"))

  # RETURN THE SINKHORN LOWER BOUND
  sink_ind <- entropic_lower_bound(y, cost, solver, solver_optns, discrete_out, scaling)

  return(sink_ind)
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
    "sinkhorn_stable" = sinkhorn_stable,
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

#' Irrelevance threshold for optimal transport sensitivity indices
#'
#' Calculate irrelevance threshold using dummy variable for Optimal Transport sensitivity indices
#'
#' @inheritParams ot_indices
#' @param y An array or a matrix containing the output values.
#' @param solver Solver for the Optimal Transport problem. Currently supported
#'   options are:
#' * `"1d"`, the one-dimensional analytic solution.
#' * `"wasserstein-bures"`, the Wasserstein-Bures solution.
#' * `"sinkhorn"` (default), the Sinkhorn's solver \insertCite{cuturi2013sinkhorn}{gsaot}.
#' * `"sinkhorn_log"`, the Sinkhorn's solver in log scale \insertCite{peyre2019computational}{gsaot}.
#' * `"transport"`, a solver of the non regularized OT problem using [transport::transport()].
#' @param dummy_optns (default `NULL`) A list containing the options on the
#'   distribution of the dummy variable. See `details` for more information.
#'
#' @details The function allows the computation of irrelevance threshold.
#'   The function samples from a distribution defined in
#'   `dummy_optns` (by default a standard normal), independent from the output
#'   `y` and then computes the indices using the algorithm specified in
#'   `solver`. Under the hood, `lower_bound` calls the other available functions
#'   in the package:
#' * [ot_indices_1d()] (for `solver="1d"`)
#' * [ot_indices_wb()] (for `solver="wasserstein-bures"`)
#' * [ot_indices()] (for `solver %in% c("sinkhorn", "sinkhorn_log", "wasserstein")`)
#'   The user can choose the distribution of the dummy variable using the
#'   argument `dummy_optns`. `dummy_optns` should be a named list with at least
#'   a term called `"distr"` defining the sampling function. The other terms in
#'   the list are used as arguments to the sampling function.
#'
#' @return An object of class `gsaot_indices`.
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
#' dummy_lb <- irrelevance_threshold(y, M)
#'
#' # Custom sampling funtion and network simplex solver
#' dummy_optns <- list(distr = "rgamma", shape = 3)
#' dummy_lb_cust <- irrelevance_threshold(y, M,
#'                                       dummy_optns = dummy_optns,
#'                                       solver = "transport")
irrelevance_threshold <- function(y,
                                  M,
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

  # Check if the cost is defined correctly
  cost_type <- identical(cost, "L2")
  if (!cost_type & !is.function(cost))
    stop("`cost` should be \"L2\" or a function")

  # Check the logical values
  if (!is.logical(discrete_out)) stop("`discrete_out` should be logical")
  if (!is.logical(scaling)) stop("`scaling` should be logical")

  # Check if the solver is present in the pool
  match.arg(solver, c("1d", "wasserstein-bures", "sinkhorn", "sinkhorn_log", "transport"))

  # RETURN THE DUMMY INDICES
  # ----------------------------------------------------------------------------
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
                      "transport" = ot_indices(x,
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

# Check that the bound_optns contain the required components, otherwise initialize them
init_dummy_optns <- function(dummy_optns) {
  if (is.null(dummy_optns)) {
    dummy_optns <- list(distr = "rnorm")
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
    "sinkhorn_stable" = sinkhorn_stable,
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

# Check that the bound_optns contain the required components, otherwise initialize them
init_dummy_optns <- function(dummy_optns) {
  if (is.null(dummy_optns)) {
    dummy_optns <- list(distr = "rnorm")
  }

  return(dummy_optns)
}
