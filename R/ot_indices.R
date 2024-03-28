#' Calculate Optimal Transport sensitivity indices for multivariate y
#'
#' @description
#' `ot_indices` calculates sensitivity indices using Optimal Transport (OT) for multivariate output data `y` with respect to input data `x`.
#' Sensitivity indices measure the influence of input variables on output variables, with values ranging between 0 and 1.
#'
#' @param x A data.frame containing the input(s) values. The values can be numeric, factors or strings.
#' @param y A matrix containing the output values. Each column represents a different output variable, and each row represents a different observation. Only numeric values are allowed.
#' @param M A scalar representing the number of partitions for continuous inputs.
#' @param solver (optional) Solver for the Optimal Transport problems: `sinkhorn` (default), `sinkhorn_log`, or `wasserstein`.
#' @param solver_optns (optional) Options for the Optimal Transport solver. See details for allowed options.
#' @param scaling (default `TRUE`) Logical that sets whether or not to scale the cost matrix.
#' @param extended_out (default `FALSE`) Logical indicating if the function should return the inner statistics and the partitions.
#'
#' @details
#' The solvers of the OT problem implemented in this package can be divided into two categories: standard and entropic. And then bla, blabla, blablabla
#'
#' @return An object containing:
#' * `indices`: sensitivity indices between 0 and 1 for each column in x, indicating the influence of each input variable on the output variables.
#' * `IS`: values of the inner statistics for the partitions defined by `partitions`. Returned only if `extended_out = TRUE`.
#' * `partitions`: the partitions built to calculate the sensitivity indices. Returned only if `extended_out = TRUE`.
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
#' x <- data.frame(x)
#'
#' M <- 25
#'
#' # Calculate sensitivity indices
#' sensitivity_indices <- ot_indices(x, y, M)
#' sensitivity_indices
#'
#' # With extended output
#' sensitivity_indices_extended <- ot_indices(x, y, M, extended_out = TRUE)
#' sensitivity_indices_extended
#'
#' @export
#'
ot_indices <- function(x,
                       y,
                       M,
                       solver = "sinkhorn",
                       solver_optns = NULL,
                       scaling = TRUE,
                       extended_out = FALSE) {
  # Input checks
  stopifnot(is.data.frame(x), is.numeric(y))
  stopifnot(dim(x)[1] == dim(y)[1])
  stopifnot(dim(x)[1] > M)

  # Remove any NA in output
  y_na <- apply(y, 1, function(row) any(is.na(row)))
  y <- y[!y_na, ]
  col_names <- colnames(x)
  x <- data.frame(x[!y_na, ])
  colnames(x) <- col_names
  if (any(y_na))
    cat("Removed", sum(y_na), "NA(s) in output\n")

  # Build cost matrix
  C <- as.matrix(stats::dist(y, method = "euclidean")) ^ 2
  scaling_param <- 1
  if (scaling) {
    scaling_param <- max(C)
    C <- C / scaling_param
  }

  # Retrieve dimensions of the inputs
  N <- dim(x)[1]
  K <- dim(x)[2]

  # Define the partitions for the estimation
  partitions <- build_partition(x, M)

  # Define the OT solver
  solver_fun <- switch (
    solver,
    "sinkhorn" = sinkhorn,
    "sinkhorn_log" = sinkhorn_log,
    "wasserstein" = transport::transport,
    default = NULL
  )

  # Check the consistency of the options provided
  solver_optns <- check_solver_optns(solver, solver_optns)

  # Build return structure
  W <- array(dim = K)
  names(W) <- colnames(x)
  if (extended_out) IS <- list()

  # Evaluate the upper bound of the indices
  V <- 2 * sum(diag(stats::cov(y)))

  # Define the histogram of the whole output
  a <- rep(1 / N, N)

  for (k in seq_len(K)) {
    # Get the current partition
    partition <- partitions[[k]]

    # Get the current number of partition elements
    M <- length(partition)

    # Build the inner return structure
    Wk <- matrix(nrow = 1, ncol = M)
    n <- matrix(nrow = M)

    for (m in seq(M)) {
      # Retrieve the elements in the partition
      partition_element <- partition[[m]]

      # Build the partition histogram
      n[m] <- length(partition_element)
      b <- rep(1 / n[m], n[m])

      # Call the solver
      ret <- do.call(solver_fun,
                     c(list(a = a, b = b, costm = C[, partition_element]),
                       solver_optns))

      # Save the estimated distance
      Wk[, m] <- ret$cost * scaling_param

      # Quick exit from the loop if NaNs are present
      if (is.na(ret$cost)) m <- M
    }

    W[k] <- ((Wk[1,] %*% n) / (V * N))[1, 1]
    if (extended_out) IS[[k]] <- Wk / V
  }

  out <- gsaot_indices(method = solver, indices = W, bound = V,
                       IS = IS, partitions = partitions,
                       x = x, y = y, extended_out)

  return(out)
}
