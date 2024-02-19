#' Title
#'
#' @param x a data.frame containing the input(s) values
#' @param y a matrix containing the outputs values
#' @param M a scalar representing the number of partitions for continuous inputs
#' @param solver maximum number of iterations of the Sinkhorn algorithm
#' @param solver_optns algorithm to solve the Optimal Transport problem
#' @param scaling Lore ipsum
#' @param extended_out a double representing the coefficient of the entropic regularization. If negative, the cost matrix is normalized.
#'
#' @return A sensitivity index between 0 and 1 for each column in x
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
#' x <- data.frame(x)
#'
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
  }

  return(W)
}
