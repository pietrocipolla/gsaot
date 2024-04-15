#' Calculate Optimal Transport sensitivity indices for multivariate y
#'
#' @description `ot_indices` calculates sensitivity indices using Optimal
#'   Transport (OT) for multivariate output data `y` with respect to input data
#'   `x`. Sensitivity indices measure the influence of input variables on output
#'   variables, with values ranging between 0 and 1.
#'
#' @param x A data.frame containing the input(s) values. The values can be
#'   numeric, factors or strings.
#' @param y A matrix containing the output values. Each column represents a
#'   different output variable, and each row represents a different observation.
#'   Only numeric values are allowed.
#' @param M A scalar representing the number of partitions for continuous
#'   inputs.
#' @param discrete_out (default `FALSE`) Logical, by default the output sample
#'   in `y` are equally weighted. If `discrete_out=TRUE`, the function tries to
#'   create an histogram of the realizations and to use the histogram as
#'   weigths. This works if the output is assumed to be discrete or mixed and
#'   the number of realizations is high. The main advantage of this option is to
#'   reduce the dimension of the cost matrix.
#' @param solver (optional) Solver for the Optimal Transport problems:
#'   `sinkhorn` (default), `sinkhorn_log`, or `wasserstein`.
#' @param solver_optns (optional) Options for the Optimal Transport solver. See
#'   details for allowed options.
#' @param scaling (default `TRUE`) Logical that sets whether or not to scale the
#'   cost matrix.
#'
#' @details The solvers of the OT problem implemented in this package can be
#'   divided into two categories: standard and entropic. And then bla, blabla,
#'   blablabla
#'
#' @return An object containing:
#' * `indices`: sensitivity indices between 0 and 1 for each column in x, indicating the influence of each input variable on the output variables.
#' * `IS`: values of the inner statistics for the partitions defined by `partitions`.
#' * `partitions`: the partitions built to calculate the sensitivity indices.
#'
#' @seealso [ot_indices_1d()], [ot_indices_wb()]
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
#' @export
#'
ot_indices <- function(x,
                       y,
                       M,
                       discrete_out = FALSE,
                       solver = "sinkhorn",
                       solver_optns = NULL,
                       scaling = TRUE) {
  # Check if x is a data.frame or a matrix
  if (!(is.data.frame(x) | is.matrix(x))) stop("`x` must be a matrix or a data.frame!")
  x <- as.data.frame(x)

  # Check if the output is a numerical
  if (!is.numeric(y) | !is.matrix(y)) stop("`y` must be a matrix of numerical values!")

  # Check if the dimensions match
  if (!(nrow(x) == nrow(y))) stop("The number of samples in `x` and `y` should be the same")

  # Check if the number of partitions is lower than the number of samples
  if (nrow(x) <= M) stop("The number of partitions should be lower than the number of samples")

  # Check the logical values
  if (!is.logical(discrete_out)) stop("`discrete_out` should be logical")
  if (!is.logical(scaling)) stop("`scaling` should be logical")

  # Check if the solver is present in the pool
  match.arg(solver, c("sinkhorn", "sinkhorn_log", "wasserstein"))

  # Remove any NA in output
  y_na <- apply(y, 1, function(row) any(is.na(row)))
  y <- y[!y_na, ]
  col_names <- colnames(x)
  x <- data.frame(x[!y_na, ])
  colnames(x) <- col_names
  if (any(y_na))
    cat("Removed", sum(y_na), "NA(s) in output\n")

  # Retrieve dimensions of the inputs
  N <- dim(x)[1]
  K <- dim(x)[2]

  # Build an histogram if requested, otherwise assume equal weights
  a <- rep(1 / N, N)
  if (discrete_out) {
    y_unique <- unique(y)

    a <- apply(y_unique, MARGIN = 1, FUN = function(u) sum(
      apply(y, 1, FUN = function(row) all(row == u))
    )) / N
  }

  if (identical(a, rep(1 / N, N)) & discrete_out)
    warning("The output is continuous, consider using `discrete_out=FALSE`")

  # Build cost matrix
  if (discrete_out) {
    C <- as.matrix(stats::dist(y_unique, method = "euclidean")) ^ 2
  } else
    C <- as.matrix(stats::dist(y, method = "euclidean")) ^ 2

  scaling_param <- 1
  if (scaling) {
    scaling_param <- max(C)
    C <- C / scaling_param
  }

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
  IS <- list()

  # Evaluate the upper bound of the indices
  V <- 2 * sum(diag(stats::cov(y)))

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
      # The idea is to use the marginal histogram b both as histogram and as selector for the cost matrix parts. Indeed, as.logical transforms everything > 0 to TRUE and everything else t FALSE
      n[m] <- length(partition_element)
      if (discrete_out) {
        b <- apply(y_unique, MARGIN = 1, FUN = function(u) sum(
          apply(y[partition_element, ], 1, FUN = function(row) all(row == u))
        )) / n[m]
      } else {
        b <- array(0, dim = N)
        b[partition_element] <- 1 / n[m]
      }

      # Call the solver
      ret <- do.call(solver_fun,
                     c(list(a = a, b = b[b > 0], costm = C[, (b > 0)]),
                       solver_optns))

      # Save the estimated distance
      Wk[, m] <- ret$cost * scaling_param

      # Quick exit from the loop if NaNs are present
      if (is.na(ret$cost)) m <- M
    }

    W[k] <- ((Wk[1,] %*% n) / (V * N))[1, 1]
    IS[[k]] <- Wk / V
  }

  out <- gsaot_indices(method = solver,
                       indices = W,
                       bound = V,
                       IS = IS,
                       partitions = partitions,
                       x = x, y = y)

  return(out)
}
