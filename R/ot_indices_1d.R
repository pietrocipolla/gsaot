#' Evaluate Optimal Transport indices on one dimensional outputs
#'
#' @inheritParams ot_indices
#' @param y An array containing the output values.
#' @param grid The spacing of the grid over which the integral of \eqn{\lVert F_Y - F_{Y|X_i \in \mathcal{X}_i} \rVert}.
#'
#' @return An Optimal Transport sensitivity index between 0 and 1 for each of
#'   the columns in x
#' @export
#'
#' @seealso [ot_indices()], [ot_indices_wb()]
#'
#' @examples
#' x <- rnorm(1000)
#' y <- 10 * x
#' ot_indices_1d(data.frame(x), y, 30)
ot_indices_1d <- function(x, y, M, grid = 1e-3) {
  # Input checks
  stopifnot(is.data.frame(x))
  stopifnot(dim(x)[1] == length(y))
  stopifnot(dim(x)[1] > M)

  # Remove any NA in output
  y_na <- is.na(y)
  y <- y[!y_na]
  x <- as.data.frame(x[!y_na, ])
  if (any(y_na))
    cat("Removed", sum(y_na), "NA(s) in output\n")

  # Build partitions for estimator
  partitions <- build_partition(x, M)

  # Sort y
  y_ecdf <- stats::ecdf(y)

  # Get inputs features
  N <- dim(x)[1]
  K <- dim(x)[2]

  # Evaluate upper bound
  V <- 2 * stats::var(y)

  # Build the return structure
  W <- array(dim = K)
  names(W) <- colnames(x)
  IS <- list()

  for (k in seq(K)) {
    # Get the partitions for the current input
    partition <- partitions[[k]]

    # Set the number of partition elements
    M <- length(partition)

    # Initialize the return structure
    Wk <- matrix(nrow = 1, ncol = M)
    n <- matrix(nrow = M)

    for (m in seq(M)) {
      partition_element <- partition[[m]]
      yc <- sort(y[partition_element])

      # Compute the OT distance in each partition
      Wk[1, m] <- optimal_transport_1d(partition_element, y_ecdf, yc, grid)
      n[m] <- length(partition_element)
    }

    W[k] <- ((Wk %*% n) / (V * N))[1, 1]
    IS[[k]] <- Wk / V
  }

  out <- gsaot_indices(method = "1-dimensional",
                       indices = W,
                       bound = V,
                       IS = IS,
                       partitions = partitions,
                       x = x, y = y)

  return(out)
}

optimal_transport_1d <- function(partition, y_ecdf, yc, grid) {
  # Compute the empirical CDF for the conditioned distribution
  yc_ecdf <- stats::ecdf(yc)

  # Fix the grid and compute quantiles
  probs <- seq(0, 1, grid)
  qy <- stats::quantile(y_ecdf, probs)
  qyc <- stats::quantile(yc_ecdf, probs)

  # Evaluate the L2 norm of the empirical CDF
  W <- mean((qy - qyc)^2)

  return(W)
}
