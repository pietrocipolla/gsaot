#' Evaluate Wasserstein-Bures approximation of the Optimal Transport solution
#'
#' @inheritParams ot_indices
#'
#' @inherit ot_indices return
#'
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
#' y <- y
#'
#' ot_indices_wb(x, y, 100)
ot_indices_wb <- function(x, y, M, extended_out = FALSE) {
  # Input checks
  stopifnot(is.data.frame(x), is.numeric(y))
  stopifnot(dim(x)[1] == dim(y)[1])
  stopifnot(dim(x)[1] > M)

  # Remove any NA in output
  y_na <- apply(y, 1, function(row) any(is.na(row)))
  y <- y[!y_na, ]
  x <- data.frame(x[!y_na, ])
  if (any(y_na))
    cat("Removed", sum(y_na), "NA(s) in output\n")

  # Compute the statistics for the unconditioned distribution
  my <- colMeans(y)
  Cy <- stats::cov(y)
  traceCy <- sum(diag(Cy))
  Ry <- sqrtm(Cy)

  # Retrieve values useful for the algorithm
  N <- dim(x)[1]
  K <- dim(x)[2]

  # Define the partitions for the estimation
  partitions <- build_partition(x, M)

  # Initialize the result matrices
  W <- array(dim = K)
  names(W) <- colnames(x)
  Adv <- array(dim = K)
  names(Adv) <- colnames(x)
  Diff <- array(dim = K)
  names(Diff) <- colnames(x)

  if (extended_out) IS <- list()

  # Evaluate the upper bound of the indices
  V <- 2 * traceCy

  # Run the algorithm for each variable
  for (k in seq_len(K)) {
    # Get the current partition
    partition <- partitions[[k]]

    # Get the current number of partition elements
    M <- length(partition)

    # Build the inner return structure
    Wk <- matrix(nrow = 3, ncol = M)
    n <- matrix(nrow = M)

    for (m in seq(M)) {
      partition_element <- partition[[m]]
      Wk[, m] <- optimal_trasport_bw(partition_element, y, my, Cy, traceCy, Ry)
      n[m] <- length(partition_element)
    }

    W[k] <- ((Wk[1, ] %*% n) / (V * N))[1, 1]
    Adv[k] <- ((Wk[2, ] %*% n) / (V * N))[1, 1]
    Diff[k] <- ((Wk[3, ] %*% n) / (V * N))[1, 1]

    if (extended_out) IS[[k]] <- Wk / V
  }

  out <- gsaot_indices(method = "wasserstein-bures", indices = W, bound = V,
                       IS = IS, partitions = partitions,
                       x = x, y = y, extended_out,
                       Adv = Adv, Diff = Diff)

  return(out)
}
