ot_indices_entropic <- function(x, y, M, eps, d = 2, num_iterations = 10000, cost = "euclidean") {
  # Input checks
  stopifnot(is.data.frame(x), is.data.frame(y))
  stopifnot(dim(x)[1] == dim(y)[1])
  stopifnot(dim(x)[1] > M)

  # Build cost matrix
  C <- as.matrix(stats::dist(y, method = cost))^d
  if (eps < 0) {
    eps <- - eps * max(C)
  }

  # Retrieve values useful for the algorithm
  N <- dim(x)[1]
  K <- dim(x)[2]

  # Define the partitions for the estimation
  partitions <- build_partition(x, M)

  # Build return structure
  W <- array(dim = K)

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
      partition_element <- partition[[m]]
      ret <- optimal_transport_sinkhorn(C[partition_element, ], num_iterations, eps)
      Wk[, m] <- ret$W22
      n[m] <- length(partition_element)
    }

    W[k] <- ((Wk[1, ] %*% n) / (V * N))[1, 1]
  }

  return(W)
}
