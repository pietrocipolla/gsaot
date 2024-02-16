#' Evaluate Wasserstein-Bures approximation of the Optimal Transport solution
#'
#' @param x a data.frame containing the input(s) values
#' @param y a data.frame containing the outputs values
#' @param M a scalar representing the number of partitions for continuous inputs
#' @param ext_out logical indicating if the function should return the inner statistics and the partitions
#'
#' @return A sensitivity index between 0 and 1 for each of the columns in x
#' @export
#'
#' @examples
#' N <- 10000
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
#' ot_indices_bw(x, y, 100)
ot_indices_bw <- function(x, y, M, ext_out = FALSE) {
  # Input checks
  stopifnot(is.data.frame(x), is.numeric(y))
  stopifnot(dim(x)[1] == dim(y)[1])
  stopifnot(dim(x)[1] > M)

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
  Adv <- array(dim = K)
  Diff <- array(dim = K)

  if (ext_out) IS <- list()

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

    if (ext_out) IS[[k]] <- Wk
  }

  if (ext_out) {
    return(list(W = W, Adv = Adv, Diff = Diff, IS = IS, partitions = partitions))
  } else {
    return(list(W = W, Adv = Adv, Diff = Diff))
  }
}
