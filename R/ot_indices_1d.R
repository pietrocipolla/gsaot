#' Evaluate optimal transport indices on one dimensional vectors
#'
#' @param x a data.frame containing the input(s) values
#' @param y an array containing the output values
#' @param M a scalar representing the number of partitions for continuous inputs
#' @param ext_out logical indicating if the function should return the inner statistics and the partitions
#'
#' @return A sensitivity index between 0 and 1 for each of the columns in x
#' @export
#'
#' @examples
#' x <- rnorm(1000)
#' y <- 10*x
#' ot_indices_1d(data.frame(x), y, 30)

ot_indices_1d <- function(x, y, M, ext_out = FALSE) {
  # Input checks
  stopifnot(is.data.frame(x))
  stopifnot(dim(x)[1] == length(y))
  stopifnot(dim(x)[1] > M)

  # Build partitions for estimator
  partitions <- build_partition(x, M)

  # Sort y
  y_sort <- sort(y)

  # Get inputs features
  N <- dim(x)[1]
  K <- dim(x)[2]

  # Evaluate upper bound
  V <- 2 * stats::var(y)

  # Build the return structure
  W <- array(dim = K)
  if (ext_out) IS <- list()

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

      Wk[1, m] <- optimal_transport_1d(partition_element, y_sort, yc)
      n[m] <- length(partition_element)
    }

    W[k] <- ((Wk %*% n) / (V * N))[1, 1]
    if (ext_out) IS[[k]] <- Wk
  }

  if (ext_out) {
    return(list(W = W, IS = IS, partitions = partitions))
  } else {
    return(list(W = W))
  }
}
