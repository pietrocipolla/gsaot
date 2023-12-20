#' Evaluate optimal transport indices on one dimensional vectors
#'
#' @param x an array containing the input values
#' @param y an array containing the output values
#' @param M a scalar representing the number of partitions
#'
#' @return A sensitivity index between 0 and 1
#' @export
#'
#' @examples
#' x <- rnorm(1000)
#' y <- 10*x
#' ot_indices_1d(x, y, 30)
ot_indices_1d <- function(x, y, M) {
  # Input checks
  stopifnot(is.double(x), is.double(y))
  stopifnot(length(x) == length(y))
  stopifnot(length(x) > M)

  # Retrieve x ranks
  ord <- rank(x)

  # Sort y
  y_sort <- sort(y)

  # Get the number of realizations
  N <- length(x)

  # Build the partitions. Each partition has ~ the same number of elements
  partitions_indices <- floor(seq(from = 1, to = N + 1, length.out = M + 1))

  # Get the number of elements in each partition
  n <- diff(partitions_indices)

  # Initialize the return structure
  W <- matrix(nrow = 1, ncol = M)

  for (m in seq(M)) {
    partition <- partitions_indices[m]:(partitions_indices[m+1] - 1)
    partition <- which(ord %in% partition)
    yc <- sort(y[partition])

    W[1, m] <- optimal_transport_1d(partition, y_sort, yc)
  }

  return(((W %*% n) / (2 * stats::var(y) * N))[1, 1])
}
