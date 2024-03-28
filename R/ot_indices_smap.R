#' Evaluate one-dimensional Optimal Transport indices on multivariate outputs
#'
#' @param x a data.frame containing the input(s) values
#' @param y a matrix containing the output(s) values
#' @param M a scalar representing the number of partitions for continuous inputs
#'
#' @return An Optimal Transport sensitivity index between 0 and 1 for each of the columns in x
#' @export
#'
#' @examples
#' x <- rnorm(1000)
#' y <- 10 * x
#' ot_indices_1d(data.frame(x), y, 30)
ot_indices_smap <- function(x, y, M) {
  # Input checks
  stopifnot(is.data.frame(x), is.matrix(y))
  stopifnot(dim(x)[1] == dim(y)[1])
  stopifnot(dim(x)[1] > M)

  # Get sizes of input and output
  K <- dim(x)[2]
  P <- dim(y)[2]

  # Build the return structure
  W <- matrix(nrow = P, ncol = K)
  colnames(W) <- colnames(x)
  rownames(W) <- colnames(y)

  # Evaluate an OT index for each univariate output separately
  for (p in seq(P)) {
    ret <- ot_indices_1d(x, y[, p], M)
    W[p, ] <- ret$indices
  }

  return(W)
}
