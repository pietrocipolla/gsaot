#' Evaluate sensitivity maps using Optimal Transport indices
#'
#' @inheritParams ot_indices_1d
#' @param y A matrix containing the output values. Each column is interpreted as
#'   a different output.
#'
#' @return A matrix where each column represents an input and each row
#'   represents an output. The values are indices between 0 and 1 computed using
#'   [ot_indices_1d()].
#' @export
#'
#' @examples
#' N <- 1000
#'
#' x1 <- rnorm(N)
#' x2 <- rnorm(N)
#' x <- cbind(x1, x2)
#'
#' y1 <- 10 * x1
#' y2 <- x1 + x2
#' y <- cbind(y1, y2)
#'
#' ot_indices_smap(data.frame(x), y, 30)
ot_indices_smap <- function(x, y, M) {
  # Check if the output is a numerical
  if (!is.numeric(y) | !is.matrix(y)) stop("`y` must be a matrix of numerical values!")

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
