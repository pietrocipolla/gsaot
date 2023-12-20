optimal_transport_1d <- function(partition, y_sort, yc) {
  # Set the scaling parameters
  N <- length(y_sort)
  Nc <- length(yc)

  # Expand the empirical CDF quantiles to match with y1 length
  yc <- yc[floor(seq(1/Nc, 1, length.out = N) * Nc + 0.5)]

  # Evaluate the L2 norm of the empirical CDF
  W <- mean((y_sort - yc)^2)

  return(W)
}
