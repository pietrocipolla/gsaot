optimal_trasport_bw <- function(partition, y, my, Cy, traceCy, Ry) {
  # Get the current partition indices for y
  yc <- y[partition, ]

  # Get the conditioned statistics
  mc <- colMeans(yc)
  Cc <- stats::cov(yc)

  Adv <- sum((my - mc)^2)
  Diff <- traceCy + sum(diag(Cc)) - 2 * tracesqrtm(Ry %*% Cc %*% Ry)

  W <- Adv + Diff

  return(cbind(W, Adv, Diff))
}
