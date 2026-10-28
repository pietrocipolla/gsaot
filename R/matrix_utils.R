sqrtm <- function(A) {
  if (ncol(A) == 1) return(sqrt(A))

  eig <- eigen(A, symmetric = TRUE)

  eig$values <- ifelse(eig$values < 0, 0, eig$values)

  RA <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)

  return(RA)
}

# Compute squared Euclidean distances without as.matrix(dist()), whose integer
# indexing overflows once the dense matrix has more than 2^31 - 1 entries.
squared_euclidean_cost <- function(y) {
  y <- as.matrix(y)
  squared_norms <- rowSums(y * y)
  cost <- tcrossprod(y)

  for (j in seq_len(ncol(cost))) {
    cost[, j] <- pmax(0, squared_norms + squared_norms[j] - 2 * cost[, j])
  }

  diag(cost) <- 0

  return(cost)
}

tracesqrtm <- function(A) {
  if (ncol(A) == 1) return(sqrt(A))

  eig <- eigen(A)
  eig$values <- ifelse(eig$values < 0, 0, eig$values)

  return(sum(sqrt(eig$values)))
}

# Function to compute the higher bound of the indices directly from the cost
# matrix
higher_bound <- function(C) {
  # Compute the bound using the U-statistic estimator
  bound <- mean(C[lower.tri(C)])

  return(bound)
}
