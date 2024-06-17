sqrtm <- function(A) {
  if (ncol(A) == 1) return(sqrt(A))

  eig <- eigen(A, symmetric = TRUE)

  eig$values <- ifelse(eig$values < 0, 0, eig$values)

  RA <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)

  return(RA)
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
