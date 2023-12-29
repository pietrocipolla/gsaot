sqrtm <- function(A) {
  eig <- eigen(A, symmetric = TRUE)

  eig$values <- ifelse(eig$values < 0, 0, eig$values)

  RA <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)

  return(RA)
}
