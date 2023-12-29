tracesqrtm <- function(A) {
  eig <- eigen(A)
  eig$values <- ifelse(eig$values < 0, 0, eig$values)
  
  return(sum(sqrt(eig$values)))
}
