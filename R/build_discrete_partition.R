build_discrete_partition <- function(x) {
  # Find the unique elements of the input
  if (is.integer(x)) {
    x_unique <- sort(unique(x))
  } else if (is.factor(x)) {
    x_unique <- levels(x)
  } else {
    x_unique <- unique(x)
  }

  # Set the number of partitions
  M <- length(x_unique)

  # Build the return structure
  partitions <- list()

  # Find the indices assigned to each partition
  for (m in seq(M)) {
    partitions[[m]] <- which(x == x_unique[m])
  }

  return(partitions)
}
