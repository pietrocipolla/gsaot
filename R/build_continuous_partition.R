build_continuous_partition <- function(x, M) {
  # Retrieve x ranks
  ord <- rank(x, ties.method = "first")

  # Get the number of realizations
  N <- length(x)

  # Build the partitions. Each partition has ~ the same number of elements
  partitions_indices <- floor(seq(from = 1, to = N + 1, length.out = M + 1))

  # Build the return structure
  partitions <- list()

  # Find the indices assigned to each partition
  for (m in seq(M)) {
    partitions[[m]] <- partitions_indices[m]:(partitions_indices[m+1] - 1)
    partitions[[m]] <- which(ord %in% partitions[[m]])

    if (length(partitions[[m]]) == 1)
      warning("Partition ", m, " has only 1 element")
  }

  return(partitions)
}
