build_partition <- function(x, M) {
  # Get number of inputs
  n <- dim(x)[1]
  K <- dim(x)[2]

  # Build the return structure
  partition_indices <- list()

  # Build the partition for each input
  for (k in seq(K)) {
    # If the variable is continuous build the partititon using quantiles
    if (is.double(x[, k])) {
      partition_indices[[k]] <- build_continuous_partition(x[!is.na(x[, k]), k], M)
    # Otherwise use the unique different elements of the inputs
    } else {
      partition_indices[[k]] <- build_discrete_partition(x[!is.na(x[, k]), k])
    }

    if (any(is.na(x[, k])))
      cat("Removed", sum(is.na(x[, k])), "NA(s) in column", k, "\n")
  }

  return(partition_indices)
}

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
    partitions[[m]] <- partitions_indices[m]:(partitions_indices[m + 1] - 1)
    partitions[[m]] <- which(ord %in% partitions[[m]])
  }

  return(partitions)
}

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

    if (length(partitions[[m]]) == 1) {
      warning("Partition ", m, " has only 1 element")
    }
  }

  return(partitions)
}
