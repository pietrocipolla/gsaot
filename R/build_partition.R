build_partition <- function(x, M) {
  # Get number of inputs
  N <- dim(x)[1]
  K <- dim(x)[2]

  # Build the return structure
  partition_indices <- matrix(nrow = N, ncol = K)

  # Build the partition for each input
  for (k in seq(K)) {
    # If the variable is continuous build the partititon using quantiles
    if (is.double(x[, k])) {
      partition_indices[, k] <- build_continuous_partition(x[, k], M)
    # Otherwise use the unique different elements of the inputs
    } else {
      partition_indices[, k] <- build_discrete_partition(x[, k])
    }

    if (any(is.na(x[, k])))
      warning("Removed", sum(is.na(x[, k])), "NA(s) in column", k, "\n")
  }

  return(partition_indices)
}

build_continuous_partition <- function(x, M) {
  # Get the number of realizations
  N <- length(x)

  # Get the number of NaNs
  N_nan <- sum(is.na(x))

  # Retrieve x ranks and the transformation from the sorted x to the original
  # ones
  ord <- rank(x, na.last = FALSE, ties.method = "first")

  # Build the partitions. Each partition has ~ the same number of elements.
  partitions_indices <- floor(seq(from = sum(is.na(x)) + 1,
                                  to = N + 1,
                                  length.out = M + 1))

  # Build the return structure
  partitions <- c(rep(NA, times = sum(is.na(x))),
    rep(seq(M), times = diff(partitions_indices)))[ord]

  return(partitions)
}

build_discrete_partition <- function(x) {
  # Find the unique elements of the input
  if (is.integer(x)) {
    x_unique <- sort(unique(x[!is.na(x)]))
  } else if (is.factor(x)) {
    x_unique <- levels(x[!is.na(x)])
  } else {
    x_unique <- unique(x[!is.na(x)])
  }

  # Set the number of partitions
  M <- length(x_unique)

  # Define the return structure
  partitions <- rep(NA, times = N)

  # Find the indices assigned to each partition
  for (m in seq(M)) {
    partitions[x == x_unique[m]] <- m

    if (sum(partitions == m, na.rm = TRUE) == 1) {
      warning("Partition ", m, " has only 1 element")
    }
  }

  return(partitions)
}
