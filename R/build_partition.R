build_partition <- function(x, M) {
  # Get number of inputs
  K <- dim(x)[2]

  # Build the return structure
  partition_indices <- list()

  # Build the partition for each input
  for (k in seq(K)) {
    # If the variable is continuous build the partititon using quantiles
    if (is.numeric(x[, k])) {
      partition_indices[[k]] <- build_continuous_partition(x[, k], M)
      # Otherwise use the unique different elements of the inputs
    } else {
      partition_indices[[k]] <- build_discrete_partition(x[, k])
    }
  }

  return(partition_indices)
}
