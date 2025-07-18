#' Estimate optimal transport indices on one dimensional outputs
#'
#' @inheritParams ot_indices
#' @param y An array containing the output values.
#' @param p A numeric representing the p-norm Lp used as ground cost in the
#'   Optimal Transport problem.
#'
#' @inherit ot_indices return
#' @export
#'
#' @seealso [ot_indices()], [ot_indices_wb()]
#'
#' @examples
#' x <- rnorm(1000)
#' y <- 10 * x
#' ot_indices_1d(data.frame(x), y, 10)
ot_indices_1d <- function(x,
                          y,
                          M,
                          p = 2,
                          boot = FALSE,
                          R = NULL,
                          parallel = "no",
                          ncpus = 1,
                          conf = 0.95,
                          type = "norm") {
  # INPUT CHECKS
  # ----------------------------------------------------------------------------
  # Check if x is a data.frame or a matrix
  if (!(is.data.frame(x) | is.matrix(x)))
    stop("`x` must be a matrix or a data.frame!")
  x <- as.data.frame(x)

  # Check if the output is a numerical
  if (!is.numeric(y) | !is.vector(y))
    stop("`y` must be a vector of numerical values!")

  # Check if the dimensions match
  if (!(nrow(x) == length(y)))
    stop("The number of samples in `x` and `y` should be the same")

  # Check if the number of partitions is lower than the number of samples
  if (nrow(x) <= M)
    stop("The number of partitions should be lower than the number of samples")

  # Check the value of the input p-norm
  if (!is.numeric(p) || p < 1)
    stop("The argument `p` should be a numeric value greater or equal than 1")

  # Check that bootstrapping is correctly set
  if ((!boot & !is.null(R)) | (boot & is.null(R))) {
    stop("Bootstrapping requires boot = TRUE and an integer in R")
  }

  # Check that the bootstrapping type is in the correct set
  match.arg(type, c("norm", "basic", "stud", "perc", "bca"))

  # REMOVE ANY NA IN OUTPUT
  # ----------------------------------------------------------------------------
  y_na <- is.na(y)
  y <- y[!y_na]
  col_names <- colnames(x)
  x <- data.frame(x[!y_na, ])
  colnames(x) <- col_names

  # Build partitions for estimator
  # ----------------------------------------------------------------------------
  partitions <- build_partition(x, M)

  # Get inputs features
  N <- dim(x)[1]
  K <- dim(x)[2]

  # NO BOOT ESTIMATIONS
  # ----------------------------------------------------------------------------
  if (!boot) {
    # Sort y
    y_sort <- sort(y)
    # Evaluate upper bound
    if (p == 2) {
      V <- 2 * stats::var(y)
    } else {
      V <- as.matrix(stats::dist(y, method = "minkowski", p = p)) ** p
      V <- higher_bound(V)
    }
  }

  # BUILD THE RETURN STRUCTURE
  # ----------------------------------------------------------------------------
  W <- array(dim = K)
  names(W) <- colnames(x)
  IS <- list()
  if (boot) {
    V <- array(dim = K)
    W_ci <- data.frame(matrix(nrow = K,
                              ncol = 5,
                              dimnames = list(NULL,
                                              c("input", "original", "bias",
                                                "low.ci", "high.ci"))))
    W_ci$input <- names(W)
    IS_ci <- list()
    V_ci <- list()

    W_boot <- list()
  }

  # ESTIMATE THE INDICES FOR EACH PARTITION
  # ----------------------------------------------------------------------------
  for (k in seq(K)) {
    # Get the partitions for the current input
    partition <- partitions[, k]

    # Set the number of partition elements
    M <- max(partition, na.rm = T)

    # Initialize the return structure
    Wk <- matrix(nrow = 1, ncol = M)
    n <- matrix(nrow = M)

    # NO BOOTSTRAP ESTIMATION
    if (!boot) {
      for (m in seq(M)) {
        partition_element <- which(partition == m)
        yc <- sort(y[partition_element])

        # Compute the OT distance in each partition
        Wk[1, m] <- optimal_transport_1d(partition_element, y_sort, yc, p)
        n[m] <- length(partition_element)
      }

      # Save the results
      W[k] <- ((Wk %*% n) / (V * N))[1, 1]
      IS[[k]] <- Wk / V
    } else {
    # BOOTSTRAP ESTIMATION
      dat <- cbind(partition, y)

      # Do boostrap estimation stratified by the partitions
      W_boot[[k]] <- boot::boot(data = dat,
                                statistic = ot_1d_boot,
                                R = R,
                                strata = partition,
                                p = p,
                                parallel = parallel,
                                ncpus = ncpus)

      # Transform the results into readable quantities
      W_stats <- bootstats(W_boot[[k]], type = type, conf = conf)

      # Save the results
      W[k] <- W_stats$index[1]
      IS[[k]] <- matrix(W_stats$index[2:(M + 1)], nrow = 1)
      V[k] <- W_stats$index[M + 2]
      W_ci[k, 2:5] <- c(W_stats$original[1], W_stats$bias[1],
                        W_stats$low.ci[1], W_stats$high.ci[1])
      IS_ci[[k]] <- cbind(W_stats$low.ci[2:(M + 1)], W_stats$high.ci[2:(M + 1)])
      V_ci[[k]] <- cbind(W_stats$low.ci[M + 2], W_stats$high.ci[M + 2])
    }
  }

  if (boot) {
    out <- gsaot_indices(method = "1-dimensional",
                         indices = W,
                         bound = V,
                         IS = IS,
                         partitions = partitions,
                         x = x, y = y,
                         indices_ci = W_ci,
                         bound_ci = V_ci,
                         IS_ci = IS_ci,
                         R = R,
                         conf = conf,
                         type = type,
                         W_boot = W_boot)

    return(out)
  }

  out <- gsaot_indices(method = paste0("1-dimensional (p=", p, ")"),
                       indices = W,
                       bound = V,
                       IS = IS,
                       partitions = partitions,
                       x = x, y = y)

  return(out)
}

# Function generating bootstrap statistics
ot_1d_boot <- function(d,
                       indices,
                       p) {
  # Retrieve the partitions
  partition <- d[indices, 1]

  # Retrieve the output
  y <- d[indices, 2:ncol(d)]

  # Get the number of realizations
  N <- length(partition)

  # Sort y
  y_sort <- sort(y)

  # Evaluate upper bound
  if (p == 2) {
    V <- 2 * stats::var(y)
  } else {
    V <- as.matrix(stats::dist(y, method = "minkowski", p = p)) ** p
    V <- higher_bound(V)
  }

  # Get the number of partitions
  M <- max(partition, na.rm = T)

  # Initialize the return structure
  Wk <- matrix(nrow = 1, ncol = M)
  n <- matrix(nrow = M)

  # Estimate the inner statistic
  for (m in seq(M)) {
    partition_element <- which(partition == m)
    yc <- sort(y[partition_element])

    # Compute the OT distance in each partition
    Wk[1, m] <- optimal_transport_1d(partition_element, y_sort, yc, p)
    n[m] <- length(partition_element)
  }

  # Calculate the 1d index
  W <- ((Wk %*% n) / (V * N))[1, 1]

  return(c(W, Wk, V))
}

# Function for computing OT using L2 cost between two numerical vectors
optimal_transport_1d <- function(partition, y_sort, yc, p) {
  # Set the scaling parameters
  N <- length(y_sort)
  Nc <- length(yc)

  # Expand the empirical CDF quantiles to match with y length
  yc <- yc[floor(seq(1 / Nc, 1, length.out = N) * Nc + 0.5)]

  # Evaluate the L2 norm of the empirical CDF
  W <- mean(abs(y_sort - yc)^p)

  return(W)
}
