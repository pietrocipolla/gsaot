#' Evaluate Wasserstein-Bures approximation of the Optimal Transport solution
#'
#' @inheritParams ot_indices
#'
#' @inherit ot_indices return
#'
#' @export
#'
#' @seealso [ot_indices()], [ot_indices_1d()]
#'
#' @examples
#' N <- 1000
#'
#' mx <- c(1, 1, 1)
#' Sigmax <- matrix(data = c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), nrow = 3)
#'
#' x1 <- rnorm(N)
#' x2 <- rnorm(N)
#' x3 <- rnorm(N)
#'
#' x <- cbind(x1, x2, x3)
#' x <- mx + x %*% chol(Sigmax)
#'
#' A <- matrix(data = c(4, -2, 1, 2, 5, -1), nrow = 2, byrow = TRUE)
#' y <- t(A %*% t(x))
#'
#' x <- data.frame(x)
#' y <- y
#'
#' ot_indices_wb(x, y, 100)
ot_indices_wb <- function(x,
                          y,
                          M,
                          boot = FALSE,
                          R = NULL,
                          parallel = "no",
                          ncpus = 1,
                          conf = 0.95,
                          type = "norm") {
  # INPUT CHECKS
  # ----------------------------------------------------------------------------
  # Check if x is a data.frame or a matrix
  if (!(is.data.frame(x) | is.matrix(x))) stop("`x` must be a matrix or a data.frame!")
  x <- as.data.frame(x)

  # Check if the output is a numerical
  if (!is.numeric(y) & !is.matrix(y)) stop("`y` must be a matrix or vector of numerical values!")

  # Conversion to matrix in case
  if (!is.matrix(y)) y <- matrix(y, ncol = 1)

  # Check if the dimensions match
  if (!(nrow(x) == nrow(y))) stop("The number of samples in `x` and `y` should be the same")

  # Check if the number of partitions is lower than the number of samples
  if (nrow(x) <= M) stop("The number of partitions should be lower than the number of samples")

  # Check that bootstrapping is correctly set
  if ((!boot & !is.null(R)) | (boot & is.null(R))) {
    stop("Bootstrapping requires boot = TRUE and an integer in R")
  }

  # Check that the bootstrapping type is in the correct set
  match.arg(type, c("norm", "basic", "stud", "perc", "bca"))

  # REMOVE ANY NA IN OUTPUT
  # ----------------------------------------------------------------------------
  y_na <- apply(y, 1, function(row) any(is.na(row)))
  y <- as.matrix(y[!y_na, ])
  col_names <- colnames(x)
  x <- data.frame(x[!y_na, ])
  colnames(x) <- col_names

  # Compute the statistics for the unconditioned distribution if no boot
  # ----------------------------------------------------------------------------
  if (!boot) {
    my <- colMeans(y)
    Cy <- stats::cov(y)
    traceCy <- sum(diag(Cy))
    Ry <- sqrtm(Cy)

    # Evaluate the upper bound of the indices
    V <- 2 * traceCy
  }

  # Retrieve values useful for the algorithm
  N <- dim(x)[1]
  K <- dim(x)[2]

  # Define the partitions for the estimation
  # ----------------------------------------------------------------------------
  partitions <- build_partition(x, M)

  # Initialize the result matrices
  # ----------------------------------------------------------------------------
  W <- array(dim = K)
  names(W) <- colnames(x)
  Adv <- array(dim = K)
  names(Adv) <- colnames(x)
  Diff <- array(dim = K)
  names(Diff) <- colnames(x)

  IS <- list()

  if (boot) {
    V <- array(dim = K)
    W_ci <- data.frame(matrix(nrow = K * 3,
                              ncol = 6,
                              dimnames = list(NULL,
                                              c("input", "component",
                                                "original", "bias",
                                                "low.ci", "high.ci"))))
    W_ci$input <- rep(names(W), times = 3)
    W_ci$component <- rep(c("wass-bures", "advective", "diffusive"), each = K)

    IS_ci <- list()
    V_ci <- list()
  }

  # ESTIMATE THE INDICES FOR EACH PARTITION
  # ----------------------------------------------------------------------------
  for (k in seq(K)) {
    # Get the partitions for the current input
    partition <- partitions[, k]

    # Set the number of partition elements
    M <- max(partition, na.rm = T)

    # NO BOOTSTRAP ESTIMATION
    if (!boot) {
      # Initialize the return structure
      Wk <- matrix(nrow = 3, ncol = M)
      n <- matrix(nrow = M)

      for (m in seq(M)) {
        partition_element <- which(partition == m)
        Wk[, m] <- optimal_trasport_bw(partition_element, y, my, Cy, traceCy, Ry)
        n[m] <- length(partition_element)
      }

      W[k] <- ((Wk[1, ] %*% n) / (V * N))[1, 1]
      Adv[k] <- ((Wk[2, ] %*% n) / (V * N))[1, 1]
      Diff[k] <- ((Wk[3, ] %*% n) / (V * N))[1, 1]

      IS[[k]] <- Wk / V
    } else {
      # BOOTSTRAP ESTIMATION
      dat <- cbind(partition, y)

      # Do boostrap estimation stratified by the partitions
      W_boot <- boot::boot(data = dat,
                           statistic = ot_wb_boot,
                           R = R,
                           strata = partition,
                           parallel = parallel,
                           ncpus = ncpus)

      # Transform the results into readable quantities
      W_stats <- bootstats(W_boot, type = type, conf = conf)

      # Save the results
      # It is a bit of a mess in this case, but it is all because everything is
      # replicated three times
      # As in the base case
      W[k] <- W_stats$index[1]
      Adv[k] <- W_stats$index[2]
      Diff[k] <- W_stats$index[3]
      IS[[k]] <- rbind(W_stats$index[4:(M + 3)],
                       W_stats$index[(M + 4):(2 * M + 3)],
                       W_stats$index[(2 * M + 4):(3 * M + 3)])

      # Bootstrap estimate of the variance
      V[k] <- W_stats$index[3 * M + 4]

      # Boostrap estimates of the indices
      W_ci[k, 3:6] <- c(W_stats$original[1], W_stats$bias[1],
                        W_stats$low.ci[1], W_stats$high.ci[1])
      W_ci[K + k, 3:6] <- c(W_stats$original[2], W_stats$bias[2],
                            W_stats$low.ci[2], W_stats$high.ci[2])
      W_ci[2 * K + k, 3:6] <- c(W_stats$original[3], W_stats$bias[3],
                                W_stats$low.ci[3], W_stats$high.ci[3])

      # Bootstrap estimates of the inner statistics
      IS_ci[[k]] <- rbind(cbind(W_stats$low.ci[4:(M + 3)],
                                W_stats$high.ci[4:(M + 3)]),
                          cbind(W_stats$low.ci[(M + 4):(2 * M + 3)],
                                W_stats$high.ci[(M + 4):(2 * M + 3)]),
                          cbind(W_stats$low.ci[(2 * M + 4):(3 * M + 3)],
                                W_stats$high.ci[(2 * M + 4):(3 * M + 3)]))

      # Bootstrap estimates for the variance
      V_ci[[k]] <- cbind(W_stats$low.ci[3 * M + 4], W_stats$high.ci[3 * M + 4])
    }
  }

  if (boot) {
    out <- gsaot_indices(method = "wass-bures",
                         indices = W,
                         bound = V,
                         IS = IS,
                         partitions = partitions,
                         x = x, y = y,
                         Adv = Adv,
                         Diff = Diff,
                         indices_ci = W_ci,
                         bound_ci = V_ci,
                         IS_ci = IS_ci,
                         R = R,
                         conf = conf,
                         type = type)

    return(out)
  }

  out <- gsaot_indices(method = "wass-bures",
                       indices = W,
                       bound = V,
                       IS = IS,
                       partitions = partitions,
                       x = x, y = y,
                       Adv = Adv,
                       Diff = Diff)

  return(out)
}

# Function generating bootstrap statistics
ot_wb_boot <- function(d,
                       indices) {
  # Retrieve the partitions
  partition <- d[indices, 1]

  # Retrieve the output
  y <- matrix(d[indices, 2:ncol(d)], ncol = ncol(d) - 1)

  # Get the number of realizations
  N <- length(partition)

  # Get the replica statistics
  my <- colMeans(y)
  Cy <- stats::cov(y)
  traceCy <- sum(diag(Cy))
  Ry <- sqrtm(Cy)

  # Evaluate the upper bound of the indices
  V <- 2 * traceCy

  # Get the number of partitions
  M <- max(partition, na.rm = T)

  # Initialize the return structure
  Wk <- matrix(nrow = 3, ncol = M)
  n <- matrix(nrow = M)

  # Estimate the inner statistic
  for (m in seq(M)) {
    partition_element <- which(partition == m)
    Wk[, m] <- optimal_trasport_bw(partition_element, y, my, Cy, traceCy, Ry)
    n[m] <- length(partition_element)
  }

  # Calculate the 1d index
  W <- ((Wk %*% n) / (V * N))[1, 1]
  Adv <- ((Wk[2, ] %*% n) / (V * N))[1, 1]
  Diff <- ((Wk[3, ] %*% n) / (V * N))[1, 1]

  return(c(W, Adv, Diff, Wk[1, ] / V, Wk[2, ] / V, Wk[3, ] / V, V))
}

optimal_trasport_bw <- function(partition, y, my, Cy, traceCy, Ry) {
  # Get the current partition indices for y
  yc <- as.matrix(y[partition, ])

  # Get the conditioned statistics
  mc <- colMeans(yc)
  Cc <- stats::cov(yc)

  # Evaluate the advective and diffusive parts separately
  Adv <- sum((my - mc)^2)
  Diff <- traceCy + sum(diag(Cc)) - 2 * tracesqrtm(Ry %*% Cc %*% Ry)

  # Evaluate the index
  W <- Adv + Diff

  return(cbind(W, Adv, Diff))
}
