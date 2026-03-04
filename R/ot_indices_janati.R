#' Estimate Janati advective component of the Sinkhorn Divergence sensitivity indices
#'
#' @inheritParams ot_indices
#'
#' @inherit ot_indices return
#'
#' @export
#'
#' @seealso [ot_indices()], [ot_indices_1d()], [ot_indices_wb()], [ot_indices_div()]
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
#' ot_indices_janati(x, y, 10)
ot_indices_janati <- function(x,
                              y,
                              M,
                              discrete_out = FALSE,
                              solver = "sinkhorn",
                              solver_optns = NULL,
                              scaling = TRUE,
                              boot = FALSE,
                              stratified_boot = TRUE,
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

  # Check if the number of partitions is lower than the number of samples
  if (nrow(x) <= M) stop("The number of partitions should be lower than the number of samples")

  # Check if the output is a numerical
  if ((!is.numeric(y) | !is.matrix(y)))
    stop("`y` must be a matrix of numerical values!")

  # Check if the dimensions match
  if (!(nrow(x) == nrow(y)))
    stop("The number of samples in `x` and `y` should be the same")

  # Check the logical values
  if (!is.logical(discrete_out)) stop("`discrete_out` should be logical")
  if (!is.logical(scaling)) stop("`scaling` should be logical")

  # Check if the solver is present in the pool
  match.arg(solver, c("sinkhorn", "sinkhorn_log"))

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

  # Retrieve dimensions of the inputs
  N <- dim(x)[1]
  K <- dim(x)[2]

  # BUILD AN HISTOGRAM IF REQUESTED (otherwise assume equal weights)
  # ----------------------------------------------------------------------------
  # Default values for continuous output
  a <- rep(1 / N, N)
  y_unique <- NULL

  if (discrete_out) {
    # y_unique is always needed, in order to avoid mismatches between the
    # pre-computed cost matrix and the resampled outputs
    y_unique <- unique(y)

    # Function to build the histogram (maybe we can put this outside to make it
    # easier to change?). If bootstrap is requested, the histogram is computed
    # inside each statistic function
    if (!boot)
      a <- apply(y_unique, MARGIN = 1, FUN = function(u) sum(
        apply(y, 1, FUN = function(row) all(row == u))
      )) / N
  }

  if (identical(a, rep(1 / N, N)) & discrete_out & !boot)
    warning("The output is continuous, consider using `discrete_out=FALSE`")

  # BUILD COST MATRIX
  # ----------------------------------------------------------------------------
  if (discrete_out) {
    C <- as.matrix(stats::dist(y_unique, method = "euclidean")) ^ 2
  } else {
    C <- as.matrix(stats::dist(y, method = "euclidean")) ^ 2
  }

  # If the scaling is requested, scale the cost matrix by the highest value
  scaling_param <- 1
  if (scaling) {
    scaling_param <- max(C)
    C <- C / scaling_param
  }

  # DEFINE THE PARTITIONS FOR THE ESTIMATION
  # ----------------------------------------------------------------------------
  partitions <- build_partition(x, M)

  # DEFINE THE OT SOLVER
  # ----------------------------------------------------------------------------
  solver_fun <- switch(
    solver,
    "sinkhorn" = sinkhorn,
    "sinkhorn_log" = sinkhorn_stable,
    default = NULL
  )

  # Check the consistency of the options provided
  solver_optns <- check_solver_optns(solver, solver_optns)

  # BUILD THE RETURN STRUCTURE
  # ----------------------------------------------------------------------------
  W <- array(dim = K)
  Sobol <- array(dim = K)
  names(W) <- colnames(x)
  names(Sobol) <- colnames(x)
  IS <- list()
  if (boot) {
    V <- array(dim = K)
    traceCy <- array(dim = K)
    W_ci <- data.frame(matrix(nrow = K * 2,
                              ncol = 6,
                              dimnames = list(NULL,
                                              c("input", "component",
                                                "original", "bias",
                                                "low.ci", "high.ci"))))
    W_ci$input <- rep(names(W), times = 2)
    W_ci$component <- rep(c("janati", "sobol"), each = K)
    IS_ci <- list()
    V_ci <- list()
    traceCy_ci <- list()

    W_boot <- list()
  }

  # Evaluate the upper bound of the indices (if no boot), otherwise it is
  # evaluated inside the bootstrap function
  # ----------------------------------------------------------------------------
  if (!boot) {
    global_correction <- do.call(solver_fun,
                                 c(list(
                                   a = a, b = a, costm = C
                                 ),
                                 solver_optns))$cost
    V <- (higher_bound(C) - 0.5 * global_correction) * scaling_param
    traceCy <- 0.5 * higher_bound(C) * scaling_param
  }

  # Calculate global mean for y (needed for Janati method)
  my <- colMeans(y)

  # ESTIMATE THE INDICES FOR EACH PARTITION
  # ----------------------------------------------------------------------------
  for (k in seq_len(K)) {
    # Get the partitions for the current input
    partition <- partitions[, k]

    # Set the number of partition elements
    M <- max(partition, na.rm = T)

    # NO BOOTSTRAP ESTIMATION
    if (!boot) {
      # Build the inner return structure
      Wk <- matrix(nrow = 1, ncol = M)
      n <- matrix(nrow = M)

      for (m in seq(M)) {
        # Retrieve the elements in the partition
        partition_element <- which(partition == m)
        Wk[, m] <- optimal_trasport_janati(partition_element, y, my)
        n[m] <- length(partition_element)
      }

      W[k] <- ((Wk[1,] %*% n) / (V * N))[1, 1]
      Sobol[k] <- ((Wk[1,] %*% n) / (traceCy * N))[1, 1]

      IS[[k]] <- rbind(Wk / V, Wk / traceCy)
    } else {
      # BOOTSTRAP ESTIMATION With stratification we pass the partitioned input.
      # Otherwise, we pass the input and we partition after the resampling
      if (stratified_boot) {
        dat <- cbind(partition, y)
        strata <- partition
      }
      else {
        dat <- cbind(x[, k], y)
        strata <- rep(1, N)
      }

      # Do boostrap estimation stratified by the partitions
      W_boot[[k]] <- boot::boot(data = dat,
                                statistic = ot_boot_janati,
                                R = R,
                                strata = strata,
                                discrete_out = discrete_out,
                                y = y,
                                M = M,
                                C = C,
                                y_unique = y_unique,
                                solver_fun = solver_fun,
                                solver_optns = solver_optns,
                                scaling_param = scaling_param,
                                parallel = parallel,
                                ncpus = ncpus,
                                stratified_boot = stratified_boot)

      # Transform the results into readable quantities
      W_stats <- bootstats(W_boot[[k]], type = type, conf = conf)

      # Save the results
      # The bootstrap returns: W, Sobol, Wk/V (M values), Wk/traceCy (M values), V, traceCy
      # Total: 2 + 2*M + 2 = 2*M + 4 elements
      W[k] <- W_stats$index[1]
      Sobol[k] <- W_stats$index[2]
      IS[[k]] <- rbind(W_stats$index[3:(M + 2)],
                       W_stats$index[(M + 3):(2 * M + 2)])

      # Bootstrap estimate of the variance bounds
      V[k] <- W_stats$index[2 * M + 3]
      traceCy[k] <- W_stats$index[2 * M + 4]

      # Boostrap estimates of the indices
      W_ci[k, 3:6] <- c(W_stats$original[1], W_stats$bias[1],
                        W_stats$low.ci[1], W_stats$high.ci[1])
      W_ci[K + k, 3:6] <- c(W_stats$original[2], W_stats$bias[2],
                            W_stats$low.ci[2], W_stats$high.ci[2])

      # Bootstrap estimates of the inner statistics
      IS_ci[[k]] <- rbind(cbind(W_stats$low.ci[3:(M + 2)],
                                W_stats$high.ci[3:(M + 2)]),
                          cbind(W_stats$low.ci[(M + 3):(2 * M + 2)],
                                W_stats$high.ci[(M + 3):(2 * M + 2)]))

      # Bootstrap estimates for the upper bound
      V_ci[[k]] <- cbind(W_stats$low.ci[2 * M + 3], W_stats$high.ci[2 * M + 3])
      traceCy_ci[[k]] <- cbind(W_stats$low.ci[2 * M + 4], W_stats$high.ci[2 * M + 4])
    }
  }

  if (boot) {
    out <- gsaot_indices(method = "janati",
                         indices = W,
                         bound = V,
                         IS = IS,
                         partitions = partitions,
                         x = x, y = y,
                         Sobol = Sobol,
                         indices_ci = W_ci,
                         bound_ci = V_ci,
                         IS_ci = IS_ci,
                         R = R,
                         conf = conf,
                         type = type,
                         W_boot = W_boot)
    
    # Store traceCy bound information
    out$traceCy <- traceCy
    out$traceCy_ci <- traceCy_ci

    return(out)
  }

  out <- gsaot_indices(method = "janati",
                       indices = W,
                       bound = V,
                       IS = IS,
                       partitions = partitions,
                       x = x, y = y,
                       Sobol = Sobol)
  
  # Store traceCy bound information
  out$traceCy <- traceCy

  return(out)
}

# Computation of bootstrap statistics
ot_boot_janati <- function(d,
                           i,
                           discrete_out,
                           y,
                           M,
                           C,
                           y_unique,
                           solver_fun,
                           solver_optns,
                           scaling_param,
                           stratified_boot) {
  # Retrieve the partitions
  if (stratified_boot)
    partition <- d[i, 1]
  else
    partition <- build_partition(as.matrix(d[i, 1]), M)

  # Retrieve the output
  y <- matrix(d[i, 2:ncol(d)], ncol = ncol(d) - 1)

  # Get the number of realizations and partitions
  N <- length(partition)
  M <- max(partition, na.rm = T)

  # Get the replica statistics
  my <- colMeans(y)

  # Compute the unconditioned histogram
  if (discrete_out) {
    a <- apply(y_unique, MARGIN = 1, FUN = function(u) sum(
      apply(y, 1, FUN = function(row) all(row == u))
    )) / N } else {
    a <- rep(1 / N, N)
    }

  if (!discrete_out)
    C <- C[i, i]

  # Compute the upper bound
  global_correction <- do.call(solver_fun,
                               c(list(
                                 a = a, b = a, costm = C
                               ),
                               solver_optns))$cost
  V <- (higher_bound(C) - 0.5 * global_correction) * scaling_param
  traceCy <- 0.5 * higher_bound(C) * scaling_param

  # Initialize the return structure
  Wk <- matrix(nrow = 1, ncol = M)
  n <- matrix(nrow = M)

  for (m in seq(M)) {
    # Retrieve the elements in the partition
    partition_element <- which(partition == m)
    n[m] <- length(partition_element)

    Wk[, m] <- optimal_trasport_janati(partition_element, y, my)
  }

  # Calculate the ot index
  W <- ((Wk %*% n) / (V * N))[1, 1]
  Sobol <- ((Wk %*% n) / (traceCy * N))[1, 1]

  return(c(W, Sobol, Wk / V, Wk / traceCy, V, traceCy))
}

# Function for computing the advective part of the OT distance between the
# unconditioned and conditioned distributions, as in Janati et al. (2019)
optimal_trasport_janati <- function(partition, y, my) {
  # Get the current partition indices for y
  yc <- as.matrix(y[partition, ])

  # Get the conditioned statistics
  mc <- colMeans(yc)

  # Evaluate the advective and diffusive parts separately
  Adv <- sum((my - mc)^2)

  return(Adv)
}
