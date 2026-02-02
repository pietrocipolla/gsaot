#' Estimate Sinkhorn divergence sensitivity indices for generic outputs
#'
#' @description `ot_indices_div` calculates sensitivity indices using Sinkhorn
#'   divergence for a multivariate output sample `y` with respect to input
#'   data `x`. Sensitivity indices measure the influence of inputs on outputs,
#'   with values ranging between 0 and 1.
#'
#' @inheritParams ot_indices
#'
#' @details ## Solvers
#'
#'   The default solver is
#'   `"sinkhorn"`, the Sinkhorn's solver introduced in
#'   \insertCite{cuturi2013sinkhorn;textual}{gsaot}. It solves the
#'   entropic-regularized version of the OT problem. The `"sinkhorn_log"` solves
#'   the same OT problem but in log scale. It is more stable for low values of
#'   the regularization parameter but slower to converge. Differently from [ot_indices()], the option
#'   `"transport"` is not available. See the next section for more information.
#'
#'   ## Solver options
#'
#'   The argument `solver_optns` should be empty (for default options) or a list
#'   with all or some of the required solver parameters. All the parameters not
#'   included in the list will be set to default values. The solvers
#'   `"sinkhorn"` and `"sinkhorn_log"` have the same options:
#'   * `numIterations` (default `1e3`): a positive integer defining the maximum number
#'   of Sinkhorn's iterations allowed. If the solver does not converge in the
#'   number of iterations set, the solver will throw an error.
#'   * `epsilon` (default `0.01`): a positive real number defining the regularization
#'   coefficient. If the value is too low, the solver may return `NA`.
#'   * `maxErr` (default `1e-9`): a positive real number defining the
#'   approximation error threshold between the marginal histogram of the
#'   partition and the one computed by the solver. The solver may fail to
#'   converge in `numIterations` if this value is too low.
#'
#' @inherit ot_indices return
#'
#' @seealso [ot_indices_1d()], [ot_indices_wb()], [ot_indices()]
#'
#' @references \insertAllCited{}
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
#'
#' M <- 25
#'
#' # Calculate sensitivity indices
#' sensitivity_indices <- ot_indices_div(x, y, M)
#' sensitivity_indices
#'
#' @export
#'
ot_indices_div <- function(x,
                           y,
                           M,
                           cost = "L2",
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

  # Check if the cost is defined correctly
  cost_type <- identical(cost, "L2")
  if (!cost_type & !is.function(cost)) stop("`cost` should be \"L2\" or a function")

  # Check if the output is a numerical
  if (cost_type & (!is.numeric(y) | !is.matrix(y)))
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
  y <- y[!y_na, ]
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
    if (cost_type) {
      C <- as.matrix(stats::dist(y_unique, method = "euclidean")) ^ 2
    } else {
      C <- cost(y_unique)
    }
  } else {
    if (cost_type) {
      C <- as.matrix(stats::dist(y, method = "euclidean")) ^ 2
    } else {
      C <- cost(y)
    }
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
  solver_fun <- switch (
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
  }

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
        n[m] <- length(partition_element)

        # Build the partition histogram. The idea is to use the marginal histogram
        # b both as histogram and as selector for the cost matrix parts. Indeed,
        # as.logical transforms everything > 0 to TRUE and everything else t FALSE
        if (discrete_out) {
          b <- apply(y_unique, MARGIN = 1, FUN = function(u) sum(
            apply(y[partition_element, ], 1, FUN = function(row) all(row == u))
          )) / n[m]
        } else {
          b <- numeric(length = N)
          b[partition_element] <- 1 / n[m]
        }

        # Call the solver
        ret <- do.call(solver_fun,
                       c(list(a = a, b = b[b > 0], costm = C[, (b > 0)]),
                         solver_optns))

        # Quick exit from the loop if NaNs are present
        if (is.na(ret$cost)) stop("NA are present, consider increasing the
                                  regularization parameter or change solver")

        # Compute the local correction term
        local_correction <- do.call(solver_fun,
                                    c(list(a = b[b > 0], b = b[b > 0],
                                           costm = C[(b > 0), (b > 0)]),
                                      solver_optns))$cost

        # Save the estimated distance
        Wk[, m] <- (ret$cost - 0.5 * local_correction - 0.5 * global_correction) * scaling_param
      }

      W[k] <- ((Wk[1,] %*% n) / (V * N))[1, 1]
      IS[[k]] <- Wk / V
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
                                statistic = ot_boot,
                                R = R,
                                strata = strata,
                                discrete_out = discrete_out,
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
    out <- gsaot_indices(method = solver,
                         solver_optns = solver_optns,
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

  out <- gsaot_indices(method = solver,
                       solver_optns = solver_optns,
                       indices = W,
                       bound = V,
                       IS = IS,
                       partitions = partitions,
                       x = x, y = y)

  return(out)
}

ot_boot <- function(d,
                    i,
                    discrete_out,
                    M,
                    C,
                    y_unique,
                    solver_fun,
                    solver_optns,
                    scaling_param,
                    stratified_boot) {
  # According to discrete_out select the correct function
  if (discrete_out) {
    ot_boot_discrete(d, i, M, C, y_unique, solver_fun, solver_optns, scaling_param, stratified_boot)
  } else {
    ot_boot_cont(d, i, M, C, solver_fun, solver_optns, scaling_param, stratified_boot)
  }
}

# Bootstrap for discrete output
ot_boot_discrete <- function(d,
                             indices,
                             M,
                             C,
                             y_unique,
                             solver_fun,
                             solver_optns,
                             scaling_param,
                             stratified_boot) {
  # Retrieve the partitions
  if (stratified_boot)
    partition <- d[indices, 1]
  else
    partition <- build_partition(as.matrix(d[indices, 1]), M)

  # Retrieve the output
  y <- d[indices, 2:ncol(d)]

  # Get the number of realizations and partitions
  N <- length(partition)
  M <- max(partition, na.rm = T)

  # Compute the unconditioned histogram
  a <- apply(y_unique, MARGIN = 1, FUN = function(u) sum(
    apply(y, 1, FUN = function(row) all(row == u))
  )) / N

  # Compute the upper bound
  global_correction <- do.call(solver_fun,
                               c(list(
                                 a = a, b = a, costm = C
                               ),
                               solver_optns))$cost
  V <- (higher_bound(C) - 0.5 * global_correction) * scaling_param

  # Initialize the return structure
  Wk <- matrix(nrow = 1, ncol = M)
  n <- matrix(nrow = M)

  for (m in seq(M)) {
    # Retrieve the elements in the partition
    partition_element <- which(partition == m)
    n[m] <- length(partition_element)

    # Build the partition histogram
    b <- apply(y_unique, MARGIN = 1, FUN = function(u) sum(
      apply(y[partition_element, ], 1, FUN = function(row) all(row == u))
    )) / n[m]

    # Call the solver
    ret <- do.call(solver_fun,
                   c(list(a = a[a > 0],
                          b = b[b > 0],
                          costm = C[(a > 0), (b > 0)]),
                     solver_optns))

    # Quick exit from the loop if NaNs are present
    if (is.na(ret$cost)) stop("NA are present, consider increasing the
                              regularization parameter or change solver")

    # Compute the local correction term
    local_correction <- do.call(solver_fun,
                                c(list(a = b[b > 0], b = b[b > 0],
                                       costm = C[(b > 0), (b > 0)]),
                                  solver_optns))$cost

    # Save the estimated distance
    Wk[, m] <- (ret$cost - 0.5 * local_correction - 0.5 * global_correction) * scaling_param
  }

  # Calculate the ot index
  W <- ((Wk %*% n) / (V * N))[1, 1]

  return(c(W, Wk / V, V))
}

# Bootstrap for continuous output
ot_boot_cont <- function(d,
                         indices,
                         M,
                         C,
                         solver_fun,
                         solver_optns,
                         scaling_param,
                         stratified_boot) {
  # Retrieve the partitions
  if (stratified_boot)
    partition <- d[indices, 1]
  else
    partition <- build_partition(as.matrix(d[indices, 1]), M)

  # Retrieve the cost matrix
  C <- C[indices, indices]

  # Get the number of realizations and partitions
  N <- length(partition)
  M <- max(partition, na.rm = T)

  # Build the histogram
  a <- rep(1 / N, N)

  # Compute the upper bound
  global_correction <- do.call(solver_fun,
                               c(list(
                                 a = a, b = a, costm = C
                               ),
                               solver_optns))$cost
  V <- (higher_bound(C) - 0.5 * global_correction) * scaling_param

  # Initialize the return structure
  Wk <- matrix(nrow = 1, ncol = M)
  n <- matrix(nrow = M)

  for (m in seq(M)) {
    # Retrieve the elements in the partition
    partition_element <- which(partition == m)
    n[m] <- length(partition_element)

    # Build the partition histogram.
    b <- numeric(length = N)
    b[partition_element] <- 1 / n[m]

    # Call the solver
    ret <- do.call(solver_fun,
                   c(list(
                     a = a, b = b[b > 0], costm = C[, (b > 0)]
                   ),
                   solver_optns))

    # Quick exit from the loop if NaNs are present
    if (is.na(ret$cost))
      stop("NA are present, consider increasing the
                                  regularization parameter or change solver")

    # Compute the local correction term
    local_correction <- do.call(solver_fun,
                                c(list(a = b[b > 0], b = b[b > 0],
                                       costm = C[(b > 0), (b > 0)]),
                                  solver_optns))$cost

    # Save the estimated distance
    Wk[, m] <- (ret$cost - 0.5 * local_correction - 0.5 * global_correction) * scaling_param
  }

  # Calculate the ot index
  W <- ((Wk %*% n) / (V * N))[1, 1]

  return(c(W, Wk / V, V))
}
