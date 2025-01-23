#' Calculate Optimal Transport sensitivity indices for multivariate y
#'
#' @description `ot_indices` calculates sensitivity indices using Optimal
#'   Transport (OT) for a multivariate output sample `y` with respect to input
#'   data `x`. Sensitivity indices measure the influence of inputs on outputs,
#'   with values ranging between 0 and 1.
#'
#' @param x A matrix or data.frame containing the input(s) values. The values
#'   can be numeric, factors, or strings. The type of data changes the
#'   partitioning. If the values are continuous (double), the function
#'   partitions the data into `M` sets. If the values are discrete (integers,
#'   strings, factors), the number of partitioning sets is data-driven.
#' @param y A matrix containing the output values. Each column represents a
#'   different output variable, and each row represents a different observation.
#'   Only numeric values are allowed.
#' @param M A scalar representing the number of partitions for continuous
#'   inputs.
#' @param cost (default `"L2"`) A string or function defining the cost function
#'   of the Optimal Transport problem. It should be "L2" or a function taking as
#'   input y and returning a cost matrix. If `cost="L2"`, `ot_indices` uses the
#'   squared Euclidean metric.
#' @param discrete_out (default `FALSE`) Logical, by default the output sample
#'   in `y` are equally weighted. If `discrete_out=TRUE`, the function tries to
#'   create an histogram of the realizations and to use the histogram as
#'   weights. It works if the output is discrete or mixed and the number of
#'   realizations is large. The advantage of this option is to reduce the
#'   dimension of the cost matrix.
#' @param solver Solver for the Optimal Transport problem. Currently supported
#'   options are:
#' * `"sinkhorn"` (default), the Sinkhorn's solver \insertCite{cuturi2013sinkhorn}{gsaot}.
#' * `"sinkhorn_log"`, the Sinkhorn's solver in log scale \insertCite{peyre2019computational}{gsaot}.
#' * `"transport"`, a solver of the non regularized OT problem using [transport::transport()].
#' @param solver_optns (optional) A list containing the options for the Optimal
#'   Transport solver. See details for allowed options and default ones.
#' @param scaling (default `TRUE`) Logical that sets whether or not to scale the
#'   cost matrix.
#' @param boot (default `FALSE`) Logical that sets whether or not to perform
#'   bootstrapping of the OT indices.
#' @param stratified_boot (default `FALSE`) Logical that sets the type of
#'   resampling performed. With `stratified_boot=FALSE`, the function resamples
#'   the dataset and then creates the partitions. Otherwise, first, it
#'   creates the partitions and then it performs stratified bootstrapping with
#'   strata being the partitions.
#' @param R (default `NULL`) Positive integer, number of bootstrap replicas.
#' @param parallel (default `"no"`) The type of parallel operation to be used
#'   (if any). If missing, the default is taken from the option `boot.parallel`
#'   (and if that is not set, `"no"`). Only considered if `boot = TRUE`. For
#'   more information, check the [boot::boot()] function.
#' @param ncpus (default `1`) Positive integer: number of processes to be used
#'   in parallel operation: typically one would chose this to the number of
#'   available CPUs. Check the `ncpus` option in the [boot::boot()] function of
#'   the boot package.
#' @param conf (default `0.95`) Number between `0` and `1` representing the
#'   confidence level. Only considered if `boot = TRUE`.
#' @param type (default `"norm"`) Method to compute the confidence interval.
#'   Only considered if `boot = TRUE`. For more information, check the `type`
#'   option of [boot::boot.ci()].
#'
#' @details ## Solvers
#'
#'   OT is a widely studied topic in Operational Research and Calculus. The
#'   reference for the OT solvers in this package is
#'   \insertCite{peyre2019computational;textual}{gsaot}. The default solver is
#'   `"sinkhorn"`, the Sinkhorn's solver introduced in
#'   \insertCite{cuturi2013sinkhorn;textual}{gsaot}. It solves the
#'   entropic-regularized version of the OT problem. The `"sinkhorn_log"` solves
#'   the same OT problem but in log scale. It is more stable for low values of
#'   the regularization parameter but slower to converge. The option
#'   `"transport"` is used to choose a solver for the non-regularized OT
#'   problem. Under the hood, the function calls [transport::transport()] from
#'   package `transport`. This option does not define the solver per se, but the
#'   solver should be defined with the argument `solver_optns`. See the next
#'   section for more information.
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
#'   The solver `"transport"` has the parameters:
#'   * `method` (default `"networkflow`): string defining the solver of the OT
#'   problem.
#'   * `control`: a named list of parameters for the chosen method or the result
#'   of a call to [transport::trcontrol()].
#'   * `threads` (default `1`): an Integer specifying the number of threads used
#'   in parallel computing.
#'
#'   For details regarding this solver, check the [transport::transport()] help
#'   page.
#'
#' @returns A `gsaot_indices` object containing:
#' * `method`: a string that identifies the type of indices computed.
#' * `indices`: a names array containing the sensitivity indices between 0 and 1
#'   for each column in x, indicating the influence of each input variable on
#'   the output variables.
#' * `bound`: a double representing the upper bound of the separation measure or
#'   an array representing the mean of the separation for each input according
#'   to the bootstrap replicas.
#' * `x`, `y`: input and output data provided as arguments of the function.
#' * `inner_statistic`: a list of matrices containing the values of the inner
#'   statistics for the partitions defined by `partitions`. If `method =
#'   wasserstein-bures`, each matrix has three rows containing the
#'   Wasserstein-Bures indices, the Advective, and the Diffusive components.
#' * `partitions`: a matrix containing the partitions built to calculate the
#'   sensitivity indices. Each column contains the partition associated to the
#'   same column in `x`. If `boot = TRUE`, the object contains also:
#' * `indices_ci`: a `data.frame` with first column the input, second and third
#'   columns the lower and upper bound of the confidence interval.
#' * `inner_statistic_ci`: a list of matrices. Each element of the list contains
#'   the lower and upper confidence bounds for the partition defined by the row.
#' * `bound_ci`: a list containing the lower and upper bounds of the confidence
#'   intervals of the separation measure bound.
#' * `type`, `conf`: type of confidence interval and confidence level, provided
#'   as arguments.
#'
#' @seealso [ot_indices_1d()], [ot_indices_wb()]
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
#' sensitivity_indices <- ot_indices(x, y, M)
#' sensitivity_indices
#'
#' @export
#'
ot_indices <- function(x,
                       y,
                       M,
                       cost = "L2",
                       discrete_out = FALSE,
                       solver = "sinkhorn",
                       solver_optns = NULL,
                       scaling = TRUE,
                       boot = FALSE,
                       stratified_boot = FALSE,
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
  if (!is.numeric(y) | !is.matrix(y)) stop("`y` must be a matrix of numerical values!")

  # Check if the dimensions match
  if (!(nrow(x) == nrow(y))) stop("The number of samples in `x` and `y` should be the same")

  # Check if the number of partitions is lower than the number of samples
  if (nrow(x) <= M) stop("The number of partitions should be lower than the number of samples")

  # Check if the cost is defined correctly
  cost_type <- identical(cost, "L2")
  if (!cost_type & !is.function(cost)) stop("`cost` should be \"L2\" or a function")

  # Check the logical values
  if (!is.logical(discrete_out)) stop("`discrete_out` should be logical")
  if (!is.logical(scaling)) stop("`scaling` should be logical")

  # Check if the solver is present in the pool
  match.arg(solver, c("sinkhorn", "sinkhorn_stable", "transport"))

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
    "sinkhorn_stable" = sinkhorn_stable,
    "transport" = transport::transport,
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
                              ncol = 3,
                              dimnames = list(NULL,
                                              c("Inputs", "low.ci", "high.ci"))))
    W_ci$Inputs <- names(W)
    IS_ci <- list()
    V_ci <- list()
  }

  # Evaluate the upper bound of the indices (if no boot), otherwise it is
  # evaluated inside the bootstrap function
  # ----------------------------------------------------------------------------
  if (!boot)
    V <- higher_bound(C) * scaling_param

  # ESTIMATE THE INDICES FOR EACH PARTITION
  # ----------------------------------------------------------------------------
  for (k in seq_len(K)) {
    # Get the partitions for the current input
    partition <- partitions[, k]

    # Set the number of partition elements
    M <- max(partition)

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

        # Save the estimated distance
        Wk[, m] <- ret$cost * scaling_param
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
      W_boot <- boot::boot(data = dat,
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
      W_stats <- bootstats(W_boot, type = type, conf = conf)

      # Save the results
      W[k] <- W_stats$original[1]
      IS[[k]] <- matrix(W_stats$original[2:(M + 1)], nrow = 1)
      V[k] <- W_stats$original[M + 2]
      W_ci[k, 2:3] <- c(W_stats$low.ci[1], W_stats$high.ci[1])
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
                         type = type)

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
  M <- max(partition)

  # Compute the unconditioned histogram
  a <- apply(y_unique, MARGIN = 1, FUN = function(u) sum(
    apply(y, 1, FUN = function(row) all(row == u))
  )) / N

  # Compute the variance
  V <- higher_bound(C) * scaling_param

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

    # Save the estimated distance
    Wk[, m] <- ret$cost * scaling_param
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
  M <- max(partition)

  # Build the histogram
  a <- rep(1 / N, N)

  # Compute the variance
  V <- higher_bound(C) * scaling_param

  # Get the number of partitions
  M <- max(partition)

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

    # Save the estimated distance
    Wk[, m] <- ret$cost * scaling_param
  }

  # Calculate the ot index
  W <- ((Wk %*% n) / (V * N))[1, 1]

  return(c(W, Wk / V, V))
}
