#' #' Title
#' #'
#' #' @param x a data.frame containing the input(s) values
#' #' @param y a data.frame containing the outputs values
#' #' @param M a scalar representing the number of partitions for continuous inputs
#' #' @param eps a double representing the coefficient of the entropic regularization. If negative the value is scaled according to the maximum of the cost matrix.
#' #' @param num_iterations maximum number of iterations of the Sinkhorn algorithm
#' #'
#' #' @return A sensitivity index between 0 and 1 for each of the columns in x
#' #' @export
#' #'
#' #' @examples
#' #' N <- 1000
#' #'
#' #' mx <- c(1, 1, 1)
#' #' Sigmax <- matrix(data = c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), nrow = 3)
#' #'
#' #' x1 <- rnorm(N)
#' #' x2 <- rnorm(N)
#' #' x3 <- rnorm(N)
#' #'
#' #' x <- cbind(x1, x2, x3)
#' #' x <- mx + x %*% chol(Sigmax)
#' #'
#' #' A <- matrix(data = c(4, -2, 1, 2, 5, -1), nrow = 2, byrow = TRUE)
#' #' y <- t(A %*% t(x))
#' #'
#' #' x <- data.frame(x)
#' #' y <- data.frame(y)
#' #'
#' ot_indices_divergence <- function(x, y, M, eps,
#'                                 num_iterations = 1e8) {
#'   # Input checks
#'   stopifnot(is.data.frame(x), is.data.frame(y))
#'   stopifnot(dim(x)[1] == dim(y)[1])
#'   stopifnot(dim(x)[1] > M)
#'
#'   # Build cost matrix
#'   C <- as.matrix(stats::dist(y, method = "euclidean"))^2
#'   scaling <- 1
#'   if (eps < 0) {
#'     scaling <- max(C)
#'     eps <- -eps
#'     C <- C / scaling
#'   }
#'
#'   # Retrieve values useful to the algorithm
#'   N <- dim(x)[1]
#'   K <- dim(x)[2]
#'
#'   # Define the partitions for the estimation
#'   partitions <- build_partition(x, M)
#'
#'   # Build return structure
#'   W <- array(dim = K)
#'   Wpart <- array(dim = K)
#'   Wtot <- optimal_transport_sinkhorn(C, num_iterations, eps)$W22 * scaling
#'
#'   # Evaluate the upper bound of the indices
#'   V <- 2 * sum(diag(stats::cov(y)))
#'
#'   for (k in seq_len(K)) {
#'     # Get the current partition
#'     partition <- partitions[[k]]
#'
#'     # Get the current number of partition elements
#'     M <- length(partition)
#'
#'     # Build the inner return structure
#'     Wk <- matrix(nrow = 1, ncol = M)
#'     n <- matrix(nrow = M)
#'
#'     for (m in seq(M)) {
#'       partition_element <- partition[[m]]
#'       ret <- optimal_transport_sinkhorn(C[partition_element, ], num_iterations, eps)
#'       ret_part <- optimal_transport_sinkhorn(C[partition_element, partition_element], num_iterations, eps)
#'       Wk[, m] <- (ret$W22 - ret_part$W22) * scaling
#'       n[m] <- length(partition_element)
#'     }
#'
#'     W[k] <- ((Wk[1, ] %*% n) / (V * N))[1, 1] - Wtot / V
#'   }
#'
#'   return(W)
#' }
