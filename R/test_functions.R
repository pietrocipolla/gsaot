#' Sobol G function evaluation
#'
#' This function evaluates the Sobol G function on a set of input samples
#' generated via crude Monte Carlos. It returns both the sampled
#' inputs and the corresponding function outputs.
#'
#' @param N Integer. Number of input samples to generate.
#' @param a (default: \code{c(0, 1, 4.5, 9, 99, 99, 99, 99)}) Numeric vector of non-negative parameters of length 8.
#' These parameters control the sensitivity of each input dimension.
#'
#' @details
#' The Sobol G function is defined as:
#' \deqn{
#'   Y = \prod_{j=1}^{8} \frac{|\ 4 X_j - 2\ | + a_j}{1 + a_j}
#' }
#' where \eqn{X_j \sim \mathcal{U}(0, 1)} independently.
#'
#' @returns A list with two elements:
#' * `x`: a numeric matrix of size \code{N x 8} containing the input samples.
#' * `y`: a numeric vector of length \code{N} with the corresponding function outputs.
#' @export
#'
#' @examples
#' result <- sobol_fun(1000)
#' head(result$x)
#' head(result$y)
#'
#' @seealso \code{\link{ishi_homma_fun}}, \code{\link{gaussian_fun}}
#'
sobol_fun <- function(N, a = c(0, 1, 4.5, 9, 99, 99, 99, 99)) {
  # Generate input samples using Latin Hypercube Sampling
  x <- sobol_fun_inputs(N)

  # Compute the Sobol function output
  y <- sobol_fun_model(x, a)

  # Return the inputs and corresponding outputs as a list
  return(list(x = x, y = y))
}

# Function to generate input samples using Latin Hypercube Sampling
sobol_fun_inputs <- function(N) {
  # Generate N samples from a uniform distribution over [0, 1]^8
  x <- matrix(stats::runif(N * 8), ncol = 8)

  # Assign column names
  colnames(x) <- paste0("X", 1:8)

  return(x)
}

# Function to compute the Sobol function model output
sobol_fun_model <- function(x, a) {
  # Initialize the output variable y
  y <- 1

  # Loop over the 8 input dimensions
  for (j in 1:8) {
    # Compute the product of each term in the Sobol function
    y <- y * (abs(4 * x[, j] - 2) + a[j]) / (1 + a[j])
  }

  # Return the computed response
  return(y)
}

#' Ishigami-Homma function evaluation
#'
#' Evaluates the Ishigami-Homma function.
#' Input samples are drawn from a uniform distribution over \eqn{[-\pi, \pi]^3}
#'
#' @param N Number of input samples to generate.
#' @param A (default: \code{2}) Numeric, amplitude of the second sine component .
#' @param B  (default: \code{1}) Numeric, coefficient of the interaction term.
#'
#' @returns A list with two elements:
#' * `x`: a numeric matrix of size \code{N x 8} containing the input samples.
#' * `y`: a numeric vector of length \code{N} with the corresponding function outputs.
#'
#' @details The Ishigami-Homma function is defined as:
#' \deqn{Y = \sin(X_1) + A \cdot \sin^2(X_2) + B \cdot X_3^4 \cdot \sin(X_1)}
#' where \eqn{X_i \sim \mathcal{U}(-\pi, \pi)}.
#'
#' @examples
#' result <- ishi_homma_fun(1000)
#' head(result$x)
#' head(result$y)
#'
#' @seealso \code{\link{sobol_fun}}, \code{\link{gaussian_fun}}
#'
#' @export
ishi_homma_fun <- function(N, A = 2, B = 1) {
  # Generate input samples using Latin Hypercube Sampling
  x <- ishi_homma_input(N)

  # Compute the Ishigami-Homma function output
  y <- ishi_homma_model(x, A, B)

  # Return the inputs and corresponding outputs as a list
  return(list(x = x, y = y))
}

# Function to generate input samples using Latin Hypercube Sampling
ishi_homma_input <- function(N) {
  # Generate N samples from a uniform distribution over [-pi, pi]^3
  x <- matrix(stats::runif(N * 8, min = -pi, max = pi), ncol = 3)

  # Assign column names
  colnames(x) <- paste0("X", 1:3)

  return(x)
}

# Function to compute the Ishigami-Homma function model output
ishi_homma_model <- function(x, A, B) {
  # Compute the Ishigami-Homma function for each sample
  y <- sin(x[, 1]) + A * sin(x[, 2]) ^ 2 + B * x[, 3] ^ 4 * sin(x[, 1])

  return(y)
}

#' Multivariate Gaussian linear model evaluation
#'
#' Generates samples from a multivariate Gaussian distribution and evaluates a simple
#' linear transformation model.
#'
#' @param N Number of input samples to generate.
#'
#' @returns A list with two elements:
#' * `x`: a numeric matrix of size \code{N x 8} containing the input samples.
#' * `y`: a numeric vector of length \code{N} with the corresponding function outputs.
#'
#' @details
#' Inputs \code{x} are sampled from:
#' \deqn{
#' \mathbf{X} \sim \mathcal{N}(\boldsymbol{\mu}, \Sigma), \quad \boldsymbol{\mu} = [1, 1, 1], \quad \Sigma = \begin{bmatrix} 1 & 0.5 & 0.5 \\ 0.5 & 1 & 0.5 \\ 0.5 & 0.5 & 1 \end{bmatrix}
#' }
#'
#' The output is given by:
#' \deqn{
#' \mathbf{Y} = A \mathbf{X}^{\top}, \quad A = \begin{bmatrix} 4 & -2 & 1 \\ 2 & 5 & -1 \end{bmatrix}
#' }
#'
#' @examples
#' result <- gaussian_fun(1000)
#' head(result$x)
#' head(result$y)
#'
#' @seealso \code{\link{sobol_fun}}, \code{\link{ishi_homma_fun}}
#'
#' @export
gaussian_fun <- function(N) {
  # Generate correlated input samples
  x <- gaussian_inputs(N)

  # Compute the model output
  y <- gaussian_model(x)

  # Return the inputs and corresponding outputs as a list
  return(list(x = x, y = y))
}

# Function to generate correlated input samples
gaussian_inputs <- function(N) {
  # Define the mean vector for the multivariate normal distribution
  mx <- c(1, 1, 1)

  # Define the covariance matrix
  Sigmax <- matrix(data = c(1, 0.5, 0.5,
                            0.5, 1, 0.5,
                            0.5, 0.5, 1), nrow = 3)

  # Sample from standard normal distributions
  x1 <- stats::rnorm(N)
  x2 <- stats::rnorm(N)
  x3 <- stats::rnorm(N)

  # Combine the sampled values into a matrix
  x <- cbind(x1, x2, x3)

  # Apply Cholesky transformation to induce correlation
  x <- mx + x %*% chol(Sigmax)

  # Assign column names
  colnames(x) <- paste0("X", 1:3)

  return(x)
}

# Function to compute the model output
gaussian_model <- function(x) {
  # Define the coefficient matrix for the linear transformation
  A <- matrix(data = c(4, -2, 1,
                       2, 5, -1), nrow = 2, byrow = TRUE)

  # Compute the output y = Ax'
  y <- t(A %*% t(x))

  # Assign column names to the output
  colnames(y) <- c("Y1", "Y2")

  return(y)
}
