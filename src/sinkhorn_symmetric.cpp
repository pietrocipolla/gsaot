#include <RcppEigen.h>
using namespace Rcpp;

// Define the logsumexp for Eigen matrices
static Eigen::VectorXd logsumexp(const Eigen::MatrixXd x,
                          int axis) {
  if (axis == 0) {
    // Compute log-sum-exp along rows
    Eigen::VectorXd maxCoeffs = x.rowwise().maxCoeff();
    return maxCoeffs.array() + (x.colwise() - maxCoeffs).array().exp().rowwise().sum().log();
  }
  if (axis == 1) {
    // Compute log-sum-exp along columns
    Eigen::VectorXd maxCoeffs = x.colwise().maxCoeff();
    return maxCoeffs.transpose().array() + (x.rowwise() - maxCoeffs.transpose()).array().exp().colwise().sum().log();
  }

  // Invalid axis value
  return Eigen::VectorXd::Zero(x.rows());
}

// Define the symmetric Sinkhorn algorithm in log-space
// Based on formula (24) and (25) from the paper
List sinkhorn_symmetric_cpp(Eigen::VectorXd a,
                            Eigen::MatrixXd costMatrix,
                            int numIterations,
                            double epsilon,
                            double maxErr) {
  const int N = costMatrix.rows();
  const int M = costMatrix.cols();

  // Initialize the K matrix
  Eigen::MatrixXd K(- costMatrix / epsilon);
  K = K.array().exp();

  // Initialize f to zero
  Eigen::VectorXd f = Eigen::VectorXd::Ones(N);

  // Precompute log(alpha)
  Eigen::VectorXd log_a = a.array().log();

  // Initialize the marginal weights and the other useful variables
  Eigen::VectorXd estimated_marginal;

  // Initialize algorithm values
  int iter = 1;
  double err = std::numeric_limits<double>::infinity();
  Eigen::MatrixXd M_lse = Eigen::MatrixXd::Constant(N, M, 0.0);

  // Iterate until convergence
  while (iter == 1 || (iter <= numIterations &&
         err > maxErr && err < std::numeric_limits<double>::infinity())) {
    // Vectorized update of M_lse using Eigen broadcasting
    // M_lse[i,k] = log_a[k] + f[i]/epsilon - costMatrix[i,k]/epsilon
    M_lse = (f / epsilon).replicate(1, M) + log_a.transpose().replicate(N, 1) - (costMatrix / epsilon);

    // Compute logsumexp along columns (axis = 1)
    Eigen::VectorXd lse_result = logsumexp(M_lse, 1);

    // Update f using formula (25): f_i ← (1/2)(f_i - epsilon * LSE(...))
    f.array() = 0.5 * (f.array() - epsilon * lse_result.array());

    // Compute convergence error (same logic as sinkhorn.cpp - check every 10 iterations)
    if (iter % 10 == 0 || iter == 1) {
      // In symmetric case: estimated_marginal = exp(f/eps) .* (K * exp(f/eps))
      Eigen::VectorXd u_tilde = (f / epsilon).array().exp();
      estimated_marginal = (u_tilde.array() * (K * u_tilde).array()).inverse();
      // Rcout << "marginal " << estimated_marginal << std::endl;
      err = (estimated_marginal - a).cwiseAbs().sum();
      // Rcout << "err " << err << std::endl;
    }

    // Update iteration
    iter++;
  }

  if (iter >= numIterations) {
    Rcerr << "Warning: Maximum iterations reached without convergence" << std::endl;
  }

  double cost = 2.0 * f.dot(a);

  // Return results
  return List::create(
    Named("iter") = iter,
    Named("cost") = cost
  );
}

// Expose the symmetric Sinkhorn function to R
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List sinkhorn_symmetric(Eigen::VectorXd a,
                        Eigen::MatrixXd costm,
                        int numIterations,
                        double epsilon,
                        double maxErr) {
  return sinkhorn_symmetric_cpp(a, costm, numIterations, epsilon, maxErr);
}

/***R
set.seed(1)
N <- 10
M <- 10
alpha <- rep(1 / N, N)
C <- as.matrix(dist(rnorm(N)))
result <- sinkhorn_symmetric(alpha, C, 1e5, 0.01, 1e-6)
*/
