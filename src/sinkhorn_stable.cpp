#include <RcppEigen.h>
using namespace Rcpp;

// Define the logsumexp for Eigen matrices
Eigen::VectorXd logsumexp(const Eigen::MatrixXd x,
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

// Define the kernel matrix computation
Eigen::MatrixXd getK(Eigen::MatrixXd costMatrix,
                     Eigen::VectorXd alpha,
                     Eigen::VectorXd beta,
                     double epsilon) {
  // Remove marginals
  costMatrix.colwise() -= alpha;
  costMatrix.rowwise() -= beta.transpose();

  // Scale and exponentiate
  costMatrix = - costMatrix / epsilon;
  costMatrix = costMatrix.array().exp();

  return costMatrix;
}

// Define the Sinkhorn algorithm function
List sinkhorn_stable_cpp(Eigen::VectorXd a,
                      Eigen::VectorXd b,
                      Eigen::MatrixXd costMatrix,
                      int numIterations,
                      double epsilon,
                      double maxErr,
                      double tau,
                      double trunc_thresh) {
  const int numRows = costMatrix.rows();
  const int numCols = costMatrix.cols();

  // Initialize the overloaded potentials to zeros
  Eigen::VectorXd alpha = Eigen::VectorXd::Zero(numRows);
  Eigen::VectorXd beta = Eigen::VectorXd::Zero(numCols);

  // Initialize the K matrix
  Eigen::MatrixXd K = getK(costMatrix, alpha, beta, epsilon);

  // Initialize all the extra dual variables to vectors of 1
  Eigen::VectorXd u_tilde = Eigen::VectorXd::Ones(numRows);
  Eigen::VectorXd v_tilde = Eigen::VectorXd::Ones(numCols);

  // Initialize the marginal weights and the other useful variables
  Eigen::VectorXd estimated_marginal;

  // Initialize algorithms values
  int iter = 1;
  double err = std::numeric_limits<double>::infinity();

  // External loop for convergence
  while (iter == 1 || (iter <= numIterations &&
         err > maxErr && err < std::numeric_limits<double>::infinity())) {
    // First step: updates
    // Compute v updates
    v_tilde = K.transpose() * u_tilde;
    v_tilde = b.array() * v_tilde.array().inverse();
    // Rcout << "v " << v_tilde << std::endl;

    // Compute u updates
    u_tilde = K * v_tilde;
    u_tilde = a.array() * u_tilde.array().inverse();
    // Rcout << "u " << u_tilde << std::endl;

    // Second step: absorption iteration
    if (u_tilde.cwiseAbs().maxCoeff() > tau ||
        v_tilde.cwiseAbs().maxCoeff() > tau ||
        iter == numIterations) {
      // Absorb the big values
      alpha = alpha.array() + epsilon * u_tilde.array().log();
      // Rcout << "alpha " << alpha << std::endl;
      beta = beta.array() + epsilon * v_tilde.array().log();
      // Rcout << "beta " << beta << std::endl;

      // Reset the potentials
      u_tilde = Eigen::VectorXd::Ones(numRows);
      v_tilde = Eigen::VectorXd::Ones(numCols);

      // Update the cost matrix
      K = getK(costMatrix, alpha, beta, epsilon);
      // Rcout << "K " << K << std::endl;
    }

    // Error control
    if (iter % 10 == 0 || iter == 1) {
      estimated_marginal = v_tilde.array() * (K.transpose() * u_tilde).array();
      // Rcout << "marginal " << estimated_marginal << std::endl;
      err = (estimated_marginal - b).cwiseAbs().sum();
    }

    // Update iteration
    iter++;
  }

  if (iter >= numIterations) {
    Rcerr << "Increase number of iterations" << std::endl;
    return List::create(Named("Error") = "Increase number of iterations");
  }

  // // Potentials (dual variables)
  // Eigen::VectorXd f(epsilon * u);
  // Eigen::VectorXd g(epsilon * v);
  //
  // // Wasserstein dual
  // double W22 = f.dot(a) + g.dot(b);

  // Optimal coupling
  Eigen::MatrixXd P((K.array().colwise() * u_tilde.array()).rowwise() * v_tilde.array().transpose());
  // Rcout << "P " << P.maxCoeff() << std::endl;

  // Wasserstein distance
  double W22_prime = (P.transpose() * costMatrix).trace();

  // Return u and v as a List
  return List::create(
    Named("iter") = iter,
    // Named("f") = f,
    // Named("g") = g,
    // Named("P") = P,
    Named("cost") = W22_prime
    // Named("cost_dual") = W22
  );
}

// Expose the Sinkhorn function to R
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List sinkhorn_stable(Eigen::VectorXd a,
                  Eigen::VectorXd b,
                  Eigen::MatrixXd costm,
                  int numIterations,
                  double epsilon,
                  double maxErr,
                  double tau) {
  return sinkhorn_stable_cpp(a, b, costm, numIterations, epsilon, maxErr, tau);
}

// # /***R
// # set.seed(1)
// # n <- 100
// # m <- 100
// # a <- rep(1 / n, n)
// # b <- rep(1 / m, m)
// # C <- as.matrix(dist(rnorm(100)))#[, 1:50]
// # sinkhorn_log(a, b, C, 1e5, 0.0000001, 1e-3, 1e5)
// # */
