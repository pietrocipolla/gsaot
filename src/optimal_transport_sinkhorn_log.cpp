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

// Define the Sinkhorn algorithm function
List sinkhorn_log_cpp(Eigen::VectorXd a,
                      Eigen::VectorXd b,
                      Eigen::MatrixXd costMatrix,
                      int numIterations,
                      double epsilon,
                      double maxErr) {
  const int numRows = costMatrix.rows();
  const int numCols = costMatrix.cols();

  // Initialize the K matrix
  Eigen::MatrixXd K(-costMatrix / epsilon);

  // Initialize all the dual variables to vectors of 1
  Eigen::VectorXd u(numRows);
  Eigen::VectorXd v(numCols);

  u.setZero();
  v.setZero();

  // Build potential matrices
  Eigen::MatrixXd U(u.array().exp().matrix().asDiagonal());
  Eigen::MatrixXd V(v.array().exp().matrix().asDiagonal());

  // Initialize the marginal weights and the other useful variables
  Eigen::VectorXd estimated_marginal;
  Eigen::VectorXd loga = a.array().log();
  Eigen::VectorXd logb = b.array().log();

  // Initialize algorithms values
  int iter = 1;
  double err = std::numeric_limits<double>::infinity();

  while (iter == 1 || (iter <= numIterations &&
         err > maxErr && err < std::numeric_limits<double>::infinity())) {
    // Compute v updates
    v = logb - logsumexp(K.colwise() + u, 1);
    // Rcout << "v " << v << std::endl;

    // Compute u updates
    u = loga - logsumexp(K.transpose().colwise() + v, 1);
    // Rcout << "u " << u << std::endl;

    if (iter % 10 == 0 || iter == 1) {
      // Update potential matrices
      U = u.array().exp().matrix().asDiagonal();
      V = v.array().exp().matrix().asDiagonal();
      // Rcout << "U " << U << std::endl;
      // Rcout << "V " << V << std::endl;

      // Update error
      estimated_marginal = (V * K.transpose().array().exp().matrix() * U).array().rowwise().sum();
      err = (estimated_marginal - b).cwiseAbs().sum();
      // Rcout << "err " << err << std::endl;
    }

    // Update iteration
    iter++;
  }

  if (iter >= numIterations) {
    Rcerr << "Increase number of iterations" << std::endl;
    return List::create(Named("Error") = "Increase number of iterations");
  }

  // Potentials (dual variables)
  Eigen::VectorXd f(epsilon * u);
  Eigen::VectorXd g(epsilon * v);

  // Wasserstein dual
  double W22 = f.dot(a) + g.dot(b);

  // Optimal coupling
  Eigen::MatrixXd P(U * K.array().exp().matrix() * V);

  // Wasserstein distance
  double W22_prime = (P.transpose() * costMatrix).trace();

  // Return u and v as a List
  return List::create(
    Named("iter") = iter,
    Named("f") = f,
    Named("g") = g,
    Named("P") = P,
    Named("cost") = W22_prime,
    Named("cost_dual") = W22
  );
}

// Expose the Sinkhorn function to R
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List sinkhorn_log(Eigen::VectorXd a,
                  Eigen::VectorXd b,
                  Eigen::MatrixXd costm,
                  int numIterations,
                  double epsilon,
                  double maxErr) {
  return sinkhorn_log_cpp(a, b, costm, numIterations, epsilon, maxErr);
}

// # /***R
// # n <- 100
// # m <- 50
// # a <- rep(1 / n, n)
// # b <- rep(1 / m, m)
// # C <- as.matrix(dist(rnorm(100)))[, 1:50]
// # sinkhorn_log(a, b, C, 1e5, 0.1, 1e-3)
// # */
