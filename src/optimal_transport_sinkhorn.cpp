#include <RcppEigen.h>
using namespace Rcpp;

// Define the Sinkhorn algorithm function
// The function is similar to Python Optimal Transport implementation
List sinkhorn_cpp(Eigen::VectorXd a,
                  Eigen::VectorXd b,
                  Eigen::MatrixXd costMatrix,
                  int numIterations,
                  double epsilon,
                  double maxErr) {
  int numRows = costMatrix.rows();
  int numCols = costMatrix.cols();

  // Initialize the K matrix
  Eigen::MatrixXd K(- costMatrix / epsilon);
  K = K.array().exp();

  // Rcout << "K " << K << std::endl;
  // Rcout << "a " << a << std::endl;
  // Rcout << "b " << b << std::endl;

  // Initialize all the dual variables to vectors of 1
  Eigen::VectorXd u(numRows);
  Eigen::VectorXd v(numCols);

  u.setOnes();

  // Initialize the marginal weights and the other useful variables
  Eigen::VectorXd estimated_marginal;

  // Initialize algorithms values
  int iter = 1;
  double err = std::numeric_limits<double>::infinity();

  while (iter == 1 || (iter <= numIterations &&
         err > maxErr && err < std::numeric_limits<double>::infinity())) {
    // Compute v updates
    v = K.transpose() * u;
    v = b.array() * v.array().inverse();

    // Rcout << "v " << v << std::endl;

    // Compute u updates
    u = K * v;
    u = a.array() * u.array().inverse();

    // Rcout << "u " << u << std::endl;

    // Update error
    if (iter % 10 == 0 || iter == 1) {
      estimated_marginal = v.array() * (K.transpose() * u).array();
      err = (estimated_marginal - b).cwiseAbs().sum();
    }

    // Update iteration
    iter++;
  }

  if (iter >= numIterations) {
    Rcerr << "Increase number of iterations" << std::endl;
    return List::create(Named("Error") = "Increase number of iterations");
  }

  // Potentials (dual variables)
  Eigen::VectorXd f(epsilon * (numRows * u).array().log());
  Eigen::VectorXd g(epsilon * (numCols * v).array().log());

  // Wasserstein dual
  double W22 = f.mean() + g.mean();

  // Optimal coupling
  Eigen::MatrixXd U(u.asDiagonal());
  Eigen::MatrixXd V(v.asDiagonal());
  Eigen::MatrixXd P(U * K * V);

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
List sinkhorn(Eigen::VectorXd a,
              Eigen::VectorXd b,
              Eigen::MatrixXd costm,
              int numIterations,
              double epsilon,
              double maxErr) {
  return sinkhorn_cpp(a, b, costm, numIterations, epsilon, maxErr);
}
