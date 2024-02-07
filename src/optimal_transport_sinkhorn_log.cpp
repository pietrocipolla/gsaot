#include <RcppEigen.h>
using namespace Rcpp;

// Define the logsumexp for Eigen matrices
Eigen::VectorXd logsumexp(const Eigen::MatrixXd x,
                          int axis) {
  if (axis == 0) {
    // Compute log-sum-exp along rows
    Eigen::VectorXd maxCoeffs = x.rowwise().maxCoeff();
    return maxCoeffs.array() + (x.colwise() - maxCoeffs).array().exp().rowwise().sum().log();
  } else if (axis == 1) {
    // Compute log-sum-exp along columns
    Eigen::VectorXd maxCoeffs = x.colwise().maxCoeff();
    return maxCoeffs.transpose().array() + (x.rowwise() - maxCoeffs.transpose()).array().exp().colwise().sum().log();
  } else {
    // Invalid axis value
    return Eigen::VectorXd::Zero(x.rows());
  }

  return Eigen::VectorXd::Zero(x.rows());
}

// Define the Sinkhorn algorithm function
List sinkhorn_log_cpp(Eigen::MatrixXd costMatrix,
                      int numIterations,
                      double epsilon,
                      double maxErr) {
  int numRows = costMatrix.rows();
  int numCols = costMatrix.cols();

  // If eps is negative, set to relative value
  if (epsilon < 0) epsilon = - epsilon * costMatrix.maxCoeff();

  // Initialize the K matrix
  Eigen::MatrixXd K(- costMatrix / epsilon);
  //std::cout << "K:\n" << K << std::endl;

  // Initialize all the dual variables to vectors of 1
  Eigen::VectorXd u(numRows);
  Eigen::VectorXd v(numCols);

  u.setZero();
  v.setZero();

  //std::cout << "u:\n" << u << std::endl;
  //std::cout << "v:\n" << v << std::endl;

  // Build potential matrices
  Eigen::MatrixXd U(u.array().exp().matrix().asDiagonal());
  Eigen::MatrixXd V(v.array().exp().matrix().asDiagonal());

  // Initialize the marginal weights and the other useful variables
  Eigen::VectorXd a(numCols);
  a.setOnes();
  Eigen::VectorXd Wm(a / numCols);

  //std::cout << "a:\n" << Wm << std::endl;

  Eigen::VectorXd Estimated_marginal;

  // Initialize algorithms values
  int iter = 1;
  double err = std::numeric_limits<double>::infinity();

  while (iter == 1 || (iter <= numIterations &&
         err > maxErr && err < std::numeric_limits<double>::infinity())) {
    //std::cout << "iter: " << iter << std::endl;

    // Compute v updates
    //std::cout << "Inner term v: " << K.transpose().colwise() + u << std::endl;
    v = - logsumexp( K.colwise() + u, 1).array() - log(numCols);

    //std::cout << "v: " << v.array().exp() << std::endl;

    // Compute u updates
    //std::cout << "Inner term u: " << K.colwise() + v << std::endl;
    u = - logsumexp(K.transpose().colwise() + v, 1).array() - log(numRows);

    //std::cout << "u: " << u.array().exp() << std::endl;

    // Update potential matrices
    U = u.array().exp().matrix().asDiagonal();
    V = v.array().exp().matrix().asDiagonal();

    // Update error
    Estimated_marginal = V * K.transpose().array().exp().matrix() * U * a;
    //std::cout << "Estimated_marginal: " << Estimated_marginal << std::endl;
    err = (Estimated_marginal - Wm).cwiseAbs().sum();

    // Update iteration
    iter++;
  }

  if (iter >= numIterations) {
    Rcerr << "Increase number of iterations" << std::endl;
    return List::create(Named("Error") = "Increase number of iterations");
  }

  // Potentials (dual variables)
  Eigen::VectorXd f(epsilon * (numRows * u.array().exp()).log());
  Eigen::VectorXd g(epsilon * (numCols * v.array().exp()).log());

  // Wasserstein dual
  double W22 = f.mean() + g.mean();

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
    Named("W22") = W22_prime,
    Named("W22_dual") = W22
  );
}

// Expose the Sinkhorn function to R
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List optimal_transport_sinkhorn_log(Eigen::MatrixXd costMatrix,
                  int numIterations,
                  double epsilon,
                  double maxErr = 1e-9) {
  return sinkhorn_log_cpp(costMatrix, numIterations, epsilon, maxErr);
}
