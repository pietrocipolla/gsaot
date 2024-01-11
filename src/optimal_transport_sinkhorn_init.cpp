#include <RcppEigen.h>
using namespace Rcpp;

// Define the Sinkhorn algorithm function with marginal inputs
List sinkhorn_cpp(Eigen::MatrixXd costMatrix,
                       int numIterations,
                       double epsilon,
                       Eigen::VectorXd u,
                       Eigen::VectorXd v,
                       double maxErr) {
  int numRows = costMatrix.rows();
  int numCols = costMatrix.cols();

  // If eps is negative, set to relative value
  if (epsilon < 0) epsilon = -epsilon * costMatrix.maxCoeff();

  // Initialize the K matrix
  Eigen::MatrixXd K(- costMatrix / epsilon);
  K = K.array().exp();

  // Initialize all the dual variables to vectors of 1
  // Eigen::VectorXd u(numRows);
  // Eigen::VectorXd v(numCols);
  //
  // u.setOnes();
  // v.setOnes();

  if (u.size() != numRows || v.size() != numCols) {
    Rcerr << "Wrong initialization" << std::endl;
    return List::create(Named("Error") = "Wrong initialization");
  }


  // Initialize the marginal weights and the other useful variables
  Eigen::VectorXd Wm(numCols);
  Wm.setOnes();
  Wm = Wm / numCols;

  // Rcout << "Empirical marginal " << Wm << std::endl;
  // Rcout << "Initialized u " << u << std::endl;
  // Rcout << "Initialized v " << v << std::endl;

  Eigen::VectorXd Estimated_marginal;

  // Initialize algorithms values
  int iter = 1;
  double err = std::numeric_limits<double>::infinity();
  std::vector<double> errors;

  while (iter == 1 || (iter <= numIterations &&
         err > maxErr && err < std::numeric_limits<double>::infinity())) {
    // Compute v updates
    v = K.transpose() * u;
    v = v.cwiseInverse() / numCols;

    //Rcout << "v: " << v << std::endl;

    // Compute u updates
    u = K * v;
    u = u.cwiseInverse() / numRows;

    //Rcout << "u: " << u << std::endl;

    // Update error
    Estimated_marginal = v.array() * (K.transpose() * u).array();
    err = (Estimated_marginal - Wm).cwiseAbs().sum();

    errors.push_back(err);

    // Rcout << "Error at iteration " << iter << ": " << err << std::endl;

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
    Named("u") = u,
    Named("f") = f,
    Named("g") = g,
    Named("P") = P,
    Named("W22") = W22_prime,
    Named("W22_dual") = W22,
    Named("Errors") = errors
  );
}

// Expose the Sinkhorn function to R
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List optimal_transport_sinkhorn_init(Eigen::MatrixXd costMatrix,
                                     int numIterations,
                                     double epsilon,
                                     Eigen::VectorXd u,
                                     Eigen::VectorXd v,
                                     double maxErr = 1e-9) {
  return sinkhorn_cpp(costMatrix, numIterations, epsilon, u, v, maxErr);
}