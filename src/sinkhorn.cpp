#include <RcppEigen.h>
using namespace Rcpp;

// Define the Sinkhorn algorithm function
// The function is similar to Python Optimal Transport implementation
List sinkhorn_cpp(Eigen::VectorXd a,
                  Eigen::VectorXd b,
                  const Eigen::MatrixXd costMatrix,
                  int numIterations,
                  double epsilon,
                  double maxErr) {
  int numRows = costMatrix.rows();
  int numCols = costMatrix.cols();

  // Initialize the K matrix
  Eigen::MatrixXd K(-costMatrix / epsilon);
  K = K.array().exp();

  // Initialize the Sinkhorn scaling vectors
  //
  // u = exp(f / epsilon)
  // v = exp(g / epsilon)
  Eigen::VectorXd u(numRows);
  Eigen::VectorXd v(numCols);

  u.setOnes();

  // Initialize the marginal weights and other useful variables
  Eigen::VectorXd estimated_marginal;

  // Initialize algorithm values
  int iter = 1;
  double err = std::numeric_limits<double>::infinity();

  while (iter == 1 ||
         (iter <= numIterations &&
         err > maxErr &&
         err < std::numeric_limits<double>::infinity())) {

    // Update v
    v = K.transpose() * u;
    v = b.array() * v.array().inverse();

    // Update u
    u = K * v;
    u = a.array() * u.array().inverse();

    // Update error
    if (iter % 10 == 0 || iter == 1) {
      estimated_marginal =
        v.array() * (K.transpose() * u).array();

      err = (estimated_marginal - b).cwiseAbs().sum();
    }

    // Update iteration
    iter++;
  }

  if (iter >= numIterations) {
    Rcerr << "Increase number of iterations" << std::endl;

    return List::create(
      Named("Error") = "Increase number of iterations"
    );
  }

  // Recover the dual potentials:
  //
  // f = epsilon * log(u)
  // g = epsilon * log(v)
  Eigen::VectorXd f = epsilon * u.array().log().matrix();
  Eigen::VectorXd g = epsilon * v.array().log().matrix();

  // Dual objective:
  //
  // <f, a> + <g, b>
  // - epsilon <exp(f / epsilon), K exp(g / epsilon)>
  //
  // Since exp(f / epsilon) = u and exp(g / epsilon) = v:
  //
  // <f, a> + <g, b> - epsilon u^T K v
  double entropicCost =
    epsilon * (
        a.dot((u.array() / a.array()).log().matrix()) +
          b.dot((v.array() / b.array()).log().matrix()) -
          u.dot(K * v) +
          1.0
    );

  return List::create(
    Named("iter") = iter,
    Named("cost") = entropicCost
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

// # /***R
// # set.seed(1)
// # n <- 100
// # m <- 100
// # a <- rep(1 / n, n)
// # b <- rep(1 / m, m)
// # C <- as.matrix(dist(rnorm(100)))#[, 1:50]
// # ret <- sinkhorn(a, b, C, 1e3, 0.1, 1e-3)
// # */
