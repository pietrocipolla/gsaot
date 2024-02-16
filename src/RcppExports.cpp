// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// sinkhorn
List sinkhorn(Eigen::VectorXd a, Eigen::VectorXd b, Eigen::MatrixXd costm, int numIterations, double epsilon, double maxErr);
RcppExport SEXP _gsaot_sinkhorn(SEXP aSEXP, SEXP bSEXP, SEXP costmSEXP, SEXP numIterationsSEXP, SEXP epsilonSEXP, SEXP maxErrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type a(aSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type b(bSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type costm(costmSEXP);
    Rcpp::traits::input_parameter< int >::type numIterations(numIterationsSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< double >::type maxErr(maxErrSEXP);
    rcpp_result_gen = Rcpp::wrap(sinkhorn(a, b, costm, numIterations, epsilon, maxErr));
    return rcpp_result_gen;
END_RCPP
}
// optimal_transport_sinkhorn_hist
List optimal_transport_sinkhorn_hist(Eigen::VectorXd a, Eigen::VectorXd b, Eigen::MatrixXd costMatrix, int numIterations, double epsilon);
RcppExport SEXP _gsaot_optimal_transport_sinkhorn_hist(SEXP aSEXP, SEXP bSEXP, SEXP costMatrixSEXP, SEXP numIterationsSEXP, SEXP epsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type a(aSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type b(bSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type costMatrix(costMatrixSEXP);
    Rcpp::traits::input_parameter< int >::type numIterations(numIterationsSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(optimal_transport_sinkhorn_hist(a, b, costMatrix, numIterations, epsilon));
    return rcpp_result_gen;
END_RCPP
}
// optimal_transport_sinkhorn_init
List optimal_transport_sinkhorn_init(Eigen::MatrixXd costMatrix, int numIterations, double epsilon, Eigen::VectorXd u, Eigen::VectorXd v, double maxErr);
RcppExport SEXP _gsaot_optimal_transport_sinkhorn_init(SEXP costMatrixSEXP, SEXP numIterationsSEXP, SEXP epsilonSEXP, SEXP uSEXP, SEXP vSEXP, SEXP maxErrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type costMatrix(costMatrixSEXP);
    Rcpp::traits::input_parameter< int >::type numIterations(numIterationsSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type u(uSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type maxErr(maxErrSEXP);
    rcpp_result_gen = Rcpp::wrap(optimal_transport_sinkhorn_init(costMatrix, numIterations, epsilon, u, v, maxErr));
    return rcpp_result_gen;
END_RCPP
}
// sinkhorn_log
List sinkhorn_log(Eigen::MatrixXd costMatrix, int numIterations, double epsilon, double maxErr);
RcppExport SEXP _gsaot_sinkhorn_log(SEXP costMatrixSEXP, SEXP numIterationsSEXP, SEXP epsilonSEXP, SEXP maxErrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type costMatrix(costMatrixSEXP);
    Rcpp::traits::input_parameter< int >::type numIterations(numIterationsSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< double >::type maxErr(maxErrSEXP);
    rcpp_result_gen = Rcpp::wrap(sinkhorn_log(costMatrix, numIterations, epsilon, maxErr));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gsaot_sinkhorn", (DL_FUNC) &_gsaot_sinkhorn, 6},
    {"_gsaot_optimal_transport_sinkhorn_hist", (DL_FUNC) &_gsaot_optimal_transport_sinkhorn_hist, 5},
    {"_gsaot_optimal_transport_sinkhorn_init", (DL_FUNC) &_gsaot_optimal_transport_sinkhorn_init, 6},
    {"_gsaot_sinkhorn_log", (DL_FUNC) &_gsaot_sinkhorn_log, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_gsaot(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
