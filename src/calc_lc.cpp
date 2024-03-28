// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::Map;
using Eigen::VectorXd;

// [[Rcpp::export]]
NumericMatrix calc_lc_cpp1(NumericVector A, NumericVector lambda, NumericVector w, NumericVector x) {
    int p = A.size();
    int n = x.size();
    NumericMatrix lc(p, n);
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < n; j++) {
            lc(i, j) = abs(A[i] * (lambda[i] / (pow(lambda[i], 2) + pow((x[j] - w[i]), 2))));
        }
    }
    return lc;
}

// [[Rcpp::export]]
NumericMatrix calc_lc_cpp2(NumericVector A, NumericVector lambda, NumericVector w, NumericVector x) {
    int p = A.size();
    int n = x.size();
    NumericMatrix lc(p, n);
    for (int i = 0; i < p; i++) {
        if (A[i] != 0) {
            lc(i, _) = abs(A[i] * (lambda[i] / (pow(lambda[i], 2) + pow((x - w[i]), 2))));;
        }
    }
    return lc;
}

// [[Rcpp::export]]
NumericMatrix calc_lc_cpp3(NumericVector A_, NumericVector lambda_, NumericVector w_, NumericVector x_) {
    Map<VectorXd> A(as<Map<VectorXd> >(A_));
    Map<VectorXd> lambda(as<Map<VectorXd> >(lambda_));
    Map<VectorXd> w(as<Map<VectorXd> >(w_));
    Map<VectorXd> x(as<Map<VectorXd> >(x_));
    
    int p = A.size();
    int n = x.size();
    NumericMatrix lc(p, n);
    
    for (int i = 0; i < p; i++) {
        if (A[i] != 0) {
            VectorXd temp = (A[i] * (lambda[i] / ((lambda[i] * lambda[i]) + ((x.array() - w[i]).array().square())))).array().abs();
            lc.row(i) = as<NumericVector>(wrap(temp));
        }
    }
    return lc;
}
