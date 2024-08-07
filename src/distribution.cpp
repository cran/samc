// Copyright (c) 2024 Andrew Marx where applicable.
// Licensed under AGPLv3.0. See LICENSE file in the project root for details.

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::export(".qpow_row")]]
Rcpp::List qpow_row(Eigen::Map<Eigen::SparseMatrix<double> > &M, const Eigen::Map<Eigen::VectorXd> &vec, Rcpp::NumericVector steps)
{
  int n = steps.size();

  Rcpp::List res = Rcpp::List::create();

  Eigen::RowVectorXd time_res = vec;

  for(int i = 1; i < n; i++) {
    for (int j = steps[i - 1]; j < steps[i]; j++) {
      if(i % 1000 == 0) Rcpp::checkUserInterrupt();
      time_res = time_res * M;
    }

    res.push_back(time_res, std::to_string((int)steps[i]));
  }

  return res;
}

// [[Rcpp::export(".qpow_col")]]
Rcpp::List qpow_col(Eigen::Map<Eigen::SparseMatrix< double> > &M, const Eigen::Map<Eigen::VectorXd> &vec, Rcpp::NumericVector steps)
{
  int n = steps.size();

  Rcpp::List res = Rcpp::List::create();

  Eigen::VectorXd time_res = M * vec;

  for(int i = 1; i < n; i++) {
    for (int j = steps[i - 1]; j < steps[i]; j++) {
      if(i % 1000 == 0) Rcpp::checkUserInterrupt();
      time_res = M * time_res;
    }

    res.push_back(time_res, std::to_string((int)steps[i]));
  }

  return res;
}
