#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
SEXP approximate(NumericVector t_obs, SEXP matData, IntegerVector knotIndex, NumericVector knotDiff, LogicalVector continuous, double r) {
  // Assume that dimensions match
  // No sanity check performed!

  NumericVector vec_data(matData);
  IntegerVector dims = vec_data.attr("dim");

  arma::vec tt(t_obs.begin(), t_obs.size(), false);
  arma::vec knotSpan(knotDiff.begin(), knotDiff.size(), false);
  
  int p = knotIndex.size()-1;    // NB: 1 less in order to simplify code
  int max_obs = knotSpan.max();
  
  arma::mat data(vec_data.begin(), dims[0], dims[1]*dims[2], false);
  arma::vec s(max_obs);
  arma::mat A(max_obs,r);
  arma::mat B(max_obs,dims[1]*dims[2]);
  arma::mat g(r*p, dims[1]*dims[2]);
  g.fill(0.0);
  
  // loop over sampling intervals
  for (int j=0; j<p; j++) {
    if (knotSpan(j)>0) {
      // there internal intervals points and an OLS to be solved
      // make interval data by removing linear interpolation
      s.rows(0,knotSpan(j)-1) = (tt.rows(knotIndex(j)+1,knotIndex(j+1)-1)-tt(knotIndex(j))) / (tt(knotIndex(j+1))-tt(knotIndex(j)));
      for (int k=0; k<knotSpan(j); k++) B.row(k) = data.row(knotIndex(j)+1+k)-(1-s(k))*data.row(knotIndex(j))-s(k)*data.row(knotIndex(j+1));
      // make design matrix
      if (continuous(j)) {
        // enforce continuity in this interval
        int q = std::min(r-1,knotSpan(j));
        for (int k=0; k<knotSpan(j); k++) for (int l=0; l<q; l++) A(k,l) = pow(s(k),l+2)-s(k);
        g.rows(j*r+1,j*r+q) = arma::solve(A.submat(0,0,knotSpan(j)-1,q-1),B.rows(0,knotSpan(j)-1));
        for (int k=1; k<r; k++) g.row(j*r) -= g.row(j*r+k);
      } else {
        // allow for discontinuity in this interval
        int q = std::min(r,knotSpan(j));
        for (int k=0; k<knotSpan(j); k++) for (int l=0; l<q; l++) A(k,l) = pow(s(k),l+1);
        g.rows(j*r,j*r+q-1) = arma::solve(A.submat(0,0,knotSpan(j)-1,q-1),B.rows(0,knotSpan(j)-1));
      }
    }
  }
  
  return wrap(arma::trans(g));
}
