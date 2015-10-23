#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
SEXP inverse(NumericVector f, SEXP g, double order, double super, LogicalVector continuous) {
  // Assume that dimensions match
  // No sanity check performed!
  
  NumericVector vec_g(g);
  IntegerVector dims = vec_g.attr("dim");

  int p = dims[2];      // NB: 1 less in order to simplify code
  int r = dims[3];
  
  arma::cube ff(f.begin(), dims[0], dims[1], p+1, false);
  arma::cube gg(vec_g.begin(), dims[0], dims[1], r*p, false);
  arma::cube inv_f(dims[0], dims[1], p+1);
  arma::mat  inv_g(std::max(1.0,order)*p, dims[0]*dims[1]); inv_g.fill(0.0);
  arma::mat  tmp(dims[0], dims[1]);
  
  // find inverse at knot points
  for (int j=0; j<p; j++) inv_f.slice(j) = arma::inv(ff.slice(j));
  tmp = ff.slice(p);
  for (int k=0; k<r; k++) tmp += gg.slice(p-1+k*p);
  inv_f.slice(p) = arma::inv(tmp);
  
  // if order==0 make piecewise constant function
  if (order==0) {
    for (int j=0; j<p; j++) inv_g.row(j) = arma::trans(arma::vectorise(inv_f.slice(j)-inv_f.slice(j+1)));
  }
  
  // make interpolation of inverse between knot points
  if ((order>1) && (super>0)) {
    // there are internal intervals points and an OLS to be solved
    // set-up matrices
    int q_cont = std::min(order-1,super);
    int q_discont = std::min(order,super);
    arma::vec s(super);
    arma::mat A_cont(super,q_cont);
    arma::mat A_discont(super,q_discont);
    arma::mat B(super,dims[0]*dims[1]);

    // precomputations
    for (int k=0; k<super; k++) s.row(k) = (k+1)/(super+1);
    for (int k=0; k<super; k++) for (int l=0; l<q_cont; l++) A_cont(k,l) = pow(s(k),l+2)-s(k);
    for (int k=0; k<super; k++) for (int l=0; l<q_discont; l++) A_discont(k,l) = pow(s(k),l+1);

    // loop over sampling intervals
    for (int j=0; j<p; j++) {
      // make interval data by removing linear interpolation
      for (int k=0; k<super; k++) {
        tmp = (1-s(k))*ff.slice(j) + s(k)*ff.slice(j+1);
        for (int l=0; l<r; l++) tmp += pow(s(k),l+1)*gg.slice(j+l*p);
        B.row(k) = arma::trans(arma::vectorise(arma::inv(tmp)-(1-s(k))*inv_f.slice(j)-s(k)*inv_f.slice(j+1)));
      }
      // make design matrix
      if (continuous(j)) {
        // enforce continuity in this interval
        inv_g.rows(j*order+1,j*order+q_cont) = arma::solve(A_cont,B);
        for (int k=1; k<=q_cont; k++) inv_g.row(j*order) -= inv_g.row(j*order+k);
      } else {
        // allow for discontinuity in this interval
        inv_g.rows(j*order,j*order+q_discont-1) = arma::solve(A_discont,B);
      }
    }
  }
  
  // return to R
  return(List::create(Rcpp::Named("f")=inv_f, Rcpp::Named("g")=arma::trans(inv_g)));
}
