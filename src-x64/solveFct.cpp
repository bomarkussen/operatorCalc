#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
SEXP solveFct(NumericVector f1, SEXP g1, NumericVector f2, SEXP g2, double order, double super, LogicalVector continuous) {
  // Assume that dimensions match
  // No sanity check performed!
  
  NumericVector vec_g1(g1);
  NumericVector vec_g2(g2);
  IntegerVector dims1 = vec_g1.attr("dim");
  IntegerVector dims2 = vec_g2.attr("dim");

  int p = dims1[2];      // NB: 1 less in order to simplify code
  int r1 = dims1[3];
  int r2 = dims2[3];
  
  arma::cube fa(f1.begin(), dims1[0], dims1[1], p+1, false);
  arma::cube fb(f2.begin(), dims2[0], dims2[1], p+1, false);
  arma::cube fx(dims1[0], dims2[1], p+1);
  arma::cube ga(vec_g1.begin(), dims1[0], dims1[1], r1*p, false);
  arma::cube gb(vec_g2.begin(), dims2[0], dims2[1], r2*p, false);
  arma::mat  gx(std::max(1.0,order)*p, dims1[0]*dims2[1]);
  gx.fill(0.0);
  arma::mat  tmp1(dims1[0], dims1[1]);
  arma::mat  tmp2(dims2[0], dims2[1]);
  
  // find inverse at knot points
  for (int j=0; j<p; j++) fx.slice(j) = arma::solve(fa.slice(j),fb.slice(j));
  tmp1 = fa.slice(p);
  tmp2 = fb.slice(p);
  for (int k=0; k<r1; k++) tmp1 += ga.slice(p-1+k*p);
  for (int k=0; k<r2; k++) tmp2 += gb.slice(p-1+k*p);
  fx.slice(p) = arma::solve(tmp1,tmp2);
  
  // Assume order>1 and super>0:
  // make interpolation of inverse between knot points
  // there are internal intervals points and an OLS to be solved
  // set-up matrices
  int q_cont = std::min(order-1,super);
  int q_discont = std::min(order,super);
  arma::vec s(super);
  arma::mat A_cont(super,q_cont);
  arma::mat A_discont(super,q_discont);
  arma::mat B(super,dims1[0]*dims2[1]);
  // precomputations
  for (int k=0; k<super; k++) s.row(k) = (k+1)/(super+1);
  for (int k=0; k<super; k++) for (int l=0; l<q_cont; l++) A_cont(k,l) = pow(s(k),l+2)-s(k);
  for (int k=0; k<super; k++) for (int l=0; l<q_discont; l++) A_discont(k,l) = pow(s(k),l+1);
  // loop over sampling intervals
  for (int j=0; j<p; j++) {
    // make interval data by removing linear interpolation
    for (int k=0; k<super; k++) {
      tmp1 = (1-s(k))*fa.slice(j) + s(k)*fa.slice(j+1);
      tmp2 = (1-s(k))*fb.slice(j) + s(k)*fb.slice(j+1);
      for (int l=0; l<r1; l++) tmp1 += pow(s(k),l+1)*ga.slice(j+l*p);
      for (int l=0; l<r2; l++) tmp2 += pow(s(k),l+1)*gb.slice(j+l*p);
      B.row(k) = arma::trans(arma::vectorise(arma::solve(tmp1,tmp2)-(1-s(k))*fx.slice(j)-s(k)*fx.slice(j+1)));
    }
    // make design matrix
    if (continuous(j)) {
      // enforce continuity in this interval
      gx.rows(j*order+1,j*order+q_cont) = arma::solve(A_cont,B);
      for (int k=1; k<=q_cont; k++) gx.row(j*order) -= gx.row(j*order+k);
    } else {
      // allow for discontinuity in this interval
      gx.rows(j*order,j*order+q_discont-1) = arma::solve(A_discont,B);
    }
  }
  
  // return to R
  return(List::create(Rcpp::Named("f")=fx, Rcpp::Named("g")=arma::trans(gx)));
}
