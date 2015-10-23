#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
SEXP triangular(NumericVector diffs, NumericVector f, SEXP g) {
  // Assume that dimensions match
  // No sanity check performed!
  NumericVector vec_g(g);
  IntegerVector dims = vec_g.attr("dim");

  int p = dims[2];     // NB: 1 less in order to simplify code
  int r = dims[3];

  arma::cube ff(f.begin(), dims[0], dims[1], p+1, false);
  arma::cube gg(vec_g.begin(), dims[0], dims[1], r*p, false);
  arma::cube integrals(dims[0], dims[1], p+1);
  integrals.fill(0.0);
  
  // find upwind integrals
  for (int j=0; j<p; j++) {
    integrals.slice(j+1) = diffs(j)/6*ff.slice(j)+diffs(j)/3*ff.slice(j+1);
    for (int k=0; k<r; k++) integrals.slice(j+1) += diffs(j)/(k+3)*gg.slice(j+k*p);
  }

  // add downwind integrals
  for (int j=0; j<p; j++) {
    integrals.slice(j) += diffs(j)/3*ff.slice(j)+diffs(j)/6*ff.slice(j+1);
    for (int k=0; k<r; k++) integrals.slice(j) += diffs(j)/(k+2)/(k+3)*gg.slice(j+k*p);
  }

  // return to R
  return wrap(integrals);
}
