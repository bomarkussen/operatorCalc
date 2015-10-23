#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
SEXP integral(NumericVector diffs, NumericVector f, SEXP g) {
  // Assume that dimensions match
  // No sanity check performed!
  NumericVector vec_g(g);
  IntegerVector dims = vec_g.attr("dim");

  int p = dims[2];     // NB: 1 less in order to simplify code
  int r = dims[3];
  
  arma::cube ff(f.begin(), dims[0], dims[1], p+1, false);
  arma::cube gg(vec_g.begin(), dims[0], dims[1], r*p, false);
  arma::mat  integral(dims[0], dims[1]);

  // find integral
  integral.fill(0.0); 
  for (int j=0; j<p; j++) {
    integral += diffs(j)/2*(ff.slice(j)+ff.slice(j+1));
    for (int k=0; k<r; k++) integral += diffs(j)/(k+2)*gg.slice(j+k*p);
  }

  // return to R
  return wrap(integral);
}
