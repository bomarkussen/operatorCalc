#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
SEXP evaluate(IntegerVector jj, NumericVector s, NumericVector f, SEXP g) {
  // Assume that dimensions match
  // No sanity check performed!
  NumericVector vec_g(g);
  IntegerVector dims = vec_g.attr("dim");

  int p = dims[2];      // NB: 1 less in order to simplify code
  int r = dims[3];
  int pp = jj.size();
  
  arma::cube ff(f.begin(), dims[0], dims[1], p+1, false);
  arma::cube gg(vec_g.begin(), dims[0], dims[1], r*p, false);
  arma::cube res(dims[0], dims[1], pp);

  // find evaluation
  for (int j=0; j<pp; j++) {
    res.slice(j) = (1-s(j))*ff.slice(jj(j))+s(j)*ff.slice(jj(j)+1);
    for (int k=0; k<r; k++) res.slice(j) += pow(s(j),k+1)*gg.slice(jj(j)+k*p);
  }
  
  // return to R
  return wrap(res);
}
