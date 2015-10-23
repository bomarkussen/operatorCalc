#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
SEXP multFct_left(SEXP x, NumericVector f, SEXP g) {
  // Multiply a matFct by a matrix from the left
  // No sanity check performed!
  NumericVector vec_x(x);
  NumericVector vec_g(g);
  IntegerVector dims1 = vec_x.attr("dim");
  IntegerVector dims2 = vec_g.attr("dim");

  int p = dims2[2];  // NB: 1 less in order to simplify code
  int r = dims2[3];
  int rp = r*p;
  
  arma::mat xx(vec_x.begin(), dims1[0], dims1[1], false);
  arma::cube yf(f.begin(), dims2[0], dims2[1], p+1, false);
  arma::cube yg(vec_g.begin(), dims2[0], dims2[1], rp, false);
  arma::cube ff(dims1[0], dims2[1], p+1);
  arma::cube gg(dims1[0], dims2[1], rp);
  
  // find f and g functions
  for (int j=0; j<=p; j++) ff.slice(j) = xx*yf.slice(j);
  for (int j=0; j<rp; j++) gg.slice(j) = xx*yg.slice(j);
            
  // return to R
  return(List::create(Rcpp::Named("f")=ff, Rcpp::Named("g")=gg));
}
