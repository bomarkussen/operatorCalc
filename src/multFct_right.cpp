#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
SEXP multFct_right(NumericVector f, SEXP g, SEXP y) {
  // Multiply a matFct by a matrix from the right
  // No sanity check performed!
  NumericVector vec_g(g);
  NumericVector vec_y(y);
  IntegerVector dims1 = vec_g.attr("dim");
  IntegerVector dims2 = vec_y.attr("dim");

  int p = dims1[2];  // NB: 1 less in order to simplify code
  int r = dims1[3];
  int rp = r*p;
  
  arma::cube xf(f.begin(), dims1[0], dims1[1], p+1, false);
  arma::cube xg(vec_g.begin(), dims1[0], dims1[1], rp, false);
  arma::mat  yy(vec_y.begin(), dims2[0], dims2[1], false);
  arma::cube ff(dims1[0], dims2[1], p+1);
  arma::cube gg(dims1[0], dims2[1], rp);
  
  // find f and g functions
  for (int j=0; j<=p; j++) ff.slice(j) = xf.slice(j)*yy;
  for (int j=0; j<rp; j++) gg.slice(j) = xg.slice(j)*yy;
            
  // return to R
  return(List::create(Rcpp::Named("f")=ff, Rcpp::Named("g")=gg));
}
