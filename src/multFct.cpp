#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
SEXP multFct(NumericVector f1, SEXP g1, NumericVector f2, SEXP g2) {
  // Assume that two factors are sampled at same mesh, and matrix dimension congruent
  // No sanity check performed!
  NumericVector vec_g1(g1);
  NumericVector vec_g2(g2);
  IntegerVector dims1 = vec_g1.attr("dim");
  IntegerVector dims2 = vec_g2.attr("dim");

  int p = dims1[2];   // NB: 1 less in order to simplify code
  int r1 = dims1[3];
  int r2 = dims2[3];
  int r = std::max(1,r1)+std::max(1,r2);
  
  arma::cube xf(f1.begin(), dims1[0], dims1[1], p+1, false);
  arma::cube yf(f2.begin(), dims2[0], dims2[1], p+1, false);
  arma::cube xg(vec_g1.begin(), dims1[0], dims1[1], r1*p, false);
  arma::cube yg(vec_g2.begin(), dims2[0], dims2[1], r2*p, false);
  arma::cube f(dims1[0], dims2[1], p+1);
  arma::cube g(dims1[0], dims2[1], r*p);
  
  // find f function
  for (int j=0; j<=p; j++) f.slice(j) = xf.slice(j)*yf.slice(j);
 
  // find g function
  g.fill(0.0);  
  for (int j=0; j<p; j++) {
    g.slice(j)   -= (xf.slice(j+1)-xf.slice(j))*(yf.slice(j+1)-yf.slice(j));
    g.slice(j+p) += (xf.slice(j+1)-xf.slice(j))*(yf.slice(j+1)-yf.slice(j));
    for (int l=0; l<r2; l++) {
      g.slice(j+l*p)     += xf.slice(j)*yg.slice(j+l*p);
      g.slice(j+(l+1)*p) += (xf.slice(j+1)-xf.slice(j))*yg.slice(j+l*p);
    }
    for (int k=0; k<r1; k++) {
      g.slice(j+k*p)     += xg.slice(j+k*p)*yf.slice(j);
      g.slice(j+(k+1)*p) += xg.slice(j+k*p)*(yf.slice(j+1)-yf.slice(j));
      for (int l=0; l<r2; l++) g.slice(j+(k+l+1)*p) += xg.slice(j+k*p)*yg.slice(j+l*p);
    }
  }
            
  // return to R
  return(List::create(Rcpp::Named("f")=f, Rcpp::Named("g")=g));
}
