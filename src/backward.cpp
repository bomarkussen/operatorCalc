#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
SEXP backward(NumericVector diffs, NumericVector f, SEXP g) {
  // Assume that dimensions match
  // No sanity check performed!
  NumericVector vec_g(g);
  IntegerVector dims = vec_g.attr("dim");

  int p = dims[2];     // NB: 1 less in order to simplify code
  int r = dims[3];
  int new_r = std::max(1,r)+1;
  
  arma::cube ff(f.begin(), dims[0], dims[1], p+1, false);
  arma::cube gg(vec_g.begin(), dims[0], dims[1], r*p, false);
  arma::cube integral_f(dims[0], dims[1], p+1);
  arma::cube integral_g(dims[0], dims[1], new_r*p);

  // find f function
  integral_f.slice(p).fill(0.0); 
  for (int j=p; j>0; j--) {
    integral_f.slice(j-1) = integral_f.slice(j)+diffs(j-1)/2*(ff.slice(j-1)+ff.slice(j));
    for (int k=0; k<r; k++) integral_f.slice(j-1) += diffs(j-1)/(k+2)*gg.slice(j-1+k*p);
  }
  
  // find g function
  for (int j=0; j<p; j++) {
    integral_g.slice(j+p) = -diffs(j)/2*(ff.slice(j+1)-ff.slice(j));
    if (r>0) integral_g.slice(j+p) -= diffs(j)/2*gg.slice(j);
    for (int k=2; k<new_r; k++) integral_g.slice(j+k*p) = -diffs(j)/(k+1)*gg.slice(j+(k-1)*p);
    integral_g.slice(j) = -integral_g.slice(j+p);
    for (int k=2; k<new_r; k++) integral_g.slice(j) -= integral_g.slice(j+k*p);
  }

  // return to R
  return(List::create(Rcpp::Named("f")=integral_f, Rcpp::Named("g")=integral_g));
}
