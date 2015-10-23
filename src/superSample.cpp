#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
SEXP superSample(SEXP choose, IntegerVector jj, NumericVector a, NumericVector b, NumericVector f, SEXP g) {
  // Assume that dimensions match
  // No sanity check performed!
  NumericVector vec_binom(choose);
  NumericVector vec_g(g);
  IntegerVector dims_binom = vec_binom.attr("dim");
  IntegerVector dims       = vec_g.attr("dim");
  
  int p = dims[2];                  // NB: 1 less in order to simplify code
  int new_p = jj.size();            // NB: 1 less in order to simplify code
  int r = dims[3];
  
  arma::mat  binom(vec_binom.begin(), dims_binom[0], dims_binom[1], false);
  arma::cube ff(f.begin(), dims[0], dims[1], p+1, false);
  arma::cube gg(vec_g.begin(), dims[0], dims[1], r*p, false);
  arma::cube new_f(dims[0], dims[1], new_p+1);
  arma::cube new_g(dims[0], dims[1], r*new_p);
  new_f.fill(0.0);

  // is there an non-trivial g-part?
  if (r > 0) {
    // find contribution from g function to new f function
    for (int j=0; j<new_p; j++) {
      for (int l=0; l<r; l++) new_f.slice(j) += pow(a(j),l+1)*gg.slice(jj(j)+l*p);
    }
    for (int l=0; l<r; l++) new_f.slice(new_p) += gg.slice(p-1+l*p);
  
    // find new g function
    new_g.fill(0.0);  
    for (int j=0; j<new_p; j++) {
      for (int k=0; k<r; k++) {
        for (int l=k; l<r; l++) {
          new_g.slice(j+k*new_p) += binom(l,k)*pow(a(j),l-k)*pow(b(j),k+1)*gg.slice(jj(j)+l*p);
        }
      }  
    }
    // subtract linear part in new f
    for (int j=0; j<new_p; j++) new_g.slice(j) -= new_f.slice(j+1)-new_f.slice(j);
  }
  
  // find contribution from f function to new f function
  for (int j=0; j<new_p; j++) new_f.slice(j) += (1-a(j))*ff.slice(jj(j))+a(j)*ff.slice(jj(j)+1);
  new_f.slice(new_p) += ff.slice(p);

  // return to R
  return(List::create(Rcpp::Named("f")=new_f, Rcpp::Named("g")=new_g));
}
