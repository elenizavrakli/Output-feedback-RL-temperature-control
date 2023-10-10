// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::export]]
arma::mat sumcrossprod(arma::mat A, int low, int up) {
  
  // Dimensions of A
  int n = A.n_rows, m = A.n_cols;
  
  arma::mat S(m, m); // Output matrix
  S.zeros(); // Fill with zeros
  
  for (int i = low; i < (up+1); i++) {
    S += A.row(i).t() * A.row(i);
  }
  
  return S;
}

