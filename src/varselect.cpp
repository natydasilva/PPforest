#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector varselect(IntegerVector ids, int sampS,
                        bool replace,
                        NumericVector prob = NumericVector::create()) {
  NumericVector sampl(sampS);
 
  sampl = Rcpp::RcppArmadillo::sample(ids, sampS,replace, prob);
  NumericVector y = clone(sampl);
  std::sort(y.begin(), y.end());
  return y;
  }

