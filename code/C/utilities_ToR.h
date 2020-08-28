#pragma once
using namespace std;
// using namespace arma;
using namespace Rcpp;  // for including R with Rcpp 

IntegerVector sort_indexes_rcpp(NumericVector data);
NumericVector empirical_cdf_rcpp(NumericVector x);
