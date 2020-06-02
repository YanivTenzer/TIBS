#include <Rcpp.h>

using namespace Rcpp;

RcppExport SEXP helloA2() {
	printf("Hello World!\n");
	return(R_NilValue);
}
