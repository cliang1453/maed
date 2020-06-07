#include <Rcpp.h>
#include <RcppEigen.h>
#include "math.h"
using namespace Rcpp;
using namespace Eigen;
using namespace std;
//[[Rcpp::depends(RcppEigen)]]


//' Average
//'
//' @param x A vector
//' @export
// [[Rcpp::export]]
double meanC(NumericVector x) {
  int n = x.size();
  double total = 0;

  for(int i = 0; i < n; ++i) {
    total += x[i];
  }
  Rcout << "the mean is: "<<total / n;
  return total / n;
}

