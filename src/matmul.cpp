#include "math.h"
#include <vector>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
using namespace Rcpp;
using namespace std;
using namespace Eigen;
//[[Rcpp::depends(RcppEigen)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//' matrix multiplication
//'
//' @param
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd matmul(Eigen::MatrixXd &A1, Eigen::MatrixXd &A2){
  return A1*A2;
}
