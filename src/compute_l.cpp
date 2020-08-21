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

//'
//'
//' @param L0, L1
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd compute_l(Eigen::VectorXd &L0, List &L1){

  Eigen::VectorXd L10 = L1[0];
  Eigen::MatrixXd L = L0 * L10.transpose();

  for(int i=1; i<L1.size(); i++){
    L.resize(L.rows()*L.cols(), 1);

    Eigen::VectorXd L1i = L1[i];
    L = L * L1i.transpose();
  }

  return L;
}
