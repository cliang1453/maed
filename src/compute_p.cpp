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

//' Construct P based on current discretization
//'
//' @param
//' @export
// [[Rcpp::export]]
List compute_p(IntegerMatrix &idx, IntegerVector &R_sizes, int G_size, int n_stage){

  List L(n_stage);
  List P(G_size);

  for(int c=0; c<n_stage; c++){

    IntegerVector v = seq_len(R_sizes[c]);
    NumericVector xc = as<NumericVector>(v); //[R_size[c]]
    Eigen::MatrixXd Pc(G_size, R_sizes[c]);

    for(int r=0; r<G_size; r++){
      double m = idx(r,c);
      NumericVector v = dnorm(xc, m);
      Pc.row(r) = as<Eigen::VectorXd>(v);
    }

    L[c] = Pc;
  }

  for(int g=0; g<G_size; g++){

    Eigen::MatrixXd P0 = L[0];
    Eigen::MatrixXd P1 = L[1];
    Eigen::MatrixXd R = P0.row(g).transpose() * P1.row(g);

    for(int c=2; c<n_stage; c++){
      R.resize(R.rows()*R.cols(), 1);

      Eigen::MatrixXd Pc = L[c];
      R = R * Pc.row(g);
    }

    P[g] = R;
  }

  return P;
}
