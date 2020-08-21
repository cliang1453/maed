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

//' Compute Discretization Grid Unit
//'
//' @param R_sizes, w, grid_len
//' @export
// [[Rcpp::export]]
NumericVector grid(NumericVector R_sizes, double w, int n_stage){
  NumericVector grid_unit(n_stage);
  for(int i=0; i<n_stage; i++){
    grid_unit[i] = w / sqrt(R_sizes[i]) + 0.0001;
  }
  return grid_unit;
}

//' element-wise floor operation on matrix
//'
//' @param
//' @export
// [[Rcpp::export]]
IntegerMatrix floormat(NumericMatrix G){

  IntegerMatrix G_(G.nrow(), G.ncol());

  for(int r=0; r<G.nrow(); r++){
    for(int c=0; c<G.ncol(); c++){
      G_(r,c) = int(floor(G(r,c)));
    }
  }

  return G_;
}

//' Discretization
//'
//' @param
//' @export
// [[Rcpp::export]]
IntegerMatrix discretize(NumericMatrix &G, NumericVector &R_sizes, int G_size, double w, int n_pop, int n_stage){

  NumericVector grid_unit = grid(R_sizes, w, n_stage);
  IntegerMatrix idx(G_size, n_stage);

  for(int i=0; i<n_stage; i++){

    IntegerMatrix idx_n = floormat(G/grid_unit[i]); //G_size, n_pop

    for(int j=0; j<n_pop-1; j++){
      idx_n(_, j+1) = idx_n(_, j) * sqrt(R_sizes[i]) + idx_n(_, j+1);
    }

    idx(_, i) = idx_n(_, n_pop-1);
  }

  return idx;
}

