#' @export
maed.p <- function(G=NULL, R_sizes=c(100, 100)){
  library(abind)

  G_size = dim(G)[1]
  n_pop = dim(G)[2]
  n_stage = length(R_sizes)
  w = max(G)-min(G)

  idx <- .Call("_maed_discretize", G-min(G), R_sizes, G_size, w, n_pop, n_stage)

  P <- .Call("_maed_compute_p", idx, R_sizes, G_size, n_stage)

  P <- abind(P, along=0)
  P <- array(P, dim=c(G_size, R_sizes))

  return(list(P, idx))
}

#' @useDynLib maed
#' @importFrom Rcpp sourceCpp
NULL
