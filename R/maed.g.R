maed#' @export
maed.g <- function(lambda=NULL, G_size=500, n_ppl=2, prop_ppl=c(0.5, 0.5)){

  library(pracma)
  library(base)

  if(lambda == "normal"){
    x = matrix(rnorm((n_ppl-1)*G_size,0,1), G_size, n_ppl-1)
    prop = repmat(matrix(prop_ppl[1:n_ppl-1], 1, n_ppl-1), G_size, 1)
    x_ = rep(0, G_size) - rowSums(x * prop)
    G = cbind(x, x_)
  }
  return(G) #[G_size, n_ppl]
}

#' @useDynLib maed
#' @importFrom Rcpp sourceCpp
NULL
