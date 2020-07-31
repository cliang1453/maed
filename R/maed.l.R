#' @export
maed.l <- function(G_idx=NULL, D=NULL, loss=NULL){
  library(abind)

  L = list()
  for (j in 1:length(loss)){
    L_table <- maed.loss_lookup_table(G_idx, D, loss[j]) #[n_stage + 1, ]
    L_j = L_table[[1]]
    for(n in 1:(length(L_table)-1)){
      L_j = L_j %o% L_table[[n+1]]
    }
    L[[j]] = L_j
  }

  L = abind(L, along=0) # [J, |G|, |D[1]|, |D[2]|]
  return(L)

}

#' @export
maed.loss_lookup_table <- function(G_idx=NULL, D=NULL, loss="default"){

  L_table = list()
  n_stage <- length(D)

  if(loss == "default"){
    # two-stage, two-sub-population
    L_table[[1]] = rowSums(G_idx)
    for(n in 1:(n_stage-1)){
      L_table[[n+1]] = D[[n]]
    }
    L_table[[n_stage+1]] = c(0, 1, 1, 1, 2, 2, 2, 3)
  }
  return(L_table)
}

#' @useDynLib maed
#' @importFrom Rcpp sourceCpp
NULL
