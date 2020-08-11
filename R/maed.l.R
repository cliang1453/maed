#' @export
maed.l <- function(G_idx=NULL, D=NULL, loss=NULL, R_sizes=NULL){
  library(abind)

  L = list()
  for (j in 1:length(loss)){
    L_table <- maed.loss_lookup_table(G_idx, D, loss[j], R_sizes) #[n_stage + 1, ]
    L_j = L_table[[1]]
    for(i in 2:length(L_table)){
      L_j = L_j %o% L_table[[i]]
    }
    L[[j]] = L_j
  }

  L = abind(L, along=0) # [J, |G|, |D[1]|, |D[2]|]
  return(L)

}

#' @export
maed.loss_lookup_table <- function(G_idx=NULL, D=NULL, loss="default", R_sizes=NULL){

  L_table = list()
  n_stage <- length(D)
  G_size = dim(G_idx)[1]
  R1_size = R_sizes[1]
  R2_size = R_sizes[2]

  if(loss == "default"){
    # two-stage, two-sub-population
    table = rnorm(R1_size*R2_size)
    L_table[[1]] = table[G_idx[,1]*R2_size + G_idx[,2]]
    for(n in 2:n_stage){
      L_table[[n]] = D[[n-1]]
    }
    L_table[[n_stage+1]] = c(-0.3, -0.1, -0.1, -0.1, 0.1, 0.1, 0.1, 0.3)
  }
  return(L_table)
}

#' @useDynLib maed
#' @importFrom Rcpp sourceCpp
NULL
