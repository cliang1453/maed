#' @export
maed.p <- function(G=NULL, R_sizes=c(100, 100)){
  # G: [G_size, n_pop]

  idx = list()
  P = list()
  n_stage = length(R_sizes)

  grid_unit = rep(max(G)-min(G), n_stage)/sqrt(R_sizes) + 0.001 #grid_unit: [n_stage]
  idx = maed.discretize(G-min(G), grid_unit, R_sizes) #idx: [G_size, n_stage + 1]
  P = array(0, dim=c(dim(G)[1], R_sizes)) #P: [G_size, D1, D2, ..., Dn_stage]
  P[idx] = 1

  return(list(P, idx[,2:3]))
}

#' @export
maed.discretize <- function(G=NULL, grid_unit=NULL, R_sizes=c(100, 100)){

  G_size = dim(G)[1]
  n_pop = dim(G)[2]
  n_stage = length(R_sizes)
  idx = zeros(G_size, n_stage + 1)
  idx[,1] = 1:G_size

  for(n in 1:n_stage){
    idx_n = G%/%grid_unit[n] #[G_size, n_pop]
    for(m in 1:(n_pop-1)){
      idx_n[,m+1] = idx_n[,m] * sqrt(R_sizes[n]) + idx_n[,m+1]
    }
    idx[,n+1] = idx_n[,n_pop]
  }

  return(idx)
}

#' @useDynLib maed
#' @importFrom Rcpp sourceCpp
NULL
