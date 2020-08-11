#' @export
maed.p <- function(G=NULL, R_sizes=c(100, 100)){
  # G: [G_size, n_pop]

  P_ = list()
  G_size = dim(G)[1]
  n_stage = length(R_sizes)
  P = array(0, dim=c(G_size, R_sizes)) #[G_size, R1_size, R2_size]

  grid_unit = rep(max(G)-min(G), n_stage)/sqrt(R_sizes) + 0.001 # grid_unit: [n_stage]
  idx = maed.discretize(G-min(G), grid_unit, R_sizes) #idx: [G_size, n_stage]

  for(r in 1:length(R_sizes)){
    x = array(1, dim=c(G_size,1)) %*% seq(1:R_sizes[r]) # [G_size, R1_size]
    P_[[r]] = dnorm(x, mean=array(idx[,r], dim=c(G_size,1)), sd=R_sizes[r]*0.2) # [G_size, R1_size]
    P_[[r]] = P_[[r]]/rowSums(P_[[r]])
  }

  for(g in 1:G_size){
    prod = P_[[1]][g,]
    for(r in 2:length(R_sizes)){
      prod = prod %o% P_[[r]][g,]
    } #[R1_size, R2_size]
    P[g,,] = prod
  }

  return(list(P, idx))
}

#' @export
maed.discretize <- function(G=NULL, grid_unit=NULL, R_sizes=c(100, 100)){

  G_size = dim(G)[1]
  n_pop = dim(G)[2]
  n_stage = length(R_sizes)
  idx = zeros(G_size, n_stage)

  for(n in 1:n_stage){
    idx_n = G%/%grid_unit[n] #[G_size, n_pop]
    for(m in 1:(n_pop-1)){
      idx_n[,m+1] = idx_n[,m] * sqrt(R_sizes[n]) + idx_n[,m+1]
    }
    idx[,n] = idx_n[,n_pop] + 1
  }

  return(idx)
}

#' @useDynLib maed
#' @importFrom Rcpp sourceCpp
NULL
