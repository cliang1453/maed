#' @export
maed.l <- function(G_idx=NULL, D=NULL, loss_tables=NULL, R_sizes=NULL, J=1){
  library(abind)

  G_size <- dim(G_idx)[1]
  R1_size <- R_sizes[1]
  R2_size <- R_sizes[2]
  D_size <- c(length(D[[1]]), length(D[[2]]))
  L <- list()

  for (j in 1:J){

    if(is.null(loss_tables)){
      loss_table <- rnorm(R1_size * R2_size)
    }else{
      loss_table <- loss_tables[j] # R1_size * R2_size
    }

    L0 <- loss_table[G_idx[,1]*R2_size + G_idx[,2] + 1] # G_size
    Lj <- .Call("_maed_compute_l", L0, D)

    L[[j]] <- array(Lj, dim=c(G_size, D_size))
  }

  L <- abind(L, along=0) # [J, |G|, |D[1]|, |D[2]|]
  return(L)
}

#' @useDynLib maed
#' @importFrom Rcpp sourceCpp
NULL
