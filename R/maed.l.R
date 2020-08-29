#-------------------------------------------------------------------------------#
# Package: Multistage Adaptive Enrichment Design                                #
# maed.l(): The function for loss matrix generation                             #
#-------------------------------------------------------------------------------#

#' The function perform L matrix generation given the user specified loss values.
#'
#' @param G_idx A matrix of discretized indices of size \code{G_size} by \code{n_stage}.
#' @param D A list of discrete decision values after each stage. Specified by users.
#' @param loss_tables A loss lookup table with user specified loss under each decision pattern. Table size should agree with with \code{R_sizes}.
#' @param R_sizes A list of numbers of discretization at each stage. Default to be \code{c(100, 100)}.
#' @param J Total number of C2 constraints. Specified by users.
#' @return
#' An matrix is returned:
#' \item{L}{
#'   A loss matrix of size \code{J} by \code{G_size} by \code{D[1]} by ... by \code{D[-1]}.
#' }

#' @seealso \code{\link{maed}}, \code{\link{maed.g}}, \code{\link{maed.p}}, \code{\link{maed.2aed.main}}, and \code{\link{maed-package}}.
#' @example
#'
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
