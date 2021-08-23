#-------------------------------------------------------------------------------#
# Package: Multistage Adaptive Enrichment Design                                #
# maed.l(): The function for loss matrix generation                             #
#-------------------------------------------------------------------------------#

#' The function perform L matrix generation given the user specified loss values.
#'
#' @param G_idx A matrix of discretized indices of size \code{G_size} by \code{n_stage}.
#' @param reward A list of discrete decision values after each stage. Specified by users.
#' @param loss_tables A loss lookup table with user specified loss under each decision pattern. Table size should agree with with \code{state_dim}.
#' @param state_dim A list of numbers of discretization at each stage. Default to be \code{c(100, 100)}.
#' @param J Total number of C2 constraints. Specified by users.
#' @return
#' An matrix is returned:
#' \item{L}{
#'   A loss matrix of size \code{J} by \code{G_size} by \code{|reward[1]|} by ... by \code{|reward[-1]|}.
#' }

#' @seealso \code{\link{maed}}, \code{\link{maed.g}}, \code{\link{maed.p}}, \code{\link{maed.2aed.main}}, and \code{\link{maed-package}}.
#' @example
maed.l <- function(G_idx=NULL, reward=NULL, loss_tables=NULL, state_dim=NULL, J=1){
  library(abind)

  G_size <- dim(G_idx)[1]
  R1_size <- state_dim[1]
  R2_size <- state_dim[2]
  D_size <- c(length(reward[[1]]), length(reward[[2]]))
  L <- list()

  for (j in 1:J){

    if(is.null(loss_tables)){
      loss_table <- rnorm(R1_size * R2_size)
    }else{
      loss_table <- loss_tables[j] # R1_size * R2_size
    }

    L0 <- loss_table[G_idx[,1]*R2_size + G_idx[,2] + 1] # G_size
    Lj <- .Call("_maed_compute_l", L0, reward)

    L[[j]] <- array(Lj, dim=c(G_size, D_size))
  }

  L <- abind(L, along=0) # [J, |G|, |reward[1]|, |reward[2]|]
  return(L)
}

#' @useDynLib maed
#' @importFrom Rcpp sourceCpp
NULL
