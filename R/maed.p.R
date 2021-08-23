#-------------------------------------------------------------------------------#
# Package: Multistage Adaptive Enrichment Design                                #
# maed.p(): The function for probability matrix generation.                     #
#-------------------------------------------------------------------------------#

#' The function performs space discretization and P matrix generation.
#'
#' @param state_dim A list of numbers of discretization at each stage. Default to be \code{c(100, 100)}.
#' @param G A subset of mean of the random treatment effects (in the LP objective), e.g. in 2-sub-population case, \code{G \subset \{(\mu_1, \mu_2,)|\mu_1=0, or \mu_2=0, or \rho_1\mu_1 + \rho_2\mu_2 = 0\}}. The size of G is \code{G_size, n_ppl}.
#' @return
#' An list is returned:
#' \item{P}{
#'   A probability matrix of size \code{G_size} by \code{state_dim[1]} by ... by \code{state_dim[-1]}.
#' }
#' \item{idx}{
#'   A matrix of discretized indices of size \code{G_size} by \code{n_stage}.
#' }

#' @seealso \code{\link{maed}}, \code{\link{maed.l}}, \code{\link{maed.g}}, \code{\link{maed.2aed.main}}, and \code{\link{maed-package}}.
#' @export
maed.p <- function(G=NULL, state_dim=c(100, 100)){
  library(abind)

  G_size = dim(G)[1]
  n_pop = dim(G)[2]
  n_stage = length(state_dim)
  w = max(G)-min(G)

  idx <- .Call("_maed_discretize", G-min(G), state_dim, G_size, w, n_pop, n_stage)

  P <- .Call("_maed_compute_p", idx, state_dim, G_size, n_stage)

  P <- abind(P, along=0)
  P <- array(P, dim=c(G_size, state_dim))

  return(list(P, idx))
}

#' @useDynLib maed
#' @importFrom Rcpp sourceCpp
NULL
