#-------------------------------------------------------------------------------#
# Package: Multistage Adaptive Enrichment Design                                #
# maed.g(): The function for sampling a subspace of the original space of       #
#           mean vectors of the random treatment effects                        #
#-------------------------------------------------------------------------------#

#' The function perform G matrix generation if user not specified.
#'
#' @param lambda The prior distribution where the mean of the random treatment effects of each sub-populations are initially drawn. Default to be \code{normal}.
#' @param G_size Dimension of G. Default to be \code{500}.
#' @param n_ppl Number of sub-populations. Default to be 2.
#' @param prop_ppl A list of proportion of each sub-population, e.g., in 2-stage case, prop_ppl = \code{(\rho_1, \rho_2)}. Default to be \code{c(0.5, 0.5)}.
#' @return
#' An matrix is returned:
#' \item{G}{
#'   A probability matrix of size \code{G_size} by \code{n_ppl}.
#' }

#' @seealso \code{\link{maed}}, \code{\link{maed.l}}, \code{\link{maed.p}}, \code{\link{maed.2aed.main}}, and \code{\link{maed-package}}.
#' @example
#'
#' @export
maed.g <- function(lambda=NULL, G_size=500, n_ppl=2, prop_ppl=c(0.5, 0.5)){

  library(pracma)
  library(base)

  if(lambda == "normal"){
    x = matrix(rnorm((n_ppl-1)*G_size,0,1), G_size, n_ppl-1)
    prop = repmat(matrix(prop_ppl[1:n_ppl-1], 1, n_ppl-1), G_size, 1)
    x_ = rep(0, G_size) - rowSums(x * prop)
    G = cbind(x, x_)
  }
  else{
    stop("Lambda type is currently not supported.")
  }
  return(G) #[G_size, n_ppl]
}

#' @useDynLib maed
#' @importFrom Rcpp sourceCpp
NULL
