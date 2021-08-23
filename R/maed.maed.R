#-------------------------------------------------------------------------------#
# Package: Multistage Adaptive Enrichment Design                                #
# maed.maed.main(): The main function for solving the formulated LP for         #
#                   multi-stage, 2-sub-population AED problem.                  #
#-------------------------------------------------------------------------------#

#' The main function for solving the formulated LP for multi-stage, 2-sub-population adaptive enrichment design problem.
#'
#' The solver is defaulted to be a pre-exsiting LP library \code{'PRIMAL'}. For high-dimensional & high-sparsity setting, we adopt the sub-space optimization methods and row generation methods. We leave these for future implementation.
#'
#' @param L A loss matrix (in the objective) of size \code{n_dim}.
#' @param P A probability matrix (in the objective) of size \code{n_dim}.
#' @param L_c2 A loss matrix (in the C2 constraints) of size \code{J} by \code{n_dim}.
#' @param P_c2 A probability matrix (in the C2 constraints) of size \code{J} by \code{n_dim}.
#' @param beta C2 constraints values that the user should specify. Its size should agree with \code{J}.
#' @param solver The solver solving the formulated LP problem. Default to be \code{PRIMAL}.
#' @param J Total number of C2 constraints that the user should specify.
#' @param n_dim Total number of possible trails.
#' @return
#' An object with S3 class \code{"maed.maed"} is returned:
#' \item{lambda}{
#'   The sequence of regularization parameters \code{lambda} obtained in the program.
#' }
#' \item{value}{
#'   The sequence of optimal value of the object function corresponded to the sequence of \code{lambda}.
#' }

#' @seealso \code{\link{maed.p}}, \code{\link{maed.l}}, \code{\link{maed.g}}, and \code{\link{maed-package}}.
#' @export
maed.maed.main <- function(L=NULL, P=NULL, L_c2=NULL, P_c2=NULL, beta=NULL, solver="primal", J=NULL, n_dim=NULL){

  # print(sum(P==0))
  # print(sum(L==0))
  # print(sum(P_c2==0))
  # print(sum(L_c2==0))

  c <- -1 * L * P # [n_dim]

  A <- L_c2 * P_c2 #[J, n_dim]

  b <- beta #[J]

  est = list()

  if(solver=="primal"){

    library(PRIMAL)

    m = dim(A)[1]
    n = dim(A)[2]

    b_bar = rep(0,m)
    c_bar = rep(0,m+n)
    B_init = seq(n,n+m-1)

    fit <- PSM_solver(A=A, b=b, b_bar=b_bar, c=c, c_bar=c_bar, B_init=B_init)
    est$lambda = fit$lambda
    est$value = fit$value
    print(fit)
    rm(fit)

  }else{
    stop("The specified solver package is not supported.")
  }

  class(est) = "maed.maed"
  return(est)
}


#' @useDynLib maed
#' @importFrom Rcpp sourceCpp
NULL
