#-------------------------------------------------------------------------------#
# Package: Multistage Adaptive Enrichment Design                                #
# maed(): The user interface for maed.g(), maed.l(), maed.p(), maed.2aed.main() #
#-------------------------------------------------------------------------------#

#' Multistage adaptive enrichment design
#'
#' The main function for multistage adaptive enrichment design.
#'
#' This function targets on the formulation of the sparse linear programming problem for 2-stage, 2-sub-population adaptive enrichment design. The solver is defaulted to be a pre-exsiting LP library \code{'PRIMAL'}.
#' For high-dimensional & high-sparsity setting, we adopt the sub-space optimization methods and row generation methods. We leave these for future implementation.
#'
#' @param n_stage Number of stages. Default to be 2.
#' @param n_ppl Number of sub-populations. Default to be 2.
#' @param prop_ppl A list of proportion of each sub-population, e.g., in 2-stage case, prop_ppl = \code{(\rho_1, \rho_2)}. Default to be \code{c(0.5, 0.5)}.
#' @param R_sizes A list of numbers of discretization at each stage. Default to be \code{c(100, 100)}.
#' @param D A list of discrete decision values after each stage. Specified by users.
#' @param lambda_0 The prior distribution where the mean of the random treatment effects (in the LP objective) of each sub-populations are initially drawn. Default to be \code{normal}.
#' @param loss_tables_0  A loss lookup table with user specified loss under each decision pattern (in the LP objective). Table size should agree with \code{R_sizes}.
#' @param lambda_c2 The prior distribution where the mean of the random treatment effects (in the LP constraints) of each sub-populations are initially drawn. Default to be the same as \code{lambda_0}.
#' @param loss_tables_c2 A loss lookup table with user specified loss under each decision pattern (in the LP constaints). Table size should agree with \code{R_sizes}. Default to be the same as \code{loss_tables_0}.
#' @param G A subspace of the original space of mean vectors of the random treatment effects (in the LP objective), e.g. in 2-sub-population case, \code{G \subset \{(\mu_1, \mu_2,)|\mu_1=0, or \mu_2=0, or \rho_1\mu_1 + \rho_2\mu_2 = 0\}}. The size of G is \code{G_size, n_ppl}.
#' @param G_c2 A subspace of the original space of mean vectors of the random treatment effects (in the LP constaints).
#' @param G_size Dimension of G and G_c2. Default to be \code{500}.
#' @param alpha C1 constraints values that the user should specify. Its size should agree with G_size.
#' @param J Total number of C2 constraints. Specified by users.
#' @param beta C2 constraints values that the user should specify. Its size should agree with J.
#' @solver The solver solving the formulated LP problem. Default to be \code{PRIMAL}.
#' @return
#' An object with S3 class \code{"maed"} is returned:
#' \item{lambda}{
#'   The sequence of regularization parameters \code{lambda} obtained in the program.
#' }
#' \item{value}{
#'   The sequence of optimal value of the object function corresponded to the sequence of \code{lambda}.
#' }
#' \item{n_stage}{
#'   The number of stages.
#' }
#' \item{n_ppl}{
#'   The number of sub-populations.
#' }
#' \item{solver}{
#'   The name of solver, e.g. \code{PRIMAL}
#' }


#' @seealso \code{\link{maed.p}}, \code{\link{maed.l}}, \code{\link{maed.g}}, \code{\link{maed.2aed.main}}, and \code{\link{maed-package}}.
#' @example
#'
#' @export
maed <- function(n_stage=2, n_ppl=2, prop_ppl=c(0.5,0.5), R_sizes=c(100,100),
                 D=list(c(-0.3,-0.1,0.1,0.3), c(-0.3, -0.1, -0.1, -0.1, 0.1, 0.1, 0.1, 0.3)),
                 lambda_0="normal", loss_tables_0=NULL,
                 lambda_c2="normal", loss_tables_c2=NULL,
                 G=NULL, G_c2=NULL, G_size=500, alpha=NULL, J=2, beta=NULL,
                 solver="primal"){

  # general
  assert(n_stage == length(R_sizes))
  assert(n_stage == length(D))
  est = list()
  est$n_stage = n_stage
  est$n_ppl = n_ppl
  est$solver = solver

  if(n_stage==2 && n_ppl==2)
  {
    if(is.null(G)){
      G <- maed.g(lambda_0, G_size, n_ppl, prop_ppl) # [|G|, n_ppl]
    }

    vars <- maed.p(G, R_sizes)
    P = vars[[1]] # [|G|, |R1|, |R2|]
    G_idx = vars[[2]] #[|G|, n_stage]

    L <- maed.l(G_idx, D, loss_tables_0, R_sizes) # [1, |G|, |D[1]|, |D[2]|]

    # C2 constraints
    if(J > 0){

      if(is.null(G_c2)){
        G_c2 <- maed.g(lambda_c2, G_size, n_ppl, prop_ppl) #[|G|, n_ppl]
      }

      vars <- maed.p(G_c2, R_sizes)
      P_c2 = vars[[1]] #[|G|, |R1|, |R2|]
      G_c2_idx = vars[[2]] #[|G|, n_stage]

      L_c2 <- maed.l(G_c2_idx, D, loss_tables_c2, R_sizes, J) #[J, |G|, |D[1]|, |D[2]|]
    }

    fit = maed.2aed.main(L, P, L_c2, P_c2, alpha, beta, solver,
                         G_size, D, R_sizes, J)

    est$lambda = fit$lambda
    est$value = fit$value
    rm(fit)
  }
  else{
    stop("Design with >=2 stages is currently not supported.")
  }

  est$data = G
  est$loss = loss_tables_c2
  class(est) = "maed"
  return(est)
}

#' @useDynLib maed
#' @importFrom Rcpp sourceCpp
NULL
