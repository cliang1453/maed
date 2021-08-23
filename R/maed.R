#-------------------------------------------------------------------------------#
# Package: Multistage Adaptive Enrichment Design                                #
# maed(): The user interface for maed.maed.main(), maed.2aed.main().            #
#-------------------------------------------------------------------------------#

#' Multistage adaptive enrichment design
#'
#' The main function for multistage adaptive enrichment design.
#'
#' This function targets on the formulation of the sparse linear programming problem for multi-stage, 2-sub-population adaptive enrichment design. The solver is defaulted to be a pre-exsiting LP library \code{'PRIMAL'}.
#' For high-dimensional & high-sparsity setting, we adopt the sub-space optimization methods and row generation methods. We leave these for future implementation.
#'
#' @param n_stage Number of stages. Default to be 3.
#' @param state_dim Dimension of state space in each stage. Default to be c(3,3,3).
#' @param decision_dim Dimension of decision space in each stage. Default to be c(4,4,4).
#' @param prob_reward_file User input probability and reward for individual trails. A sample format is provided in "reward_prob_s4_a3.txt". Now the format only support G_size=1 case.
#' @param solver The solver solving the formulated LP problem. Default to be \code{PRIMAL}.
#' @param J Total number of C2 constraints. Default to be 1.
#' @param beta C2 constraints values that the user should specify. Its size should agree with J.
#' @param n_ppl Number of sub-populations. No need to specify if \code{prob_reward_file} is specified.
#' @param prop_ppl A list of proportion of each sub-population, e.g., in 2-stage case, prop_ppl = \code{(\rho_1, \rho_2)}. No need to specify if \code{prob_reward_file} is specified.
#' @param reward A list of discrete rewards after each stage. No need to specify if \code{prob_reward_file} is specified.
#' @param lambda_0 The prior distribution where the mean of the random treatment effects (in the LP objective) of each sub-populations are initially drawn. No need to specify if \code{prob_reward_file} is specified.
#' @param loss_tables_0  A loss lookup table with user specified loss under each decision pattern (in the LP objective). Table size should agree with \code{state_dim}. No need to specify if \code{prob_reward_file} is specified.
#' @param lambda_c2 The prior distribution where the mean of the random treatment effects (in the LP constraints) of each sub-populations are initially drawn. No need to specify if \code{prob_reward_file} is specified.
#' @param loss_tables_c2 A loss lookup table with user specified loss under each decision pattern (in the LP constaints). Table size should agree with \code{state_dim}. No need to specify if \code{prob_reward_file} is specified.
#' @param G A subspace of the original space of mean vectors of the random treatment effects (in the LP objective), e.g. in 2-sub-population case, \code{G \subset \{(\mu_1, \mu_2,)|\mu_1=0, or \mu_2=0, or \rho_1\mu_1 + \rho_2\mu_2 = 0\}}. The size of G is \code{G_size, n_ppl}. No need to specify if \code{prob_reward_file} is specified.
#' @param G_c2 A subspace of the original space of mean vectors of the random treatment effects (in the LP constaints). No need to specify if \code{prob_reward_file} is specified.
#' @param G_size Dimension of G and G_c2. No need to specify if \code{prob_reward_file} is specified.
#' @param alpha C1 constraints values that the user should specify. Its size should agree with \code{G_size}. No need to specify if \code{n_stages} > 2.
#' @return
#' An object with S3 class \code{"maed"} is returned:
#' \item{lambda}{
#'   The sequence of regularization parameters \code{lambda} obtained in the program.
#' }
#' \item{value}{
#'   The sequence of optimal value of the object function corresponded to the sequence of \code{lambda}.
#' }
#' }


#' @seealso \code{\link{maed.p}}, \code{\link{maed.l}}, \code{\link{maed.g}}, \code{\link{maed.2aed.main}}, \code{\link{maed.maed.main}}, and \code{\link{maed-package}}.
#' @export
maed <- function(n_stage=3,
                 state_dim=c(4,4,4),
                 decision_dim=c(3,3,3),
                 prob_reward_file=NULL,
                 solver="primal",
                 J=1,
                 beta=c(1),
                 n_ppl=2, prop_ppl=c(0.5,0.5),
                 reward=list(c(-0.3,-0.1,0.1,0.3), c(-0.3,-0.1,0.1,0.3), c(-0.3,-0.1,0.1,0.3)),
                 lambda_0="normal", loss_tables_0=NULL,
                 lambda_c2="normal", loss_tables_c2=NULL,
                 G=NULL, G_c2=NULL, G_size=500, alpha=NULL){


  if(!is.null(prob_reward_file)){

    n_dim = prod(state_dim)*prod(decision_dim)

    # P =  # [state_dim[0]*decision_dim[0], ..., state_dim[-1]*decision_dim[-1]]
    # P_c2 # [J, state_dim[0]*decision_dim[0], ..., state_dim[-1]*decision_dim[-1]]
    # L # [state_dim[0]*decision_dim[0], ..., state_dim[-1]*decision_dim[-1]]
    # L_c2 # [J, state_dim[0]*decision_dim[0], ..., state_dim[-1]*decision_dim[-1]]

    P <- rep(0, n_dim)
    P_c2 <- array(rep(0, J*n_dim), append(J, n_dim))
    L <- rep(0, n_dim)
    L_c2 <- array(rep(0, J*n_dim), append(J, n_dim))

    data <- as.matrix(read.delim(prob_reward_file, header=FALSE, sep ="\t"))
    dim <- 1
    for (r in 1:nrow(data)){
      P[dim] <- data[r,n_stage*2+1]
      L[dim] <- data[r,n_stage*2+2]
      for (j in 1:J){
        P_c2[j, dim] <- data[r,n_stage*2+1]
        L_c2[j, dim] <- data[r,n_stage*2+1+j]
      }
      dim <- dim + 1
    }

    fit = maed.maed.main(L, P, L_c2, P_c2, beta, solver, J, n_dim)
    est = list()
    est$lambda = fit$lambda
    est$value = fit$value
    rm(fit)

    class(est) = "maed"
    return(est)
  }
  else{
    stop("No user input specified.")
  }

  ##################################################################################################
  # If you wish to integrate the P and L generation codes into this package,                       #
  # we provide below a 2AED framework for further development.                                     #
  # Please also check out maed.g, maed.p, maed.l and maed.p.                                       #
  ##################################################################################################

  if(is.null(G)){
    G <- maed.g(lambda_0, G_size, n_ppl, prop_ppl) # [|G|, n_ppl]
  }

  vars <- maed.p(G, state_dim)
  P = vars[[1]] # [|G|, |R1|, |R2|]
  G_idx = vars[[2]] #[|G|, n_stage]

  L <- maed.l(G_idx, reward, loss_tables_0, state_dim) # [1, |G|, |D[1]|, |D[2]|]

  # C2 constraints
  if(J > 0){

    if(is.null(G_c2)){
      G_c2 <- maed.g(lambda_c2, G_size, n_ppl, prop_ppl) #[|G|, n_ppl]
    }

    vars <- maed.p(G_c2, state_dim)
    P_c2 = vars[[1]] #[|G|, |R1|, |R2|]
    G_c2_idx = vars[[2]] #[|G|, n_stage]

    L_c2 <- maed.l(G_c2_idx, reward, loss_tables_c2, state_dim, J) #[J, |G|, |D[1]|, |D[2]|]
  }

  fit = maed.2aed.main(L, P, L_c2, P_c2,
                       alpha, beta, solver,
                       G_size, reward, state_dim, J)
  est = list()
  est$lambda = fit$lambda
  est$value = fit$value
  rm(fit)

  class(est) = "maed"
  return(est)

}

#' @useDynLib maed
#' @importFrom Rcpp sourceCpp
NULL
