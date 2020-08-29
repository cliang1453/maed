#-------------------------------------------------------------------------------#
# Package: Multistage Adaptive Enrichment Design                                #
# maed.2aed.main(): The interface for maed.2aed.construct_param,                #
#                   maed.2aed.construct_variables, maed.2aed.alg                #
#-------------------------------------------------------------------------------#

#' The main function for solving the formulated LP for 2-stage, 2-sub-population adaptive enrichment design problem.
#'
#' The solver is defaulted to be a pre-exsiting LP library \code{'PRIMAL'}. For high-dimensional & high-sparsity setting, we adopt the sub-space optimization methods and row generation methods. We leave these for future implementation.
#'
#' @param L A loss matrix (in the objective) of size \code{1} by \code{G_size} by \code{D[1]} by ... by \code{D[-1]}.
#' @param P A probability matrix (in the objective) of size \code{G_size} by \code{R_sizes[1]} by ... by \code{R_sizes[-1]}.
#' @param L_c2 A loss matrix (in the C2 constraints) of size \code{J} by \code{G_size} by \code{D[1]} by ... by \code{D[-1]}.
#' @param P_c2 A probability matrix (in the C2 constraints) of size \code{G_size} by \code{R_sizes[1]} by ... by \code{R_sizes[-1]}.
#' @param alpha C1 constraints values that the user should specify. Its size should agree with G_size.
#' @param beta C2 constraints values that the user should specify. Its size should agree with J.
#' @param solver The solver solving the formulated LP problem. Default to be \code{PRIMAL}.
#' @param G_size Dimension of G and G_c2. Default to be \code{500}.
#' @param D A list of discrete decision values after each stage. Specified by users.
#' @param R_sizes A list of numbers of discretization at each stage. Default to be \code{c(100, 100)}.
#' @param J Total number of C2 constraints that the user should specify.
#' @return
#' An object with S3 class \code{"maed.2aed"} is returned:
#' \item{lambda}{
#'   The sequence of regularization parameters \code{lambda} obtained in the program.
#' }
#' \item{value}{
#'   The sequence of optimal value of the object function corresponded to the sequence of \code{lambda}.
#' }
#' \item{solver}{
#'   The solver solving the formulated LP problem.
#' }

#' @seealso \code{\link{maed.p}}, \code{\link{maed.l}}, \code{\link{maed.g}}, and \code{\link{maed-package}}.
#' @example
#'
#' @export
maed.2aed.main <- function(L=NULL, P=NULL, L_c2=NULL, P_c2=NULL,
                           alpha=NULL, beta=NULL, solver="primal",
                           G_size=500, D=NULL, R_sizes=c(100,100), J=NULL){

  R1_size <- R_sizes[1]
  R2_size <- R_sizes[2]
  D1_size <- length(D[[1]])
  D2_size <- length(D[[2]])

  vars <- maed.2aed.construct_param(L, P, L_c2, P_c2, alpha, beta,
                                     G_size, R1_size, R2_size, D1_size, D2_size, J)
  c <- vars[[1]]
  A <- vars[[2]]
  b <- vars[[3]]

  est = list()

  if(solver=="primal"){

    library(PRIMAL)

    m = dim(A)[1]
    n = dim(A)[2]

    A = cbind(A,diag(rep(1,m)))
    c = c(c,rep(0,m))
    b_bar = rep(1,m)
    c_bar = rep(0,m+n)
    B_init = seq(n,n+m-1)

    fit <- PSM_solver(A, b, b_bar, c, c_bar, B_init)

    est$lambda = fit$lambda
    est$value = fit$value
    rm(fit)

  }else{

    vars <- maed.2aed.construct_variables(R1_size, R2_size, D1_size, D2_size, G_size, J)
    v_0 <- vars[[1]]
    pi_0 <- vars[[2]]

    fit <- maed.2aed.alg(c, A, b, v_0, pi_0, J, rowgen_mode="direct_sep")

    est$lambda = fit$lambda
    est$value = fit$value
    rm(fit)

  }

  est$solver = solver
  class(est) = "maed.2aed"
  return(est)
}


#-------------------------------------------------------------------------------#
# Package: Multistage Adaptive Enrichment Design                                #
# maed.2aed.construct_param(): construct data and constraint matrices for the   #
#                              LP problem given the component matrices          #
#-------------------------------------------------------------------------------#

#' A helper function for construct data and constraint matrices for the LP problem under 2-stage 2-sub-population adaptive enrichment design.
#'
#' @param L A loss matrix (in the objective) of size \code{1} by \code{G_size} by \code{D[1]} by ... by \code{D[-1]}.
#' @param P A probability matrix (in the objective) of size \code{G_size} by \code{R_sizes[1]} by ... by \code{R_sizes[-1]}.
#' @param L_c2 A loss matrix (in the C2 constraints) of size \code{J} by \code{G_size} by \code{D[1]} by ... by \code{D[-1]}.
#' @param P_c2 A probability matrix (in the C2 constraints) of size \code{G_size} by \code{R_sizes[1]} by ... by \code{R_sizes[-1]}.
#' @param alpha C1 constraints values that the user should specify. Its size should agree with G_size.
#' @param beta C2 constraints values that the user should specify. Its size should agree with J.
#' @param G_size Dimension of G and G_c2. Default to be \code{500}.
#' @param R1_size The number of discretization at first stage. Default to be \code{100}.
#' @param R2_size The number of discretization at second stage. Default to be \code{100}.
#' @param D1_size The number of discrete decision values after the first stage. Specified by users.
#' @param D2_size The number of discrete decision values after the second stage. Specified by users.
#' @param J Total number of C2 constraints. Specified by users.
#' @return
#' A list is returned:
#' \item{c}{
#'   The objective vector of length \code{R1_size}x\code{R2_size}x\code{D1_size}x\code{D2_size}.
#' }
#' \item{A}{
#'   The constraint matrix of size \code{G_size+J} by \code{R1_size}x\code{R2_size}x\code{D1_size}x\code{D2_size}.
#' }
#' \item{b}{
#'   The value matrix of length \code{G_size+J}.
#' }

#' @seealso \code{\link{maed.p}}, \code{\link{maed.l}}, \code{\link{maed.g}}, and \code{\link{maed-package}}.
#' @example
#'
#' @export
maed.2aed.construct_param <- function(L=NULL, P=NULL, L_c2=NULL, P_c2=NULL,
                                      alpha=NULL, beta=c(0.5,0.5),
                                      G_size=500, R1_size=100, R2_size=100, D1_size=NULL, D2_size=NULL, J=NULL){

  # input:
  # L, L_c2: matrix, [1, |G|, |D[1]|, |D[2]|], [J, |G|, |D[1]|, |D[2]|]
  # P, P_c2: matrix, [|G|, |R1|, |R2|], [|G|, |R1|, |R2|]

  n <- D1_size*D2_size*R1_size*R2_size
  A_c1 <- rep(rep(P, D1_size), D2_size)
  dim(A_c1) <- c(G_size, n) #[|G|, |R1|*|R2|*|D[1]|*|D[2]|]

  L <- array(aperm(L, c(1,3,4,2)), dim=c(D1_size*D2_size, G_size)) #[|D1|*|D2|, G]
  L_c2 <- array(aperm(L_c2, c(1,3,4,2)), dim=c(J*D1_size*D2_size, G_size)) #[J*|D1|*|D2|, G]
  P <- array(P, dim=c(G_size,R1_size*R2_size)) #[G, R1 * R2]
  P_c2 <- array(P_c2, dim=c(G_size,R1_size*R2_size)) #[G, R1 * R2]


  c <- .Call("_maed_matmul", L, P) #[|R1|*|R2|*|D[1]|*|D[2]|]
  A_c2 <- .Call("_maed_matmul", L_c2, P_c2)
  dim(A_c2) <- c(J, n) #[J, |R1|*|R2|*|D[1]|*|D[2]|]

  A <- rbind(A_c1, A_c2) #[G_size + J, |R1|*|R2|*|D[1]|*|D[2]|]

  b <- abind(c(alpha, beta), along=2) #[G_size + J, 1]
  return(list(c, A, b))
}


#-------------------------------------------------------------------------------#
# Package: Multistage Adaptive Enrichment Design                                #
# maed.2aed.construct_variables(): initialize decision variables for sub-space  #
#                                  optimization and row generation methods      #
#-------------------------------------------------------------------------------#

#' A helper function for initialize decision variables for the LP problem under 2-stage 2-sub-population adaptive enrichment design.
#'
#' @param R1_size The number of discretization at first stage. Default to be \code{100}.
#' @param R2_size The number of discretization at second stage. Default to be \code{100}.
#' @param D1_size The number of discrete decision values after the first stage. Specified by users.
#' @param D2_size The number of discrete decision values after the second stage. Specified by users.
#' @param G_size Dimension of G and G_c2. Default to be \code{500}.
#' @param J Total number of C2 constraints. Specified by users.
#' @return
#' A list is returned:
#' \item{v}{
#'   The decision variable vector of length \code{R1_size}x\code{R2_size}x\code{D1_size}x\code{D2_size}.
#' }
#' \item{pi}{
#'   The decision variable pi of length \code{G_size+J}.
#' }

#' @seealso \code{\link{maed.p}}, \code{\link{maed.l}}, \code{\link{maed.g}}, and \code{\link{maed-package}}.
#' @example
#'
#' @export
maed.2aed.construct_variables <- function(R1_size=100, R2_size=100, D1_size=NULL, D2_size=NULL,
                                          G_size=500, J=NULL){

  n <- D1_size*D2_size*R1_size*R2_size
  # X [|R1|, |D[1]|]
  x <- matrix(runif(R1_size*D1_size), R1_size, D1_size)
  x <- x/rowSums(x)

  # V [|R1|, |D[1]|, |R2|, |D[2]|]
  v <- rep(rep(x, R2_size), D2_size)
  dim(v) <- c(R1_size*D1_size*R2_size, D2_size)
  v_ <- matrix(runif(n), R1_size*D1_size*R2_size, D2_size)
  v <- v * v_/rowSums(v_)
  dim(v) <- c(n) # [|R1|*|D[1]|*|R2|*|D[2]|]

  # pi
  pi <- -1 * matrix(runif(G_size+J)) #[G_size + J]

  return(list(v, pi))
}

#' @useDynLib maed
#' @importFrom Rcpp sourceCpp
NULL
