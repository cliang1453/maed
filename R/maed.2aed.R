#' @export
maed.2aed.main <- function(L=NULL, P=NULL, L_c2=NULL, P_c2=NULL,
                           alpha=NULL, beta=NULL, solver="primal",
                           rowgen_mode="direct_sep",
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

  vars <- maed.2aed.construct_variables(R1_size, R2_size, D1_size, D2_size, G_size, J)
  v_0 <- vars[[1]]
  pi_0 <- vars[[2]]

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
    print(fit$lambda)
    print(fit$value)

  }else{
    vars <- maed.2aed.alg(c, A, b, v_0, pi_0, J, rowgen_mode)
    v <- vars[1]
    pi <- vars[2]
  }

}

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
