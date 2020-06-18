#' @export
maed.2aed.main <- function(L=NULL, P=NULL, L_c2=NULL, P_c2=NULL,
                           alpha=0.8, beta=c(0.5,0.5), rowgen_mode="direct_sep",
                           G_size=500, D=NULL, R_sizes=c(100,100), J=NULL){

  R1_size <- R_sizes[1]
  R2_size <- R_sizes[2]
  D1_size <- length(D[1])
  D2_size <- length(D[2])

  param <- maed.2aed.construct_param(L, P, L_c2, P_c2, alpha, beta,
                                     G_size, R1_size, R2_size, D1_size, D2_size, J)
  C <- param[1]
  A <- param[2]
  b <- param[3]

  vars <- maed.2aed.construct_variables(R1_size, R2_size, D1_size, D2_size, G_size, J)
  v_0 <- vars[1]
  pi_0 <- vars[2]

  vars <- maed.2aed.alg(C, A, b, v_0, pi_0, J, rowgen_mode)
  v <- vars[1]
  pi <- vars[2]

  return(list(v, C %*% v))
}

#' @export
maed.2aed.construct_param <- function(L=NULL, P=NULL, L_c2=NULL, P_c2=NULL,
                                      alpha=0.8, beta=c(0.5,0.5),
                                      G_size=500, R1_size=100, R2_size=100, D1_size=NULL, D2_size=NULL, J=NULL){

  # input:
  # L, L_c2: matrix, [|G|, |R1|, |R2|]
  # P, P_c2: matrix, [|J|, |G|, |D[1]|, |D[2]|]

  n = D1_size * D2_size * R1_size * R2_size
  A_c1 <- rep(rep(P, R1_size), R2_size)
  dim(A_c1) <- c(G_size, n) #[|G|, |R1|*|R2|*|D[1]|*|D[2]|]

  dim(L) <- c(R1_size * R2_size, G_size)
  dim(L_c2) <- c(R1_size * R2_size, G_size)
  dim(P) <- c(G_size, D1_size * D2_size)
  dim(P_c2) <- c(G_size, J * D1_size * D2_size)

  C <- L %*% P # rcpp
  dim(C) <- c(1, n) #[1, |R1|*|R2|*|D[1]|*|D[2]|]

  A_c2 <- L_c2 %*% P_c2 # rcpp
  dim(A_c2) <- c(J, n) #[J, |R1|*|R2|*|D[1]|*|D[2]|]

  A <- dim(rbind(A_c1, A_c2)) #[G_size + J, |R1|*|R2|*|D[1]|*|D[2]|]
  b <- c(rep(alpha, G_size), beta) #[G_size + J]

  return(list(C, A, b))
}

#' @export
maed.2aed.construct_variables <- function(R1_size=100, R2_size=100, D1_size=NULL, D2_size=NULL,
                                          G_size=500, J=NULL){

  # X [|R1|, |D[1]|]
  x <- matrix(runif(R1_size*D1_size), R1_size, D1_size)
  x <- x/rowSums(x)

  # V [|R1|, |D[1]|, |R2|, |D[2]|]
  v <- rep(rep(x, R2_size), D2_size)
  dim(v) <- c(R1_size*R2_size*D1_size, D2_size)
  v_ <- matrix(runif(R1_size*R2_size*D1_size*D2_size), R1_size*D1_size*R2_size, D2_size)
  v <- v * v_/rowSums(v_)
  dim(v) <- c(R1_size*R2_size*D1_size*D2_size) # [|R1|*|D[1]|*|R2|*|D[2]|]

  # pi
  pi <- -1 * matrix(runif(G_size*J)) #[G_size + J]

  return(list(v, pi))
}

# out <- .Call("_maed_meanC", x)

#' @useDynLib maed
#' @importFrom Rcpp sourceCpp
NULL
