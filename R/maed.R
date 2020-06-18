#' @export
maed <- function(n_stage=2, n_ppl=2, prop_ppl=c(0.5, 0.5), R_sizes=c(100,100), D=c(c(0,1,2,3)),
                 lambda_0="normal", loss_0=NULL, lambda_c2=c("normal", "normal"), loss_c2=NULL,
                 G_size=500, alpha=0.8, J=2, beta=c(0.5,0.5), rowgen_mode="direct_sep"){

  # general
  assert(n_stage == length(R_sizes))
  assert(n_stage == length(D)+1)
  H <- n_ppl + 1
  D <- c(D, seq(0,2**H))
  J <- length(lambda_c2)

  if(n_stage==2 && n_ppl==2)
  {

    # objective
    G <- maed.g(lambda_0, G_size, H) # \delta \in G
    P <- maed.p(G, R_sizes) # [|G|, |R1|, |R2|]
    L <- maed.l(G, D, loss_0) # [1, |G|, |D[1]|, |D[2]|]

    # C2 constraints
    if(J > 0){
      G_c2 <- maed.g(lambda_c2, G_size, H)
      P_c2 <- maed.p(G_c2, R_sizes) #[|G|, |R1|, |R2|]
      L_c2 <- maed.l(G_c2, D, loss_c2) #[|J|, |G|, |D[1]|, |D[2]|]
    }

    ret <- maed.2aed.main(L, P, L_c2, P_c2, alpha, beta,
                          rowgen_mode, G_size, D, R_sizes, J)
  }
}

#' @useDynLib maed
#' @importFrom Rcpp sourceCpp
NULL
