#' @export
maed <- function(n_stage=2, n_ppl=2, prop_ppl=c(0.5,0.5), R_sizes=c(100,100),
                 D=list(c(-0.3,-0.1,0.1,0.3)),
                 lambda_0="normal", loss_0=list("default"),
                 lambda_c2="normal", loss_c2=list("default", "default"),
                 G_size=500, alpha=NULL, J=2, beta=NULL,
                 solver="primal", rowgen_mode="direct_sep"){

  library(testit)

  # general
  assert(n_stage == length(R_sizes))
  assert(n_stage == length(D)+1)
  H <- n_ppl+1
  D[[length(D)+1]] <- c(1:2**H)
  J <- length(loss_c2)

  if(n_stage==2 && n_ppl==2)
  {

    G <- maed.g(lambda_0, G_size, n_ppl, prop_ppl) # [|G|, n_ppl]
    vars <- maed.p(G, R_sizes)
    P = vars[[1]] # [|G|, |R1|, |R2|]
    G_idx = vars[[2]] #[|G|, n_stage]
    L <- maed.l(G_idx, D, loss_0, R_sizes) # [1, |G|, |D[1]|, |D[2]|]

    # C2 constraints
    if(J > 0){
      G_c2 <- maed.g(lambda_c2, G_size, n_ppl, prop_ppl) #[|G|, n_ppl]
      vars <- maed.p(G_c2, R_sizes)
      P_c2 = vars[[1]] #[|G|, |R1|, |R2|]
      G_c2_idx = vars[[2]] #[|G|, n_stage]
      L_c2 <- maed.l(G_c2_idx, D, loss_c2, R_sizes) #[J, |G|, |D[1]|, |D[2]|]
    }

    maed.2aed.main(L, P, L_c2, P_c2, alpha, beta, solver,
                  rowgen_mode, G_size, D, R_sizes, J)
  }
}

#' @useDynLib maed
#' @importFrom Rcpp sourceCpp
NULL
