#' @export
maed.2aed.alg <- function(C=NULL, A=NULL, b=NULL, v_0=NULL, pi_0=NULL, J=NULL,
                          T_max=1000, rowgen_mode="direct_sep"){

  t <- 1
  I <- c(dim(A)[1]-J:J)


  while (t < T_max){

    # subspace optimization
    vars <- maed.2aed.subopt(C[I], A[I,], b, v_0, pi_0[I])
    v <- vars[1]
    pi_I <- vars[2]

    # row generation
    I_new <- maed.2aed.rowgen(C, A, b, v, pi_I, I, rowgen_mode)
    if (length(I_new) == length(I)){
      break
    }

    I <- I_new
    t = t+1
  }

  return(list(v, pi_I))

}

#' @export
maed.2aed.subopt <- function(C=NULL, A=NULL, b=NULL, v_0=NULL, pi_0=NULL){

  # L-infinity bundle trust region method

}

#' @export
maed.2aed.rowgen <- function(C=NULL, A=NULL, b=NULL, v=NULL, pi_I=NULL,
                             I=NULL, rowgen_mode="direct_sep"){

  # direct separation method

  # adaptive discretization method
}
