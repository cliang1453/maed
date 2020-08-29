#-------------------------------------------------------------------------------#
# Package: Multistage Adaptive Enrichment Design                                #
# maed.2aed.alg(): The interface for maed.2aed.subopt, maed.2aed.rowgen,        #
#                  maed.2aed.fitwrapper                                         #
#-------------------------------------------------------------------------------#

#' The function for using sub-space optimization and row-generation method to solve the formulated LP for 2-stage, 2-sub-population adaptive enrichment design problem.
#'
#' @param C The objective vector of length \code{R1_size}x\code{R2_size}x\code{D1_size}x\code{D2_size}.
#' @param A The constraint matrix of size \code{G_size+J} by \code{R1_size}x\code{R2_size}x\code{D1_size}x\code{D2_size}.
#' @param b The value matrix of length \code{G_size+J}.
#' @param v_0 The decision variable vector of length \code{R1_size}x\code{R2_size}x\code{D1_size}x\code{D2_size}.
#' @param pi_0 The decision variable pi of length \code{G_size+J}.
#' @param J Total number of C2 constraints that the user should specify.
#' @param T_max Max number of interation. Default to be \code{1000}.
#' @param rowgen_mode Mode of performing row generation. choices: \code{direct_sep}(direct seperation) or \code{ada_disc}adaptive discretization. Default to be the former.
#' @return
#' An object with S3 class \code{"maed.2aed.fitwrapper"} is returned.

#' @seealso \code{\link{maed.p}}, \code{\link{maed.l}}, \code{\link{maed.g}}, and \code{\link{maed-package}}.
#' @example
#'
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

  fit <- maed.2aed.fitwrapper(v, pi_I)
  return(fit)
}

#' @export
maed.2aed.subopt <- function(C=NULL, A=NULL, b=NULL, v_0=NULL, pi_0=NULL){

  stop("Subspace optimization method is currently under-construction.")
  # L-infinity bundle trust region method
}

#' @export
maed.2aed.rowgen <- function(C=NULL, A=NULL, b=NULL, v=NULL, pi_I=NULL,
                             I=NULL, rowgen_mode="direct_sep"){

  stop("Methods of row generation is currently under-construction.")
  # direct separation method
  # adaptive discretization method
}

#' @export
maed.2aed.fitwrapper <- function(v=NULL, pi_I=NULL){

  stop("Class object wrapper is currently under-construction.")
  # wrap the return object
}
