#' @export
hello <- function() {
  print("Hello, world!")
}

#' @export
add <- function(x, y){
  print(x + y)
}

#' @export
take_mean <- function(x){
  out <- .Call("_maed_meanC", x)
}

#' @useDynLib maed
#' @importFrom Rcpp sourceCpp
NULL
