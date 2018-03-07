#' Distribution Moments.
#'
#' General function to compute the k-th moment of the Zipf-PE distribution for any integer value \eqn{k \geq 1},
#' when it exists. The k-th moment exists if and only if  \eqn{\alpha > k + 1}.
#' For k = 1, this function returns the same value as the \link{zipfpeMean} function.
#'
#' @param k Order of the moment to compute.
#' @param alpha Value of the \eqn{\alpha} parameter (\eqn{\alpha > k + 1}).
#' @param beta Value of the \eqn{\beta} parameter (\eqn{\beta \in (-\infty, +\infty)}).
#' @param tolerance Tolerance used in the calculations (default = \eqn{10^{-4}}).
#'
#' @return A positive real value corresponding to the k-th moment of the distribution.
#'
#' @details
#' The k-th moment of the Zipf-PE distribution is finite for \eqn{\alpha} values strictly greater than \eqn{k + 1}.
#' It is computed by calculating the partial sums of the serie, and stopping when two
#' consecutive partial sums differ less than the \code{tolerance} value.
#' The value of the last partial sum is returned.
#'
#' @examples
#' zipfpeMoments(3, 4.5, 1.3)
#' zipfpeMoments(3, 4.5, 1.3,  1*10^(-3))
#' @export
zipfpeMoments <- function(k, alpha, beta, tolerance = 10^(-4)){
  if(!is.numeric(k) || !is.numeric(alpha) || !is.numeric(beta) || !is.numeric(tolerance)){
    stop("Wrong input parameters!!")
  }

  if(alpha < k + 1){
    stop(sprintf('Alpha value must be greater than %s.', k + 1))
  }

  if(!k%%1 == 0 || k < 1){
    stop('Wrong moment value!!. You have to provide a possitive and integer value.')
  }

  aux <- 1
  x <- 1
  result <- 0

  while(aux > tolerance) {
    pk <- dzipfpe(x, alpha, beta)
    aux <- x^k * pk
    result <- result + aux
    x <- x + 1
  }

  return(result)
}

