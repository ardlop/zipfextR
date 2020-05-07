#' Moments of the Zipf-Polylog Distribution.
#'
#' General function to compute the k-th moment of the ZipfPolylog distribution for any integer value \eqn{k \geq 1},
#' when it exists. #'
#' For k = 1, this function returns the same value as the \link{zipfpoylogMean} function.
#'
#' @param k Order of the moment to compute.
#' @param alpha Value of the \eqn{\alpha} parameter (\eqn{\alpha > k + 1}).
#' @param beta Value of the \eqn{\beta} parameter (\eqn{\beta \in (-\infty, +\infty)}).
#' @param tolerance Tolerance used in the calculations (default = \eqn{10^{-4}}).
#' @param nSum The number of terms used for computing the Polylogarithm function (default = 1000).
#'
#' @return A positive real value corresponding to the k-th moment of the distribution.
#'
#' @details
#' The k-th moment of the Zipf-Polylog distribution is always finite, but,
#' for \eqn{\alpha >1} and \eqn{\beta = 0} the k-th moment is only finite for all \eqn{\alpha > k + 1}.
#' It is computed by calculating the partial sums of the serie, and stopping when two
#' consecutive partial sums differ less than the \code{tolerance} value.
#' The value of the last partial sum is returned.
#'
#' @examples
#' zipfpolylogMoments(1, 0.2, 0.90)
#' zipfpolylogMoments(3, 4.5, 0.90,  1*10^(-3))
#' @export
zipfpolylogMoments <- function(k, alpha, beta, tolerance = 10^(-4), nSum = 1000){
  if(!is.numeric(k) || !is.numeric(alpha) || !is.numeric(beta) || !is.numeric(tolerance)){
    stop("Wrong input parameters!!")
  }

  if(alpha > 1 && beta == 0 && alpha < k + 1){
    stop(sprintf('Alpha value must be greater than %s.', k + 1))
  }

  if(!k%%1 == 0 || k < 1){
    stop('Wrong moment value!!. You have to provide a positive and integer value.')
  }

  aux <- 1
  x <- 1
  result <- 0
  while(aux > tolerance) {
    pk <- dzipfpolylog(x, alpha, beta, nSum = nSum)
    aux <- x^k * pk
    #aux <- copula::polylog(z = beta, s = (alpha-k), method = 'sum', n.sum = nSum)/ copula::polylog(z = beta, s = (alpha), method = 'sum', n.sum = nSum)
    result <- result + aux
    # print(c(x, pk, aux, result))
    x <- x + 1
  }
  return(result)
}

