#' Distribution Moments.
#'
#' General function to compute the k-th moment of the Zipf-PSS distribution for any integer value \eqn{k \geq 1},
#' when it exists. The k-th moment exists if and only if  \eqn{\alpha > k + 1}.
#'
#' @param k Order of the moment to compute.
#' @param alpha Value of the \eqn{\alpha} parameter (\eqn{\alpha > k + 1}).
#' @param lambda Value of the \eqn{\lambda} parameter (\eqn{\lambda > 0}).
#' @param isTruncated Logical; if TRUE, the truncated version of the distribution is returned.
#' @param tolerance Tolerance used in the calculations (default = \eqn{10^{-4}}).
#'
#' @return A positive real value corresponding to the k-th moment of the distribution.
#'
#' @details
#' The k-th moment of the Zipf-PSS distribution is finite for \eqn{\alpha} values
#' strictly greater than \eqn{k + 1}.
#' It is computed by calculating the partial sums of the serie, and stopping when two
#' consecutive partial sums differ less than the \code{tolerance} value.
#' The value of the last partial sum is returned.
#'
#' @examples
#' zipfpssMoments(1, 2.5, 2.3)
#' zipfpssMoments(1, 2.5, 2.3, TRUE)
#' @export
zipfpssMoments <- function(k, alpha, lambda, isTruncated = FALSE, tolerance = 10^(-4)){
  if(!is.numeric(k) || !is.numeric(alpha) || !is.numeric(lambda) || !is.numeric(tolerance) || !is.logical(isTruncated)){
    stop("Wrong input parameters!!")
  }

  if(alpha < k + 1){
    stop(sprintf('Alpha value must be greater than %s.', k + 1))
  }

  if(!k%%1 == 0 || k < 1){
    stop('Wrong moment value!!. You have to provide a possitive and integer value.')
  }

  probs <- dzipfpss(1:1000, alpha, lambda, isTruncated = isTruncated)
  result <- sum(((1:1000)^k)*probs)
  aux <- 1000*probs[length(probs)]
  x <- 1001
  while(aux > tolerance) {
      probs <- .getProbs(x, alpha, lambda, isTruncated=isTruncated, probs)
      aux <- x^k * probs[if(isTruncated) x else (x+1)]
      result <- result + aux
      x <- x + 1
  }
  return(result)
}
