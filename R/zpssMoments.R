#' Distribution Moments.
#'
#' General function to compute the k-th moment of the Z-PSS distribution. Note that the k-th moment
#' exists if and only if  \eqn{\alpha > k + 1}.
#'
#' @param k Order of the moment to compute.
#' @param alpha Value of the \eqn{\alpha} parameter (\eqn{\alpha > k + 1}).
#' @param lambda Value of the \eqn{\lambda} parameter (\eqn{\lambda \geq 0}).
#' @param isTruncated Logical; if TRUE, the truncated version of the distribution is returned.
#' @param tolerance Tolerance used in the calculations (default = \eqn{10^{-4}}).
#'
#' @return A positive real value corresponding to the k-th moment of the distribution.
#'
#' @details
#' The k-th moment of the Z-PSS distribution is finite for \eqn{\alpha} values strictly greater than \eqn{k + 1}.
#' It is computed calculating the partial sums of the serie, and it stops when two consecutive
#' partial sums differs less than the \code{tolerance} value. The last partial sum is returned.
#'
#' @examples
#' zpssMoments(3, 4.5, 1.3)
#' zpssMoments(3, 4.5, 1.3,  1*10^(-3))
#' @export
zpssMoments <- function(k, alpha, lambda, isTruncated = FALSE, tolerance = 10^(-4)){
  if(!is.numeric(k) || !is.numeric(alpha) || !is.numeric(lambda) || !is.numeric(tolerance)){
    stop("Wrong input parameters!!")
  }

  if(alpha < k + 1){
    stop(sprintf('Alpha value must be greater than %s.', k + 1))
  }

  if(!k%%1 == 0 || k < 1){
    stop('Wrong moment value!!. You have to provide a possitive and integer value.')
  }

  aux <- 1
  # x <- if(isTruncated) 1 else 0
  x <- 1
  result <- 0

  # while(aux > tolerance) {
  #   px <- dzpss(x, alpha, lambda, isTruncated = isTruncated)
  #   print(px)
  #   aux <- x^k * px
  #   print(aux)
  #   result <- result + aux
  #   x <- x + 1
  # }

  while(aux > tolerance) {
      px <- dzipfpss(x, alpha, lambda, isTruncated = isTruncated)
      # print(px)
      aux <- x^k * px
      # print(aux)
      result <- result + aux
      x <- x +1
  }

  return(result)
}
