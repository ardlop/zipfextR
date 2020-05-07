#' Variance of the ZipfPolylog distribution.
#'
#' Computes the variance of the ZipfPolylog distribution for given values of \eqn{\alpha} and \eqn{\beta}.
#' @param alpha Value of the \eqn{\alpha} parameter (\eqn{\alpha > 3}).
#' @param beta Value of the \eqn{\beta} parameter (\eqn{\beta \in (-\infty, +\infty)}).
#' @param tolerance Tolerance used in the calculations. (default = \eqn{10^{-4}})
#' @return A positive real value corresponding to the variance of the distribution.
#'
#' @details
#' The variance of the distribution only exists for \eqn{\alpha} strictly greater than 3.
#'
#' @examples
#' zipfpoylogVariance(0.5, 0.75)
#' @seealso \code{\link{zipfpolylogMoments}}, \code{\link{zipfpoylogMean}}.
#' @export
zipfpoylogVariance <- function(alpha, beta, tolerance = 10^(-4)){
  moment1 <- zipfpolylogMoments(1, alpha, beta, tolerance)
  moment2 <- zipfpolylogMoments(2, alpha, beta, tolerance)

  return(moment2 - (moment1^2))
}
