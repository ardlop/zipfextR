#' Expected value of the ZipfPolylog distribution.
#'
#' Computes the expected value of the ZipfPolylog distribution for given values of parameters
#' \eqn{\alpha} and \eqn{\beta}.
#'
#' @param alpha Value of the \eqn{\alpha} parameter (\eqn{\alpha > 2}).
#' @param beta Value of the \eqn{\beta} parameter (\eqn{\beta \in (-\infty, +\infty)}).
#' @param tolerance Tolerance used in the calculations (default = \eqn{10^{-4}}).
#'
#' @return A positive real value corresponding to the mean value of the ZipfPolylog distribution.
#'
#' @examples
#' zipfpolylogMean(0.5, 0.8)
#' zipfpolylogMean(2.5, 0.8, 10^(-3))
#' @export
zipfpolylogMean <- function(alpha, beta, tolerance = 10^(-4)){
  return(zipfpolylogMoments(1, alpha, beta, tolerance = tolerance))
}
