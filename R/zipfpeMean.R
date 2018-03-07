#' Expected value of the Zipf-PE distribution.
#'
#' Computes the expected value of the Zipf-PE distribution for given values of parameters
#' \eqn{\alpha} and \eqn{\beta}.
#'
#' @param alpha Value of the \eqn{\alpha} parameter (\eqn{\alpha > 2}).
#' @param beta Value of the \eqn{\beta} parameter (\eqn{\beta \in (-\infty, +\infty)}).
#' @param tolerance Tolerance used in the calculations (default = \eqn{10^{-4}}).
#'
#' @return A positive real value corresponding to the mean value of the Zipf-PE distribution.
#'
#' @details
#' The mean of the distribution only exists for \eqn{\alpha} strictly greater than 2.
#' It is computed by calculating the partial sums of the serie, and stopping when two
#' consecutive partial sums differ less than the \code{tolerance} value.
#' The value of the last partial sum is returned.
#'
#' @examples
#' zipfpeMean(2.5, 1.3)
#' zipfpeMean(2.5, 1.3, 10^(-3))
#' @export
zipfpeMean <- function(alpha, beta, tolerance = 10^(-4)){
  return(zipfpeMoments(1, alpha, beta, tolerance = tolerance))
}
