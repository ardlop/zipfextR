#' Expected value.
#'
#' Computes the expected value of the MOEZipf distribution for given values of parameters
#' \eqn{\alpha} and \eqn{\beta}.
#'
#' @param alpha Value of the \eqn{\alpha} parameter (\eqn{\alpha > 2}).
#' @param beta Value of the \eqn{\beta} parameter (\eqn{\beta > 0}).
#' @param tolerance Tolerance used in the calculations (default = \eqn{10^{-4}}).
#'
#' @return A positive real value corresponding to the mean value of the distribution.
#'
#' @details
#' The mean of the distribution only exists for \eqn{\alpha} strictly greater than 2.
#' It is computed by calculating the partial sums of the serie, and stopping when two
#' consecutive partial sums differ less than the \code{tolerance} value.
#' The value of the last partial sum is returned.
#'
#' @examples
#' moezipfMean(2.5, 1.3)
#' moezipfMean(2.5, 1.3, 10^(-3))
#' @export
moezipfMean <- function(alpha, beta, tolerance = 10^(-4)){
  return(moezipfMoments(1, alpha, beta, tolerance = tolerance))
}

