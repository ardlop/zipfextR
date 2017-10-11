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
#' The expected value of the MOEZipf distribution only exists for \eqn{\alpha} values strictly greater than 2.
#' In this case, if Y is a random variable that follows a MOEZipf distribution with parameters \eqn{\alpha}
#' and \eqn{\beta}, the expected value is computed as:
#' \deqn{E(Y) = \sum_{x = 1} ^\infty \frac{\beta \zeta(\alpha) x^{-\alpha + 1}}{[\zeta(\alpha) - \bar{\beta}\zeta(\alpha, x)][\zeta(\alpha) - \bar{\beta}\zeta(\alpha, x + 1)]}\,, \alpha > 2\,,  \beta > 0}
#'
#' The mean is computed calculating the partial sums of the serie, and it stops when two
#' consecutive partial sums differs less than the \code{tolerance} value.
#' The last partial sum is returned.
#'
#' @examples
#' moezipfMean(2.5, 1.3)
#' moezipfMean(2.5, 1.3, 10^(-3))
#' @export
moezipfMean <- function(alpha, beta, tolerance = 10^(-4)){
  return(moezipfMoments(1, alpha, beta, tolerance = tolerance))
}

