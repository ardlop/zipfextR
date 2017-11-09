#' Expected value of the ZPE distribution.
#'
#' Computes the expected value of the ZPE distribution for given values of parameters
#' \eqn{\alpha} and \eqn{\beta}.
#'
#' @param alpha Value of the \eqn{\alpha} parameter (\eqn{\alpha > 2}).
#' @param beta Value of the \eqn{\beta} parameter (\eqn{\beta \in [-\infty, +\infty]}).
#' @param tolerance Tolerance used in the calculations (default = \eqn{10^{-4}}).
#'
#' @return A positive real value corresponding to the mean value of the ZPE distribution.
#'
#' @details
#' The mean of the distribution only exists for \eqn{\alpha} strictly greater than 2. It is computed
#' calculating the partial sums of the serie, and it stops when two consecutive partial
#' sums differs less than the \code{tolerance} value. The last partial sum is returned.
#'
#' @examples
#' zpeMean(2.5, 1.3)
#' zpeMean(2.5, 1.3, 10^(-3))
#' @export
zpeMean <- function(alpha, beta, tolerance = 10^(-4)){
  return(zpeMoments(1, alpha, beta, tolerance = tolerance))
}
