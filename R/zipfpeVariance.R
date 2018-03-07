#' Variance of the Zipf-PE distribution.
#'
#' Computes the variance of the Zipf-PE distribution for given values of \eqn{\alpha} and \eqn{\beta}.
#' @param alpha Value of the \eqn{\alpha} parameter (\eqn{\alpha > 3}).
#' @param beta Value of the \eqn{\beta} parameter (\eqn{\beta \in (-\infty, +\infty)}).
#' @param tolerance Tolerance used in the calculations. (default = \eqn{10^{-4}})
#' @return A positive real value corresponding to the variance of the distribution.
#'
#' @details
#' The variance of the distribution only exists for \eqn{\alpha} strictly greater than 3.
#'
#' @examples
#' zipfpeVariance(3.5, 1.3)
#' @seealso \code{\link{zipfpeMoments}}, \code{\link{zipfpeMean}}.
#' @export
zipfpeVariance <- function(alpha, beta, tolerance = 10^(-4)){
  moment1 <- zipfpeMoments(1, alpha, beta, tolerance)
  moment2 <- zipfpeMoments(2, alpha, beta, tolerance)

  return(moment2 - (moment1^2))
}
