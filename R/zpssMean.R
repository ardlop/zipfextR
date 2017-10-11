#' Expected value of the Z-PSS distribution.
#'
#' Computes the expected value of the Z-PSS distribution for given values of parameters
#' \eqn{\alpha} and \eqn{\lambda}.
#'
#' @param alpha Value of the \eqn{\alpha} parameter (\eqn{\alpha > 2}).
#' @param lambda Value of the \eqn{\lambda} parameter (\eqn{\lambda > 0}).
#' @param isTruncated Logical; if TRUE Use the zero-truncated version of the distribution to calculate the expected value (default = FALSE).
#'
#' @return A positive real value corresponding to the mean value of the distribution.
#'
#' @details
#' The expected value of the Z-PSS distribution only exists for \eqn{\alpha} values strictly greater than 2.
#' The value is derive from \eqn{E[Y] = E[N]\, E[X]} where E[X] is the mean value of the Zipf distribution and
#' E[N] is the expected value of a Poisson one. The resulting expression is set to be equal to:
#'\deqn{E[Y] = \lambda\, \frac{\zeta(\alpha - 1)}{\zeta(\alpha)}}.
#'Particularlly, if one is dealing with the zero-truncated version of the Z-PSS distribution. This values es calculated as:
#'\deqn{E[Y^{ZT}] = \frac{\lambda\, \zeta(\alpha - 1)}{\zeta(\alpha)\, (1 - e^{-\lambda})}}
#'
#' @references {
#' Sarabia Alegría, JM. and Gómez Déniz, E. and Vázquez Polo, F. Estadística actuarial: teoría y aplicaciones. Pearson Prentice Hall.
#' }
#' @examples
#' zpssMean(2.5, 1.3)
#' zpssMean(2.5, 1.3, TRUE)
#' @export
zpssMean <- function(alpha, lambda, isTruncated = FALSE){
  if(!is.numeric(alpha) || !is.numeric(lambda)){
    stop("Wrong input parameters!!")
  }

  if(alpha < 2){
    stop('The alpha parameter must be greater than 2.')
  }
  zeta_a <- VGAM::zeta(alpha)
  if(!isTruncated){
    return(lambda * (VGAM::zeta(alpha - 1)/zeta_a))
  }
  lambda*VGAM::zeta(alpha - 1)/zeta_a*(1 - exp(-lambda))
}
