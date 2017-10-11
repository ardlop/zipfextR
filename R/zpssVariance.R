#' Variance of the Z-PSS distribution.
#'
#' Computes the variance of the Z-PSS distribution for given values of parameters
#' \eqn{\alpha} and \eqn{\lambda}.
#'
#' @param alpha Value of the \eqn{\alpha} parameter (\eqn{\alpha > 3}).
#' @param lambda Value of the \eqn{\lambda} parameter (\eqn{\lambda \geq 0}).
#' @param isTruncated Logical; if TRUE Use the zero-truncated version of the distribution to calculate the expected value (default = FALSE).
#'
#' @return A positive real value corresponding to the variance of the distribution.
#'
#' @details
#' The variance of the Z-PSS distribution only exists for \eqn{\alpha} values strictly greater than 3.
#' The value is derive from \eqn{Var[Y] = E[N]\, Var[X] + E[X]^2 \, Var[N]} where E[X] and E[N] ares the expected value of the
#' Zipf and the Poisson distributions respectively. In the same way the values of Var[X] and Var[N] stand for the variances of the Zipf and the Poisson
#' distributions. The resulting expression is set to be equal to:
#' \deqn{Var[Y] = \lambda\, \frac{\zeta(\alpha - 2)}{\zeta(\alpha)}}.
#' Particularlly, the variance of the zero-truncated version of the Z-PSS distribution is calculated as:
#' \deqn{Var[Y^{ZT}] = \frac{\lambda\, \zeta(\alpha)\, \zeta(\alpha - 2)\, (1 - e^{-\lambda}) - \lambda^2 \, \zeta(\alpha - 1)^2 \, e^{-\lambda}}{\zeta(\alpha)^2 \, (1 - e^{-\lambda})^2}}
#'
#' @references {
#' Sarabia Alegría, JM. and Gómez Déniz, E. and Vázquez Polo, F. Estadística actuarial: teoría y aplicaciones. Pearson Prentice Hall.
#' }
#' @examples
#' zpssVariance(4.5, 2.3)
#' zpssVariance(4.5, 2.3, TRUE)
#' @export
zpssVariance <- function(alpha, lambda, isTruncated = FALSE){
  if(!is.numeric(alpha) || !is.numeric(lambda)){
    stop("Wrong input parameters!!")
  }

  if(alpha < 3){
    stop('The alpha parameter must be greater than 2.')
  }
  zeta_a <- VGAM::zeta(alpha)
  if(!isTruncated){
    return(lambda*VGAM::zeta(alpha - 2)/zeta_a)
  }
  (lambda*zeta_a*VGAM::zeta(alpha - 2)*(1-exp(-lambda)) - lambda^2 * VGAM::zeta(alpha - 1)^2 * exp(-lambda))/(zeta_a^2 * (1 - exp(-lambda))^2)
}
