#' Variance of the Zipf-PSS distribution.
#'
#' Computes the variance of the Zipf-PSS distribution for given values of parameters
#' \eqn{\alpha} and \eqn{\lambda}.
#'
#' @param alpha Value of the \eqn{\alpha} parameter (\eqn{\alpha > 3}).
#' @param lambda Value of the \eqn{\lambda} parameter (\eqn{\lambda > 0}).
#' @param isTruncated Logical; if TRUE Use the zero-truncated version of the distribution to calculate the expected value (default = FALSE).
#'
#' @return A positive real value corresponding to the variance of the distribution.
#'
#' @details
#' The variance of the Zipf-PSS distribution only exists for \eqn{\alpha} values strictly greater than 3.
#' The value is obtained from the \emph{law of total variance} that says that: \deqn{Var[Y] = E[N]\, Var[X] + E[X]^2 \, Var[N],}
#' where X follows a Zipf distribution with parameter \eqn{\alpha}, and N follows a Poisson distribution with
#' parameter \eqn{\lambda}. From where one has that:
#'
#' \deqn{Var[Y] = \lambda\, \frac{\zeta(\alpha - 2)}{\zeta(\alpha)}}
#' Particularlly, if one is working with the zero-truncated version of the Zipf-PSS distribution.
#' This values is computed as:
#' \deqn{Var[Y^{ZT}] = \frac{\lambda\, \zeta(\alpha)\, \zeta(\alpha - 2)\, (1 - e^{-\lambda}) - \lambda^2 \, \zeta(\alpha - 1)^2 \, e^{-\lambda}}{\zeta(\alpha)^2 \, (1 - e^{-\lambda})^2}}
#'
#' @references {
#' Sarabia Alegría, JM. and Gómez Déniz, E. and Vázquez Polo, F. Estadística actuarial: teoría y aplicaciones. Pearson Prentice Hall.
#' }
#' @examples
#' zipfpssVariance(4.5, 2.3)
#' zipfpssVariance(4.5, 2.3, TRUE)
#' @export
zipfpssVariance <- function(alpha, lambda, isTruncated = FALSE){
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

  fc <- lambda/(zeta_a*(1 - exp(-lambda)))
  term1 <- (lambda * (VGAM::zeta(alpha - 1)^2) * exp(-lambda))/(zeta_a * (1 - exp(-lambda)))

  (fc * (VGAM::zeta(alpha - 2) - term1))


  # (lambda*zeta_a*VGAM::zeta(alpha - 2)*(1-exp(-lambda)) - lambda^2 * VGAM::zeta(alpha - 1)^2 * exp(-lambda))/(zeta_a^2 * (1 - exp(-lambda))^2)
}
