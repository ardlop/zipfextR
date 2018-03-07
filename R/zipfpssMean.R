#' Expected value of the Zipf-PSS distribution.
#'
#' Computes the expected value of the Zipf-PSS distribution for given values of parameters
#' \eqn{\alpha} and \eqn{\lambda}.
#'
#' @param alpha Value of the \eqn{\alpha} parameter (\eqn{\alpha > 2}).
#' @param lambda Value of the \eqn{\lambda} parameter (\eqn{\lambda > 0}).
#' @param isTruncated Logical; if TRUE Use the zero-truncated version of the distribution to calculate the expected value (default = FALSE).
#'
#' @return A positive real value corresponding to the mean value of the distribution.
#'
#' @details
#' The expected value of the Zipf-PSS distribution only exists for \eqn{\alpha} values strictly
#' greater than 2. The value is obtained from the \emph{law of total expectation} that says that: \deqn{E[Y] = E[N]\, E[X],}
#' where E[X] is the mean value of the Zipf distribution and E[N] is the expected value of a Poisson one.
#' From where one has that:
#' \deqn{E[Y] = \lambda\, \frac{\zeta(\alpha - 1)}{\zeta(\alpha)}}
#'
#' Particularlly, if one is working with the zero-truncated version of the Zipf-PSS distribution.
#' This values is computed as:
#' \deqn{E[Y^{ZT}] = \frac{\lambda\, \zeta(\alpha - 1)}{\zeta(\alpha)\, (1 - e^{-\lambda})}}
#'
#' @references {
#' Sarabia Alegría, J. M., Gómez Déniz, E. M. I. L. I. O., & Vázquez Polo, F. (2007).
#' Estadística actuarial: teoría y aplicaciones. Pearson Prentice Hall.
#' }
#' @examples
#' zipfpssMean(2.5, 1.3)
#' zipfpssMean(2.5, 1.3, TRUE)
#' @export
zipfpssMean <- function(alpha, lambda, isTruncated = FALSE){
  if(!is.numeric(alpha) || !is.numeric(lambda)){
    stop("Wrong input parameters!!")
  }

  if(alpha < 2){
    stop('The alpha parameter must be greater than 2.')
  }

  meanVal <-lambda * VGAM::zeta(alpha - 1)/VGAM::zeta(alpha)

  if(!isTruncated){
    return(meanVal)
  }
  return(meanVal/(1 - exp(-lambda)))
}
