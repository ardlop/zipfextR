.moment <- function(x, k, zeta_alpha, alpha, beta){
  (zeta_alpha * beta * x^(-alpha + k))/((zeta_alpha - (1 - beta) * .zeta_x(alpha, x)) * (zeta_alpha - (1 - beta) * .zeta_x(alpha, x + 1)))
}


#' Distribution Moments.
#'
#' General function to compute the k-th moment of the MOEZipf distribution, for any integer value \eqn{k \geq 1}
#' when it exists. Note that the k-th moment exists if and only if  \eqn{\alpha > k + 1}.
#' When k = 1, this function returns the same value as the
#' \link{moezipfMean} function.
#'
#' @param k Order of the moment to compute.
#' @param alpha Value of the \eqn{\alpha} parameter (\eqn{\alpha > k + 1}).
#' @param beta Value of the \eqn{\beta} parameter (\eqn{\beta > 0}).
#' @param tolerance Tolerance used in the calculations (default = \eqn{10^{-4}}).
#'
#' @return A positive real value corresponding to the k-th moment of the distribution.
#'
#' @details
#' The k-th moment of the MOEZipf distribution is finite for \eqn{\alpha} values strictly greater than \eqn{k + 1}.
#' For a random variable Y that follows a MOEZipf distribution with parameters \eqn{\alpha} and \eqn{\beta},
#' the k-th moment is equal to:
#'
#' \deqn{E(Y^k) = \sum_{x = 1} ^\infty \frac{\beta \zeta(\alpha) x^{-\alpha + k}}{[\zeta(\alpha) - \bar{\beta}\zeta(\alpha, x)][\zeta(\alpha) - \bar{\beta}\zeta(\alpha, x + 1)]}\,, \alpha \geq k + 1\,, \beta > 0}
#'
#' The k-th moment is computed calculating the partial sums of the serie, and stopping when two
#' consecutive partial sums differs less than the \code{tolerance} value.
#' The value of the last partial sum is returned.
#'
#' @examples
#' moezipfMoments(3, 4.5, 1.3)
#' moezipfMoments(3, 4.5, 1.3,  1*10^(-3))
#' @export
moezipfMoments <- function(k, alpha, beta, tolerance = 10^(-4)){
  if(!is.numeric(k) || !is.numeric(alpha) || !is.numeric(beta) || !is.numeric(tolerance)){
    stop("Wrong input parameters!!")
  }

  if(alpha < k + 1){
    stop(sprintf('Alpha value must be greater than %s.', k + 1))
  }

  if(!k%%1 == 0 || k < 1){
    stop('Wrong moment value!!. You have to provide a possitive and integer value.')
  }

  aux <- 1
  x <- 1
  result <- 0
  zeta_alpha <- VGAM::zeta(alpha)

  while(aux > tolerance) {
    aux <- sapply(x, .moment, k = k, zeta_alpha = zeta_alpha, alpha = alpha, beta = beta)
    result <- result + aux
    x <- x + 1
  }

  return(result)
}

