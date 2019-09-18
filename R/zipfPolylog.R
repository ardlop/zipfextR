#' The Zipf-Polylog Distribution (Zipf-Polylog).
#'
#' Probability mass function of the Zipf-Polylog distribution with parameters \eqn{\alpha} and \eqn{\beta}.
#' The support of the Zipf-Polylog distribution are the strictly positive integer numbers large or equal
#' than one.
#'
#' @name zipfPolylog
#' @aliases dzipfpolylog
#'
#' @param x Vector of positive integer values.
#' @param alpha Value of the \eqn{\alpha} parameter (\eqn{\alpha > 1} ).
#' @param beta Value of the \eqn{\beta} parameter (\eqn{\beta > 0} ).
#' @param nSum The number of terms used for computing the Polilogarithm function [Default = 1000].
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @details The \emph{probability mass function} at a positive integer value \eqn{x} of the Zipf-Polylog distribution with
#' parameters \eqn{\alpha} and \eqn{\beta} is computed as follows:
#'
#' @return {
#' \code{dzipfpolylog} gives the probability mass function
#' }
#' @importFrom copula polylog
#' @examples
#' dzipfpolylog(1:10, 1.61, 0.98)
#'
NULL
#> NULL

.prec.zipfpolylog.checkXvalue <- function(x){
  if(!is.numeric(x) || any(x < 1) || any(x%%1 != 0)) {
    stop('The x value is not included into the support of the distribution.')
  }
}

.prec.zipfpolylog.checkparams <- function(alpha, beta){
  if(!is.numeric(beta) | beta < 0 | beta > 1){
    stop('Incorrect beta parameter. You should provide a numeric value.')
  }

  # if(!is.numeric(alpha) || (beta == 1 && alpha < 1) || (beta < 1 && alpha < 0)){
  #   stop('Incorrect alpha parameter. This parameter should be greater than one.')
  # }
}

.dzipfPolylog.default <- function(x, alpha, beta, nSum){
  return((beta^x * x^(-alpha))/copula::polylog(z = beta, s = alpha, method = 'sum', n.sum = nSum))
}

#' @rdname zipfPolylog
#' @export
dzipfpolylog <- function(x, alpha, beta, log = FALSE, nSum = 1000){
  .prec.zipfpolylog.checkXvalue(x)
  .prec.zipfpolylog.checkparams(alpha, beta)

  values <- sapply(x, .dzipfPolylog.default, alpha = alpha, beta = beta, nSum = nSum)
  if(log){
    return(log(values))
  }
  return(values)
}
