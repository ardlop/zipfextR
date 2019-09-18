
#' The Zero Inflated Zipf-Poisson Stop Sum Distribution (ZI Zipf-PSS).
#'
#' Probability mass function for the zero inflated Zipf-PSS distribution with parameters \eqn{\alpha}, \eqn{\lambda} and \eqn{w}.
#' The support of thezero inflated Zipf-PSS distribution are the positive integer numbers including the zero value.
#'
#' @name zi_zipfpss
#' @aliases d_zi_zipfpss
#'
#' @param x Vector of positive integer values.
#' @param alpha Value of the \eqn{\alpha} parameter (\eqn{\alpha > 1} ).
#' @param lambda Value of the \eqn{\lambda} parameter (\eqn{\lambda > 0} ).
#' @param w Value of the \eqn{w} parameter (0 < \eqn{w < 1} ).
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#'
#' @details
#' The support of the \eqn{\lambda} parameter increases when the distribution is truncated at zero being
#' \eqn{\lambda \geq 0}. It has been proved that when \eqn{\lambda = 0} one has the degenerated version of the distribution at one.
#'
#' @references {
#' Panjer, H. H. (1981). Recursive evaluation of a family of compound
#' distributions. ASTIN Bulletin: The Journal of the IAA, 12(1), 22-26.
#'
#' Sundt, B., & Jewell, W. S. (1981). Further results on recursive evaluation of
#' compound distributions. ASTIN Bulletin: The Journal of the IAA, 12(1), 27-39.
#' }
NULL
#> NULL

.prec.zi_zipfpss.checkparams <- function(alpha, lambda, w){
  if(!is.numeric(alpha) | alpha <= 1){
    stop('Incorrect alpha parameter. This parameter should be greater than one.')
  }

  if(!is.numeric(lambda) | lambda < 0){
    stop('Incorrect lambda parameter. You should provide a numeric value.')
  }

  if(!is.numeric(w) | any(w <= 0) | any(w > 1)){
    stop('Incorrect w parameter. You should provide a numeric value.')
  }
}

#' @rdname zi_zipfpss
#' @export
d_zi_zipfpss <- function(x, alpha, lambda, w, log = FALSE){
  .prec.zipfpss.checkXvalue(x)
  .prec.zi_zipfpss.checkparams(alpha, lambda, w)

   values <- sapply(x, function(i, alpha, lambda, w, log){
    if(i == 0){
      return(w + (1 - w)*dzipfpss(i, alpha, lambda, log))
    } else {
      return((1-w)*dzipfpss(i, alpha, lambda, log))
    }
  }, alpha = alpha, lambda = lambda, w = w, log = log)

  return(values)
}
