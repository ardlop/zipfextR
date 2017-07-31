#' The Zipf-Poisson Extreme Distribution (ZPE).
#'
#' Probability Mass Function, Cumulative Function of the ZPE distribution
#' with parameter \eqn{\alpha} and \eqn{\beta}.
#'
#' @name zpe
#' @aliases dzpe
#' @aliases pzpe
#'
#' @return {
#' \code{dzpe} gives the probability mass function,
#' \code{pzpe} gives the cumulative function.
#'
NULL
#> NULL

#' @rdname zpe
#' @export
dzpe <-  function(x, alpha, lambda){
  # comprobar parametros.
  # unificar ecuaciones.
  # renombrar lambda as beta.

  if(x < 1) {
    stop('Wrong x value!!!')
  }
  z_alpha <- zeta(alpha)
  if(x == 1){
    return((exp(lambda/z_alpha) - 1)/(exp(lambda) - 1))
  }

  if(x >= 2){
    zeta_x1 <-  zeta_x(alpha, x+1)
    zx <- zeta_x(alpha, x)
    value <- (exp(lambda * ((z_alpha -zeta_x1)/z_alpha)) - exp(lambda*((z_alpha - zx)/z_alpha)))/(exp(lambda) - 1)
    return(value)
  }
}

#' @rdname zpe
#' @export
pzpe <- function(q, alpha, beta){
  return(NULL)
}
