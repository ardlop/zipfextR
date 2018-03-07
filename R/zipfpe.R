#' The Zipf-Poisson Extreme Distribution (Zipf-PE).
#'
#' Probability mass function, cumulative distribution function, quantile function and random number
#' generation for the Zipf-PE distribution with parameters \eqn{\alpha} and \eqn{\beta}. The support of the Zipf-PE
#' distribution are the strictly positive integer numbers large or equal than one.
#'
#' @name zipfpe
#' @aliases dzipfpe
#' @aliases pzipfpe
#' @aliases qzipfpe
#' @aliases rzipfpe
#'
#' @param x,q Vector of positive integer values.
#' @param p Vector of probabilities.
#' @param n Number of random values to return.
#' @param alpha Value of the \eqn{\alpha} parameter (\eqn{\alpha > 1} ).
#' @param beta Value of the \eqn{\beta} parameter (\eqn{\beta\in (-\infty, +\infty)} ).
#' @param log,log.p Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE (default), probabilities are \eqn{P[X \leq x]}, otherwise, \eqn{P[X > x]}.
#' @details The \emph{probability mass function} of the Zipf-PE distribution with parameters \eqn{\alpha} and \eqn{\beta}
#' at a positive integer value \eqn{x} is computed as follows:
#'
#' \deqn{p(x | \alpha, \beta) = \frac{e^{\beta (1 - \frac{\zeta(\alpha, x)}{\zeta(\alpha)})} (e^{\beta \frac{x^{-\alpha}}{\zeta(\alpha)}} - 1)}
#' {e^{\beta} - 1},\, x= 1,2,...,\, \alpha > 1,\, -\infty < \beta < +\infty,}
#'
#' where \eqn{\zeta(\alpha)} is the Riemann-zeta function at \eqn{\alpha}, and \eqn{\zeta(\alpha, x)}
#' is the Hurtwitz zeta function with arguments \eqn{\alpha} and x.
#'
#' The \emph{cumulative distribution function} at a given positive
#' integer value \eqn{x}, \eqn{F(x)}, is equal to:
#' \deqn{F(x) = \frac{e^{\beta (1 - \frac{\zeta(\alpha, x + 1)}{\zeta(\alpha)})} - 1}{e^{\beta} -1}}
#'
#' The quantile of the Zipf-PE\eqn{(\alpha, \beta)} distribution of a given probability value p
#' is equal to the quantile of the Zipf\eqn{(\alpha)} distribution at the value:
#'
#' \deqn{p\prime = \frac{log(p\, (e^{\beta} - 1) + 1)}{\beta}}
#' The quantiles of the Zipf\eqn{(\alpha)} distribution are computed by means of the \emph{tolerance}
#' package.
#'
#' To generate random data from a Zipf-PE one applies the \emph{quantile} function over \emph{n} values randomly generated
#' from an Uniform distribution in the interval (0, 1).
#'
#' @return {
#' \code{dzipfpe} gives the probability mass function,
#' \code{pzipfpe} gives the cumulative function,
#' \code{qzipfpe} gives the quantile function, and
#' \code{rzipfpe} generates random values from a Zipf-PE distribution.  }
#'
#' @references {
#' Young, D. S. (2010). \emph{Tolerance: an R package for estimating tolerance intervals}. Journal of Statistical Software, 36(5), 1-39.
#' }
#'
#' @examples
#' dzipfpe(1:10, 2.5, -1.5)
#' pzipfpe(1:10, 2.5, -1.5)
#' qzipfpe(0.56, 2.5, 1.3)
#' rzipfpe(10, 2.5, 1.3)
#'
NULL
#> NULL

.prec.zipfpe.checkXvalue <- function(x){
  if(!is.numeric(x) || x < 1 || x%%1 != 0) {
    stop('The x value is not included into the support of the distribution.')
  }
}

.prec.zipfpe.checkparams <- function(alpha, beta){
  if(!is.numeric(alpha) | alpha <= 1){
    stop('Incorrect alpha parameter. This parameter should be greater than one.')
  }

  if(!is.numeric(beta)){
    stop('Incorrect beta parameter. You should provide a numeric value.')
  }
}

.dzpe.default <- function(x, alpha, beta, z){
  .prec.zipfpe.checkXvalue(x)
  zetaX <-.zeta_x(alpha, x)
  return((exp(beta*(1 - (zetaX/z)))*(exp(beta*(x^(-alpha)/z)) - 1))/(exp(beta) -1))
}

#' @rdname zipfpe
#' @export
dzipfpe <-  function(x, alpha, beta, log = FALSE){
  .prec.zipfpe.checkparams(alpha, beta)

  z <- VGAM::zeta(alpha)
  probs <- sapply(x, .dzpe.default, alpha = alpha, beta = beta, z = z)

  if(log) {
    return(log(probs))
  }

  return(probs)
}

.pzpe.default <- function(v, alpha, beta, z){
  .prec.zipfpe.checkXvalue(v)
  zetaX <- .zeta_x(alpha, v + 1)
  return((exp(beta * (1 - (zetaX/z))) - 1)/(exp(beta) - 1))
}

#' @rdname zipfpe
#' @export
pzipfpe <- function(q, alpha, beta, log.p = FALSE, lower.tail = TRUE){
  .prec.zipfpe.checkparams(alpha, beta)

  z <- VGAM::zeta(alpha)
  probs <- sapply(q, .pzpe.default, alpha = alpha, beta = beta, z = z)

  if(!log.p & lower.tail){
    return(probs)
  } else{
    if(!log.p & !lower.tail){
      return(1 - probs)
    } else{
      if(log.p & !lower.tail){
        return(log(1 - probs))
      }
      return(log(probs))
    }
  }
}

.getUprime <- function(u, beta){
  p <- log(u*(exp(beta) - 1) + 1)/beta
  return(p)
}

#' @rdname zipfpe
#' @export
qzipfpe <- function(p, alpha, beta, log.p = FALSE, lower.tail = TRUE){
  .prec.zipfpe.checkparams(alpha, beta)

  if(length(p) < 1){
    stop('Wrong value(s) for the p parameter.')
  }

  if(log.p && lower.tail){
    p <- exp(p)
  } else{
    if(log.p && !lower.tail){
      p <- 1-exp(p)
    } else{
      if(!log.p && !lower.tail){
        p <- 1-p
      }
    }
  }

  if(length(which(p > 1 || p < 0 )) > 0){
    stop('There is a wrong value(s) in the p parameter.')
  }

  u <- sapply(p, .getUprime, beta = beta)
  data <- tolerance::qzipfman(u, s = alpha, b = NULL, N = Inf)
  return(data)
}

#' @rdname zipfpe
#' @export
rzipfpe <- function(n, alpha, beta){
  .prec.zipfpe.checkXvalue(n)
  .prec.zipfpe.checkparams(alpha, beta)

  uValues <- stats::runif(n, 0, 1)
  data <- qzipfpe(uValues, alpha, beta)
  return(data)
}


