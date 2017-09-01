#' The Zipf-Poisson Extreme Distribution (ZPE).
#'
#' Probability Mass Function, Cumulative Function of the ZPE distribution
#' with parameter \eqn{\alpha} and \eqn{\beta}.
#'
#' @name zpe
#' @aliases dzpe
#' @aliases pzpe
#' @aliases qzpe
#' @aliases rzpe
#'
#' @param x,q Vector of positive integer values.
#' @param p Vector of probabilities.
#' @param n Number of random numbers to return.
#' @param alpha Value of the \eqn{\alpha} parameter (\eqn{\alpha > 1} ).
#' @param beta Value of the \eqn{\beta} parameter (\eqn{\beta > 0} ).
#' @param log,log.p Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE (default), probabilities are \eqn{P[X \leq x]}, otherwise, \eqn{P[X > x]}.
#' @details The \emph{probability mass function} at a positive integer value \eqn{x} of the ZPE distribution with
#' parameters \eqn{\alpha} and \eqn{\beta} is computed as follows:
#'
#' \deqn{p(x | \alpha, \beta) = \frac{e^{\beta (1 - \frac{\zeta(\alpha, x)}{\zeta(\alpha)})} (e^{\beta \frac{x^{-alpha}}{\zeta(\alpha)}} - 1)}
#' {e^{beta} - 1}, \alpha > 1, -\infty < \beta < +\infty,}
#'
#' where \eqn{\zeta(\alpha)} is the Riemann-zeta function at \eqn{\alpha}, \eqn{\zeta(\alpha, x)}
#' is the Hurtwitz zeta function with arguments \eqn{\alpha} and x.
#'
#' The \emph{cumulative distribution function}, \eqn{F_{\alpha, \beta}(x)}, at a given positive integer value \eqn{x},
#'  is calcuted as:
#' \deqn{F(x) = \frac{e^{beta (1 - \frac{\zeta(\alpha, x + 1)}{\zeta(\alpha)})} - 1}{e^{beta} -1}}
#'
#' The \emph{quantiles} of a ZPE distribution for a given probability
#' vector \code{p}, are obtained by computing the quantiles associated to a Zipf distribution with
#' the same parameter \eqn{\alpha}, and probability vector equal to:
#' \deqn{p\prime = \frac{log(p\, (e^{\beta} - 1) + 1)}{\beta}}
#'
#' @return {
#' \code{dzpe} gives the probability mass function,
#' \code{pzpe} gives the cumulative function.
#' \code{qzpe} gives the quantile function, and
#' \code{rzpe} generates random deviates.  }
#'
#' @references {
#' Young, D. S. (2010). \emph{Tolerance: an R package for estimating tolerance intervals}. Journal of Statistical Software, 36(5), 1-39.
#' }
#'
#' @examples
#' dzpe(1:10, 2.5, -1.5)
#' pzpe(1:10, 2.5, -1.5)
#' qzpe(0.56, 2.5, 1.3)
#' rzpe(10, 2.5, 1.3)
#'
NULL
#> NULL

.prec.zpe.checkXvalue <- function(x){
  if(!is.numeric(x) || x < 1 || x%%1 != 0) {
    stop('The x value is not included into the support of the distribution.')
  }
}

.prec.zpe.checkparams <- function(alpha, beta){
  if(!is.numeric(alpha) | alpha <= 1){
    stop('Incorrect alpha parameter. This parameter should be greater than one.')
  }

  if(!is.numeric(beta)){
    stop('Incorrect beta parameter. You should provide a numeric value.')
  }
}

.dzpe.default <- function(x, alpha, beta, z){
  .prec.zpe.checkXvalue(x)
  zetaX <-.zeta_x(alpha, x)
  return((exp(beta*(1 - (zetaX/z)))*(exp(beta*(x^(-alpha)/z)) - 1))/(exp(beta) -1))
}

#' @rdname zpe
#' @export
dzpe <-  function(x, alpha, beta, log = FALSE){
  .prec.zpe.checkparams(alpha, beta)

  z <- VGAM::zeta(alpha)
  probs <- sapply(x, .dzpe.default, alpha = alpha, beta = beta, z = z)

  if(log) {
    return(log(probs))
  }

  return(probs)
}

.pzpe.default <- function(v, alpha, beta, z){
  .prec.zpe.checkXvalue(v)
  zetaX <- .zeta_x(alpha, v + 1)
  return((exp(beta * (1 - (zetaX/z))) - 1)/(exp(beta) - 1))
}

#' @rdname zpe
#' @export
pzpe <- function(q, alpha, beta, log.p = FALSE, lower.tail = TRUE){
  .prec.zpe.checkparams(alpha, beta)

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

#' @rdname zpe
#' @export
qzpe <- function(p, alpha, beta, log.p = FALSE, lower.tail = TRUE){
  .prec.zpe.checkparams(alpha, beta)

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


#' @rdname zpe
#' @export
rzpe <- function(n, alpha, beta){
  .prec.moezipf.checkXvalue(n)
  .prec.moezipf.checkparams(alpha, beta)

  uValues <- stats::runif(n, 0, 1)
  data <- qzpe(uValues, alpha, beta)
  return(data)
}


