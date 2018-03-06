#' The Marshal-Olkin Extended Zipf Distribution (MOEZipf).
#'
#' Probability mass function, cumulative distribution function, quantile function and random number
#' generation for the MOEZipf distribution with parameters \eqn{\alpha} and \eqn{\beta}. The support of the MOEZipf
#' distribution are the strictly positive integer numbers large or equal than one.
#'
#' @name moezipf
#' @aliases dmoezipf
#' @aliases pmoezipf
#' @aliases qmoezipf
#' @aliases rmoezipf
#'
#' @param x,q Vector of positive integer values.
#' @param p Vector of probabilities.
#' @param n Number of random values to return.
#' @param alpha Value of the \eqn{\alpha} parameter (\eqn{\alpha > 1} ).
#' @param beta Value of the \eqn{\beta} parameter (\eqn{\beta > 0} ).
#' @param log,log.p Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE (default), probabilities are \eqn{P[X \leq x]}, otherwise, \eqn{P[X > x]}.
#' @details The \emph{probability mass function} at a positive integer value \eqn{x} of the MOEZipf distribution with
#' parameters \eqn{\alpha} and \eqn{\beta} is computed as follows:
#'
#' \deqn{p(x | \alpha, \beta) = \frac{x^{-\alpha} \beta \zeta(\alpha) }{[\zeta(\alpha) - \bar{\beta} \zeta (\alpha, x)] [\zeta (\alpha) - \bar{\beta} \zeta (\alpha, x + 1)]},\, x = 1,2,...,\, \alpha > 1, \beta > 0, }
#'
#' where \eqn{\zeta(\alpha)} is the Riemann-zeta function at \eqn{\alpha}, \eqn{\zeta(\alpha, x)}
#' is the Hurtwitz zeta function with arguments \eqn{\alpha} and x, and \eqn{\bar{\beta} = 1 - \beta}.
#'
#' The \emph{cumulative distribution function}, at a given positive integer value \eqn{x},
#'  is computed as \eqn{F(x) = 1 - S(x)}, where the survival function \eqn{S(x)} is equal to:
#' \deqn{S(x) = \frac{\beta\, \zeta(\alpha, x + 1)}{\zeta(\alpha) - \bar{\beta}\,\zeta(\alpha, x + 1)},\, x = 1, 2, .. }
#'
#' The quantile of the MOEZipf\eqn{(\alpha, \beta)} distribution of a given probability value p
#' is equal to the quantile of the Zipf\eqn{(\alpha)} distribution at the value:
#' \deqn{p\prime = \frac{p\,\beta}{1 + p\,(\beta - 1)}}
#'
#' The quantiles of the Zipf\eqn{(\alpha)} distribution are computed by means of the \emph{tolerance}
#' package.
#'
#' To generate random data from a MOEZipf one applies the \emph{quantile} function over \emph{n} values randomly generated
#' from an Uniform distribution in the interval (0, 1).
#'
#' @return {
#' \code{dmoezipf} gives the probability mass function,
#' \code{pmoezipf} gives the cumulative distribution function,
#' \code{qmoezipf} gives the quantile function, and
#' \code{rmoezipf} generates random values from a MOEZipf distribution.}
#'
#' @references {
#' Casellas, A. (2013) \emph{La distribució Zipf Estesa segons la transformació Marshall-Olkin}. Universitat Politécnica de Catalunya.
#'
#' Devroye L. (1986) Non-Uniform Random Variate Generation. Springer, New York, NY.
#'
#' Duarte-López, A., Prat-Pérez, A., & Pérez-Casany, M. (2015). \emph{Using the Marshall-Olkin Extended Zipf Distribution in Graph Generation}. European Conference on Parallel Processing, pp. 493-502, Springer International Publishing.
#'
#' Pérez-Casany, M. and Casellas, A. (2013) \emph{Marshall-Olkin Extended Zipf Distribution}. arXiv preprint arXiv:1304.4540.
#'
#' Young, D. S. (2010). \emph{Tolerance: an R package for estimating tolerance intervals}. Journal of Statistical Software, 36(5), 1-39.
#' }
#'
#' @examples
#' dmoezipf(1:10, 2.5, 1.3)
#' pmoezipf(1:10, 2.5, 1.3)
#' qmoezipf(0.56, 2.5, 1.3)
#' rmoezipf(10, 2.5, 1.3)
#'
NULL
#> NULL

.prec.moezipf.checkXvalue <- function(x){
  if(!is.numeric(x) || x < 1 || x%%1 != 0) {
    stop('The x value is not included into the support of the distribution.')
  }
}

.prec.moezipf.checkparams <- function(alpha, beta){
  if(!is.numeric(alpha) || alpha <= 1){
    stop('Incorrect alpha parameter. This parameter should be greater than one.')
  }

  if(!is.numeric(beta) || beta < 0){
    stop('Incorrect beta parameter. You should provide a numeric value.')
  }
}

.dmoezipf.default <- function(x, alpha, beta, z){
  .prec.moezipf.checkXvalue(x)
  num <- beta * z * x^(-alpha)
  den <- (z - (1 - beta)*.zeta_x(alpha, x))*(z - (1 - beta) * .zeta_x(alpha, x + 1))
  return(num/den)
}

#' @rdname moezipf
#' @export
dmoezipf <- function(x, alpha, beta, log = FALSE){
  .prec.moezipf.checkparams(alpha, beta)
  z <- VGAM::zeta(alpha)
  values <- sapply(x, .dmoezipf.default, alpha = alpha, beta = beta, z = z)
  if(log){
    return(log(values))
  }
  return(values)
}

.survival.default <- function(x, alpha, beta, z){
  .prec.moezipf.checkXvalue(x)
  zetaX <- .zeta_x(alpha, x + 1)
  num <- beta * (zetaX)
  den <- (z - (1 - beta)*(zetaX))
  return(num/den)
}

#' @rdname moezipf
#' @export
pmoezipf <- function(q, alpha, beta, log.p = FALSE, lower.tail = TRUE){
  .prec.moezipf.checkparams(alpha, beta)
  z <- VGAM::zeta(alpha)
  srvvl <- sapply(q, .survival.default, alpha = alpha, beta = beta, z = z, simplify = TRUE)

  if(!log.p && lower.tail){
    return(1 - srvvl)
  } else{
    if(!log.p && !lower.tail){
      return(srvvl)
    } else{
      if(log.p && !lower.tail){
        return(log(srvvl))
      }
      return(log(1-srvvl))
    }
  }
}

.qmoezipf.default <- function(x, beta){
  if(x > 1 || x < 0){
    stop('Wrong values for the p parameter.')
  }
  return((x*beta)/(1 + x*(beta-1)))
}

#' @rdname moezipf
#' @export
qmoezipf <- function(p, alpha, beta, log.p = FALSE, lower.tail = TRUE){
  .prec.moezipf.checkparams(alpha, beta)

  if(length(p) < 1){
    stop('Wrong values for the p parameter.')
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

  u <- sapply(p, .qmoezipf.default, beta = beta)
  data <- tolerance::qzipfman(u, s = alpha, b = NULL, N = Inf)
  return(data)
}

#' @rdname moezipf
#' @export
rmoezipf <- function(n, alpha, beta){
  .prec.moezipf.checkXvalue(n)
  .prec.moezipf.checkparams(alpha, beta)

  uValues <- stats::runif(n, 0, 1)
  data <- qmoezipf(uValues, alpha, beta)
  # print(data)
  # lvl <- names(sort(table(data), decreasing = TRUE))
  # uniqVal <- 1:length(unique(data))
  # data <- uniqVal[match(data, lvl)]
  return(data)
}

