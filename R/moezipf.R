#' The Marshal-Olkin Extended Zipf Distribution (MOEZipf).
#'
#' Probability Mass Function, Cumulative Function, Quantile Function and Random
#' Generation for the MOEZipf distribution with parameter \eqn{\alpha} and \eqn{\beta}.
#'
#' @name moezipf
#' @aliases dmoezipf
#' @aliases pmoezipf
#' @aliases qmoezipf
#' @aliases rmoezipf
#'
#' @param x,q Vector of positive integer values.
#' @param p Vector of probabilities.
#' @param n Number of random numbers to return.
#' @param alpha Value of the \eqn{\alpha} parameter (\eqn{\alpha > 1} ).
#' @param beta Value of the \eqn{\beta} parameter (\eqn{\beta > 0} ).
#' @param log,log.p Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE (default), probabilities are \eqn{P[X \leq x]}, otherwise, \eqn{P[X > x]}.
#' @details The \emph{probability mass function} at a positive integer value \eqn{x} of the MOEZipf distribution with
#' parameters \eqn{\alpha} and \eqn{\beta} is computed as follows:
#'
#' \deqn{p(x | \alpha, \beta) = \frac{x^{-\alpha} \beta \zeta(\alpha) }{[\zeta(\alpha) - \bar{\beta} \zeta (\alpha, x)] [\zeta (\alpha) - \bar{\beta} \zeta (\alpha, x + 1)]}, \alpha > 1, \beta > 0, }
#'
#' where \eqn{\zeta(\alpha)} is the Riemann-zeta function at \eqn{\alpha}, \eqn{\zeta(\alpha, x)}
#' is the Hurtwitz zeta function with arguments \eqn{\alpha} and x, and \eqn{\bar{\beta} = 1 - \beta}.
#'
#' The \emph{cumulative distribution function}, \eqn{F_{\alpha}(x)}, at a given positive real value \eqn{x},
#'  is calcuted from the survival function \eqn{S(x)} as:
#' \deqn{F(x) = 1 - S(x), }
#' the survival function \eqn{S(x)} is equal to:
#' \deqn{S(x) = \frac{\beta \zeta(\alpha, x + 1)}{\zeta(\alpha) - \bar{\beta}\zeta(\alpha, x + 1)}, \forall x > 0}
#'
#' The \emph{quantiles} of a MOEZipf distribution for a given probability
#' vector \code{p}, are obtained by computing the quantiles associated to a Zipf distribution with
#' the same parameter \eqn{\alpha}, and probability vector equal to:
#' \deqn{p\prime = \frac{p\,\beta}{1 + p\,(\beta - 1)}}
#'
#' @return {
#' \code{dmoezipf} gives the probability mass function,
#' \code{pmoezipf} gives the cumulative function,
#' \code{qmoezipf} gives the quantile function, and
#' \code{rmoezipf} generates random deviates.}
#'
#' @references {
#' Young, D. S. (2010). \emph{Tolerance: an R package for estimating tolerance intervals}. Journal of Statistical Software, 36(5), 1-39.
#'
#' Casellas, A. (2013) \emph{La distribució Zipf Estesa segons la transformació Marshall-Olkin}. Universitat Politécnica de Catalunya.
#'
#' Pérez-Casany, M. and Casellas, A. (2013) \emph{Marshall-Olkin Extended Zipf Distribution}. arXiv preprint arXiv:1304.4540.
#'
#' Duarte-López, A., Prat-Pérez, A., & Pérez-Casany, M. (2015, August). \emph{Using the Marshall-Olkin Extended Zipf Distribution in Graph Generation}. In European Conference on Parallel Processing (pp. 493-502). Springer International Publishing.
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

.dmoezipf.default <- function(x, alpha, beta){
  num <- beta * VGAM::zeta(alpha) * x^(-alpha)
  den <- (VGAM::zeta(alpha) - (1 - beta)*.zeta_x(alpha, x))*(VGAM::zeta(alpha) - (1 - beta) * .zeta_x(alpha, x + 1))
  values <- num/den

  return(values)
}

#' @rdname moezipf
#' @export
dmoezipf <- function(x, alpha, beta, log = FALSE){
  values <- sapply(x, .dmoezipf.default, alpha = alpha, beta = beta)
  if(log){
    return(log(values))
  }
  return(values)
}

.survival.default <- function(x, alpha, beta){
  num <- beta * (.zeta_x(alpha, x + 1))
  den <- (VGAM::zeta(alpha) - (1 - beta)*(.zeta_x(alpha, x + 1)))
  return(num/den)
}

#' @rdname moezipf
#' @export
pmoezipf <- function(q, alpha, beta, log.p = FALSE, lower.tail = TRUE){
  srvvl <- sapply(q, .survival.default, alpha = alpha, beta = beta, simplify = T)

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
  p <- (x*beta)/(1 + x*(beta-1))
  return(p)
}

#' @rdname moezipf
#' @export
qmoezipf <- function(p, alpha, beta, log.p = FALSE, lower.tail = TRUE){
  if(!is.numeric(alpha) || !is.numeric(beta)){
    stop('Wrong values for the parameters.')
  }

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

  if(length(which(p > 1 || p < 0 )) > 0){
    stop('There is a wrong value(s) in the p parameter.')
  }

  u <- sapply(p, .qmoezipf.default, beta = beta)
  data <- tolerance::qzipfman(u, s = alpha, b = NULL, N = Inf)
  return(data)
}

#' @rdname moezipf
#' @export
rmoezipf <- function(n, alpha, beta){
  if(!is.numeric(alpha) || !is.numeric(beta)){
    stop('Wrong values for the parameters.')
  }

  uValues <- stats::runif(n, 0, 1)
  data <- qmoezipf(uValues, alpha, beta)
  print(data)
  lvl <- names(sort(table(data), decreasing = TRUE))
  uniqVal <- 1:length(unique(data))
  data <- uniqVal[match(data, lvl)]
  return(data)
}
