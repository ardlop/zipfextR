
#' The Zipf-Poisson Stop Sum Distribution (Z-PSS).
#'
#' Probability Mass function, Cumulative function, Quantile function and Random generation for the Z-PSS distribution
#' with parameters \eqn{\alpha} and \eqn{\lambda}.
#'
#' @name zpss
#' @aliases dzpss
#' @aliases pzpss
#' @aliases rzpss
#'
#' @param x,q Vector of positive integer values.
#' @param p Vector of probabilities.
#' @param alpha Value of the \eqn{\alpha} parameter (\eqn{\alpha > 1} ).
#' @param lambda Value of the \eqn{\lambda} parameter (\eqn{\lambda \geq 0} ).
#' @param n Number of random values to return.
#' @param log,log.p Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE (default), probabilities are \eqn{P[X \leq x]}, otherwise, \eqn{P[X > x]}.
#' @param isTruncated Logical; if TRUE, the truncated version of the distribution is returned.
NULL
#> NULL

.prec.zpss.checkXvalue <- function(x){
  if(!is.numeric(x) || x < 0 || x%%1 != 0) {
    stop('The x value is not included into the support of the distribution.')
  }
}

.prec.zpss.checkparams <- function(alpha, lambda){
  if(!is.numeric(alpha) || alpha <= 1){
    stop('Incorrect alpha parameter. This parameter should be greater than one.')
  }

  if(!is.numeric(lambda) || lambda < 0){
    stop('Incorrect lambda parameter. You should provide a numeric value.')
  }
}

.panjerRecursion <- function(k, alpha, lambda){
  p0 <- exp(-lambda)

  if(k == 0){
    return(p0)
  }

  probs <- array(0, k + 1)
  probs[1] <- p0
  z_a <- VGAM::zeta(alpha)

  for(i in 1:k){
    probs[i+1] <- (lambda/(i*z_a)) * sum((1:i)^(-alpha + 1) * probs[i:1])
  }
  return(probs)
}

.getProbs <- function(x, alpha, lambda, isTruncated = FALSE){
  k <- max(x)
  .prec.zpss.checkXvalue(k)
  .prec.zpss.checkparams(alpha, lambda)

  probs <- .panjerRecursion(k, alpha, lambda)

  if(isTruncated){
    probs <- (probs)/(1 - probs[1])
    probs <- probs[-1]
  }
  return(probs)
}

#' @rdname zpss
#' @export
dzpss <- function(x, alpha, lambda, log = FALSE, isTruncated = FALSE){
  # k <- max(x)
  # .prec.zpss.checkXvalue(k)
  # .prec.zpss.checkparams(alpha, lambda)
  #
  # probs <- .panjerRecursion(k, alpha, lambda)
  #
  # if(isTruncated){
  #   probs <- (probs)/(1 - probs[1])
  #   probs <- probs[-1]
  # }
  probs <- .getProbs(x, alpha, lambda)

  finalProbs <- probs[if(isTruncated) x else (x+1)]#probs[x]#
  if(log){
    return(log(finalProbs))
  }

  return(finalProbs)
}

#' @rdname zpss
#' @export
pzpss <- function(q, alpha, lambda, log.p = FALSE, lower.tail = TRUE, isTruncated = FALSE){
  # k <- max(q)
  # .prec.zpss.checkXvalue(k)
  # .prec.zpss.checkparams(alpha, lambda)
  #
  # probs <- .panjerRecursion(k, alpha, lambda)
  # if(isTruncated){
  #   probs <- (probs)/(1 - probs[1])
  #   probs <- probs[-1]
  # }
  probs <- .getProbs(q, alpha, lambda)
  finalProbs <- array(0, length(q))
  # finalProbs <- sapply(1:length(q), function(i, q, probs, isTruncated){
  #   index <- if(isTruncated) i else (i+1)
  #   print(c(index, sum(probs[1:q[index]])))
  #   finalProbs[index] <- sum(probs[1:q[index]])
  # }, q = q, probs = probs, isTruncated = isTruncated)

  for(i in 1:length(q)){
    # print(if(isTruncated) 1 else 1:(q[i]+1))
    # print(probs)
    # print(probs[if(isTruncated) 1 else 0:q[i]])
    finalProbs[i] <- sum(probs[1:if(isTruncated) q[i] else (q[i]+1)])
  }
  #finalProbs <- probs[q]
  if(!log.p & lower.tail){
    return(finalProbs)
  } else{
    if(!log.p & !lower.tail){
      return(1 - finalProbs)
    } else{
      if(log.p & !lower.tail){
        return(log(1 - finalProbs))
      }
      return(log(finalProbs))
    }
  }
}


#' #' @rdname zpss
#' #' @export
#' rzpss <- function(n, alpha, lambda){
#'   .prec.zpe.checkparams(alpha, lambda)
#'
#'   data <- array(0, n)
#'   for(i in 1:n){
#'     nPois <- stats::rpois(1, lambda = lambda)
#'     if(nPois == 0){
#'       data[i] <- 0
#'       print(c(i, 0))
#'     } else{
#'       xZipfs <- tolerance::rzipfman(nPois, s = alpha, b = NULL, N = Inf)
#'       data[i] <- sum(xZipfs)
#'       # print(c(i, nPois, sum(xZipfs)))
#'     }
#'   }
#'    return(data)
#' }


#' @rdname zpss
#' @export
rzpss <- function(n, alpha, lambda, log.p = FALSE, lower.tail = TRUE, isTruncated = FALSE){
  .prec.zpe.checkparams(alpha, lambda)

  # data <- array(0, n)
  # for(i in 1:n){
  #   nPois <- stats::rpois(1, lambda = lambda)
  #   if(nPois == 0){
  #     data[i] <- 0
  #     print(c(i, 0))
  #   } else{
  #     xZipfs <- tolerance::rzipfman(nPois, s = alpha, b = NULL, N = Inf)
  #     data[i] <- sum(xZipfs)
  #     print(c(i, nPois, sum(xZipfs)))
  #   }
  # }

  u <- stats::runif(n)
  data <- sapply(u, qzpss, alpha, lambda, log.p = FALSE, lower.tail = TRUE, isTruncated = FALSE)
  #data <- tolerance::qzipfman(u, s = alpha, b = NULL, N = Inf)
  return(data)
}

.invertMethod <- function(p, alpha, lambda, log.p, lower.tail, isTruncated) {
  i <- if(isTruncated) 1 else 0
  p_i <- pzpss(i, alpha = alpha, lambda = lambda, log.p, lower.tail, isTruncated)
  repeat{
    if(p <= p_i){
      return(i)
    }
    i <- i + 1
    p_i = pzpss(i, alpha = alpha, lambda = lambda, log.p, lower.tail, isTruncated)
  }
  # return(i)
}


#' @rdname zpss
#' @export
qzpss <- function(p, alpha, lambda, log.p = FALSE, lower.tail = TRUE, isTruncated = FALSE){
  .prec.zpe.checkparams(alpha, lambda)

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

  data <- sapply(p, .invertMethod, alpha, lambda, log.p, lower.tail, isTruncated)

  #u <- sapply(p, .getUprime, beta = beta)
  #data <- tolerance::qzipfman(u, s = alpha, b = NULL, N = Inf)
  return(data)
}


