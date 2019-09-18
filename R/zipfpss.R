
#' The Zipf-Poisson Stop Sum Distribution (Zipf-PSS).
#'
#' Probability mass function, cumulative distribution function, quantile function and random number
#' generation for the Zipf-PSS distribution with parameters \eqn{\alpha} and \eqn{\lambda}. The support of the Zipf-PSS
#' distribution are the positive integer numbers including the zero value. In order to work with its zero-truncated version
#' the parameter \code{isTruncated} should be equal to True.
#'
#' @name zipfpss
#' @aliases dzipfpss
#' @aliases pzipfpss
#' @aliases rzipfpss
#'
#' @param x,q Vector of positive integer values.
#' @param p Vector of probabilities.
#' @param alpha Value of the \eqn{\alpha} parameter (\eqn{\alpha > 1} ).
#' @param lambda Value of the \eqn{\lambda} parameter (\eqn{\lambda > 0} ).
#' @param n Number of random values to return.
#' @param log,log.p Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE (default), probabilities are \eqn{P[X \leq x]}, otherwise, \eqn{P[X > x]}.
#' @param isTruncated Logical; if TRUE, the zero truncated version of the distribution is returned.
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

.prec.zipfpss.checkXvalue <- function(x){
  if(!is.numeric(x) | any(x < 0) | any(x%%1 != 0)) {
    stop('The x value is not included into the support of the distribution.')
  }
}

.prec.zipfpss.checkparams <- function(alpha, lambda){
  if(!is.numeric(alpha) | alpha <= 1){
    stop('Incorrect alpha parameter. This parameter should be greater than one.')
  }

  if(!is.numeric(lambda) | lambda < 0){
    stop('Incorrect lambda parameter. You should provide a numeric value.')
  }
}

# .panjerRecursion <- function(k, alpha, lambda, previousProb = NULL){
#   p0 <- exp(-lambda)
#   if(k == 0){
#     return(p0)
#   }
#
#   z_a <- VGAM::zeta(alpha)
#   probs <- NULL
#
#   if(is.null(previousProb)) {
#     probs <- array(0, k + 1)
#     probs[1] <- p0
#
#     for(i in 1:k){
#       probs[i+1] <- (lambda/(i*z_a)) * sum((1:i)^(-alpha + 1) * probs[i:1])
#     }
#   } else{
#     probs <- array(0, k + 1)
#     lastProbLength <- length(previousProb)
#     probs[1:lastProbLength] <- previousProb[1:lastProbLength]
#     for(i in lastProbLength:k){
#       probs[i+1] <- (lambda/(i*z_a)) * sum((1:i)^(-alpha + 1) * probs[i:1])
#     }
#   }
#   return(probs)
# }

.getPanjerProbs <- function(startIndex, finalIndex, probs, alpha, lambda){
  z_a <- VGAM::zeta(alpha)
  for(i in startIndex:finalIndex){
    probs[i+1] <- (lambda/(i*z_a)) * sum((1:i)^(-alpha + 1) * probs[i:1])
  }
  return(probs)
}

.panjerRecursion <- function(k, alpha, lambda, previousProb = NULL){
  p0 <- exp(-lambda)
  if(k == 0){
    return(p0)
  }

  probs <- NULL
  if(is.null(previousProb)) {
    probs <- array(0, k + 1)
    probs[1] <- p0
    probs <- .getPanjerProbs(1, k, probs, alpha, lambda)
  } else{
    probs <- array(0, k + 1)
    lastProbLength <- length(previousProb)
    probs[1:lastProbLength] <- previousProb[1:lastProbLength]
    probs <- .getPanjerProbs(lastProbLength, k, probs, alpha, lambda)
  }
  return(probs)
}




.getProbs <- function(x, alpha, lambda, isTruncated = FALSE, previousProb = NULL){
  k <- max(x)
  .prec.zipfpss.checkXvalue(k)
  .prec.zipfpss.checkparams(alpha, lambda)

  probs <- .panjerRecursion(k, alpha, lambda, previousProb)

  if(isTruncated){
    probs <- (probs)/(1 - probs[1])
    probs <- probs[-1]
  }
  return(probs)
}

#' @rdname zipfpss
#' @export
dzipfpss <- function(x, alpha, lambda, log = FALSE, isTruncated = FALSE){
  if(lambda == 0){
    return(.999999)
  } else {
    probs <- .getProbs(x, alpha, lambda, isTruncated = isTruncated)

    finalProbs <- probs[if(isTruncated) x else (x+1)]#probs[x]#
    if(log){
      return(log(finalProbs))
    }
    return(finalProbs)
  }
}

#' @rdname zipfpss
#' @export
pzipfpss <- function(q, alpha, lambda, log.p = FALSE, lower.tail = TRUE, isTruncated = FALSE){
  probs <- .getProbs(q, alpha, lambda)
  finalProbs <- array(0, length(q))

  for(i in seq_along(q)){
    finalProbs[i] <- sum(probs[1:if(isTruncated) q[i] else (q[i]+1)])
  }

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

#' @rdname zipfpss
#' @export
rzipfpss <- function(n, alpha, lambda, log.p = FALSE, lower.tail = TRUE, isTruncated = FALSE){
  .prec.zipfpss.checkparams(alpha, lambda)

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
  data <- sapply(u, qzipfpss, alpha, lambda, log.p, lower.tail, isTruncated)
  #data <- tolerance::qzipfman(u, s = alpha, b = NULL, N = Inf)
  return(data)
}

.invertMethod <- function(p, alpha, lambda, log.p, lower.tail, isTruncated) {
  i <- if(isTruncated) 1 else 0
  p_i <- pzipfpss(i, alpha = alpha, lambda = lambda, log.p, lower.tail, isTruncated)
  repeat{
    if(p <= p_i){
      return(i)
    }
    i <- i + 1
    p_i <- pzipfpss(i, alpha = alpha, lambda = lambda, log.p, lower.tail, isTruncated)
  }
}


#' @rdname zipfpss
#' @export
qzipfpss <- function(p, alpha, lambda, log.p = FALSE, lower.tail = TRUE, isTruncated = FALSE){
  .prec.zipfpss.checkparams(alpha, lambda)

  if(length(p) < 1){
    stop('Wrong value(s) for the p parameter.')
  }

  if(log.p & lower.tail){
    p <- exp(p)
  } else{
    if(log.p & !lower.tail){
      p <- 1-exp(p)
    } else{
      if(!log.p & !lower.tail){
        p <- 1-p
      }
    }
  }

  if(any(p > 1) | any(p < 0 )){
    stop('There is a wrong value(s) in the p parameter.')
  }

  data <- sapply(p, .invertMethod, alpha, lambda, log.p, lower.tail, isTruncated)
  return(data)
}


