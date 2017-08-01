
#' The Zipf-Poisson Stop Sum Distribution (Z-PSS).
#'
#' Probability Mass Function, Cumulative Function of the Z-PSS distribution
#' with parameter \eqn{\alpha} and \eqn{\lambda}.
#'
#' @name zpss
#' @aliases dzpss
#' @aliases pzpss
#'
#' @param x,q Vector of positive integer values.
#' @param alpha Value of the \eqn{\alpha} parameter (\eqn{\alpha > 1} ).
#' @param lambda Value of the \eqn{\lambda} parameter (\eqn{\lambda >= 0} ).
#' @param log,log.p Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE (default), probabilities are \eqn{P[X \leq x]}, otherwise, \eqn{P[X > x]}.
#' @param isTruncated Logical; if TRUE, the truncated version of the distribution is returned.
NULL
#> NULL

.prec.zpss.checkXvalue <- function(x){
  if(!is.numeric(x) || x < 1 || x%%1 != 0) {
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

#' @rdname zpss
#' @export
dzpss <- function(x, alpha, lambda, log = FALSE, isTruncated = FALSE){
  k <- max(x)
  .prec.zpss.checkXvalue(k)
  .prec.zpss.checkparams(alpha, lambda)

  probs <- .panjerRecursion(k, alpha, lambda)

  if(isTruncated){
    probs <- (probs)/(1 - probs[1])
    probs <- probs[-1]
  }

  finalProbs <- probs[x]
  if(log){
    return(log(finalProbs))
  }

  return(finalProbs)
}

#' @rdname zpss
#' @export
pzpss <- function(q, alpha, lambda, log.p = FALSE, lower.tail = TRUE, isTruncated = FALSE){
  k <- max(q)
  .prec.zpss.checkXvalue(k)
  .prec.zpss.checkparams(alpha, lambda)

  probs <- .panjerRecursion(k, alpha, lambda)
  if(isTruncated){
    probs <- (probs)/(1 - probs[1])
    probs <- probs[-1]
  }

  finalProbs <- array(0, length(q))
  finalProbs <- sapply(1:length(q), function(i, q, probs){
    finalProbs[i] <- sum(probs[1:q[i]])
  }, q = q, probs = probs)
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





