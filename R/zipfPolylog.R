#' The Zipf-Polylog Distribution (Zipf-Polylog).
#'
#' Probability mass function of the Zipf-Polylog distribution with parameters \eqn{\alpha} and \eqn{\beta}.
#' The support of the Zipf-Polylog distribution are the strictly positive integer numbers large or equal
#' than one.
#'
#' @name zipfPolylog
#' @aliases dzipfpolylog
#' @aliases pzipfpolylog
#' @aliases qzipfpolylog
#' @aliases rzipfpolylog
#'
#' @param x Vector of positive integer values.
#' @param p Vector of probabilities.
#' @param n Number of random values to return.
#' @param alpha Value of the \eqn{\alpha} parameter (\eqn{\alpha > 1} ).
#' @param beta Value of the \eqn{\beta} parameter (\eqn{\beta > 0} ).
#' @param nSum The number of terms used for computing the Polylogarithm function (Default = 1000).
#' @param log,log.p Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE (default), probabilities are \eqn{P[X \leq x]}, otherwise, \eqn{P[X > x]}.
#' @details The \emph{probability mass function} at a positive integer value \eqn{x} of the Zipf-Polylog distribution with
#' parameters \eqn{\alpha} and \eqn{\beta} is computed as follows:
#'
#' @return {
#' \code{dzipfpolylog} gives the probability mass function
#' }
#' @importFrom copula polylog
#' @examples
#' dzipfpolylog(1:10, 1.61, 0.98)
#' pzipfpolylog(1:10, 1.61, 0.98)
#' qzipfpolylog(0.8, 1.61, 0.98)
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

.pzipfpoly <- function(x, alpha, beta, nSum = 1000){
  liValue <- copula::polylog(z = beta, s = alpha, method = 'sum', n.sum = nSum)
  sumValues <- sum(sapply(1:x, function(i, alpha, beta){
    return(beta^(i) * i^(-alpha))
  }, alpha = alpha, beta = beta))
  cumProb <- (1/liValue) * sumValues
  return(cumProb)
}

#' @rdname zipfPolylog
#' @export
pzipfpolylog <- function(x, alpha, beta, log.p = FALSE, lower.tail = TRUE, nSum = 1000){
  .prec.zipfpolylog.checkXvalue(x)
  .prec.zipfpolylog.checkparams(alpha, beta)

  finalProbs <-sapply(x, .pzipfpoly, alpha = alpha, beta = beta, nSum = nSum)

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

#' @rdname zipfPolylog
#' @export
qzipfpolylog <- function(p, alpha, beta, log.p = FALSE, lower.tail = TRUE, nSum = 1000){
  .prec.zipfpolylog.checkparams(alpha, beta)

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

  # if(length(which(p > 1 || p < 0 )) > 0){
  if(any(p > 1) | any(p < 0 )){
    stop('There is a wrong value(s) in the p parameter.')
  }

  # liValue <- copula::polylog(z = -log(beta), s = alpha, method = 'sum', n.sum = nSum)
  # uprime <- p*liValue
  i <- 1
  p_i <- pzipfpolylog(i, alpha = alpha, beta = beta)

  repeat{
    if(p <= p_i){
      # print(i)
      return(i)
    }
    i <- i + 1
    p_i <- pzipfpolylog(i, alpha = alpha, beta = beta)
  }


  # x <- 1
  # # print((uprime <= ((beta^x)*(x^(-alpha)))))
  # cumValue <- exp(-alpha*log(x)-(-log(beta)*x)) #((beta^x)*(x^(-alpha)))
  # while(uprime <= cumValue){
  #   print(c(sprintf('x = %s, uprime = %s, cumulative = %s', x, uprime,cumValue)))
  #   x <- x + 1
  #   # print(cumValue)
  #   # print(((beta^x)*(x^(-alpha))))
  #   # print(beta^x)
  #   # print(beta)
  #   # print(x^-alpha)
  #   print(exp(-alpha*log(x)-(-log(beta)*x)))
  #   cumValue <- cumValue + exp(-alpha*log(x)-(-log(beta)*x))#((beta^x)*(x^(-alpha)))
  #   print(cumValue)
  # }
  return(x)
}

#' @rdname zipfPolylog
#' @export
rzipfpolylog <- function(n, alpha, beta, nSum = 1000){
  .prec.zipfpolylog.checkXvalue(n)
  .prec.zipfpolylog.checkparams(alpha, beta)

  uValues <- stats::runif(n, 0, 1)
  # print(uValues)
  data <- sapply(uValues, qzipfpolylog, alpha = alpha, beta = beta, nSum = nSum)
  # data <- qzipfpolylog(uValues, alpha, beta, nSum = nSum)
  return(data)
}
#
# data <- rzipfpolylog(100, -0.05, 0.95)
# data1 <- table(data)
# data1 <- data.frame(data1)
# data1[,1] <- as.numeric(data1[,1])
# data1[,2] <- as.numeric(data1[,2])
# plot(data1[,1], data1[,2], log = 'xy')
# zipfPolylogFit(data1, init_alpha = -0.08, init_beta = -0.02)

