.zpss_mle <- function(par, values, freq, truncated = TRUE) {
  alpha <- par[1]
  lambda <- par[2]
  probs <- .panjerRecursion(max(values), alpha, lambda)
  # probs <- getProbsFaster(alpha, lambda, max(values))
  if (truncated) {
    probs <- (probs) / (1 - probs[1])
    probs <- probs[-1]
  }
  - (sum(freq * log(probs[values+1])))
}

#' ZPSS parameters estimation.
#'
#' For a given count data set,  usually of the type of ranking data or frequencies of frequencies
#' data, estimates the parameters of the Z-PSS distribution by means of the maximum likelihood
#' method.
#' @param data Matrix of count data.
#' @param init_alpha Initial value of \eqn{\alpha} parameter (\eqn{\alpha > 1}).
#' @param init_lambda Initial value of \eqn{\lambda} parameter (\eqn{\lambda \geq 0}).
#' @param level Confidence level used to calculate the confidence intervals (default 0.95).
#' @param isTruncated Logical; if TRUE, the truncated version of the distribution is returned.(default = FALSE)
#' @param object An object from class "zpssR" (output of \emph{zpssFit} function).
#' @param x An object from class "zpssR" (output of \emph{zpssFit} function).
#' @param ... Further arguments to the generic functions. The extra arguments are passing
#' to the \emph{\link{optim}} function.
#'
#' @details
#' The argument \code{data} is a matrix where, for each row, the first column contains a count,
#' and the second column contains its corresponding frequency.
#'
#' The log-likelihood function is computed by means of the following equation:
#' \deqn{l(\alpha, \lambda, x) = \sum_{i =1} ^{m} f_a(x_i)\, log(P(Y = x_i))},
#' where \eqn{m} is the number of different values \eqn{x_i} in the sample, and \eqn{f_a(x_i)}
#' is the absolute frequency of \eqn{x_i}. The probabilities are calculated applying the Panjer recursion.
#'
#' The function \emph{\link{optim}} is used to estimate the parameters.
#' @return Returns a \emph{zpssR} object composed by the maximum likelihood parameter estimations,
#' their standard deviation, their confidence intervals and the log-likelihood value.
#' @examples
#' data <- rmoezipf(100, 2.5, 1.3)
#' data <- zipfExtR_getDataMatrix(data)
#' obj <- zpssFit(data, 1.001, 0.001)
#' @seealso \code{\link{zipfExtR_getDataMatrix}}, \code{\link{moezipf_getInitialValues}}.
#' @export
zpssFit <- function(data, init_alpha, init_lambda, level=0.95, isTruncated = FALSE, ...){
  Call <- match.call()
  if(!is.numeric(init_alpha) || !is.numeric(init_lambda)){
    stop('Wrong intial values for the parameters.')
  }

  tryCatch({
    res <- stats::optim(par = c(init_alpha, init_lambda), .zpss_mle, values = data[, 1], freq = data[, 2],
                  truncated = isTruncated, hessian = T, ...)
    estAlpha <- as.numeric(res$par[1])
    estLambda <- as.numeric(res$par[2])
    paramSD <- sqrt(diag(solve(res$hessian)))
    paramsCI <- .getConfidenceIntervals(paramSD, estAlpha, estLambda, level)

    structure(class = "zpssR", list(alphaHat = estAlpha,
                                   betaHat = estLambda,
                                   alphaSD = paramSD[1],
                                   lambdaSD = paramSD[2],
                                   alphaCI = c(paramsCI[1,1],paramsCI[1,2]),
                                   lambdaCI = c(paramsCI[2,1],paramsCI[2,2]),
                                   logLikelihood = -res$value,
                                   call = Call))

  }, error = function(e){
    print(c('Error', e))
  })
}


#' @rdname zpssFit
#' @export
residuals.zpssR <- function(object, ...){
  dataMatrix <- get(as.character(object[['call']]$data))
  fitted.values <- fitted(object)
  residual.values <- dataMatrix[, 2] - fitted.values
  return(residual.values)
}

#' @rdname zpssFit
#' @export
fitted.zpssR <- function(object, ...) {
  dataMatrix <- get(as.character(object[['call']]$data))
  N <- sum(dataMatrix[, 2])
  fitted.values <- N*sapply(dataMatrix[,1], dmoezipf, alpha = object[['alphaHat']],
                            beta = object[['betaHat']])
  return(fitted.values)
}

#' @rdname zpssFit
#' @export
coef.zpssR <- function(object, ...){
  estimation <- matrix(nrow = 2, ncol = 4)
  estimation[1, ] <- c(object[['alphaHat']], object[['alphaSD']], object[['alphaCI']][1], object[['alphaCI']][2])
  estimation[2, ] <- c(object[['betaHat']], object[['betaSD']], object[['betaCI']][1], object[['betaCI']][2])
  colnames(estimation) <- c("MLE", "Std. Dev.", paste0("Inf. ", "95% CI"),
                            paste0("Sup. ", "95% CI"))
  rownames(estimation) <- c("alpha", "lambda")
  estimation
}

#' @rdname zpssFit
#' @export
plot.zpssR <- function(x, ...){
  dataMatrix <- get(as.character(x[['call']]$data))
  graphics::plot(dataMatrix[,1], dataMatrix[,2], log="xy",
                 xlab="Observation", ylab="Frequency",
                 main="Fitting Z-PSS Distribution", ...)

  graphics::lines(dataMatrix[,1], fitted(x), col="blue")

  graphics::legend("topright",  legend = c('Observations', 'Z-PSS Distribution'),
                   col=c('black', 'blue'), pch=c(21,NA),
                   lty=c(NA, 1), lwd=c(NA, 2))
}

#' @rdname zpssFit
#' @export
print.zpssR <- function(x, ...){
  cat('Call:\n')
  print(x[['call']])
  cat('\n')
  cat('Initial Values:\n')
  cat(sprintf('Alpha: %s\n', format(eval(x[['call']]$init_alpha), digits = 3)))
  cat(sprintf('Lambda: %s\n', format(eval(x[['call']]$init_beta), digits = 3)))
  cat('\n')
  cat('Coefficients:\n')
  print(coef(x))
  cat('\n')
  cat('Metrics:\n')
  cat(sprintf('Log-likelihood: %s\n', logLik(x)))
  cat(sprintf('AIC: %s\n', AIC(x)))
  cat(sprintf('BIC: %s\n', BIC(x)))
}

#' @rdname zpssFit
#' @export
summary.zpssR <- function(object, ...){
  print(object)
  cat('\n')
  cat('Fitted values:\n')
  print(fitted(object))
}

#' @rdname zpssFit
#' @export
logLik.zpssR <- function(object, ...){
  if(!is.na(object[['logLikelihood']]) || !is.null(object[['logLikelihood']])){
    return(object[['logLikelihood']])
  }
  return(NA)
}

#' @rdname zpssFit
#' @export
AIC.zpssR <- function(object, ...){
  aic <- .get_AIC(object[['logLikelihood']], 2)
  return(aic)
}

#' @rdname zpssFit
#' @export
BIC.zpssR <- function(object, ...){
  dataMatrix <- get(as.character(object[['call']]$data))
  bic <- .get_BIC(object[['logLikelihood']], 2, sum(dataMatrix[, 2]))
  return(bic)
}

