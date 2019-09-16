.zpolyloglikelihood <- function(alpha, betaNew, values, freq, nSize){
  poly <- copula::polylog(z = exp(betaNew), s = alpha, method = 'sum', n.sum = 1000)
  alpha*sum(freq*log(values)) +nSize*log(poly) - betaNew*sum(freq*values)
}

.zPoly <- function(params, values, freq, nSize){
  beta <- min(0, params[2])
  alpha <- params[1]
  # print(c(alpha, beta))
  .zpolyloglikelihood(alpha, beta, values, freq, nSize)
}

#' ZipfPolylog parameters estimation.
#'
#' For a given sample of strictly positive integer numbers,  usually of the type of ranking data or
#' frequencies of frequencies data, estimates the parameters of the ZipfPolylog distribution by means of
#' the maximum likelihood method. The input data should be provided as a frequency matrix.
#'
#' @param data Matrix of count data in form of a table of frequencies.
#' @param init_alpha Initial value of \eqn{\alpha} parameter (\eqn{\alpha > 1}).
#' @param init_beta Initial value of \eqn{\beta} parameter (\eqn{\beta > 0}).
#' @param level Confidence level used to calculate the confidence intervals (default 0.95).
#' @param object An object from class "zipfPolyR" (output of \emph{zipfPolylogFit} function).
#' @param x An object from class "zipfPolyR" (output of \emph{zipfPolylogFit} function).
#' @param ... Further arguments to the generic functions. The extra arguments are passing to the \emph{\link{optim}} function.
#' @details
#' The argument \code{data} is a two column matrix with the first column containing the observations and
#' the second column containing their frequencies.
#'
#' The log-likelihood function is equal to:
#'
#' The function \emph{\link{optim}} is used to estimate the parameters.
#' @return Returns a \emph{zipfPolyR} object composed by the maximum likelihood parameter estimations
#' jointly with their standard deviation and confidence intervals. It also contains
#' the value of the log-likelihood at the maximum likelihood estimator.
#'
#' @importFrom stats AIC BIC coef fitted logLik
#' @importFrom copula polylog
#' @export
zipfPolylogFit <- function(data, init_alpha, init_beta, level = 0.95, ...){
  Call <- match.call()

  if(!is.numeric(init_alpha) || !is.numeric(init_beta)){
    stop('Wrong intial values for the parameters.')
  }

  tryCatch(
    {
      estResults <- .paramEstimationBase(data, c(init_alpha, init_beta), .zPoly, ...)

      estAlpha <- as.numeric(estResults$results$par[1])
      estBeta <- exp(as.numeric(estResults$results$par[2]))
      paramSD <- sqrt(diag(solve(estResults$results$hessian)))
      paramsCI <- .getConfidenceIntervals(paramSD, estAlpha, estBeta, level)

      structure(class = "zipfPolyR", list(alphaHat = estAlpha,
                                         betaHat = estBeta,
                                         alphaSD = paramSD[1],
                                         betaSD = paramSD[2],
                                         alphaCI = c(paramsCI[1,1],paramsCI[1,2]),
                                         betaCI = c(paramsCI[2,1],paramsCI[2,2]),
                                         logLikelihood = -estResults$results$value,
                                         hessian = estResults$results$hessian,
                                         call = Call))
    },
    error=function(cond) {
      print(cond)
      return(NA)
    })
}

#' @rdname zipfPolylogFit
#' @export
residuals.zipfPolyR <- function(object, ...){
  dataMatrix <- get(as.character(object[['call']]$data))
  dataMatrix[,1] <- as.numeric(as.character(dataMatrix[,1]))
  dataMatrix[,2] <-as.numeric(as.character(dataMatrix[,2]))

  residual.values <- dataMatrix[, 2] - fitted(object)
  return(residual.values)
}

#' @rdname zipfPolylogFit
#' @export
fitted.zipfPolyR <- function(object, ...) {
  dataMatrix <- get(as.character(object[['call']]$data))
  dataMatrix[,1] <- as.numeric(as.character(dataMatrix[,1]))
  dataMatrix[,2] <-as.numeric(as.character(dataMatrix[,2]))

  N <- sum(dataMatrix[, 2])
  fitted.values <- N*sapply(dataMatrix[,1], dzipfpolylog, alpha = object[['alphaHat']],
                            beta = object[['betaHat']])
  return(fitted.values)
}

#' @rdname zipfPolylogFit
#' @export
coef.zipfPolyR <- function(object, ...){
  estimation <- matrix(nrow = 2, ncol = 4)
  estimation[1, ] <- c(object[['alphaHat']], object[['alphaSD']], object[['alphaCI']][1], object[['alphaCI']][2])
  estimation[2, ] <- c(object[['betaHat']], object[['betaSD']], object[['betaCI']][1], object[['betaCI']][2])
  colnames(estimation) <- c("MLE", "Std. Dev.", paste0("Inf. ", "95% CI"),
                            paste0("Sup. ", "95% CI"))
  rownames(estimation) <- c("alpha", "beta")
  estimation
}

#' @rdname zipfPolylogFit
#' @export
plot.zipfPolyR <- function(x, ...){
  dataMatrix <- get(as.character(x[['call']]$data))
  dataMatrix[,1] <- as.numeric(as.character(dataMatrix[,1]))
  dataMatrix[,2] <-as.numeric(as.character(dataMatrix[,2]))

  graphics::plot(dataMatrix[,1], dataMatrix[,2], log="xy",
                 xlab="Observation", ylab="Frequency",
                 main="Fitting ZipfPolylog Distribution", ...)

  graphics::lines(dataMatrix[,1], fitted(x), col="blue")

  graphics::legend("topright",  legend = c('Observations', 'ZipfPolylog Distribution'),
                   col=c('black', 'blue'), pch=c(21,NA),
                   lty=c(NA, 1), lwd=c(NA, 2))
}

#' @rdname zipfPolylogFit
#' @export
print.zipfPolyR <- function(x, ...){
  cat('Call:\n')
  print(x[['call']])
  cat('\n')
  cat('Initial Values:\n')
  cat(sprintf('Alpha: %s\n', format(eval(x[['call']]$init_alpha), digits = 3)))
  cat(sprintf('Beta: %s\n', format(eval(x[['call']]$init_beta), digits = 3)))
  cat('\n')
  cat('Coefficients:\n')
  print(coef(x))
  cat('\n')
  cat('Metrics:\n')
  cat(sprintf('Log-likelihood: %s\n', logLik(x)))
  cat(sprintf('AIC: %s\n', AIC(x)))
  cat(sprintf('BIC: %s\n', BIC(x)))
}

#' @rdname zipfPolylogFit
#' @export
summary.zipfPolyR <- function(object, ...){
  print(object)
  cat('\n')
  cat('Fitted values:\n')
  print(fitted(object))
}

#' @rdname zipfPolylogFit
#' @export
logLik.zipfPolyR <- function(object, ...){
  if(!is.na(object[['logLikelihood']]) || !is.null(object[['logLikelihood']])){
    return(object[['logLikelihood']])
  }
  return(NA)
}

#' @rdname zipfPolylogFit
#' @export
AIC.zipfPolyR <- function(object, ...){
  aic <- .get_AIC(object[['logLikelihood']], 2)
  return(aic)
}

#' @rdname zipfPolylogFit
#' @export
BIC.zipfPolyR <- function(object, ...){
  dataMatrix <- get(as.character(object[['call']]$data))
  dataMatrix[,1] <- as.numeric(as.character(dataMatrix[,1]))
  dataMatrix[,2] <-as.numeric(as.character(dataMatrix[,2]))

  bic <- .get_BIC(object[['logLikelihood']], 2, sum(dataMatrix[, 2]))
  return(bic)
}

