#' MOEZipf parameters estimation.
#'
#' For a given sample of strictly positive integer numbers,  usually of the type of ranking data or
#' frequencies of frequencies data, estimates the parameters of the MOEZipf distribution by means of
#' the maximum likelihood method. The input data should be provided as a frequency matrix.
#'
#' @param data Matrix of count data in form of a table of frequencies.
#' @param init_alpha Initial value of \eqn{\alpha} parameter (\eqn{\alpha > 1}).
#' @param init_beta Initial value of \eqn{\beta} parameter (\eqn{\beta > 0}).
#' @param level Confidence level used to calculate the confidence intervals (default 0.95).
#' @param object An object from class "moezipfR" (output of \emph{moezipfFit} function).
#' @param x An object from class "moezipfR" (output of \emph{moezipfFit} function).
#' @param ... Further arguments to the generic functions. The extra arguments are passing to the \emph{\link{optim}} function.
#' @details
#' The argument \code{data} is a two column matrix with the first column containing the observations and
#' the second column containing their frequencies.
#'
#' The log-likelihood function is equal to:
#'
#' \deqn{l(\alpha, \beta; x) = -\alpha \sum_{i = 1} ^m f_{a}(x_{i}) log(x_{i}) + N (log(\beta) + \log(\zeta(\alpha)))}
#' \deqn{ - \sum_{i = 1} ^m f_a(x_i) log[(\zeta(\alpha) - \bar{\beta}\zeta(\alpha, x_i)(\zeta(\alpha) - \bar{\beta}\zeta(\alpha, x_i + 1)))], }
#' where \eqn{f_{a}(x_i)} is the absolute frequency of \eqn{x_i}, \eqn{m} is the number of different values in the sample and \eqn{N} is the sample size,
#' i.e.  \eqn{N = \sum_{i = 1} ^m x_i f_a(x_i)}.
#'
#' By default the initial values of the parameters are computed using the function \code{getInitialValues}.
#'
#' The function \emph{\link{optim}} is used to estimate the parameters.
#' @return Returns a \emph{moezipfR} object composed by the maximum likelihood parameter estimations
#' jointly with their standard deviation and confidence intervals. It also contains
#' the value of the log-likelihood at the maximum likelihood estimator.
#'
#' @examples
#' data <- rmoezipf(100, 2.5, 1.3)
#' data <- as.data.frame(table(data))
#' data[,1] <- as.numeric(levels(data[,1])[data[,1]])
#' initValues <- getInitialValues(data, model='moezipf')
#' obj <- moezipfFit(data, init_alpha = initValues$init_alpha, init_beta = initValues$init_beta)
#' @seealso \code{\link{getInitialValues}}.
#' @importFrom stats AIC BIC coef fitted logLik
#' @export
moezipfFit <- function(data, init_alpha = NULL, init_beta = NULL, level = 0.95, ...){
  Call <- match.call()

  if(is.null(init_alpha) || is.null(init_beta)){
    initValues <- getInitialValues(data, model = 'moezipf')
    init_alpha <- initValues$init_alpha
    init_beta <- initValues$init_beta
  }

  if(!is.numeric(init_alpha) || !is.numeric(init_beta)){
    stop('Wrong intial values for the parameters.')
  }

  tryCatch(
    {
      estResults <- .paramEstimationBase(data, c(init_alpha, init_beta), .mloglikelihood, ...)
      estAlpha <- as.numeric(estResults$results$par[1])
      estBeta <- as.numeric(estResults$results$par[2])
      paramSD <- sqrt(diag(solve(estResults$results$hessian)))
      paramsCI <- .getConfidenceIntervals(paramSD, estAlpha, estBeta, level)

      structure(class = "moezipfR", list(alphaHat = estAlpha,
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

#' @rdname moezipfFit
#' @export
residuals.moezipfR <- function(object, ...){
  dataMatrix <- get(as.character(object[['call']]$data))
  fitted.values <- fitted(object)
  residual.values <- as.numeric(dataMatrix[, 2]) - fitted.values
  return(residual.values)
}

#' @rdname moezipfFit
#' @export
fitted.moezipfR <- function(object, ...) {
  dataMatrix <- get(as.character(object[['call']]$data))
  N <- sum(as.numeric(dataMatrix[, 2]))
  fitted.values <- N*sapply(as.numeric(as.character(dataMatrix[,1])), dmoezipf, alpha = object[['alphaHat']],
                            beta = object[['betaHat']])
  return(fitted.values)
}

#' @rdname moezipfFit
#' @export
coef.moezipfR <- function(object, ...){
  estimation <- matrix(nrow = 2, ncol = 4)
  estimation[1, ] <- c(object[['alphaHat']], object[['alphaSD']], object[['alphaCI']][1], object[['alphaCI']][2])
  estimation[2, ] <- c(object[['betaHat']], object[['betaSD']], object[['betaCI']][1], object[['betaCI']][2])
  colnames(estimation) <- c("MLE", "Std. Dev.", paste0("Inf. ", "95% CI"),
                            paste0("Sup. ", "95% CI"))
  rownames(estimation) <- c("alpha", "beta")
  estimation
}

#' @rdname moezipfFit
#' @export
plot.moezipfR <- function(x, ...){
  dataMatrix <- get(as.character(x[['call']]$data))
  graphics::plot(dataMatrix[,1], dataMatrix[,2], log="xy",
                 xlab="Observation", ylab="Frequency",
                 main="Fitting MOEZipf Distribution", ...)

  graphics::lines(as.numeric(dataMatrix[,1]), fitted(x), col="blue")

  graphics::legend("topright",  legend = c('Observations', 'MOEZipf Distribution'),
                   col=c('black', 'blue'), pch=c(21,NA),
                   lty=c(NA, 1), lwd=c(NA, 2))
}

#' @rdname moezipfFit
#' @export
print.moezipfR <- function(x, ...){
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

#' @rdname moezipfFit
#' @export
summary.moezipfR <- function(object, ...){
  print(object)
  cat('\n')
  cat('Fitted values:\n')
  print(fitted(object))
}

#' @rdname moezipfFit
#' @export
logLik.moezipfR <- function(object, ...){
  if(!is.na(object[['logLikelihood']]) || !is.null(object[['logLikelihood']])){
    return(object[['logLikelihood']])
  }
  return(NA)
}

#' @rdname moezipfFit
#' @export
AIC.moezipfR <- function(object, ...){
  aic <- .get_AIC(object[['logLikelihood']], 2)
  return(aic)
}

#' @rdname moezipfFit
#' @export
BIC.moezipfR <- function(object, ...){
  dataMatrix <- get(as.character(object[['call']]$data))
  bic <- .get_BIC(object[['logLikelihood']], 2, sum(as.numeric(dataMatrix[, 2])))
  return(bic)
}
