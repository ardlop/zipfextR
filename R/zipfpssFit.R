.zpss_mle <- function(par, values, freq, truncated = TRUE) {
  alpha <- par[1]
  lambda <- par[2]
  probs <- .panjerRecursion(max(values), alpha, lambda)
  # probs <- getProbsFaster(alpha, lambda, max(values))
  if (truncated) {
    probs <- (probs) / (1 - probs[1])
    probs <- probs[-1]
    return(- (sum(freq * log(probs[values]))))
  }
  return(- (sum(freq * log(probs[values+1]))))
}

#' Zipf-PSS parameters estimation.
#'
#' For a given sample of strictly positive integer numbers,  usually of the type of ranking data or
#' frequencies of frequencies data, estimates the parameters of the Zipf-PSS distribution by means of
#' the maximum likelihood method. The input data should be provided as a frequency matrix.
#'
#' @param data Matrix of count data in form of table of frequencies.
#' @param init_alpha Initial value of \eqn{\alpha} parameter (\eqn{\alpha > 1}).
#' @param init_lambda Initial value of \eqn{\lambda} parameter (\eqn{\lambda > 0}).
#' @param level Confidence level used to calculate the confidence intervals (default 0.95).
#' @param isTruncated Logical; if TRUE, the truncated version of the distribution is returned.(default = FALSE)
#' @param object An object from class "zpssR" (output of \emph{zipfpssFit} function).
#' @param x An object from class "zpssR" (output of \emph{zipfpssFit} function).
#' @param ... Further arguments to the generic functions. The extra arguments are passing
#' to the \emph{\link{optim}} function.
#'
#' @details
#' The argument \code{data} is a two column matrix with the first column containing the observations and
#' the second column containing their frequencies.
#'
#' The log-likelihood function is equal to:
#' \deqn{l(\alpha, \lambda, x) = \sum_{i =1} ^{m} f_a(x_i)\, log(P(Y = x_i)),}
#' where \eqn{m} is the number of different values in the sample, being \eqn{f_{a}(x_i)} is the absolute
#' frequency of \eqn{x_i}.The probabilities are calculated applying the Panjer recursion.
#' By default the initial values of the parameters are computed using the function \code{getInitialValues}.
#' The function \emph{\link{optim}} is used to estimate the parameters.
#' @return Returns a \emph{zpssR} object composed by the maximum likelihood parameter estimations jointly
#' with their standard deviation and confidence intervals and the value of the log-likelihood at the
#' maximum likelihood estimator.
#' @references {
#' Panjer, H. H. (1981). Recursive evaluation of a family of compound
#' distributions. ASTIN Bulletin: The Journal of the IAA, 12(1), 22-26.
#'
#' Sundt, B., & Jewell, W. S. (1981). Further results on recursive evaluation of
#' compound distributions. ASTIN Bulletin: The Journal of the IAA, 12(1), 27-39.
#'  }
#'
#' @examples
#' data <- rzipfpss(100, 2.5, 1.3)
#' data <- as.data.frame(table(data))
#' data[,1] <- as.numeric(levels(data[,1])[data[,1]])
#' initValues <- getInitialValues(data, model='zipfpss')
#' obj <- zipfpssFit(data, init_alpha = initValues$init_alpha, init_lambda = initValues$init_lambda)
#' @seealso \code{\link{getInitialValues}}.
#' @export
zipfpssFit <- function(data, init_alpha = NULL, init_lambda = NULL, level=0.95, isTruncated = FALSE, ...){
  Call <- match.call()

  if(is.null(init_alpha) || is.null(init_lambda)){
    if(isTruncated){
      model <- 'zt_zipfpss'
    } else {
      model <- 'zipfpss'
    }

    initValues <- getInitialValues(data, model = model)
    init_alpha <- initValues$init_alpha
    init_lambda <- initValues$init_lambda
  }

  if(!is.numeric(init_alpha) || !is.numeric(init_lambda)){
    stop('Wrong intial values for the parameters.')
  }

  tryCatch({
    res <- stats::optim(par = c(init_alpha, init_lambda), .zpss_mle,
                        values = as.numeric(as.character(data[, 1])),
                        freq = as.numeric(data[, 2]),
                        truncated = isTruncated, hessian = TRUE, ...)
    estAlpha <- as.numeric(res$par[1])
    estLambda <- as.numeric(res$par[2])
    paramSD <- sqrt(diag(solve(res$hessian)))
    paramsCI <- .getConfidenceIntervals(paramSD, estAlpha, estLambda, level)

    structure(class = "zipfpssR", list(alphaHat = estAlpha,
                                   lambdaHat = estLambda,
                                   alphaSD = paramSD[1],
                                   lambdaSD = paramSD[2],
                                   alphaCI = c(paramsCI[1,1],paramsCI[1,2]),
                                   lambdaCI = c(paramsCI[2,1],paramsCI[2,2]),
                                   logLikelihood = -res$value,
                                   hessian = res$hessian,
                                   call = Call))

  }, error = function(e){
    print(c('Error', e))
  })
}


#' @rdname zipfpssFit
#' @export
residuals.zipfpssR <- function(object, ...){
  dataMatrix <- get(as.character(object[['call']]$data))
  fitted.values <- fitted(object)
  residual.values <- as.numeric(dataMatrix[, 2]) - fitted.values
  return(residual.values)
}

#' @rdname zipfpssFit
#' @export
fitted.zipfpssR <- function(object, ...) {
  dataMatrix <- get(as.character(object[['call']]$data))
  N <- sum(as.numeric(dataMatrix[, 2]))
  fitted.values <- N*sapply(as.numeric(as.character(dataMatrix[,1])), dzipfpss,
                            alpha = object[['alphaHat']], lambda = object[['lambdaHat']])
  return(fitted.values)
}

#' @rdname zipfpssFit
#' @export
coef.zipfpssR <- function(object, ...){
  estimation <- matrix(nrow = 2, ncol = 4)
  estimation[1, ] <- c(object[['alphaHat']], object[['alphaSD']], object[['alphaCI']][1], object[['alphaCI']][2])
  estimation[2, ] <- c(object[['lambdaHat']], object[['lambdaSD']], object[['lambdaCI']][1], object[['lambdaCI']][2])
  colnames(estimation) <- c("MLE", "Std. Dev.", paste0("Inf. ", "95% CI"),
                            paste0("Sup. ", "95% CI"))
  rownames(estimation) <- c("alpha", "lambda")
  estimation
}

#' @rdname zipfpssFit
#' @export
plot.zipfpssR <- function(x, ...){
  dataMatrix <- get(as.character(x[['call']]$data))
  graphics::plot(as.numeric(as.character(dataMatrix[,1])), as.numeric(dataMatrix[,2]), log="xy",
                 xlab="Observation", ylab="Frequency",
                 main="Fitting Zipf-PSS Distribution", ...)

  graphics::lines(as.numeric(as.character(dataMatrix[,1])), fitted(x), col="blue")

  graphics::legend("topright",  legend = c('Observations', 'Zipf-PSS Distribution'),
                   col=c('black', 'blue'), pch=c(21,NA),
                   lty=c(NA, 1), lwd=c(NA, 2))
}

#' @rdname zipfpssFit
#' @export
print.zipfpssR <- function(x, ...){
  cat('Call:\n')
  print(x[['call']])
  cat('\n')
  cat('Initial Values:\n')
  cat(sprintf('Alpha: %s\n', format(eval(x[['call']]$init_alpha), digits = 3)))
  cat(sprintf('Lambda: %s\n', format(eval(x[['call']]$init_lambda), digits = 3)))
  cat('\n')
  cat('Coefficients:\n')
  print(coef(x))
  cat('\n')
  cat('Metrics:\n')
  cat(sprintf('Log-likelihood: %s\n', logLik(x)))
  cat(sprintf('AIC: %s\n', AIC(x)))
  cat(sprintf('BIC: %s\n', BIC(x)))
}

#' @rdname zipfpssFit
#' @export
summary.zipfpssR <- function(object, ...){
  print(object)
  cat('\n')
  cat('Fitted values:\n')
  print(fitted(object))
}

#' @rdname zipfpssFit
#' @export
logLik.zipfpssR <- function(object, ...){
  if(!is.na(object[['logLikelihood']]) || !is.null(object[['logLikelihood']])){
    return(object[['logLikelihood']])
  }
  return(NA)
}

#' @rdname zipfpssFit
#' @export
AIC.zipfpssR <- function(object, ...){
  aic <- .get_AIC(object[['logLikelihood']], 2)
  return(aic)
}

#' @rdname zipfpssFit
#' @export
BIC.zipfpssR <- function(object, ...){
  dataMatrix <- get(as.character(object[['call']]$data))
  bic <- .get_BIC(object[['logLikelihood']], 2, sum(as.numeric(dataMatrix[, 2])))
  return(bic)
}

