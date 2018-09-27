.log_lik_zi_ZipfPSS <-function(values, freq, w, alpha, lambda){
  f_0 <- freq[values == 0]
  t1 <- f_0 * log(w + (1 - w) * exp(-lambda))
  t2 <- log(1 - w) * (sum(freq) - f_0)
  t3 <- sum(freq * dzipfpss(values, alpha, lambda, log = TRUE))
  -(t1 + t2 + t3)
}

.zi_zpss_mle <- function(par, values, freq) {
  alpha <- as.numeric(par[1])
  lambda <- as.numeric(par[2])
  w <- as.numeric(par[3])
  # print(c(w, alpha, lambda))
  return(.log_lik_zi_ZipfPSS(values, freq, w, alpha, lambda))
}

#' Zero Inflated Zipf-PSS parameters estimation.
#'
#' For a given sample of strictly positive integer numbers,  usually of the type of ranking data or
#' frequencies of frequencies data, estimates the parameters of the zero inflated Zipf-PSS distribution by means of
#' the maximum likelihood method. The input data should be provided as a frequency matrix.
#'
#' @param data Matrix of count data in form of table of frequencies.
#' @param init_alpha Initial value of \eqn{\alpha} parameter (\eqn{\alpha > 1}).
#' @param init_lambda Initial value of \eqn{\lambda} parameter (\eqn{\lambda > 0}).
#' @param init_w Initial value of \eqn{w} parameter (\eqn{0 < w < 1}).
#' @param level Confidence level used to calculate the confidence intervals (default 0.95).
#' @param object An object from class "zpssR" (output of \emph{zipfpssFit} function).
#' @param x An object from class "zpssR" (output of \emph{zipfpssFit} function).
#' @param ... Further arguments to the generic functions. The extra arguments are passing
#' to the \emph{\link{optim}} function.
#'
#' @details
#' The argument \code{data} is a two column matrix with the first column containing the observations and
#' the second column containing their frequencies.
#'
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
#' data[,1] <- as.numeric(as.character(data[,1]))
#' data[,2] <- as.numeric(as.character(data[,2]))
#' obj <- zipfpssFit(data, init_alpha = 1.5, init_lambda = 1.5)
#' @seealso \code{\link{getInitialValues}}.
#' @export
zi_zipfpssFit <- function(data, init_alpha = 1.5, init_lambda = 1.5, init_w = 0.1, level=0.95, ...){
  Call <- match.call()

  if(!is.numeric(init_alpha) || !is.numeric(init_lambda) || !is.numeric(init_w)){
    stop('Wrong intial values for the parameters.')
  }

  tryCatch({
    res <- stats::optim(par = c(init_alpha, init_lambda, init_w), .zi_zpss_mle,
                        values = data[, 1],
                        freq = data[, 2],
                        method = "L-BFGS-B",
                        lower = c(1, 1, 0.0001),
                        upper = c(Inf, Inf, 0.9999),
                        hessian = TRUE)

    estAlpha <- as.numeric(res$par[1])
    estLambda <- as.numeric(res$par[2])
    estW <- as.numeric(res$par[3])
    paramSD <- sqrt(diag(solve(res$hessian)))

    paramsCI <- .getConfidenceIntervalsForZeroInflated(paramSD, estAlpha, estLambda, estW, level)

    structure(class = "zi_zipfpssR", list(alphaHat = estAlpha,
                                       lambdaHat = estLambda,
                                       wHat = estW,
                                       alphaSD = paramSD[1],
                                       lambdaSD = paramSD[2],
                                       wSD = paramSD[3],
                                       alphaCI = c(paramsCI[1,1], paramsCI[1,2]),
                                       lambdaCI = c(paramsCI[2,1], paramsCI[2,2]),
                                       wCI = c(paramsCI[3,1], paramsCI[3,2]),
                                       logLikelihood = -res$value,
                                       hessian = res$hessian,
                                       call = Call))

  }, error = function(e){
    print(c('Error', e))
  })
}


#' @rdname zi_zipfpssFit
#' @export
residuals.zi_zipfpssR <- function(object, ...){
  dataMatrix <- get(as.character(object[['call']]$data))
  dataMatrix[,1] <- as.numeric(as.character(dataMatrix[,1]))
  dataMatrix[,2] <-as.numeric(as.character(dataMatrix[,2]))

  fitted.values <- fitted(object)
  residual.values <- dataMatrix[, 2] - fitted.values
  return(residual.values)
}

#' @rdname zi_zipfpssFit
#' @export
fitted.zi_zipfpssR <- function(object,  ...) {
  dataMatrix <- get(as.character(object[['call']]$data))
  dataMatrix[,1] <- as.numeric(as.character(dataMatrix[,1]))
  dataMatrix[,2] <-as.numeric(as.character(dataMatrix[,2]))

  N <- sum(dataMatrix[, 2])

  fitted.values <- N*sapply(dataMatrix[,1], d_zi_zipfpss,
                            alpha = object[['alphaHat']], lambda = object[['lambdaHat']],
                            w = object[['wHat']])
  return(fitted.values)
}

#' @rdname zi_zipfpssFit
#' @export
coef.zi_zipfpssR <- function(object, ...){
  estimation <- matrix(nrow = 3, ncol = 4)
  estimation[1, ] <- c(object[['alphaHat']], object[['alphaSD']], object[['alphaCI']][1], object[['alphaCI']][2])
  estimation[2, ] <- c(object[['lambdaHat']], object[['lambdaSD']], object[['lambdaCI']][1], object[['lambdaCI']][2])
  estimation[3, ] <- c(object[['wHat']], object[['wSD']], object[['wCI']][1], object[['wCI']][2])
  colnames(estimation) <- c("MLE", "Std. Dev.", paste0("Inf. ", "95% CI"),
                            paste0("Sup. ", "95% CI"))
  rownames(estimation) <- c("alpha", "lambda", "w")
  estimation
}

#' @rdname zi_zipfpssFit
#' @export
plot.zi_zipfpssR <- function(x, ...){
  dataMatrix <- get(as.character(x[['call']]$data))
  dataMatrix[,1] <- as.numeric(as.character(dataMatrix[,1]))
  dataMatrix[,2] <-as.numeric(as.character(dataMatrix[,2]))

  graphics::plot(dataMatrix[,1], dataMatrix[,2], log="xy",
                 xlab="Observation", ylab="Frequency",
                 main="Fitting Zipf-PSS Distribution", ...)

  graphics::lines(dataMatrix[,1], fitted(x), col="blue")

  graphics::legend("topright",  legend = c('Observations', 'Zipf-PSS Distribution'),
                   col=c('black', 'blue'), pch=c(21,NA),
                   lty=c(NA, 1), lwd=c(NA, 2))
}

#' @rdname zi_zipfpssFit
#' @export
print.zi_zipfpssR <- function(x, ...){
  cat('Call:\n')
  print(x[['call']])
  cat('\n')
  cat('Initial Values:\n')
  cat(sprintf('alpha: %s\n', format(eval(x[['call']]$init_alpha), digits = 3)))
  cat(sprintf('lambda: %s\n', format(eval(x[['call']]$init_lambda), digits = 3)))
  cat(sprintf('w: %s\n', format(eval(x[['call']]$init_w), digits = 3)))
  cat('\n')
  cat('Coefficients:\n')
  print(coef(x))
  cat('\n')
  cat('Metrics:\n')
  cat(sprintf('Log-likelihood: %s\n', logLik(x)))
  cat(sprintf('AIC: %s\n', AIC(x)))
  cat(sprintf('BIC: %s\n', BIC(x)))
}

#' @rdname zi_zipfpssFit
#' @export
summary.zi_zipfpssR <- function(object, ...){
  print(object)
  cat('\n')
  cat('Fitted values:\n')
  print(fitted(object))
}

#' @rdname zi_zipfpssFit
#' @export
logLik.zi_zipfpssR <- function(object, ...){
  if(!is.na(object[['logLikelihood']]) || !is.null(object[['logLikelihood']])){
    return(object[['logLikelihood']])
  }
  return(NA)
}

#' @rdname zi_zipfpssFit
#' @export
AIC.zi_zipfpssR <- function(object, ...){
  aic <- .get_AIC(object[['logLikelihood']], 3)
  return(aic)
}

#' @rdname zi_zipfpssFit
#' @export
BIC.zi_zipfpssR <- function(object, ...){
  dataMatrix <- get(as.character(object[['call']]$data))
  dataMatrix[,1] <- as.numeric(as.character(dataMatrix[,1]))
  dataMatrix[,2] <-as.numeric(as.character(dataMatrix[,2]))
  bic <- .get_BIC(object[['logLikelihood']], 3, sum(dataMatrix[, 2]))
  return(bic)
}

