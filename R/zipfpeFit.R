.loglikexzp <- function(param, N, values, freq) {
  alpha <- param[1]
  beta <- param[2]
  # print(c(alpha, beta))
  zeta_a <- VGAM::zeta(alpha)
  val1 <-
    sum(sapply(base::seq_along(values), function(i, alpha, values, freq) {
    # sum(sapply(1:length(values), function(i, alpha, values, freq) {
      freq[i] * .zeta_x(alpha, values[i])
    }, alpha = alpha, values = values, freq = freq))

  # val2 <- sum(sapply(1:length(values),
  val2 <- sum(sapply(base::seq_along(values),
                     function(i, alpha, beta, values, freq) {
                       freq[i] * log(exp((beta * values[i] ^ (-alpha)) / (zeta_a)) - 1)
                     },
                     alpha = alpha, beta = beta, values = values, freq = freq))
  - ((beta * (N - zeta_a ^ (-1) * val1) + val2) - (N * log(exp(beta) - 1)))
}

.loglik3 <- function(param, N, values, freq){
  alpha <- param[1]
  beta <- param[2]

  -sum(sapply(base::seq_along(values),
              function(i, alpha, beta, values, freq) {
                freq[i] * log(dzipfpe(values[i], alpha, beta))
              }, alpha = alpha, beta = beta, values = values, freq = freq))
}

.loglik4 <- function(param, N, values, freq) {
  alpha <- param[1]
  beta <- param[2]
  val <- sum(sapply(base::seq_along(values),
                    function(i, alpha, beta, values, freq){
                      freq[i] * log((exp((beta*values[i]^(-alpha))/VGAM::zeta(alpha)) -1)/(exp(beta) - 1))
                    }, alpha = alpha, beta = beta, values = values, freq = freq))

  val1 <- sum(sapply(base::seq_along(values),function(i, alpha, values, freq){
    freq[i] * .zeta_x(alpha, values[i])
  }, alpha = alpha, values = values, freq = freq))

  -(beta*(N - VGAM::zeta(alpha)^(-1)*val1) + val)

}

#' Zipf-PE parameters estimation.
#'
#' For a given sample of strictly positive integer values,  usually of the type of ranking data or
#' frequencies of frequencies data, estimates the parameters of the Zipf-PE
#' distribution by means of the maximum likelihood method. The input data should be provided as a frequency matrix.
#'
#' @param data Matrix of count data in form of table of frequencies.
#' @param init_alpha Initial value of \eqn{\alpha} parameter (\eqn{\alpha > 1}).
#' @param init_beta Initial value of \eqn{\beta} parameter (\eqn{\beta \in (-\infty, +\infty)}).
#' @param level Confidence level used to calculate the confidence intervals (default 0.95).
#' @param object An object from class "zpeR" (output of \emph{zipfpeFit} function).
#' @param x An object from class "zpeR" (output of \emph{zipfpeFit} function).
#' @param ... Further arguments to the generic functions.The extra arguments are passing
#' to the \emph{\link{optim}} function.
#' @details
#' The argument \code{data} is a two column matrix with the first column containing the observations and
#' the second column containing their frequencies.
#'
#' The log-likelihood function is equal to:
#'
#' \deqn{l(\alpha, \beta; x) = \beta\, (N - \zeta(\alpha)^{-1}\, \sum_{i = 1} ^m  f_{a}(x_{i})\, \zeta(\alpha, x_i)) +
#' \sum_{i = 1} ^m f_{a}(x_{i})\,  log \left( \frac{e^{\frac{\beta\, x_{i}^{-\alpha}}{\zeta(\alpha)}} - 1}{e^{\beta} - 1} \right), }
#' where \eqn{f_{a}(x_i)} is the absolute frequency of \eqn{x_i}, \eqn{m} is the number of different values in the sample and \eqn{N} is the sample size,
#' i.e.  \eqn{N = \sum_{i = 1} ^m x_i f_a(x_i)}.
#'
#' By default the initial values of the parameters are computed using the function \code{getInitialValues}.
#'
#' The function \emph{\link{optim}} is used to estimate the parameters.
#' @return Returns an object composed by the maximum likelihood parameter estimations
#' jointly with their standard deviation and confidence intervals. It also contains
#' the value of the log-likelihood at the maximum likelihood estimator.
#' @examples
#' data <- rzipfpe(100, 2.5, 1.3)
#' data <- as.data.frame(table(data))
#' data[,1] <- as.numeric(as.character(data[,1]))
#' data[,2] <- as.numeric(as.character(data[,2]))
#' initValues <- getInitialValues(data, model='zipfpe')
#' obj <- zipfpeFit(data, init_alpha = initValues$init_alpha, init_beta = initValues$init_beta)
#' @seealso \code{\link{getInitialValues}}.
#' @export
zipfpeFit <- function(data, init_alpha = NULL, init_beta = NULL, level = 0.95, ...){
  Call <- match.call()

  if(is.null(init_alpha) || is.null(init_beta)){
    initValues <- getInitialValues(data, model = 'zipfpe')
    init_alpha <- initValues$init_alpha
    init_beta <- initValues$init_beta
  }

  if(!is.numeric(init_alpha) || !is.numeric(init_beta)){
    stop('Wrong intial values for the parameters.')
  }

  tryCatch({
    res <- stats::optim(par = c(init_alpha, init_beta), .loglik4, N = sum(as.numeric(data[,2])),
                        values = as.numeric(as.character(data[, 1])), freq = data[, 2],
                        hessian = TRUE, ...)

    estAlpha <- as.numeric(res$par[1])
    estBeta <- as.numeric(res$par[2])
    paramSD <- sqrt(diag(solve(res$hessian)))
    paramsCI <- .getConfidenceIntervals(paramSD, estAlpha, estBeta, level)

        structure(class = "zipfpeR", list(alphaHat = estAlpha,
                                   betaHat = estBeta,
                                   alphaSD = paramSD[1],
                                   betaSD = paramSD[2],
                                   alphaCI = c(paramsCI[1,1],paramsCI[1,2]),
                                   betaCI = c(paramsCI[2,1],paramsCI[2,2]),
                                   logLikelihood = -res$value,
                                   hessian=res$hessian,
                                   call = Call))
  },  error = function(e) {
    print(c("Error", e))
  })
}

#' @rdname zipfpeFit
#' @export
residuals.zipfpeR <- function(object, ...){
  dataMatrix <- get(as.character(object[['call']]$data))
  dataMatrix[,1] <- as.numeric(as.character(dataMatrix[,1]))
  dataMatrix[,2] <-as.numeric(as.character(dataMatrix[,2]))

  residual.values <- dataMatrix[, 2] - fitted(object)
  return(residual.values)
}

#' @rdname zipfpeFit
#' @export
fitted.zipfpeR <- function(object, ...) {
  dataMatrix <- get(as.character(object[['call']]$data))
  dataMatrix[,1] <- as.numeric(as.character(dataMatrix[,1]))
  dataMatrix[,2] <-as.numeric(as.character(dataMatrix[,2]))

  N <- sum(dataMatrix[, 2])
  fitted.values <- N*sapply(dataMatrix[,1], dzipfpe, alpha = object[['alphaHat']],
                            beta = object[['betaHat']])
  return(fitted.values)
}

#' @rdname zipfpeFit
#' @export
coef.zipfpeR <- function(object, ...){
  estimation <- matrix(nrow = 2, ncol = 4)
  estimation[1, ] <- c(object[['alphaHat']], object[['alphaSD']], object[['alphaCI']][1], object[['alphaCI']][2])
  estimation[2, ] <- c(object[['betaHat']], object[['betaSD']], object[['betaCI']][1], object[['betaCI']][2])
  colnames(estimation) <- c("MLE", "Std. Dev.", paste0("Inf. ", "95% CI"),
                            paste0("Sup. ", "95% CI"))
  rownames(estimation) <- c("alpha", "beta")
  estimation
}

#' @rdname zipfpeFit
#' @export
plot.zipfpeR <- function(x, ...){
  dataMatrix <- get(as.character(x[['call']]$data))
  dataMatrix[,1] <- as.numeric(as.character(dataMatrix[,1]))
  dataMatrix[,2] <-as.numeric(as.character(dataMatrix[,2]))

  graphics::plot(dataMatrix[,1], dataMatrix[,2], log="xy",
                 xlab="Observation", ylab="Frequency",
                 main="Fitting Zipf-PE Distribution", ...)

  graphics::lines(dataMatrix[,1], fitted(x), col="blue")

  graphics::legend("topright",  legend = c('Observations', 'Zipf-PE Distribution'),
                   col=c('black', 'blue'), pch=c(21,NA),
                   lty=c(NA, 1), lwd=c(NA, 2))
}

#' @rdname zipfpeFit
#' @export
print.zipfpeR <- function(x, ...){
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

#' @rdname zipfpeFit
#' @export
summary.zipfpeR <- function(object, ...){
  print(object)
  cat('\n')
  cat('Fitted values:\n')
  print(fitted(object))
}

#' @rdname zipfpeFit
#' @export
logLik.zipfpeR <- function(object, ...){
  if(!is.na(object[['logLikelihood']]) || !is.null(object[['logLikelihood']])){
    return(object[['logLikelihood']])
  }
  return(NA)
}

#' @rdname zipfpeFit
#' @export
AIC.zipfpeR <- function(object, ...){
  aic <- .get_AIC(object[['logLikelihood']], 2)
  return(aic)
}

#' @rdname zipfpeFit
#' @export
BIC.zipfpeR <- function(object, ...){
  dataMatrix <- get(as.character(object[['call']]$data))
  dataMatrix[,1] <- as.numeric(as.character(dataMatrix[,1]))
  dataMatrix[,2] <-as.numeric(as.character(dataMatrix[,2]))

  bic <- .get_BIC(object[['logLikelihood']], 2, sum(dataMatrix[, 2]))
  return(bic)
}
