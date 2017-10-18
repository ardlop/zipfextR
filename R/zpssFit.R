.zpss_mle <- function(par, values, freq, truncated = TRUE) {
  alpha <- par[1]
  lambda <- par[2]
  probs <- .panjerRecursion(max(values), alpha, lambda)
  # probs <- getProbsFaster(alpha, lambda, max(values))
  if (truncated) {
    probs <- (probs) / (1 - probs[1])
    probs <- probs[-1]
  }
  - (sum(freq * log(probs[values])))
}

#' ZPSS parameters estimation.
#'
#' For a given count data set,  usually of the type of ranking data or frequencies of frequencies data, estimates the parameters of the MOEZipf distribution by means of the maximum likelihood
#' method.
#' @param data Matrix of count data.
#' @param init_alpha Initial value of \eqn{\alpha} parameter (\eqn{\alpha > 1}).
#' @param init_beta Initial value of \eqn{\beta} parameter (\eqn{\beta \geq 0}).
#' @param isTruncated Logical; if TRUE, the truncated version of the distribution is returned.
#' @param ... Further arguments to the generic functions. In case of the function \emph{moezipfR.fit}
#' the extra arguments are passing to the \link{optim} function.
#' @details
#' The argument \code{data} is a matrix where, for each row, the first column contains a count,
#' and the second column contains its corresponding frequency.
#'
#' The log-likelihood function is computed by means of the following equation:
#'
#'
#' The function \emph{\link{optim}} is used to estimate the parameters.
#' @return Returns a \emph{moezipfR} object composed by the maximum likelihood parameter estimations,
#' their standard deviation, their confidence intervals and the log-likelihood value.
#' @examples
#' data <- rmoezipf(100, 2.5, 1.3)
#' data <- zipfExtR_getDataMatrix(data)
#' obj <- zpssFit(data, 1.001, 0.001)
#' @seealso \code{\link{zipfExtR_getDataMatrix}}, \code{\link{moezipf_getInitialValues}}.
#' @export
zpssFit <- function(data, init_alpha, init_lambda, isTruncated = FALSE, ...){
  tryCatch({
    stats::optim(par = c(init_alpha, init_lambda), .zpss_mle, values = data[, 1], freq = data[, 2],
                  truncated = isTruncated, hessian = T, ...)
  }, error = function(e){
    print(c('Error, e'))
  })
}
