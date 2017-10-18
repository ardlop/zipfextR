.loglikexzp <- function(param, N, values, freq) {
  alpha <- param[1]
  beta <- param[2]
  # print(c(alpha, beta))
  zeta_a <- VGAM::zeta(alpha)
  val1 <-
    sum(sapply(1:length(values), function(i, alpha, values, freq) {
      freq[i] * .zeta_x(alpha, values[i])
    }, alpha = alpha, values = values, freq = freq))

  val2 <- sum(sapply(1:length(values),
                     function(i, alpha, beta, values, freq) {
                       freq[i] * log(exp((beta * values[i] ^ (-alpha)) / (zeta_a)) - 1)
                     },
                     alpha = alpha, beta = beta, values = values, freq = freq))
  - ((beta * (N - zeta_a ^ (-1) * val1) + val2) - (N * log(exp(beta) - 1)))
}



#' ZPE parameters estimation.
#'
#' For a given count data set,  usually of the type of ranking data or frequencies of frequencies data, estimates the parameters of the MOEZipf distribution by means of the maximum likelihood
#' method.
#' @param data Matrix of count data.
#' @param init_alpha Initial value of \eqn{\alpha} parameter (\eqn{\alpha > 1}).
#' @param init_beta Initial value of \eqn{\beta} parameter (\eqn{\beta \in [ -\infity, +\infity]}).
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
#' obj <- zpeFit(data, 1.001, 0.001)
#' @seealso \code{\link{zipfExtR_getDataMatrix}}, \code{\link{moezipf_getInitialValues}}.
#' @export
zpeFit <- function(data, init_alpha, init_beta, ...){
  tryCatch({
    stats::optim(par = c(init_alpha, init_beta), .loglikexzp, N = sum(data[,2]),
                              values = data[, 1], freq = data[, 2], hessian = T, ...)
      },  error = function(e) {
    print(c("Error", e))
  })
}
