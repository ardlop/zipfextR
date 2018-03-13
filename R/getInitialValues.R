.zpssAlphaHat <- function(x, value){
  (VGAM::zeta(x)-value)^2
}


#' Calculates initial values for the parameters of the models.
#'
#' The selection of appropiate initial values to compute the maximum likelihood estimations
#' reduces the number of iterations which in turn, reduces the computation time.
#' The initial values proposed by this function are computed using the first two empirical
#' frequencies.
#'
#' @param data Matrix of count data.
#' @param model Specify the model that requests the initial values (default='zipf').
#'
#' @details
#'
#' The argument \code{data} is a two column matrix with the first column containing the observations and
#' the second column containing their frequencies. The argument \code{model} refers to the selected model of those
#' implemented in the package. The possible values are: \emph{zipf}, \emph{moezipf}, \emph{zipfpe},
#' \emph{zipfpss} or its zero truncated version \emph{zt_zipfpss}. By default, the selected model is the Zipf one.
#'
#' For the MOEZipf, the Zipf-PE and the zero truncated Zipf-PSS models that contain the Zipf model as
#' a particular case, the \eqn{\beta} value will correspond to the one of the Zipf model (i.e. \eqn{\beta = 1} for the MOEZipf,
#' \eqn{\beta = 0} for the Zipf-PE and \eqn{\lambda = 0} for the zero truncated Zipf-PSS model) and the initial value for \eqn{\alpha}
#' is set to be equal to:
#' \deqn{\alpha_0 = log_2 \big (\frac{f_r(1)}{f_r(2)} \big),}
#' where \eqn{f_r(1)} and \eqn{f_r(2)} are the empirical relative frequencies of one and two.
#' This value is obtained equating the two empirical probabilities to their theoritical ones.
#'
#' For the case of the Zipf-PSS the proposed initial values are obtained equating the empirical probability of zero
#' to the theoretical one which gives:
#' \deqn{\lambda_0 = -log(f_r(0)),}
#' where \eqn{f_r(0)} is the empirical relative frequency of zero. The initial value of \eqn{\alpha} is obtained
#' equating the ratio of the theoretical probabilities at zero and one to the empirical ones. This gives place to:
#' \deqn{\alpha_0 = \zeta^{-1}(\lambda_0 * f_r(0)/f_r(1)),}
#' where \eqn{f_r(0)} and \eqn{f_r(1)} are the empirical relative frequencies associated to the values 0 and 1 respectively.
#' The inverse of the Riemman Zeta function is obtained using the \code{optim} routine.
#'
#' @return Returns the initial values of the parameters for a given distribution.
#' @examples
#' data <- rmoezipf(100, 2.5, 1.3)
#' data <- as.data.frame(table(data))
#' data[,1] <- as.numeric(levels(data[,1])[data[,1]])
#' initials <- getInitialValues(data, model='zipf')
#' @references{ Güney, Y., Tuaç, Y., & Arslan, O. (2016). Marshall–Olkin distribution: parameter estimation and
#' application to cancer data. Journal of Applied Statistics, 1-13.}
#' @export
getInitialValues <- function(data, model='zipf'){
  freq1 <- data[which(data[,1] == 1),][2][[1]]
  freq2 <- data[which(data[,1] == 2),][2][[1]]

  if(length(freq1) == 0 || length(freq2) == 0){
    alpha0 <- 1.001
  } else{
    alpha0 <- max(log2(freq1/freq2), 1.001, na.rm = TRUE)
  }

  if(model=='zipf'){
    return(list(init_alpha = round(alpha0, 4)))
  } else if(model=='moezipf'){
    return(list(init_alpha = round(alpha0, 4), init_beta = 1.001))
  } else if(model == 'zipfpe'){
    return(list(init_alpha = round(alpha0, 4), init_beta = 0.001))
  } else if(model == 'zt_zipfpss'){
    return(list(init_alpha = round(alpha0, 4), init_lambda = 0.001))
  } else if(model=='zipfpss'){
    if(length(freq1) == 0 || length(freq2) == 0){
      lambda0 <- 0.001
      alpha0 <- 1.001
    } else{
      lambda0 <- max(-log(freq1/sum(data[,2])), 0.001)
      value <- lambda0*freq1/freq2
      est <- stats::optim(1.01,
            .zpssAlphaHat,
            value = value,
            method = 'L-BFGS-B',
            lower=1.001,
            upper = 30)
      alpha0 <- max(1.001, round(est$par, 4))
    }
    return(list(init_alpha = alpha0, init_lambda = round(lambda0, 4)))
  } else{
    stop('You should introduced a valid model')
  }
}
