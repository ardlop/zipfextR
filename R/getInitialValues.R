.zpssAlphaHat <- function(x, value){
  (VGAM::zeta(x)-value)^2
}


#' Calculates initial values for the \eqn{\alpha} and \eqn{\beta} parameters.
#'
#' The selection of robust initial values allows to reduce the number
#' of iterations which in turn, reduces the computation time. The initial values proposed
#' by this function are computed using the empirical absolute frequencies
#' of the first two consecutive values (i.e. 0 and 1, or, 1 and 2 when the distribution does not
#' contains the zero on its support). In the case where one of these the two positive integer values
#' does not appear in the data set, the default values are set to be equal to the minimum values in the support
#' of the parameters of the model.
#'
#'
#' @param data Matrix of count data.
#' @param model Specify the model that requests the initial values (default='zipf').
#'
#' @details
#'
#' The argument \code{data} is a two column matrix with the first column containing the observations and
#' the second column containing their frequencies. On the other hand, the argument \code{model} refers to
#' all the models implemented in the package. The possible values are \emph{moezipf}, \emph{zipfpe},
#' \emph{zipfpss} or its zero truncated version \emph{zt_zipfpss}.
#'
#' For those models that allow to obtain the Zipf distribution as a particular case, the value of the \eqn{\alpha}
#' parameter is set to be equal to:
#' \deqn{\alpha_0 = log_2 \big (\frac{f_a(1)}{f_a(2)} \big)}
#' where \eqn{f_a(1)} and \eqn{f_a(2)} are the absolute frequencies of one and two in the sample.
#' On the other hand, the second parameter of the model is set to be equal to the value that give place to have the same
#' probabilities as the Zipf distribution.
#'
#' In the particular case of the Zipf-PSS which does not have included the Zipf distribution the initial value proposed for
#' the \eqn{\lambda} parameter is:
#' \deqn{\lambda_0 = -log(f_r(o)),}
#' where \eqn{f_r(0)} is the relative frequency associated to the zero probability.
#' The value of the \eqn{\alpha_0} is obtaining by means of an optimization process and it is set to be equal to the solution
#' of the following equation:
#' \deqn{\alpha_0 = \zeta^{-1}(\lambda_0 * f_a(0)/f_a(1)),}
#' where \eqn{f_a(0)} and \eqn{f_a(1)} are the absolute frequencies associated to the values 0 and 1 respectively.
#'
#' Remark: The user is in charge of considerer to use or not the proposed initial values. In any case it will change the estimation
#' provided for estimating the likelihood parameters of the distribution.
#'
#' @return Returns the initial values of the parameters for a given distribution.
#' @examples
#' data <- rmoezipf(100, 2.5, 1.3)
#' data <- as.data.frame(table(data))
#' data[,1] <- as.numeric(levels(data[,1])[data[,1]])
#' initials <- getInitialValues(data)
#' @references{ Güney, Y., Tuaç, Y., & Arslan, O. (2016). Marshall–Olkin distribution: parameter estimation and
#' application to cancer data. Journal of Applied Statistics, 1-13.}
#' @export
getInitialValues <- function(data, model='zipf'){
  freq1 <- data[which(data[,1] == 1),][2][[1]]
  freq2 <- data[which(data[,1] == 2),][2][[1]]

  alpha0 <- max(log2(freq1/freq2), 1.0001, na.rm = TRUE)
  if(model=='zipf'){
    print('zipf')
    return(list(init_alpha = round(alpha0, 4), init_p2 = NULL))
  } else if(model=='moezipf'){
    print('moezipf')
    return(list(init_alpha = round(alpha0, 4), init_p2 = 1))
  } else if(model == 'zipfpe'){
    print('zipfpe')
    return(list(init_alpha = round(alpha0, 4), init_p2 = 0))
  } else if(model == 'zt_zipfpss'){
    print('zt_zipfpss')
    return(list(init_alpha = round(alpha0, 4), init_p2 = round(0, 4)))
  } else if(model=='zipfpss'){
    lambda0 <- -log(freq1/sum(data[,2]))
    value <- lambda0*freq1/freq2
    est <- stats::optim(1.01,
          .zpssAlphaHat,
          value = value,
          method = 'L-BFGS-B',
          lower=1.001,
          upper = 30)

    return(list(init_alpha = round(est$par, 4), init_p2 = round(lambda0, 4)))

  } else{
    stop('You should introduced a valid model')
  }
  # return(list(init_alpha = round(1.56565, 4), init_beta = 1))

  # if((is.na(freq1) || is.na(freq2)) &&){
  #   return(list(init_alpha = 1.0001, init_beta = 0.0001))
  # }
  #
  # if(is.na(freq1) || is.na(freq2)){
  #   return(list(init_alpha = 1.0001, init_beta = 0.0001))
  # }
  # alpha0 <- max(log2(freq1/freq2), 1.0001, na.rm = TRUE)
  # return(list(init_alpha = round(alpha0, 4), init_beta = 1))
}
