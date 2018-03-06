.zpssAlphaHat <- function(x, value){
  (VGAM::zeta(x)-value)^2
}


#' Calculates initial values for the \eqn{\alpha} and \eqn{\beta} parameters.
#'
#' Initial values of the parameters useful to find the maximum likelihood estimators and required
#' in the \emph{moezipfFit} procedure. They may be computed using the empirical absolute frequencies
#' of values one and two. The selection of robust initial values allows to reduce the number
#' of iterations which in turn, reduces the computation time.
#' In the case where one of the two first positive integer values does not appear in the data
#' set, the default values are set to be equal to \eqn{\alpha} = 1.0001 and \eqn{\beta} = 0.0001.
#'
#' @param data Matrix of count data.
#' @param model ....
#'
#' @details
#'
#' The argument \code{data} is a two column matrix such that the first column of each row contains a
#' count, while the corresponding second column contains its frequency.
#'
#' To obtain the initial value of \eqn{\alpha} and \eqn{\beta}, it is assumed that
#' the data come from a Zipf(\eqn{\alpha}) distribution. The initial value for \eqn{\beta}
#' is set to be equal to one, and the inital value for \eqn{\alpha}, denoted by \eqn{\alpha_0}, is obtained
#' equating the ratio of the theoretical probabilities at one and two to the corresponding empirical
#' ratio. this gives:
#'
#' \deqn{\alpha_0 = log_2 \big (\frac{f_a(1)}{f_a(2)} \big)}
#' where \eqn{f_1} and \eqn{f_2} are the absolute frequencies of one and two in the sample.
#' @return Returns the initial value for parameters \eqn{\alpha} and \eqn{\beta}.
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
    lambda0 <- 2*VGAM::zeta(alpha0)
    return(list(init_alpha = round(alpha0, 4), init_p2 = round(lambda0, 4)))
  } else if(model=='zipfpss'){
    lambda0 <- -log(freq1/sum(data[,2]))
    value <- lambda0*freq1/freq2
    est <- stats::optim(1.01,
          .zpssAlphaHat,
          value = value,
          method = 'L-BFGS-B',
          lower=1.001)

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
