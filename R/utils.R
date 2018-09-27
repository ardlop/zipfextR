# Hurwitz Zeta function ------------------------------
.zeta_x<-function(alpha, x) {
  if(!is.numeric(alpha) || !is.numeric(x)){
    stop("Wrong input values!")
  }

  if(x < 1){
    stop("Error!! The 'x' have to be greater or iqual than one.")
  }

  aux <- 0
  if(x == 1) {
    VGAM::zeta(alpha)
  }
  else {
    VGAM::zeta(alpha) - sum((1:(x-1))^(-alpha))
  }
}

# Likelihood functions -----------------------------

.QValue <- function(alpha, beta, x){
  log(VGAM::zeta(alpha) - (1-beta)*.zeta_x(alpha, x))
}

.mloglikelihood <- function(param, nSize, freq, values){
  a <- param[1]
  b <- param[2]
  -(nSize*log(b) + nSize * log(VGAM::zeta(a)) - a * sum(freq * log(values))
    - sum(freq * sapply(values, .QValue, alpha = a, beta = b))
    - sum(freq * sapply(values + 1, .QValue, alpha = a, beta = b)))
}

.zipf_pmf <- function(k, alpha){
  (k^(-alpha))/VGAM::zeta(alpha)
}

.zeta_Distribution <- function(alpha, nSize, freq, values){
  -( -alpha*sum(freq * log(values)) - nSize*log(VGAM::zeta(alpha)))
}

.getSpectrumValues <- function(data){
  frequencies <- as.numeric(data[, 2])
  values <- as.numeric(data[, 1])
  nSize <- sum(frequencies)
  return(list(values = values, frequencies = frequencies,
              nSize = nSize))
}

.paramEstimationBase <- function(x, initValues, likelihoodFunc, ...){
  result <- NULL

  tryCatch({
    statistics <- list(nSize = sum(as.numeric(x[, 2])),
                       values = as.numeric(x[, 1]), frequencies = as.numeric(x[, 2]))
    result <- stats::optim(par = initValues, likelihoodFunc,
                           nSize = statistics$nSize, freq = statistics$frequencies,
                           values = statistics$values, hessian = TRUE, ...)

    return(list(results = result, stats=statistics))
  },
  error = function(e) {
    print(e$message)
  })
}

# Metrics --------------------

.get_AIC <- function(loglike, K) {
  -2*loglike + 2*K
}

.get_BIC <- function(loglike, K, N) {
  -2*loglike + K*log(N)
}

# Kolmogorov - Smirnov Test ---------------------

# Utils ----------------

.getConfidenceIntervals <- function(paramSD, alpha, beta, level){
  result <- matrix(nrow=2, ncol=2)
  levelCoef <- round(stats::qnorm(1-((1-level)/2)), 2)
  offset <- levelCoef * paramSD
  result[1, ] <- c(alpha - offset[1], alpha + offset[1])
  result[2, ] <- c(beta - offset[2], beta + offset[2])
  colnames(result) <- c('Inf. CI', 'Sup. CI')
  rownames(result) <- c('alpha', 'beta')
  return(result)
}

.getConfidenceIntervalsForZeroInflated <- function(paramSD, alpha, beta, w, level){
  result <- matrix(nrow=3, ncol=2)
  levelCoef <- round(stats::qnorm(1-((1-level)/2)), 2)
  offset <- levelCoef * paramSD
  result[1, ] <- c(alpha - offset[1], alpha + offset[1])
  result[2, ] <- c(beta - offset[2], beta + offset[2])
  result[3, ] <- c(w - offset[3], w + offset[3])
  colnames(result) <- c('Inf. CI', 'Sup. CI')
  rownames(result) <- c('alpha', 'beta', 'w')
  return(result)
}

