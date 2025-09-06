#' Zipf-PSS INAR(1) Innovation parameters estimation.
#'
#' TO COMPLETE
#' @importFrom stats optim uniroot qnorm var
#' @importFrom utils str tail
#' @importFrom VGAM zeta
#' @importFrom graphics lines

#' @export


# --- Constructor ---
ZipfPSS_INAR1Innov <- function(params = NULL, data = NULL, loglik = NULL, fit_info = NULL, estimation_method = NULL, estimated_values = NULL) {
  obj <- list(
    model_name = "Zipf-PSS INAR(1) Innovation Model",
    estimation_method = estimation_method,
    estimated_values =estimated_values,
    params     = params,     # named vector or list
    data       = data,       # input data
    loglik     = loglik,     # log-likelihood
    fit_info   = fit_info    # convergence info, iterations, etc.
  )
  class(obj) <- "ZipfPSS_INAR1Innov"
  return(obj)
}

# --- Fitting function ---
zipfPSS_inar1_innovation_fit <- function(df, method = c("MLE", "CLS"), start_params = NULL, ...) {
  method <- match.arg(method)

  if (method == "MLE") {

    .getTransitionsProbs <- function(k_value, l_value, beta, alpha, lambda) {
      tryCatch({
        j <- 0:min(k_value, l_value)

        binom_part <- choose(l_value, j)
        beta_part <- beta^j * (1 - beta)^(l_value - j)
        dzip_part <- dzipfpss(k_value - j, alpha, lambda)

        sum(binom_part * beta_part * dzip_part)
      }, error = function(e) {
        message("Error: ", e$message)
        return(NA)
      })
    }

    .conditional_max_likelihood <- function(params, data) {
      alpha <- params[1]
      beta <- params[2]
      lambda <- params[3]

      k_vals <- data$value
      l_vals <- data$previous

      log_probs <- mapply(function(k, l) {
        prob <- .getTransitionsProbs(k, l, beta, alpha, lambda)
        if (prob > 0) log(prob) else -Inf
      }, k_vals, l_vals)

      -sum(log_probs, na.rm = TRUE)
    }

    res <- optim(
      par = start_params,
      fn = .conditional_max_likelihood,
      data = df,
      method = "L-BFGS-B",
      lower = c(1.2, 0.01, 0.01), #parametro
      upper = c(Inf, 0.9999, Inf) #parametro
    )
    E_zipfpss <- zipfpssMean(res$par[['alpha']], res$par[['lambda']])
    estimations <- res$par[['beta']]*df[, 2] + E_zipfpss

    model <- ZipfPSS_INAR1Innov(
      params   = res$par,
      data     = df,
      estimation_method = method,
      estimated_values = estimations,
      loglik   = -res$value,
      fit_info = res
    )

  } else if (method == "CLS") {
    T <- nrow(df)
    X_t <- df[2:T, ]
    X_t1 <- df[1:(T - 1), ]

    .step1 <- function(params, X_t, X_t1) {
      beta <- params[1]
      mu_e <- params[2]
      sum((X_t - beta * X_t1 - mu_e)^2)
    }

    # init_eta1 <- c(0.5, 1.0)
    result_eta1 <- optim(
      par = start_params[1:2], #se necesitan tres paramatros para este metodo
      fn = function(params) .step1(params, X_t, X_t1)
    )
    beta_hat <- result_eta1$par[1]
    mu_hat <- result_eta1$par[2]

    .step2 <- function(sigma2, beta, mu) {
      first_term <- (X_t - beta * X_t1 - mu)^2
      second_term <-  -beta*(1-beta)*X_t1 - sigma2
      sum((first_term + second_term)^2)
    }

    result_sigma2 <- optim(
      par = start_params[3],
      fn = function(sigma2) .step2(sigma2, beta = beta_hat, mu = mu_hat),
      method = "Brent",
      lower = 0,
      upper = 100)

    f <- function(alpha, value_ratio) {
      (zeta(alpha - 2) / zeta(alpha - 1)) - value_ratio
    }

    sigma_hat <- result_sigma2$par[1]

    result <- uniroot(f, value_ratio = sigma_hat/mu_hat, interval = c(3, 10))
    alpha_hat <- result$root
    lambda_hat <- mu_hat*zeta(alpha_hat)/zeta(alpha_hat - 1)

    parameters <- c(alpha_hat, beta_hat, lambda_hat, sigma_hat, mu_hat)
    names(parameters) <- c('alpha', 'beta', 'lambda', 'sigma', 'mu')

    E_zipfpss <- zipfpssMean(alpha_hat, lambda_hat)
    estimations <- beta_hat*df[, 2] + E_zipfpss

    model <- ZipfPSS_INAR1Innov(
      params   = parameters,
      data     = df,
      estimation_method = method,
      estimated_values = estimations,
      loglik   = NULL,
      fit_info = NULL
    )
  }

  return(model)
}

#' @rdname zipfPSS_INAR_Innovation
#' @export
residuals.ZipfPSS_INAR1Innov <- function(object, standarized = FALSE, ...){
  res_values <- object$data[, 1] - object$estimated_values

  if (standarized){
    conditined_variances <- sapply(object$data[,2], function(i, alpha, beta, lambda){
      .get_conditioned_variance(i, alpha, beta, lambda)
    }, alpha= object$params[['alpha']],
    beta = object$params[['beta']],
    lambda = object$params[['lambda']])

    std_residuals <- res_values/sqrt(var(conditined_variances))
    return(std_residuals)
  }
  return(res_values)
}


# --- Prediction method ---
#' @rdname zipfPSS_INAR_Innovation
#' @export
predict.ZipfPSS_INAR1Innov <- function(object, newdata = NULL, ...) {
  alpha  <- object$params[["alpha"]]
  beta   <- object$params[["beta"]]
  lambda <- object$params[["lambda"]]

  # last observed value in the series
  last_val <- tail(object$data$value, 1)

  # Forecasting function (Weiss, p.37)
  .get_forecasting <- function(step, last_val, alpha, beta, lambda) {
    E_zipfpss <- zipfpssMean(alpha, lambda)
    beta_h  <- beta^step
    beta_h * last_val + E_zipfpss * ((1 - beta_h) / (1 - beta))
  }

  if (is.null(newdata)) {
    stop("Please provide forecast horizons in `newdata` (e.g., 1:5).")
  }

  predicted <- sapply(newdata, function(h) {
    .get_forecasting(h, last_val, alpha, beta, lambda)
  })

  names(predicted) <- paste0("t+", newdata)
  return(predicted)
}

.get_conditioned_variance <- function(l_value, alpha, beta, lambda){
  return(beta*(beta - 1) * l_value + lambda*(zeta(alpha - 2)/zeta(alpha)))
}



# --- Plot method ---
#' @rdname zipfPSS_INAR_Innovation
#' @export
plot.ZipfPSS_INAR1Innov <- function(x, predict_range = NULL, alpha_CI = 0.05, ...) {
  alpha <- x$params[['alpha']]
  beta <- x$params[['beta']]
  lambda <- x$params[['lambda']]

  #por que no se puede utilizar la función de predict aqui
  # estimations <- predict.ZipfPSS_INAR1Innov(object, object$data[, 1])

  plot(x$data[, 1], type = 'l', ...) # expected real value in the first column
  # E_zipfpss <- zipfpssMean(alpha, lambda)
  # estimations <- beta*x$data[, 2] + E_zipfpss #(se mueve a la clase) #previous value expected in column 2
  lines(x$estimated_values[-1], type = 'l', col = 'red') #estimated values are leaded in 1 position

  if (x$estimation_method == 'MLE'){
    conditined_variances <- sapply(x$data[,2], function(i, alpha, beta, lambda){
      .get_conditioned_variance(i, alpha, beta, lambda)
    }, alpha= alpha, beta = beta, lambda = lambda)
    lowerIC <- x$estimated_values - qnorm(1 - alpha_CI)*sqrt(conditined_variances)
    upperIC <- x$estimated_values + qnorm(1 - alpha_CI)*sqrt(conditined_variances)

    lines(lowerIC[-1], type = 'l', col = 'green')
    lines(upperIC[-1], type = 'l', col = 'green')
  }
}


# --- Summary method ---
#' @rdname zipfPSS_INAR_Innovation
#' @export
summary.ZipfPSS_INAR1Innov <- function(object, ...) {
  # build a structured list instead of printing directly
  s <- list(
    model = object$model_name,
    coefficients = object$params,
    logLik = object$loglik,
    AIC = AIC(object),
    BIC = BIC(object),
    fit_info = object$fit_info
  )
  class(s) <- "summary.ZipfPSS_INAR1Innov"
  return(s)
}

# --- Print method for the summary ---
#' @rdname zipfPSS_INAR_Innovation
#' @export
print.summary.ZipfPSS_INAR1Innov <- function(x, ...) {
  cat("Model:", x$model, "\n\n")

  cat("Coefficients:\n")
  print(x$coefficients)
  cat("\n")

  if (!is.null(x$logLik)) {
    cat("Log-likelihood:", x$logLik, "\n")
    cat("AIC:", x$AIC, "\n")
    cat("BIC:", x$BIC, "\n")
  }

  if (!is.null(x$fit_info)) {
    cat("\nFit information:\n")
    if (is.list(x$fit_info)) {
      str(x$fit_info)
    } else {
      print(x$fit_info)
    }
  }

  invisible(x)
}

#' @rdname zipfPSS_INAR_Innovation
#' @export
print.ZipfPSS_INAR1Innov <- function(x, ...) {
  s <- summary(x)
  print(s)
  invisible(x)
}

#' @rdname zipfPSS_INAR_Innovation
#' @export
AIC.ZipfPSS_INAR1Innov <- function(object, ...){
  if (object$estimation_method == 'MLE') {
    k <- length(object$params)
    loglik <- object$loglik
    aic <- .get_AIC(loglik, k)
    return(aic)
  }
  NULL
}

#' @rdname zipfPSS_INAR_Innovation
#' @export
BIC.ZipfPSS_INAR1Innov <- function(object, ...){
  if (object$estimation_method == 'MLE') {
    k <- length(object$params)
    loglik <- object$loglik
    bic <- .get_BIC(loglik, k, sum(object$data[, 2]))
    return(bic)
  }
  NULL
}

