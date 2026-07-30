#' Zipf-PSS INAR(1) Marginal parameters estimation.
#'
#' TO COMPLETE
#' @importFrom stats optim uniroot qnorm var
#' @importFrom utils str tail
#' @importFrom VGAM zeta
#' @importFrom copula polylog
#' @importFrom graphics lines
#'
#' @export


#' @param params Vector of model parameters.
#' @param data The observed time series data.
#' @param loglik Logical; whether to compute the log-likelihood.
#' @param fit_info List with additional fitting information.
#' @param estimation_method Character; method used for estimation.
#' @param estimated_values Optional; previously estimated values for initialization.
ZipfPSS_INAR1Marginal <- function(params = NULL, data = NULL, loglik = NULL, fit_info = NULL, estimation_method = NULL, estimated_values = NULL) {
  obj <- list(
    model_name = "Zipf-PSS INAR(1) Marginal Model",
    estimation_method = estimation_method,
    estimated_values =estimated_values,
    params     = params,     # named vector or list
    data       = data,       # input data
    loglik     = loglik,     # log-likelihood
    fit_info   = fit_info    # convergence info, iterations, etc.
  )
  class(obj) <- "ZipfPSS_INAR1Marginal"
  return(obj)
}

# ---------------------------------------------------------------------------
# Innovations PMF via recursion (from recursion_teo11.R)
# ---------------------------------------------------------------------------

# lambda* = lambda * Li_alpha(1 - beta) / zeta(alpha)
.lambda_star <- function(lambda, beta, alpha) {
  li_value <- copula::polylog(1 - beta, alpha, method = "sum", n.sum = 100)
  lambda * li_value / VGAM::zeta(alpha)
}

# Inner infinite sum S(j, beta, alpha) = sum_{s=j}^{inf} C(s,j)*(1-beta)^(s-j)/s^alpha
.inner_sum <- function(j, beta, alpha, tol = 1e-3) {
  total <- 0
  prev  <- -Inf
  s     <- j
  while (TRUE) {
    term  <- choose(s, j) * (1 - beta)^(s - j) / s^alpha
    total <- total + term
    if (abs(total - prev) < tol) break
    prev <- total
    s    <- s + 1
  }
  total
}

# Compute the full innovations PMF vector P(epsilon = 0), ..., P(epsilon = max_k)
# Returns a numeric vector of length max_k + 1 (0-indexed via +1 offset)
.precompute_prob_eps <- function(alpha, beta, lambda, max_k, tol_inner = 1e-3) {
  lambda_star_val <- .lambda_star(lambda, beta, alpha)
  li_val          <- copula::polylog(1 - beta, alpha, method = "sum", n.sum = 100)

  pmf    <- numeric(max_k + 1)
  pmf[1] <- exp(-lambda_star_val)          # P(epsilon = 0)

  for (kk in seq_len(max_k)) {
    total <- 0
    for (j in seq_len(kk)) {
      outer_coef <- (lambda_star_val * j) / (kk * li_val)
      bracket    <- (1 / j^alpha) - beta^j * .inner_sum(j, beta, alpha, tol = tol_inner)
      total      <- total + outer_coef * bracket * pmf[kk - j + 1]
    }
    pmf[kk + 1] <- total
  }

  pmf
}

.getTransitionsProbs_vectorized <- function(k_values, l_values,
                                            beta, alpha, lambda,
                                            prob_eps_vector) {
  mapply(function(k, l) {
    j <- 0:min(k, l)
    sum(choose(l, j) * beta^j * (1 - beta)^(l - j) *
          prob_eps_vector[k - j + 1])
  }, k_values, l_values)
}


.conditional_max_likelihood <- function(params, data) {

  alpha  <- params[1]
  beta   <- params[2]
  lambda <- params[3]

  max_k <- max(data$value)
  prob_eps_vector <- .precompute_prob_eps(alpha, beta, lambda, max_k)

  log_probs <- log(
    .getTransitionsProbs_vectorized(
      data$value, data$previous,
      beta, alpha, lambda,
      prob_eps_vector
    )
  )

  log_probs[!is.finite(log_probs)] <- -1e10
  -sum(log_probs)
}


# --- Fitting function ---
#' @export
zipfPSS_inar1_marginal_fit <- function(df, method = c("MLE", "CLS"), start_params = NULL, ...) {
  method <- match.arg(method)

  if (method == "MLE") {
      res <- optim(
      par = start_params,
      fn = .conditional_max_likelihood,
      data = df,
      method = "L-BFGS-B",
      lower = c(1.2, 0.01, 0.01),
      upper = c(Inf, 0.9999, Inf),
      hessian = TRUE
    )

    E_zipfpss <- zipfpssMean(res$par[['alpha']], res$par[['lambda']])
    estimations <- res$par[['beta']]*df[, 2] + (1-res$par[['beta']]) * E_zipfpss

    model <- ZipfPSS_INAR1Marginal(
      params = res$par,
      data = df,
      estimation_method = method,
      estimated_values = estimations,
      loglik = -res$value,
      fit_info = res
    )

  } else if (method == "CLS") {
    T <- nrow(df)
    X_t <- df[2:T, ]
    X_t1 <- df[1:(T - 1), ]

    step1 <- function(params, X_t, X_t1) {
      beta <- params[1]
      mu_e <- params[2]
      sum((X_t - beta * X_t1 - mu_e)^2)
    }

    init_eta1 <- c(0.5, 1.0)
    result_eta1 <- optim(
      par = init_eta1,
      fn = function(params) step1(params, X_t, X_t1)
    )
    beta_cls <- result_eta1$par[1]
    mu_elpsilon_cls <- result_eta1$par[2]

    step2 <- function(sigma2, beta, mu) {
      first_term <- (X_t - beta * X_t1 - mu)^2
      second_term <-  -beta*(1-beta)*X_t1 - sigma2
      sum((first_term + second_term)^2)
    }

    result_sigma2 <- optim(
      par = c(1),
      fn = function(sigma2) step2(sigma2, beta = beta_cls, mu = mu_elpsilon_cls),
      method = "Brent",
      lower = 0,
      upper = 100)

    sigma2_epsilon_cls <- result_sigma2$par[1]

    f <- function(alpha, value_ratio) {
      zeta(alpha - 2) / zeta(alpha - 1) - value_ratio
    }

    # This method cannot be used if sigma2_epsilon_cls/mu_elpsilon_cls  < 1
    result <- uniroot(f, value_ratio = sigma2_epsilon_cls/mu_elpsilon_cls, interval = c(3, 10))
    alpha_cls <- result$root

    lambda_cls <- mu_elpsilon_cls*zeta(alpha_cls)/zeta(alpha_cls - 1)

    parameters <- round(c(alpha_cls, beta_cls, lambda_cls, sigma2_epsilon_cls, mu_elpsilon_cls), 4)
    names(parameters) <- c('alpha', 'beta', 'lambda', 'sigma', 'mu')

    E_zipfpss <- zipfpssMean(alpha_cls, lambda_cls)
    estimations <- beta_cls*df[, 2] + (1-beta_cls) * lambda_cls * E_zipfpss

    model <- ZipfPSS_INAR1Marginal(
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

# --- Prediction method ---
#' @rdname zipfPSS_INAR_Marginal
#' @method predict ZipfPSS_INAR1Marginal
#' @export
predict.ZipfPSS_INAR1Marginal <- function(object, newdata = NULL, ...) {
  alpha  <- object$params[["alpha"]]
  beta   <- object$params[["beta"]]
  lambda <- object$params[["lambda"]]

  # last observed value in the series
  last_val <- tail(object$data$value, 1)

  .get_expected_value <- function(alpha, beta, lambda){ #Marta confirmar equation
    (1 - beta) * lambda * (zeta(alpha - 1)/zeta(alpha))
  }

  # Forecasting function (Weiss, p.37)
  .get_forecasting <- function(step, last_val, alpha, beta, lambda) {
    E_value <- .get_expected_value(alpha, beta, lambda)
    beta_h  <- beta^step
    beta_h * last_val + E_value * ((1 - beta_h) / (1 - beta))
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
  if(alpha <= 3){
    stop("Conditioned variance cannot be computed: parameter 'alpha' must be greater than 3 for MLE estimation.")
  }
  return(beta*(1 - beta) * l_value + (lambda*(1-beta)/zeta(alpha))*((1+beta)*zeta(alpha-2) - beta*zeta(alpha-1)))
}

#' @rdname zipfPSS_INAR_Marginal
#' @method residuals ZipfPSS_INAR1Marginal
#' @export
residuals.ZipfPSS_INAR1Marginal <- function(object, standarized = FALSE, ...){
  res_values <- object$data[, 2] - object$estimated_values

  if (standarized){
    conditined_variances <- sapply(object$data[,2], function(i, alpha, beta, lambda){
      .get_conditioned_variance(i, alpha, beta, lambda)
    }, alpha= object$params[['alpha']],
    beta = object$params[['beta']],
    lambda = object$params[['lambda']])

    std_residuals <- res_values/sqrt(conditined_variances)
    return(std_residuals)
  }
  return(res_values)
}

# --- Summary method ---
#' @rdname zipfPSS_INAR_Marginal
#' @method summary ZipfPSS_INAR1Marginal
#' @export
summary.ZipfPSS_INAR1Marginal <- function(object, ...) {
  # build a structured list instead of printing directly
  s <- list(
    model = object$model_name,
    coefficients = object$params,
    logLik = object$loglik,
    AIC = AIC(object),
    BIC = BIC(object),
    fit_info = object$fit_info
  )
  class(s) <- "summary.ZipfPSS_INAR1Marginal"
  return(s)
}

# --- Print method for the summary ---
#' @rdname zipfPSS_INAR_Marginal
#' @export
print.summary.ZipfPSS_INAR1Marginal <- function(x, ...) {
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

#' @rdname zipfPSS_INAR_Marginal
#' @method print ZipfPSS_INAR1Marginal
#' @export
print.ZipfPSS_INAR1Marginal <- function(x, ...) {
  s <- summary(x)
  print(s)
  invisible(x)
}

#' @rdname zipfPSS_INAR_Marginal
#' @method AIC ZipfPSS_INAR1Marginal
#' @export
AIC.ZipfPSS_INAR1Marginal <- function(object, ...){
  if (object$estimation_method == 'MLE') {
    k <- length(object$params)
    loglik <- object$loglik
    aic <- .get_AIC(loglik, k)
    return(aic)
  }
  NULL
}

#' @rdname zipfPSS_INAR_Marginal
#' @method BIC ZipfPSS_INAR1Marginal
#' @export
BIC.ZipfPSS_INAR1Marginal <- function(object, ...){
  if (object$estimation_method == 'MLE') {
    k <- length(object$params)
    loglik <- object$loglik
    bic <- .get_BIC(loglik, k, sum(object$data[, 2]))
    return(bic)
  }
  NULL
}

# --- Plot method ---
#' @rdname zipfPSS_INAR_Marginal
#' @method plot ZipfPSS_INAR1Marginal
#' @export
plot.ZipfPSS_INAR1Marginal <- function(x, predict_range = NULL, alpha_CI = 0.05, ...) {
  alpha <- x$params[['alpha']]
  beta <- x$params[['beta']]
  lambda <- x$params[['lambda']]

  plot(x$data[, 1], type = 'l', ...) # expected real value in the first column
  lines(x$estimated_values, type = 'l', col = 'red') #estimated values are leaded in 1 position

  # x$estimation_method == 'MLE' &&
  if (alpha > 3){
    conditined_variances <- sapply(x$data[,2], function(i, alpha, beta, lambda){
      .get_conditioned_variance(i, alpha, beta, lambda)
    }, alpha= alpha, beta = beta, lambda = lambda)

    lowerIC <- x$estimated_values - qnorm(1 - alpha_CI)*sqrt(conditined_variances)
    upperIC <- x$estimated_values + qnorm(1 - alpha_CI)*sqrt(conditined_variances)

    lines(lowerIC, type = 'l', col = 'green')
    lines(upperIC, type = 'l', col = 'green')
  }
}


