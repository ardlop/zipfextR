#' Zipf-PSS INAR(1) Marginal parameters estimation.
#'
#' TO COMPLETE
#' @importFrom stats optim uniroot qnorm var
#' @importFrom utils str tail
#' @importFrom VGAM zeta
#' @importFrom graphics lines
#' @importFrom gmp asNumeric Stirling1

#' @export

# .h_kth_derivative_zero_cache <- new.env(parent = emptyenv())
# .prob_epsilon_cache <- new.env(parent = emptyenv())
# .reset_caches <- function() {
#   if (exists(".prob_epsilon_cache", inherits = FALSE)) {
#     vars_pe <- ls(envir = .prob_epsilon_cache)
#     if (length(vars_pe) > 0) {
#       rm(list = vars_pe, envir = .prob_epsilon_cache)
#     }
#   }
#
#   if (exists(".h_kth_derivative_zero_cache", inherits = FALSE)) {
#     vars_hk <- ls(envir = .h_kth_derivative_zero_cache)
#     if (length(vars_hk) > 0) {
#       rm(list = vars_hk, envir = .h_kth_derivative_zero_cache)
#     }
#   }
# }
#
# .prob_H_new <- function(k, alpha, beta, m_max = 100, tol = 1e-12) {
#   Z_alpha <- zeta(alpha)
#
#   # First term: (1 - β^k) / k^α
#   term1 <- (1 - beta^k) / (k^alpha)
#
#   # Vector of m values
#   m <- 1:m_max
#
#   # Binomial coefficient: choose(m+k, k)
#   binom_part <- choose(m + k, k)
#
#   # (1 - β)^m
#   beta_part <- (1 - beta)^m
#
#   # (m+k)^(-α)
#   power_part <- (m + k)^(-alpha)
#
#   # Summand
#   terms <- binom_part * beta_part * power_part
#
#   # Truncate early if tail becomes negligible
#   cs <- cumsum(terms)
#   last <- cs[length(cs)]
#   # If tail contribution is below tolerance → approximate convergence
#   if (terms[length(terms)] < tol) {
#     S <- last
#   } else {
#     warning("Series may not have converged. Increase m_max.")
#     S <- last
#   }
#
#   # Full expression:
#   # (1 / Z(α)) * [ term1 - β^k * S ]
#   out <- (1 / Z_alpha) * (term1 - (beta^k) * S)
#
#   out
# }


# --- Constructor ---
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

# --- Fitting function ---
zipfPSS_inar1_marginal_fit <- function(df, method = c("MLE", "CLS"), start_params = NULL, ...) {
  method <- match.arg(method)

  if (method == "MLE") {
    # .reset_caches()

    # # -----------------------------
    # # Memoized h_kth_derivative_zero
    # # -----------------------------
    # .h_kth_derivative_zero <- function(k, alpha, beta) {
    #   key <- paste0(k, "_", alpha, "_", beta)
    #   if (exists(key, envir = .h_kth_derivative_zero_cache, inherits = FALSE)) {
    #     return(.h_kth_derivative_zero_cache[[key]])
    #   }
    #   zeta_alpha <- VGAM::zeta(alpha)
    #   k_factorial <- factorial(k)
    #   j <- 0:k
    #   S <- vapply(j, function(x) gmp::asNumeric(gmp::Stirling1(k, x)), numeric(1))
    #   lerch1 <- VGAM::lerch(0, alpha - j, k)
    #   lerch2 <- VGAM::lerch(1 - beta, alpha - j, k) * beta^k
    #   result <- sum(S * (lerch1 - lerch2)) / (k_factorial * zeta_alpha)
    #   .h_kth_derivative_zero_cache[[key]] <- result
    #   return(result)
    # }

    # # -----------------------------
    # # Memoized prob_epsilon
    # # -----------------------------
    # .prob_epsilon <- function(k, alpha, beta, lambda) {
    #   key <- paste0(k, "_", alpha, "_", beta, "_", lambda)
    #   if (exists(key, envir = .prob_epsilon_cache, inherits = FALSE)) {
    #     return(.prob_epsilon_cache[[key]])
    #   }
    #
    #   zeta_alpha <- VGAM::zeta(alpha)
    #   li_value <- copula::polylog(1 - beta, alpha, method = "sum", n.sum = 100)
    #   prob_zero <- exp(-lambda * li_value / zeta_alpha)
    #
    #   if (k == 0) {
    #     .prob_epsilon_cache[[key]] <- prob_zero
    #     return(prob_zero)
    #   }
    #
    #   cache <- numeric(k + 1)
    #   cache[1] <- prob_zero
    #
    #   # h <- vapply(1:k, function(j) .h_kth_derivative_zero(j, alpha, beta), numeric(1))
    #   h <- vapply(1:k, function(j) .prob_H_new(j, alpha, beta), numeric(1))
    #
    #
    #   for (i in 1:k) {
    #     j_vec <- 1:i
    #     cache[i + 1] <- (lambda / k) * sum(j_vec * h[j_vec] * cache[i - j_vec + 1])
    #   }
    #
    #   result <- cache[k + 1]
    #   .prob_epsilon_cache[[key]] <- result
    #   return(result)
    # }

    # # -----------------------------
    # # Precompute prob_epsilon vector
    # # -----------------------------
    # .precompute_prob_eps <- function(alpha, beta, lambda, max_k) {
    #   vapply(0:max_k, function(k) .prob_epsilon(k, alpha, beta, lambda), numeric(1))
    # }


    .prob_H <- function(x, alpha, beta, m_max = 1000, tol = 1e-03) {
      zeta_alpha <- VGAM::zeta(alpha)

      if(x == 0){
        li_value <- copula::polylog(1 - beta, alpha, method = "sum", n.sum = 100)
        return(1 - li_value/zeta_alpha)
      }

      k <- x:m_max
      binom_part <- choose(k, x)
      beta_part <- (1 - beta)^(k-x)
      power_part <- k^(-alpha)
      terms <- binom_part * beta_part * power_part

      # sanity check tolerancia
      cs <- cumsum(terms)
      last <- cs[length(cs)]
      if (terms[length(terms)] < tol) {
        S <- last
      } else {
        warning("Series may not have converged. Increase m_max.")
        S <- last
      }

      return((1/zeta_alpha) * ((1/x^alpha) - (beta^x)*S))
    }

    .prob_epsilon_new <- function(x, alpha, beta, lambda, eps_cache = NULL) {
      if (is.null(eps_cache)) {
        eps_cache <- new.env(parent = emptyenv())
      }

      key <- as.character(x)
      if (exists(key, envir = eps_cache, inherits = FALSE)) {
        return(eps_cache[[key]])
      }

      zeta_alpha <- VGAM::zeta(alpha)
      if (x == 0) {
        li_value <- copula::polylog(1 - beta, alpha, method = "sum", n.sum = 100)
        p0 <- exp(-lambda * li_value / zeta_alpha)
        eps_cache[[key]] <- p0
        return(p0)
      }

      j <- 0:(x - 1)
      x_minus_j <- x - j

      h_vals <- vapply(
        x_minus_j,
        function(v) {
          .prob_H(v, alpha, beta)
        },
        numeric(1)
      )

      eps_vals <- vapply(
        j,
        function(k) {
          .prob_epsilon_new(k, alpha, beta, lambda, eps_cache)
        },
        numeric(1)
      )

      S <- sum((x - j) * h_vals * eps_vals)
      result <- (lambda / x) * S
      eps_cache[[key]] <- result
      return(result)
    }

    # -----------------------------
    # Precompute prob_epsilon vector
    # -----------------------------
    .precompute_prob_eps <- function(alpha, beta, lambda, max_k) {
      vapply(0:max_k, function(k) .prob_epsilon_new(k, alpha, beta, lambda), numeric(1))
    }

    # -----------------------------
    # Vectorized transition probability
    # -----------------------------
    .getTransitionsProbs_vectorized <- function(k_values, l_values, beta, alpha, lambda, prob_eps_vector) {
      mapply(function(k, l) {
        j <- 0:min(k, l)
        binom_part <- choose(l, j)
        beta_part <- beta^j * (1 - beta)^(l - j)
        prob_eps <- prob_eps_vector[k - j + 1]
        sum(binom_part * beta_part * prob_eps)
      }, k_values, l_values)
    }

    # -----------------------------
    # Conditional log-likelihood
    # -----------------------------
    .conditional_max_likelihood <- function(params, data) {
      alpha <- params[1]
      beta <- params[2]
      lambda <- params[3]

      max_k <- max(data$value)
      prob_eps_vector <- .precompute_prob_eps(alpha, beta, lambda, max_k)

      log_probs <- log(.getTransitionsProbs_vectorized(
        k_values = data$value,
        l_values = data$previous,
        beta = beta, alpha = alpha, lambda = lambda,
        prob_eps_vector = prob_eps_vector
      ))

      # replace -Inf if any probability is zero
      log_probs[!is.finite(log_probs)] <- -1e10
      -sum(log_probs)
    }

    # -----------------------------
    # Optimize
    # -----------------------------
    res <- optim(
      par = start_params,
      fn = .conditional_max_likelihood,
      data = df,
      method = "L-BFGS-B",
      lower = c(1.2, 0.01, 0.01),
      upper = c(Inf, 0.9999, Inf)
    )

    print(res)

    E_zipfpss <- zipfpssMean(res$par[['alpha']], res$par[['lambda']])
    estimations <- res$par[['beta']]*df[, 2] + (1-res$par[['beta']]) * E_zipfpss
    print(c(res$par[['beta']], res$par[['lambda']], E_zipfpss))

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
  return(beta*(beta - 1) * l_value + (lambda*(1-beta)/zeta(alpha))*((1+beta)*zeta(alpha-2) - beta*zeta(alpha-1)))
}

#' @rdname zipfPSS_INAR_Marginal
#' @export
residuals.ZipfPSS_INAR1Marginal <- function(object, standarized = FALSE, ...){
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


# --- Summary method ---
#' @rdname zipfPSS_INAR_Marginal
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
#' @export
print.ZipfPSS_INAR1Marginal <- function(x, ...) {
  s <- summary(x)
  print(s)
  invisible(x)
}

#' @rdname zipfPSS_INAR_Marginal
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
#' @export
plot.ZipfPSS_INAR1Marginal <- function(x, predict_range = NULL, alpha_CI = 0.05, ...) {
  alpha <- x$params[['alpha']]
  beta <- x$params[['beta']]
  lambda <- x$params[['lambda']]

  plot(x$data[, 1], type = 'l', ...) # expected real value in the first column
  lines(x$estimated_values, type = 'l', col = 'red') #estimated values are leaded in 1 position

  if (x$estimation_method == 'MLE' && alpha > 3){
    conditined_variances <- sapply(x$data[,2], function(i, alpha, beta, lambda){
      .get_conditioned_variance(i, alpha, beta, lambda)
    }, alpha= alpha, beta = beta, lambda = lambda)

    lowerIC <- x$estimated_values - qnorm(1 - alpha_CI)*sqrt(conditined_variances)
    upperIC <- x$estimated_values + qnorm(1 - alpha_CI)*sqrt(conditined_variances)

    lines(lowerIC, type = 'l', col = 'green')
    lines(upperIC, type = 'l', col = 'green')
  }
}


