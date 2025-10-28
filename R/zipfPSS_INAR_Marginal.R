#' Zipf-PSS INAR(1) Marginal parameters estimation.
#'
#' TO COMPLETE
#' @importFrom stats optim uniroot qnorm var
#' @importFrom utils str tail
#' @importFrom VGAM zeta
#' @importFrom graphics lines
#' @importFrom gmp asNumeric Stirling1

#' @export


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

    .h_kth_derivative_zero <- function(k, alpha, beta){
      k_factorial <- factorial(k)
      term1 <- 1/(k_factorial*VGAM::zeta(alpha))

      term2 <- sapply(0:k, function(j, alpha, beta, k) {
        S_value <- gmp::asNumeric(gmp::Stirling1(k, j))
        lerch_value1 <- VGAM::lerch(0, alpha - j, k)
        lerch_value2 <- VGAM::lerch(1-beta, alpha - j, k)*beta^k
        S_value * (lerch_value1 - lerch_value2)
      }, alpha = alpha, beta = beta, k = k)

      term1 * sum(term2)
    }

    .prob_epsilon <- function(k, alpha, beta, lambda) {
      zeta_alpha <- VGAM::zeta(alpha)
      li_value <- copula::polylog(z = 1 - beta, s = alpha, method = 'sum', n.sum = 100)

      #(ε = 0)
      prob_zero <- exp(-lambda * li_value / zeta_alpha)
      if(k == 0){
        return(prob_zero)
      }

      # Initialize cache for values P(ε = 0) to P(ε = k)
      cache <- numeric(k + 1)
      cache[1] <- prob_zero

      for (i in 1:k) {
        sum_term <- 0
        for (j in 1:i) {
          # P(X = j) via h_function
          pxj <- .h_kth_derivative_zero(j, alpha, beta)
          peps <- cache[i - j + 1]  # P(ε = i - j)
          sum_term <- sum_term +  j * pxj * peps
        }
        cache[i + 1] <- lambda/k * sum_term
      }
      return(cache[k + 1])
    }

    #Equation 10 paper
    .getTransitionsProbs <- function(k_value, l_value, beta, alpha, lambda) {
      tryCatch({
        j <- 0:min(k_value, l_value)
        print('j')
        print(j)

        binom_part <- choose(l_value, j)
        print("bin")
        print(binom_part)
        beta_part <- beta^j * (1 - beta)^(l_value - j)
        print("beta")
        print(beta_part)

        #eq antes section 6
        # print(c(k_value-j, alpha, beta, lambda))
        # prob_epsilon_value <- prob_epsilon(k = k_value-j, alpha = alpha, beta = beta, lambda = lambda)
        prob_epsilon_value <- sapply(j, function(i, k_value, alpha, beta, lambda){
          .prob_epsilon(k = k_value-i, alpha = alpha, beta = beta, lambda = lambda)
        }, k_value=k_value, alpha = alpha, beta = beta, lambda = lambda)
        # print("prob_eps")
        # print(prob_epsilon_value)

        sum(binom_part * beta_part*prob_epsilon_value)
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
        # print(c(k, l, beta, alpha, lambda))
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
      lower = c(1.2, 0.01, 0.01),
      upper = c(Inf, 0.9999, Inf)
    )

    model <- ZipfPSS_INAR1Marginal(
      params   = res$par,
      data     = df,
      estimation_method = method,
      estimated_values = NULL, #To comment with Marta
      loglik   = -res$value,
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

    model <- ZipfPSS_INAR1Innov(
      params   = parameters,
      data     = df,
      estimation_method = method,
      estimated_values = NULL, #to comment with Marta
      loglik   = NULL,
      fit_info = NULL
    )

  }

  return(model)
}
