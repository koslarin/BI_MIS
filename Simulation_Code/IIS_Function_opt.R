
# R code organization and documentation

### Ridge version
IIS_normal_optimized2 = function(mcmc_samples = 4000, burn_in = 2000, y, x, zeta_ppd, 
                                 alpha_g_pr_y = 0.01, beta_g_pr_y = 0.01, intercept = TRUE, 
                                 tau = NULL, eta_beta = 1, A = 1000) {
  n <- length(y)
  S <- mcmc_samples
  N <- nrow(zeta_ppd) # 1000 draws
  
  if(is.null(tau)) tau <- 1/runif(1, 0.05, 3)^2
  
  # 1. Pre-computation (Outside loop)
  mu_part <- colMeans(zeta_ppd)
  zeta_sam <- mu_part
  X_fixed <- if(intercept) cbind(rep(1,n), x) else as.matrix(x)
  beta_params <- ncol(X_fixed) + 1
  
  # Pre-allocate output
  matrix_beta_out <- matrix(NA_real_, S, beta_params - 1)
  matrix_theta_out <- matrix(NA_real_, S, 1)
  matrix_zeta_post <- matrix(NA_real_, S, n)
  matrix_sigma_y_sq <- matrix(NA_real_, S, 1)
  
  # Pre-create y-matrix for dnorm to avoid re-allocating
  Y_mat <- matrix(y, nrow = N, ncol = n, byrow = TRUE)
  
  for (j in 1:(burn_in + S)) {
    # --- Regression Update ---
    X <- cbind(X_fixed, zeta_sam)
    Xar <- crossprod(X)
    Q_B <- tau * Xar + diag(eta_beta, beta_params)
    l_B <- tau * crossprod(X, y)
    ch_Q <- chol(Q_B)
    beta_prev <- backsolve(ch_Q, forwardsolve(t(ch_Q), l_B) + rnorm(beta_params))
    
    y_hat_full <- X %*% beta_prev
    tau <- rgamma(1, alpha_g_pr_y + n/2, beta_g_pr_y + 0.5 * sum((y - y_hat_full)^2))
    sigma_y_sq <- 1 / tau
    eta_beta <- truncdist::rtrunc(1, 'gamma', a = 1/A^2, b = Inf,
                                  shape = (beta_params - 1)/2, rate = sum(beta_prev^2)/2)
    
    # --- LAYER 1: Vectorized Likelihood ---
    y_hat_no_zeta <- as.vector(X_fixed %*% beta_prev[1:(beta_params-1)])
    theta <- beta_prev[beta_params]
    
    means_matrix <- sweep(zeta_ppd * theta, 2, y_hat_no_zeta, "+")
    log_lik <- dnorm(means_matrix, mean = Y_mat, sd = sqrt(sigma_y_sq), log = TRUE)
    
    # --- LAYER 1: Robust Discrete Sampling ---
    col_maxes <- colMaxs(log_lik)
    prob_zeta <- exp(sweep(log_lik, 2, col_maxes, "-"))
    cum_probs <- colCumsums(prob_zeta)
    
    # Matrix to store R draws for each of the n locations
    zeta_iis_indices <- vector("numeric", length = n) #matrix(NA_integer_, nrow = R, ncol = n)
    for (h in 1:n) {
      # Sample R times for location h
      u <- runif(1) * cum_probs[N, h]
      # findInterval is extremely fast binary search
      zeta_iis_indices[h] <- findInterval(u, cum_probs[, h]) + 1
      zeta_sam[h] = zeta_ppd[zeta_iis_indices[h], h]
    }
    
    # Save results
    if (j > burn_in) {
      idx <- j - burn_in
      matrix_beta_out[idx, ] <- beta_prev[1:(beta_params-1)]
      matrix_theta_out[idx, ] <- theta
      matrix_zeta_post[idx, ] <- zeta_sam
      matrix_sigma_y_sq[idx, ] <- sigma_y_sq
    }
  }
  
  return(list(beta = matrix_beta_out, theta = matrix_theta_out, 
              zeta = matrix_zeta_post, sigma_epsilon_sq = matrix_sigma_y_sq))
}


## Non-Ridge Version
  
IIS_normal_opt2_nonridge =  function(mcmc_samples = 4000, burn_in = 2000, y, x, zeta_ppd,
                       sigma_sq_beta_prior = 1000^2, alpha_g_pr_y = 0.01, beta_g_pr_y = 0.01, intercept = TRUE, 
                       tau = 1/runif(1, 0.05,3)^2){
  n <- length(y)
  S <- mcmc_samples
  N <- nrow(zeta_ppd) # 1000 draws
  
  if(is.null(tau)) tau <- 1/runif(1, 0.05, 3)^2
  
  # 1. Pre-computation (Outside loop)
  mu_part <- colMeans(zeta_ppd)
  zeta_sam <- mu_part
  X_fixed <- if(intercept) cbind(rep(1,n), x) else as.matrix(x)
  beta_params <- ncol(X_fixed) + 1
  
  # Pre-allocate output
  matrix_beta_out <- matrix(NA_real_, S, beta_params - 1)
  matrix_theta_out <- matrix(NA_real_, S, 1)
  matrix_zeta_post <- matrix(NA_real_, S, n)
  matrix_sigma_y_sq <- matrix(NA_real_, S, 1)
  
  # Pre-create y-matrix for dnorm to avoid re-allocating
  Y_mat <- matrix(y, nrow = N, ncol = n, byrow = TRUE)
  
  for (j in 1:(burn_in + S)) {
    # --- Regression Update ---
    X <- cbind(X_fixed, zeta_sam)
    Xar <- crossprod(X)
    Q_B = Xar + sigma_y_sq/sigma_sq_beta_prior*diag(ncol(X))
    inv_Q_B = solve(Q_B)
    l_B = t(X) %*%y
    ch_Q = chol(Q_B)
    
    beta_prev = backsolve(ch_Q,    forwardsolve(t(ch_Q), l_B) + rnorm(beta_params))
    
    y_hat_full <- X %*% beta_prev
    tau <- rgamma(1, alpha_g_pr_y + n/2, beta_g_pr_y + 0.5 * sum((y - y_hat_full)^2))
    sigma_y_sq <- 1 / tau

    # --- LAYER 1: Vectorized Likelihood ---
    y_hat_no_zeta <- as.vector(X_fixed %*% beta_prev[1:(beta_params-1)])
    theta <- beta_prev[beta_params]
    
    means_matrix <- sweep(zeta_ppd * theta, 2, y_hat_no_zeta, "+")
    log_lik <- dnorm(means_matrix, mean = Y_mat, sd = sqrt(sigma_y_sq), log = TRUE)
    
    # --- LAYER 1: Robust Discrete Sampling ---
    col_maxes <- colMaxs(log_lik)
    prob_zeta <- exp(sweep(log_lik, 2, col_maxes, "-"))
    cum_probs <- colCumsums(prob_zeta)
    
    # Matrix to store R draws for each of the n locations
    zeta_iis_indices <- vector("numeric", length = n) #matrix(NA_integer_, nrow = R, ncol = n)
    for (h in 1:n) {
      # Sample R times for location h
      u <- runif(1) * cum_probs[N, h]
      # findInterval is extremely fast binary search
      zeta_iis_indices[h] <- findInterval(u, cum_probs[, h]) + 1
      zeta_sam[h] = zeta_ppd[zeta_iis_indices[h], h]
    }
    
    # Save results
    if (j > burn_in) {
      idx <- j - burn_in
      matrix_beta_out[idx, ] <- beta_prev[1:(beta_params-1)]
      matrix_theta_out[idx, ] <- theta
      matrix_zeta_post[idx, ] <- zeta_sam
      matrix_sigma_y_sq[idx, ] <- sigma_y_sq
    }
  }
  
  return(list(beta = matrix_beta_out, theta = matrix_theta_out, 
              zeta = matrix_zeta_post, sigma_epsilon_sq = matrix_sigma_y_sq))
}
