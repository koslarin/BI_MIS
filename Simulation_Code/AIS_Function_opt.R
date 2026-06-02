
library(matrixStats)
library(truncdist)

AIS_normal_optimized2 = function(mcmc_samples = 4000, burn_in = 2000, R = 100, y, x, zeta_ppd, 
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
  
  Sigma_part <- cov(zeta_ppd)
  Sigma_part_0_diag <- diag(Sigma_part)
  Sigma_part_inv <- solve(Sigma_part)
  Sigma_part_0_inv_diag <- 1 / Sigma_part_0_diag
  
  # Determinant ratio (log scale for stability)
  log_det_0 <- sum(log(Sigma_part_0_diag))
  log_det_full <- as.numeric(determinant(Sigma_part, logarithm = TRUE)$modulus)
  log_sqrt_det_ratio <- 0.5 * (log_det_0 - log_det_full)
  
  # Kernel for AIS weights
  Sigma_inv_dif <- Sigma_part_inv - diag(Sigma_part_0_inv_diag)
  
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
    zeta_iis_indices <- matrix(NA_integer_, nrow = R, ncol = n)
    for (h in 1:n) {
      # Sample R times for location h
      u <- runif(R) * cum_probs[N, h]
      # findInterval is extremely fast binary search
      zeta_iis_indices[, h] <- findInterval(u, cum_probs[, h]) + 1
    }
    
    # Mapping indices to actual values from zeta_ppd
    # zeta_iis_all is R x n
    zeta_iis_all <- matrix(zeta_ppd[cbind(as.vector(zeta_iis_indices), 
                                          rep(1:n, each = R))], nrow = R, ncol = n)
    
    # --- LAYER 2: Vectorized AIS Weights ---
    v <- sweep(zeta_iis_all, 2, mu_part, "-")
    # Matrix trick: diag(v %*% M %*% v') == rowSums((v %*% M) * v)
    quad_forms <- rowSums((v %*% Sigma_inv_dif) * v)
    log_ratio_AIS_vec <- log_sqrt_det_ratio - 0.5 * quad_forms
    
    # Final Importance Weighting Selection
    log_ratio_mod <- log_ratio_AIS_vec - max(log_ratio_AIS_vec)
    adj_weight_idx <- sample.int(R, 1, prob = exp(log_ratio_mod))
    zeta_sam <- zeta_iis_all[adj_weight_idx, ]
    
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

AIS_normal_optimized2_nr = function(mcmc_samples = 4000, burn_in = 2000, R = 100, y, x, zeta_ppd, 
                                    sigma_sq_beta_prior = 1000^2, alpha_g_pr_y = 0.01, beta_g_pr_y = 0.01, intercept = TRUE, 
                                 tau = NULL) {
  n <- length(y)
  S <- mcmc_samples
  N <- nrow(zeta_ppd) # 1000 draws
  
  if(is.null(tau)) tau <- 1/runif(1, 0.05, 3)^2
  
  # 1. Pre-computation (Outside loop)
  mu_part <- colMeans(zeta_ppd)
  zeta_sam <- mu_part
  X_fixed <- if(intercept) cbind(rep(1,n), x) else as.matrix(x)
  beta_params <- ncol(X_fixed) + 1
  
  Sigma_part <- cov(zeta_ppd)
  Sigma_part_0_diag <- diag(Sigma_part)
  Sigma_part_inv <- solve(Sigma_part)
  Sigma_part_0_inv_diag <- 1 / Sigma_part_0_diag
  
  # Determinant ratio (log scale for stability)
  log_det_0 <- sum(log(Sigma_part_0_diag))
  log_det_full <- as.numeric(determinant(Sigma_part, logarithm = TRUE)$modulus)
  log_sqrt_det_ratio <- 0.5 * (log_det_0 - log_det_full)
  
  # Kernel for AIS weights
  Sigma_inv_dif <- Sigma_part_inv - diag(Sigma_part_0_inv_diag)
  
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
    Q_B <- Xar + sigma_y_sq/sigma_sq_beta_prior*diag(ncol(X))
    inv_Q_B = solve(Q_B)
    l_B <- crossprod(X, y)
    ch_Q <- chol(Q_B)
    beta_prev <- backsolve(ch_Q, forwardsolve(t(ch_Q), l_B) + rnorm(beta_params))
    
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
    zeta_iis_indices <- matrix(NA_integer_, nrow = R, ncol = n)
    for (h in 1:n) {
      # Sample R times for location h
      u <- runif(R) * cum_probs[N, h]
      # findInterval is extremely fast binary search
      zeta_iis_indices[, h] <- findInterval(u, cum_probs[, h]) + 1
    }
    
    # Mapping indices to actual values from zeta_ppd
    # zeta_iis_all is R x n
    zeta_iis_all <- matrix(zeta_ppd[cbind(as.vector(zeta_iis_indices), 
                                          rep(1:n, each = R))], nrow = R, ncol = n)
    
    # --- LAYER 2: Vectorized AIS Weights ---
    v <- sweep(zeta_iis_all, 2, mu_part, "-")
    # Matrix trick: diag(v %*% M %*% v') == rowSums((v %*% M) * v)
    quad_forms <- rowSums((v %*% Sigma_inv_dif) * v)
    log_ratio_AIS_vec <- log_sqrt_det_ratio - 0.5 * quad_forms
    
    # Final Importance Weighting Selection
    log_ratio_mod <- log_ratio_AIS_vec - max(log_ratio_AIS_vec)
    adj_weight_idx <- sample.int(R, 1, prob = exp(log_ratio_mod))
    zeta_sam <- zeta_iis_all[adj_weight_idx, ]
    
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



#For zeta_ppd, rows are draws and columns locations
AIS_normal_optimized = function(mcmc_samples = 4000, burn_in = 2000, R = 500, y, x, zeta_ppd, 
                                alpha_g_pr_y = 0.01, beta_g_pr_y = 0.01, intercept = TRUE, 
                                tau = NULL, eta_beta = 1, A = 1000) {
  n <- length(y)
  S <- mcmc_samples
  N <- nrow(zeta_ppd)
  
  if(is.null(tau)) tau <- 1/runif(1, 0.05, 3)^2
  
  # 1. Pre-computation
  mu_part <- colMeans(zeta_ppd)
  zeta_sam <- mu_part
  X_fixed <- if(intercept) cbind(rep(1,n), x) else as.matrix(x)
  beta_params <- ncol(X_fixed) + 1
  
  Sigma_part <- cov(zeta_ppd)
  Sigma_part_0_diag <- diag(Sigma_part)
  
  Sigma_part_inv <- solve(Sigma_part)
  # Sigma_part_0 is diagonal, so inverse is just 1/diag
  Sigma_part_0_inv_diag <- 1 / Sigma_part_0_diag
  
  # Log Determinant Ratio
  # log(sqrt(det0/det)) = 0.5 * (log(det0) - log(det))
  log_det_0 <- sum(log(Sigma_part_0_diag))
  log_det_full <- as.numeric(determinant(Sigma_part, logarithm = TRUE)$modulus)
  log_sqrt_det_ratio <- 0.5 * (log_det_0 - log_det_full)
  
  # Difference in precision matrices
  Sigma_inv_dif <- Sigma_part_inv - diag(Sigma_part_0_inv_diag)
  
  # Pre-allocate output
  matrix_beta_out <- matrix(NA_real_, S, beta_params - 1)
  matrix_theta_out <- matrix(NA_real_, S, 1)
  matrix_zeta_post <- matrix(NA_real_, S, n)
  matrix_sigma_y_sq <- matrix(NA_real_, S, 1)
  
  # Sequence for matrix indexing
  loc_indices <- 1:n
  
  for (j in 1:(burn_in + S)) {
    # --- Standard Regression Update ---
    X <- cbind(X_fixed, zeta_sam)
    Xar <- crossprod(X)
    Q_B <- tau * Xar + base::diag(eta_beta, beta_params)
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
    log_lik <- dnorm(means_matrix, mean = matrix(y, N, n, byrow=TRUE), 
                     sd = sqrt(sigma_y_sq), log = TRUE)
    
    # --- LAYER 1: Vectorized Sampling (R candidates per location) ---
    col_maxes <- colMaxs(log_lik)
    prob_zeta <- exp(sweep(log_lik, 2, col_maxes, "-"))
    cum_probs <- colCumsums(prob_zeta)
    
    # Generate R indices for each of the n locations
    # Use a vectorized runif and findInterval for speed
    u_mat <- matrix(runif(R * n), nrow = R, ncol = n)
    u_mat <- sweep(u_mat, 2, cum_probs[N, ], "*")
    
    # This matrix contains the indices in zeta_ppd for all R candidates across n locations
    zeta_iis_indices <- apply(u_mat, 1, function(row_u) colSums(row_u > cum_probs) + 1)
    zeta_iis_indices <- t(zeta_iis_indices) # Dimension: R x n
    
    # Extract the actual values
    # Each row is a candidate vector for the spatial field
    zeta_iis_all <- matrix(zeta_ppd[cbind(as.vector(zeta_iis_indices), 
                                          rep(1:n, each = R))], nrow = R, ncol = n, byrow = FALSE)
    
    # --- LAYER 2: Vectorized AIS Weights  --- goal is to minimize memory churn here
    # Goal: Calculate 0.5 * (v %*% Sigma_inv_dif %*% t(v)) for all rows in v
    v <- sweep(zeta_iis_all, 2, mu_part, "-")
    
    # Matrix algebra trick: diag(V M V') = rowSums((V %*% M) * V)
    # This avoids the 'for (i in 1:R)' loop and the matrix transpose operations
    quad_forms <- rowSums((v %*% Sigma_inv_dif) * v)
    log_ratio_AIS_vec <- log_sqrt_det_ratio - 0.5 * quad_forms
    
    # Final Selection
    log_ratio_mod <- log_ratio_AIS_vec - max(log_ratio_AIS_vec)
    adj_weight_idx <- sample.int(R, size = 1, prob = exp(log_ratio_mod))
    zeta_sam <- zeta_iis_all[adj_weight_idx, ]
    
    # --- Store Results ---
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

