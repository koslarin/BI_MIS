
library(matrixStats)
library(truncdist)


#ridge =1 corresponds to a ridge prior on beta with a shrinkage parameter σθ ∼U(0,1000), 
#ridge = 0 to a normal prior on beta ~ N(0,sigma_sq_theta_prior)

#If you believe the standard empirical covariance matrix to be unstable, we recommend to set unstable_covariance_mat to 1,
#That makes it possible to calculate the inverse matrices outside the loop using a different method, 
#e.g using e.g. clime() function in the clime package
#and input precision matrices (Sigma_inv_dif) directly into the function

AIS_fun = function(mcmc_samples = 4000, burn_in = 2000, R = 100, y, x, zeta_ppd, 
                                alpha_g_pr_y = 0.01, beta_g_pr_y = 0.01,
                                tau = NULL, ridge = 1, 
                                eta_beta = 1, A = 1000, sigma_sq_theta_prior = 1000^2,
                                unstable_covariance_mat = 0, 
                               Sigma_part_inv = NULL, Sigma_part_0_inv = NULL) {
  n <- length(y)
  S <- mcmc_samples
  N <- nrow(zeta_ppd) # Number of available draws
  
  if(is.null(tau)) tau <- 1/runif(1, 0.05, 3)^2
  sigma_y_sq = 1/tau
  
  # 1. Pre-computation (Outside loop)
  mu_part <- colMeans(zeta_ppd)
  zeta_post <- mu_part  # to initialize posterior of zeta, use partial posterior mean
  X_fixed <-  as.matrix(x)
  beta_params <- ncol(X_fixed) + 1
  
  Sigma_part <- cov(zeta_ppd)
  Sigma_part_0_diag <- base::diag(Sigma_part)
  
  if (unstable_covariance_mat == 0) { # regular case
    chol_Sigma  <- chol(Sigma_part)          
    Sigma_part_inv  <- chol2inv(chol_Sigma)
    
    Sigma_part_0_inv_diag <- 1 / Sigma_part_0_diag
    
    # Determinant ratio (log scale for stability)
    log_det_0 <- sum(log(Sigma_part_0_diag))
    log_det_full <- 2 * sum(log(base::diag(chol_Sigma)))
    log_sqrt_det_ratio <- 0.5 * (log_det_0 - log_det_full)
    
    # Kernel for AIS weights
    Sigma_inv_dif <- Sigma_part_inv - base::diag(Sigma_part_0_inv_diag)
  } else { 
    # in the p > n case, the standard empirical covariance matrix becomes unstable
    # In those cases, we recommend we recommend to set unstable_covariance_mat to 1.
    #That makes it possible to calculate the inverse matrices outside the loop using a different method, 
    #e.g using e.g. clime() function in the clime package,
    # then input them directly into the function
    if(is.null(Sigma_part_inv)) stop("If you believe the standard empirical covariance matrix to be unstable, you need to input partial posterior precision matrix manually!")
    if(is.null(Sigma_part_0_inv)) stop("When unstable_covariance_mat = 1, you need to input precision matrix of the diagonal covariance matrix manually!")
    
    log_sqrt_det_ratio = log(sqrt(det(Sigma_part_inv)/det(Sigma_part_0_inv)))
    Sigma_inv_dif <- Sigma_part_inv - Sigma_part_0_inv
  }
  
  
  # Pre-allocate output
  matrix_beta_out <- matrix(NA_real_, S, beta_params - 1)
  matrix_theta_out <- matrix(NA_real_, S, 1)
  matrix_zeta_post <- matrix(NA_real_, S, n)
  matrix_sigma_y_sq <- matrix(NA_real_, S, 1)
  
  # Pre-create y-matrix for dnorm to avoid re-allocating
  Y_mat <- matrix(y, nrow = N, ncol = n, byrow = TRUE)
  
  # initialize beta
  beta_prev = rep(0.001, beta_params)
  
  for (j in 1:(burn_in + S)) {
    
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
    zeta_post <- zeta_iis_all[adj_weight_idx, ]
    
    # --- Regression Update ---
    X <- cbind(X_fixed, zeta_post)
    Xar <- crossprod(X)
    
    # --- Theta and Beta Update ---
    if (ridge == 1) {
      # sigma_beta update
      eta_beta <- truncdist::rtrunc(1, 'gamma', a = 1/A^2, b = Inf,
                                    shape = (beta_params - 1)/2, rate = sum(beta_prev^2)/2)
      Q_B <- tau * Xar + base::diag(eta_beta, beta_params)
    } else {
      Q_B = tau *Xar + base::diag(1/sigma_sq_beta_prior,beta_params)
    }
    l_B = tau *crossprod(X, y)
    ch_Q = chol(Q_B)
    beta_prev <- backsolve(ch_Q, forwardsolve(t(ch_Q), l_B) + rnorm(beta_params))
    
    # --- Tau Update ---
    y_hat_full <- X %*% beta_prev
    tau <- rgamma(1, alpha_g_pr_y + n/2, beta_g_pr_y + 0.5 * sum((y - y_hat_full)^2))
    sigma_y_sq <- 1 / tau
    
    
    # Save results
    if (j > burn_in) {
      idx <- j - burn_in
      matrix_beta_out[idx, ] <- beta_prev[1:(beta_params-1)]
      matrix_theta_out[idx, ] <- theta
      matrix_zeta_post[idx, ] <- zeta_post
      matrix_sigma_y_sq[idx, ] <- sigma_y_sq
    }
  }
  
  return(list(beta = matrix_beta_out, theta = matrix_theta_out, 
              zeta = matrix_zeta_post, sigma_epsilon_sq = matrix_sigma_y_sq))
}

