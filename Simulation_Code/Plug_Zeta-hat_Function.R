
# R code organization and documentation

### Ridge version

Plug_fun =  function(mcmc_samples = 4000, burn_in = 2000, y, x, plug_in_zeta,
               alpha_g_pr_y = 0.01, beta_g_pr_y = 0.01, intercept = TRUE, 
                tau = 1/runif(1, 0.05,3)^2, ridge = 1, 
               eta_beta = 1, A = 1000, sigma_sq_beta_prior = 1000^2) {
  n = length(y)
  S = mcmc_samples
  
  X = cbind(x, plug_in_zeta) 
  
  beta_params = ncol(X)
  
  matrix_beta_out <- matrix(data = NA, nrow = S, ncol = beta_params-1)
  matrix_theta_out <- matrix(data = NA, nrow = S, ncol = 1)
  matrix_sigma_y_sq <- matrix(data = NA, nrow = S, ncol = 1)
  
  for (j in 1:(burn_in+S)) { 
    
    Xar <- t(X) %*% X
    if (ridge == 1) {
      Q_B <- tau * Xar + base::diag(eta_beta, beta_params)
      # Getting sigma-beta
      eta_beta <- truncdist::rtrunc(1, 'gamma', a = 1/A^2, b = Inf,
                                    shape = (beta_params - 1)/2, rate = sum(beta_prev^2)/2)
    } else {
      Q_B = tau *Xar + base::diag(1/sigma_sq_beta_prior,beta_params)
    }
    l_B = tau *crossprod(X, y)
    ch_Q = chol(Q_B)
    beta_prev <- backsolve(ch_Q, forwardsolve(t(ch_Q), l_B) + rnorm(beta_params))
  
    
    part = (y - X %*% beta_prev)
    
    tau = rgamma(n=1, shape = alpha_g_pr_y + n/2, rate = beta_g_pr_y + 0.5* sum(part^2) )
    sigma_y_sq = 1/tau

    # save the results.
    if (j > burn_in) {
      matrix_beta_out[j-burn_in, ] <- beta_prev[1:(ncol(X)-1)]
      matrix_theta_out[j-burn_in, ] <- beta_prev[ncol(X)]
      matrix_sigma_y_sq[j - burn_in, ] <- sigma_y_sq
    }
  }
  final_val =  list(beta = matrix_beta_out, theta = matrix_theta_out, sigma_epsilon_sq = matrix_sigma_y_sq)
  return(final_val)
}
