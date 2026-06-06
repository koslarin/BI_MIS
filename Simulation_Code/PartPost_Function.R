
# R code organization and documentation

### Ridge version

PartPost_fun =  function(mcmc_samples = 4000, burn_in = 2000, y, x, zeta_ppd,
                         alpha_g_pr_y = 0.01, beta_g_pr_y = 0.01, intercept = TRUE, 
                         tau = 1/runif(1, 0.05,3)^2,ridge = 1, sigma_sq_theta_prior = 1000^2) {
  n = length(y)
  S = mcmc_samples
  N = nrow(zeta_ppd)
  zeta_sam = colMeans(zeta_ppd)
  
  if ( intercept == TRUE) { #intercept vs no intercept case 
    X = cbind(rep(1, n), x, zeta_sam) 
  } else  {
    X = cbind(x, zeta_sam) 
  } 
  
  beta_params = ncol(X)
  
  matrix_beta_out <- matrix(data = NA, nrow = S, ncol = beta_params-1)
  matrix_theta_out <- matrix(data = NA, nrow = S, ncol = 1)
  matrix_zeta_post <- matrix(data = NA, nrow = S, ncol = n)
  matrix_sigma_y_sq <- matrix(data = NA, nrow = S, ncol = 1)
  
  for (j in 1:(burn_in+S)) { 
    
    X[, ncol(X)] <- zeta_sam
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
    
    draw = runif(1, 1, N)
    zeta_sam = zeta_ppd[draw,]
    
    # save the results.
    if (j > burn_in) {
      matrix_beta_out[j-burn_in, ] <- beta_prev[1:(ncol(X)-1)]
      matrix_theta_out[j-burn_in, ] <- beta_prev[ncol(X)]
      matrix_zeta_post[j-burn_in, ] <- zeta_sam
      matrix_sigma_y_sq[j - burn_in, ] <- sigma_y_sq
    }
  }
  final_val =  list(beta = matrix_beta_out, theta = matrix_theta_out, zeta = matrix_zeta_post, sigma_epsilon_sq = matrix_sigma_y_sq)
  return(final_val)
}

