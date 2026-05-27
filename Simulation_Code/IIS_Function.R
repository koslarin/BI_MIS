
# R code organization and documentation

### Ridge version

IIS_normal =  function(mcmc_samples = 4000, burn_in = 2000, y, x, zeta_ppd, likelihood_indicator = 1, 
               alpha_g_pr_y = 0.01, beta_g_pr_y = 0.01, intercept = TRUE, 
                  tau = 1/runif(1, 0.05,3)^2, eta_beta = 1, A= 1000 ){
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
  log_lik = matrix(data = NA, nrow =N, ncol = n)
  
  for (j in 1:(burn_in+S)) { 
    
    X[, ncol(X)] <- zeta_sam

    
    Xar <- t(X) %*% X
    Q_B = tau*Xar + base::diag(eta_beta, beta_params)
    
    l_B = tau*t(X) %*%y
    ch_Q = chol(Q_B)
    
    beta_prev = backsolve(ch_Q,
                          forwardsolve(t(ch_Q), l_B) +
                            rnorm(beta_params))
  
    
    part = (y - X %*% beta_prev)
    
    tau = rgamma(n=1, shape = alpha_g_pr_y + n/2, rate = beta_g_pr_y + 0.5* sum(part^2) )
    sigma_y_sq = 1/tau
    # Getting sigma-beta
    eta_beta = rtrunc(n = 1,
                      'gamma',   # Family of distribution
                      a = 1/A^2, # Lower interval
                      b = Inf,   # Upper interval
                      shape = beta_params/2 - 1/2,
                      rate =  sum(beta_prev^2)/2)
   
    ##
    
    for (i in 1:N) {
      #sum up log densities to get log-likelihood for the weights
      if (beta_params <= 2) {
        log_lik[i,] = dnorm(X[,1:(beta_params-1)]*beta_prev[1:(beta_params-1)] + zeta_ppd[i,]*beta_prev[beta_params], mean = y, 
                                                            sd= sqrt(sigma_y_sq), log = TRUE) 
      } else {
        log_lik[i,] = dnorm(X[,1:(beta_params-1)]%*%beta_prev[1:(beta_params-1)] + zeta_ppd[i,]*beta_prev[beta_params], mean = y, 
                            sd= sqrt(sigma_y_sq), log = TRUE)
      }
      
      
    }
    
    log_lik_modded = log_lik - 
      apply(log_lik, 2, max) # max by column (i.e. S different values)
    
    
    prob_zeta= exp(log_lik_modded)
    
    for (h in 1:length(y)) {
      zeta_n <- sample(c(1:N), size = 1, prob = prob_zeta[,h])
      zeta_sam[h] <- zeta_ppd[zeta_n, h]
    }
    
    
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

## Non-Ridge Version

IIS_normal_nonridge =  function(mcmc_samples = 4000, burn_in = 2000, y, x, zeta_ppd, likelihood_indicator = 1, 
                       sigma_sq_beta_prior = 1000^2, alpha_g_pr_y = 0.01, beta_g_pr_y = 0.01, intercept = TRUE, 
                       tau = 1/runif(1, 0.05,3)^2){
  n = length(y)
  S = mcmc_samples
  N = nrow(zeta_ppd)
  zeta_sam = colMeans(zeta_ppd)
  sigma_y_sq = 1/tau
  
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
  log_lik = matrix(data = NA, nrow =N, ncol = n)
  
  for (j in 1:(burn_in+S)) { 
    
    X[, ncol(X)] <- zeta_sam
    
    
    Xar <- t(X) %*% X
    
    Q_B = Xar + sigma_y_sq/sigma_sq_beta_prior*diag(ncol(X))
    inv_Q_B = solve(Q_B)
    l_B = tau*t(X) %*%y
    ch_Q = chol(Q_B)
    
    beta_prev = mvrnorm(n = 1, inv_Q_B %*% l_B, Sigma = sigma_y_sq*inv_Q_B)
    
    
    part = (y - X %*% beta_prev)
    
    tau = rgamma(n=1, shape = alpha_g_pr_y + n/2, rate = beta_g_pr_y + 0.5* sum(part^2) )
    sigma_y_sq = 1/tau
    # Getting sigma-beta
    
    for (i in 1:N) {
      #sum up log densities to get log-likelihood for the weights
      if (beta_params <= 2) {
        log_lik[i,] = dnorm(X[,1:(beta_params-1)]*beta_prev[1:(beta_params-1)] + zeta_ppd[i,]*beta_prev[beta_params], mean = y, 
                            sd= sqrt(sigma_y_sq), log = TRUE) 
      } else {
        log_lik[i,] = dnorm(X[,1:(beta_params-1)]%*%beta_prev[1:(beta_params-1)] + zeta_ppd[i,]*beta_prev[beta_params], mean = y, 
                            sd= sqrt(sigma_y_sq), log = TRUE)
      }
      
      
    }
    
    log_lik_modded = log_lik - 
      apply(log_lik, 2, max) # max by column (i.e. S different values)
    
    
    prob_zeta= exp(log_lik_modded)
    
    for (h in 1:length(y)) {
      zeta_n <- sample(c(1:N), size = 1, prob = prob_zeta[,h])
      zeta_sam[h] <- zeta_ppd[zeta_n, h]
    }
    
    
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
