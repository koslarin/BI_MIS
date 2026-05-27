
library(clime)
library(truncdist)
# AIS version that uses solve to get inverse covariances

AIS_normal =  function(mcmc_samples = 4000, burn_in = 2000, R = 500, y, x, zeta_ppd, 
                             alpha_g_pr_y = 0.01, beta_g_pr_y = 0.01, intercept = TRUE, 
                             tau = 1/runif(1, 0.05,3)^2, eta_beta = 1, A = 1000 ){
  n = length(y)
  S = mcmc_samples
  N = nrow(zeta_ppd)
  mu_part = colMeans(zeta_ppd)
  zeta_sam = mu_part
  if ( intercept == TRUE) { #intercept vs no intercept case 
    X = cbind(rep(1, n), x, zeta_sam) 
  } else  {
    X = cbind(x, zeta_sam) 
  } 
  
  beta_params = ncol(X)
  
  Sigma_part = cov(zeta_ppd)
  Sigma_part_0 = matrix(data=0, nrow=n, ncol=n)
  diag(Sigma_part_0) = diag(Sigma_part)
  Sigma_part_inv = solve(Sigma_part)
  Sigma_part_0_inv = solve(Sigma_part_0)
  det_Sigma_part_0 = det(Sigma_part_0)
  det_Sigma_part = det(Sigma_part)
  
  
  log_sqrt_det_ratio = log(sqrt(det_Sigma_part_0/det_Sigma_part))  # ratio calculation for the second layer of weights
  Sigma_inv_dif = Sigma_part_inv - Sigma_part_0_inv
  
  matrix_beta_out <- matrix(data = NA, nrow = S, ncol = beta_params-1)
  matrix_theta_out <- matrix(data = NA, nrow = S, ncol = 1)
  matrix_zeta_post <- matrix(data = NA, nrow = S, ncol = n)
  matrix_sigma_y_sq <- matrix(data = NA, nrow = S, ncol = 1)
  log_lik = matrix(data = NA, nrow =N, ncol = n)
  
  log_ratio_AIS_vec = vector("numeric", length=R)
  zeta_iis_all= matrix(NA, nrow = R, ncol = n)
  
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
    sigma_beta = 1/sqrt(eta_beta)
    ##
    
    for (i in 1:N) {
      #sum up log densities to get log-likelihood
      if (beta_params <= 2) {
        log_lik[i,] = dnorm(X[,1:(beta_params-1)]*beta_prev[1:(beta_params-1)] + zeta_ppd[i,]*beta_prev[beta_params], mean = y, 
                            sd= sqrt(sigma_y_sq), log = TRUE) # ALT 
      } else {
        log_lik[i,] = dnorm(X[,1:(beta_params-1)]%*%beta_prev[1:(beta_params-1)] + zeta_ppd[i,]*beta_prev[beta_params], mean = y, 
                            sd= sqrt(sigma_y_sq), log = TRUE) # ALT 
      }
      
      
    }
    
    log_lik_modded = log_lik - 
      apply(log_lik, 2, max) # max by column (i.e. S different values)
    
    
    prob_zeta= exp(log_lik_modded)
    
    
    for (h in 1:length(y)) {
      zeta_iis <- sample(c(1:N), size = R, replace = TRUE, prob = prob_zeta[,h])
      zeta_iis_all[,h] <- zeta_ppd[zeta_iis, h]
    }
    # Second layer of weights
    for (i in 1:R) {
      log_ratio_AIS_vec[i] = log_sqrt_det_ratio -0.5*t((zeta_iis_all[i,] - 
                                                          mu_part))%*%Sigma_inv_dif%*%(zeta_iis_all[i,] - mu_part)
    }
    
    log_ratio_mod = log_ratio_AIS_vec - max(log_ratio_AIS_vec)
    adj_weight_zeta_scale = sample(c(1:R), size=1, prob = exp(log_ratio_mod))
    zeta_sam <- zeta_iis_all[adj_weight_zeta_scale, ]
    
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


# AIS version that uses clime to estimate inverse covariance of the first stage model, 
#for cases where # of available draws S from the partial posterior is smaller than n 

AIS_normal_clime =  function(mcmc_samples = 4000, burn_in = 2000, R = 500, y, x, zeta_ppd, 
                alpha_g_pr_y = NULL, beta_g_pr_y = NULL, intercept = TRUE, 
                 tau = 1/runif(1, 0.05,3)^2, eta_beta = 1, A = 1000 ){
  n = length(y)
  S = mcmc_samples
  N = nrow(zeta_ppd)
  mu_part = colMeans(zeta_ppd)
  zeta_sam = mu_part
  if ( intercept == TRUE) { #intercept vs no intercept case 
    X = cbind(rep(1, n), x, zeta_sam) 
  } else  {
    X = cbind(x, zeta_sam) 
  } 
  
  beta_params = ncol(X)
  
  Sigma_part = cov(zeta_ppd)
  Sigma_part_0 = matrix(data=0, nrow=n, ncol=n)
  diag(Sigma_part_0) = diag(Sigma_part)
  Sigma_part_inv = clime::clime(zeta_ppd,lambda = 0.05, standardize = FALSE)
  Sigma_part_0_inv <- clime::clime(Sigma_part_0, lambda = .05, sigma=TRUE,standardize = FALSE)
  det_Sigma_part_0 = 1/det(Sigma_part_0_inv$Omegalist[[1]])
  det_Sigma_part = 1/det(Sigma_part_inv$Omegalist[[1]])
  
  
  log_sqrt_det_ratio = log(sqrt(det_Sigma_part_0/det_Sigma_part)) # ratio calculation for the second layer of weights
  Sigma_inv_dif = Sigma_part_inv$Omegalist[[1]] - Sigma_part_0_inv$Omegalist[[1]]
  
  matrix_beta_out <- matrix(data = NA, nrow = S, ncol = beta_params-1)
  matrix_theta_out <- matrix(data = NA, nrow = S, ncol = 1)
  matrix_zeta_post <- matrix(data = NA, nrow = S, ncol = n)
  matrix_sigma_y_sq <- matrix(data = NA, nrow = S, ncol = 1)
  log_lik = matrix(data = NA, nrow =N, ncol = n)
  
  log_ratio_AIS_vec = vector("numeric", length=R)
  zeta_iis_all= matrix(NA, nrow = R, ncol = n)
  
  for (j in 1:(burn_in+S)) { 
    
    X[, ncol(X)] <- zeta_sam
    
    # Using Partial Posterior to Get Zetas
    
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
    sigma_beta = 1/sqrt(eta_beta)
    ##
    
    for (i in 1:N) {
      #sum up log densities to get log-likelihood for first layer of weights
      if (beta_params <= 2) {
        log_lik[i,] = dnorm(X[,1:(beta_params-1)]*beta_prev[1:(beta_params-1)] + zeta_ppd[i,]*beta_prev[beta_params], mean = y, 
                            sd= sqrt(sigma_y_sq), log = TRUE) # ALT 
      } else {
        log_lik[i,] = dnorm(X[,1:(beta_params-1)]%*%beta_prev[1:(beta_params-1)] + zeta_ppd[i,]*beta_prev[beta_params], mean = y, 
                            sd= sqrt(sigma_y_sq), log = TRUE) # ALT 
      }
      
      
    }
    
    log_lik_modded = log_lik - 
      apply(log_lik, 2, max) # max by column (i.e. S different values)
    
    
    prob_zeta= exp(log_lik_modded)
    
    
    for (h in 1:length(y)) {
      zeta_iis <- sample(c(1:N), size = R, replace = TRUE, prob = prob_zeta[,h])
      zeta_iis_all[,h] <- zeta_ppd[zeta_iis, h]
    }
    # second layer of weights
    for (i in 1:R) {
      log_ratio_AIS_vec[i] = log_sqrt_det_ratio -0.5*t((zeta_iis_all[i,] - 
                                                                  mu_part))%*%Sigma_inv_dif%*%(zeta_iis_all[i,] - mu_part)
    }
    
    log_ratio_mod = log_ratio_AIS_vec - max(log_ratio_AIS_vec)
    adj_weight_zeta_scale = sample(c(1:R), size=1, prob = exp(log_ratio_mod))
    zeta_sam <- zeta_iis_all[adj_weight_zeta_scale, ]
    
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

