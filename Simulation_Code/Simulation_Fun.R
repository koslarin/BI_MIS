#### Simulating Data for a basic model for simulations, where stage one is a measurement error model, and stage two is a simple linear Gaussian model

meas_err_model_sim =  function(n = 200, n_draws = 1000, zeta_mean = 0, sd_zeta = 1, corr_zeta = 0, sd_u = 1, corr_u = 0, 
                                beta0 = 0, theta_zeta = 1,  sd_y = 1) {
  
  sigma_zeta = sd_zeta^2
  sigma_u = sd_u^2
    
  zeta_cov_matrix_true = matrix(data = corr_zeta*sigma_zeta, nrow = n, ncol = n)# 
  diag(zeta_cov_matrix_true) = sigma_zeta
  
  z_cov_matrix = matrix(data= corr_u*sigma_u, nrow=n, ncol=n)
  diag(z_cov_matrix) = sigma_u
  
  zeta = zeta_mean +mvrnorm(n=1, mu = rep(0, times=n), Sigma =zeta_cov_matrix_true)
  y <- beta0  + theta_zeta*zeta + rnorm(n, sd = sd_y^2)
  z <- zeta + mvrnorm(n=1, mu = rep(0, times=n), Sigma =z_cov_matrix)
  
  zeta_cov_matrix_inv = solve(zeta_cov_matrix_true)
  z_cov_matrix_inv = solve(z_cov_matrix)
  
  zeta_ppd_draws <- matrix(data = NA, nrow = n_draws, ncol = n)
  Sigma_part = solve(zeta_cov_matrix_inv+z_cov_matrix_inv)
  mu_part = Sigma_part%*%( z_cov_matrix_inv%*%z + zeta_cov_matrix_inv%*%rep(zeta_mean,n)))
  
  zeta_ppd_draws = mvrnorm(n=n_draws, mu = mu_part, Sigma = Sigma_part)
  
  final_val =  list(y = y, z = z, zeta= zeta, zeta_ppd_draws = zeta_ppd_draws)
  return(final_val)
}
