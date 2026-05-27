#### Simulating Code for Ex1: Independent Scenario

set.seed(2023)
Runs=100

n <- 200
a_0 = 0
sigma_zeta = 1
cor_ex = 0.3
sigma_y = sqrt(2)
b0 = 0
b_z = 4
sigma_u <- 1
sigma_u_sq = sigma_u^2
cor_u = 0.3

zeta_cov_matrix_true = matrix(data = cor_ex*sigma_zeta^2, nrow = n, ncol = n)# 
diag(zeta_cov_matrix_true) = sigma_zeta^2

y_vec_Cor0_merg_mod = matrix(NA, nrow = Runs, ncol = n)
z1_vec_Cor0_merg_mod = matrix(NA, nrow = Runs, ncol = n)

z_cov_matrix = matrix(data= cor_u*sigma_u_sq, nrow=n, ncol=n)
diag(z_cov_matrix) = sigma_u_sq

for (i in 1:Runs) {
  zeta = a_0 +mvrnorm(n=1, mu = rep(0, times=n), Sigma =zeta_cov_matrix_true)
  y_vec_Cor0_merg_mod[i,] <- b0  + b_z*zeta + rnorm(n, sd = sigma_y)
  z1_vec_Cor0_merg_mod[i,] <- zeta + mvrnorm(n=1, mu = rep(0, times=n), Sigma =z_cov_matrix)
}

Stage1_Cor0_zeta_ppd  = array(0, c(Runs,N,n))

## Obtaining Partial Posterior Draws
for (t in 1:Runs) {
  y = y_vec_Cor0_merg_mod[t,]
  z1 = z1_vec_Cor0_merg_mod[t, ]
  zeta_cov_matrix_inv = solve(zeta_cov_matrix_true)
  z_cov_matrix_inv = solve(z_cov_matrix)
  y = y_vec_Cor0_merg_mod[t,]
  z1 = z1_vec_Cor0_merg_mod[t, ]
  
  N=1000
  
  zeta_ppd_draws <- matrix(data = NA, nrow = N, ncol = n)
  Sigma_part = solve(zeta_cov_matrix_inv+z_cov_matrix_inv)
  mu_part = Sigma_part%*%( z_cov_matrix_inv%*%z1 + zeta_cov_matrix_inv%*%mu_zeta)
  
  for (k in 1:N) {
    zeta_part = mvrnorm(n=1, mu = mu_part, Sigma = Sigma_part)
    Stage1_Cor0_zeta_ppd[t,k,] <- zeta_part
  }
}
