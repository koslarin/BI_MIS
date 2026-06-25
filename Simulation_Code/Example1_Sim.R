
#### Simulating Code for Ex1: Correlated Scenario, rho = 0.3

set.seed(2023)
Runs=100

n <- 200
N = 1000
sd_zeta = 1
corr_sim = 0
sd_u = 1
intercept = 0
theta_zeta = 4
sd_y = sqrt(2)


# matrix of true zeta
zeta_vec_Cor0_merg_mod = matrix(NA, nrow = Runs, ncol = n)
# matrix of y
y_vec_Cor0_merg_mod = matrix(NA, nrow = Runs, ncol = n)
# matrix of z (zeta measured with error)
z_vec_Cor0_merg_mod = matrix(NA, nrow = Runs, ncol = n)
# matrix of ppd draws of zeta
Stage1_Cor0_zeta_ppd  = array(0, c(Runs,N,n))

#obtaining the draws
set.seed(2024)
for (t in 1:Runs) {
  ds = meas_err_model_sim(n = n, n_draws = N, zeta_mean = 0, sd_zeta = sd_zeta, corr_zeta = corr_sim, 
                     sd_u = sd_u, corr_u = corr_sim, beta0= intercept, theta_zeta = theta_zeta, sd_y = sd_y)
  y_vec_Cor0_merg_mod[t,] = ds$y
  z_vec_Cor0_merg_mod[t,] = ds$z
  zeta_vec_Cor0_merg_mod[t,] = ds$zeta
  Stage1_Cor0_zeta_ppd[t,,] = ds$zeta_ppd_draws
  
}

# AIS
AIS_100_Ex1_beta = vector("numeric", length=Runs)
AIS_100_Ex1_theta = vector("numeric", length=Runs)
AIS_100_Ex1_sigma_epsilon_sq = vector("numeric", length=Runs)

for (t in 1:Runs) {
  AIS_Ex1 = AIS_fun(mcmc_samples = 4000, burn_in = 2000, R = 500, y = y_vec_Cor0_merg_mod[t,], x = rep(1,n), 
                           zeta_ppd  = Stage1_Cor0_zeta_ppd[t,,] , 
                           alpha_g_pr_y = 3, beta_g_pr_y = 6,
                           tau = 1/runif(1, 0.05,3)^2, ridge = 0)
  AIS_100_Ex1_beta[t] = mean(AIS_Ex1$beta)
  AIS_100_Ex1_theta[t] = mean(AIS_Ex1$theta)
  AIS_100_Ex1_sigma_epsilon_sq[t] = mean(AIS_Ex1$sigma_epsilon_sq)
}

# IIS
IIS_100_Ex1_beta = vector("numeric", length=Runs)
IIS_100_Ex1_theta = vector("numeric", length=Runs)
IIS_100_Ex1_sigma_epsilon_sq = vector("numeric", length=Runs)

for (t in 1:Runs) {
  IIS_Ex1 = IIS_fun(mcmc_samples = 4000, burn_in = 2000, y = y_vec_Cor0_merg_mod[t,], x = rep(1,n), 
                              zeta_ppd  = Stage1_Cor0_zeta_ppd[t,,] , 
                              alpha_g_pr_y = 3, beta_g_pr_y = 6,
                              tau = 1/runif(1, 0.05,3)^2, ridge = 0)
  IIS_100_Ex1_beta[t] = mean(IIS_Ex1$beta)
  IIS_100_Ex1_theta[t] = mean(IIS_Ex1$theta)
  IIS_100_Ex1_sigma_epsilon_sq[t] = mean(IIS_Ex1$sigma_epsilon_sq)
}

# Plug-in (zeta-hat)

Plug_100_Ex1_beta = vector("numeric", length=Runs)
Plug_100_Ex1_theta = vector("numeric", length=Runs)
Plug_100_Ex1_sigma_epsilon_sq = vector("numeric", length=Runs)

for (t in 1:Runs) {
  Plug_Ex1 = Plug_fun(mcmc_samples = 4000, burn_in = 2000, y = y_vec_Cor0_merg_mod[t,], x = rep(1,n), 
                    zeta_ppd  = colMeans(Stage1_n400_Cor0_matrix[t,,]) , 
                    alpha_g_pr_y = 3, beta_g_pr_y = 6,
                    tau = 1/runif(1, 0.05,3)^2, ridge = 0)
  Plug_100_Ex1_beta[t] = mean(Plug_Ex1$beta)
  Plug_100_Ex1_theta[t] = mean(Plug_Ex1$theta)
  Plug_100_Ex1_sigma_epsilon_sq[t] = mean(Plug_Ex1$sigma_epsilon_sq)
}

# Plug-in (z)

Plug_Z_100_Ex1_beta = vector("numeric", length=Runs)
Plug_Z_100_Ex1_theta = vector("numeric", length=Runs)
Plug_Z_100_Ex1_sigma_epsilon_sq = vector("numeric", length=Runs)

for (t in 1:Runs) {
  Plug_Z_Ex1 = Plug_fun(mcmc_samples = 4000, burn_in = 2000, y = y_vec_Cor0_merg_mod[t,], x = rep(1,n), 
                    zeta_ppd  = z_vec_Cor0_merg_mod[t,] , 
                    alpha_g_pr_y = 3, beta_g_pr_y = 6,
                    tau = 1/runif(1, 0.05,3)^2, ridge = 0)
  Plug_Z_100_Ex1_beta[t] = mean(Plug_Z_Ex1$beta)
  Plug_Z_100_Ex1_theta[t] = mean(Plug_Z_Ex1$theta)
  Plug_Z_100_Ex1_sigma_epsilon_sq[t] = mean(Plug_Z_Ex1$sigma_epsilon_sq)
}


# Partial Posterior Method

PPost_100_Ex1_beta = vector("numeric", length=Runs)
PPost_100_Ex1_theta = vector("numeric", length=Runs)
PPost_100_Ex1_sigma_epsilon_sq = vector("numeric", length=Runs)

for (t in 1:Runs) {
  PPost_Ex1 = PartPost_fun(mcmc_samples = 4000, burn_in = 2000, y = y_vec_Cor0_merg_mod[t,], x = rep(1,n), 
                    zeta_ppd  = Stage1_Cor0_zeta_ppd[t,,] , 
                    alpha_g_pr_y = 3, beta_g_pr_y = 6,
                    tau = 1/runif(1, 0.05,3)^2, ridge = 0)
  PPost_100_Ex1_beta[t] = mean(PPost_Ex1$beta)
  PPost_100_Ex1_theta[t] = mean(PPost_Ex1$theta)
  PPost_100_Ex1_sigma_epsilon_sq[t] = mean(PPost_Ex1$sigma_epsilon_sq)
}



