library(BayesLogit)

# --- 1. Initialization Function ---
init_params_negbin = function(y, x = NULL, zeta_ppd) {
  p_x <- if (is.null(x)) 0 else ncol(as.matrix(x))
  list(
    beta = if (p_x > 0) rep(0, p_x) else numeric(0),
    theta = 1.0,
    r = 2.0
  )
}

# --- 2. Likelihood Function ---
lik_fun_negbin = function(y, x = NULL, zeta_ppd, params) {
  N <- nrow(zeta_ppd)
  n <- ncol(zeta_ppd)
  
  eta_cov <- if (is.null(x)) 0 else as.vector(as.matrix(x) %*% params$beta)
  eta_mat <- sweep(params$theta * zeta_ppd, 2, eta_cov, "+")
  
  mu_mat <- params$r * exp(eta_mat)
  y_mat <- matrix(y, nrow = N, ncol = n, byrow = TRUE)
  
  matrix(
    stats::dnbinom(y_mat, size = params$r, mu = mu_mat, log = TRUE),
    nrow = N, 
    ncol = n
  )
}

# --- Helper: Slice Sampler for Dispersion Parameter r ---
slice_sample_r = function(r_curr, y, eta, e0 = 0.01, a_r = 0.001, b_r = 100) {
  log_post = function(r) {
    if (r <= a_r || r >= b_r) return(-Inf)
    sum(stats::dnbinom(y, size = r, mu = exp(eta), log = TRUE)) - log(1 + (r^2) * e0)
  }
  
  target_y <- log_post(r_curr) - stats::rexp(1)
  
  # Set initial bounds around current draw
  L <- max(a_r, r_curr - 1)
  R <- min(b_r, r_curr + 1)
  
  # Stepping-out / Shrinkage loop
  repeat {
    r_cand <- stats::runif(1, L, R)
    if (log_post(r_cand) > target_y) return(r_cand)
    if (r_cand < r_curr) L <- r_cand else R <- r_cand
  }
}

# --- 3. Gibbs Parameter Sampler ---
sample_params_negbin = function(y, x = NULL, zeta_sam, params, 
                                prior_prec = 0.001, e0 = 0.01, 
                                a_r = 0.001, b_r = 100) {
  n <- length(y)
  X <- if (is.null(x)) as.matrix(zeta_sam) else cbind(as.matrix(x), zeta_sam)
  p <- ncol(X)
  
  # Reconstruct current full coefficient vector
  beta_full_curr <- if (is.null(x)) params$theta else c(params$beta, params$theta)
  eta <- as.vector(X %*% beta_full_curr)
  
  # --- Step 1: Sample Polya-Gamma Latent Auxiliary Variables (omega) ---
  h_vec <- params$r + y
  omega <- BayesLogit::rpg(n, h_vec, eta)
  
  # --- Step 2: Sample Regression Coefficients (beta and theta) ---
  kappa <- (y - params$r) / 2
  Q_B <- crossprod(X, omega * X) + base::diag(prior_prec, p)
  l_B <- crossprod(X, kappa)
  
  ch_Q <- chol(Q_B)
  beta_full <- backsolve(ch_Q, forwardsolve(t(ch_Q), l_B) + stats::rnorm(p))
  
  if (is.null(x)) {
    beta_out <- numeric(0)
    theta_out <- beta_full[1]
  } else {
    beta_out <- beta_full[1:(p - 1)]
    theta_out <- beta_full[p]
  }
  
  # --- Step 3: Sample Size/Dispersion Parameter (r) ---
  eta_new <- as.vector(X %*% beta_full)
  r_new <- slice_sample_r(params$r, y, eta_new, e0 = e0, a_r = a_r, b_r = b_r)
  
  list(
    beta = beta_out,
    theta = theta_out,
    r = r_new
  )
}




