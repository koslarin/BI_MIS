
library(GpGp)
library(matrixStats)


# MAIN FUNCTION

gen_IS_fun = function(y, x = NULL, sample_params, init_params, lik_fun = NULL, 
                      model = c("AIS", "IIS", "plugin", "AIS_copula"), 
                      zeta_ppd, locs = NULL, 
                      R = 300, m_vecchia = 20, mcmc_samples = 1000, burn_in = 1000) {
  
  model <- match.arg(model)
  n <- length(y)
  M <- mcmc_samples
  
  # --- 1. Handle Plug-in vs. Sampling Modes for Zeta ---
  if (model == "plugin") {
    if (is.matrix(zeta_ppd)) {
      #posterior mean as the default if given a matrix of partial posterior draws for plugin
      zeta_sam <- colMeans(zeta_ppd) 
    } else {
      zeta_sam <- as.vector(zeta_ppd)
    }
    if (length(zeta_sam) != n) {
      stop(paste0("For 'plugin' model, zeta_ppd must be length ", n, 
                  " (matching length of y). Got length ", length(zeta_sam)))
    }
  } else {
    if (!is.matrix(zeta_ppd)) {
      stop("For 'AIS' and 'IIS', zeta_ppd must be an (N x n) matrix of posterior draws.")
    }
    if (ncol(zeta_ppd) != n) {
      stop(paste0("zeta_ppd must have ", n, " columns (matching length of y)."))
    }
    S <- nrow(zeta_ppd)
    zeta_sam <- colMeans(zeta_ppd)
  }
  
  # Initialize parameters
  params <- init_params(y, x, zeta_ppd)
  
  # --- 2. AIS Pre-computations (Skipped for IIS & Plug-in) ---
  ais_cache <- list()
  copula_cache <- list()
  if (model == "AIS") {
    ais_cache$mu_part <- colMeans(zeta_ppd)
    Sigma_part <- cov(zeta_ppd)
    Sigma_part_0_diag <- diag(Sigma_part)
    
    Sigma_part_inv <- solve(Sigma_part)
    Sigma_part_0_inv_diag <- 1 / Sigma_part_0_diag
    
    log_det_0 <- sum(log(Sigma_part_0_diag))
    log_det_full <- as.numeric(determinant(Sigma_part, logarithm = TRUE)$modulus)
    ais_cache$log_sqrt_det_ratio <- 0.5 * (log_det_0 - log_det_full)
    ais_cache$Sigma_inv_dif <- Sigma_part_inv - diag(Sigma_part_0_inv_diag)
  } else if (model == "AIS_copula") { # Pre-computations for AIS Copula
    if (is.null(locs)) stop("Spatial coordinates 'locs' required for AIS_copula.")
    
    locs <- as.matrix(locs)
    storage.mode(locs) <- "double"
    
    copula_cache <- setup_copula_cache(
      zeta_ppd = zeta_ppd, 
      locs = locs, 
      m_vecchia = m_vecchia, 
      nugget_prop = 0.10  # Enforces a 10% nugget floor
    )
  }
  
  # --- 3. Storage Allocation ---
  param_names <- names(params)
  storage <- list()
  for (p in param_names) {
    val <- params[[p]]
    if (is.matrix(val)) {
      storage[[p]] <- array(NA_real_, dim = c(M, dim(val)))
    } else if (length(val) == 0 || is.null(val)) {
      storage[[p]] <- NULL
    } else if (length(val) == 1) {
      storage[[p]] <- matrix(NA_real_, nrow = M, ncol = 1)
    } else {
      storage[[p]] <- matrix(NA_real_, nrow = M, ncol = length(val))
    }
  }
  storage[["zeta"]] <- matrix(NA_real_, nrow = M, ncol = n)
  
  # --- 4. Main MCMC Loop ---
  for (j in 1:(burn_in + M)) {
    
    # Step A: Update Stage-2 Parameters
    params <- sample_params(y, x, zeta_sam, params)
    
    # Step B: Update Zeta (Only for IIS and AIS)
    
    if (model == "AIS_copula") {
      log_lik <- lik_fun(y, x, zeta_ppd, params)
      zeta_sam <- sample_zeta_ais_copula_vecchia(
        log_lik = log_lik, 
        zeta_ppd = zeta_ppd, 
        sorted_zeta_ppd = copula_cache$sorted_zeta_ppd, #sorted_zeta_ppd,
        locs = locs, 
        covparms = copula_cache$covparms, 
        NNarray = copula_cache$NNarray,
        R = R, 
        m = m_vecchia
      )
    } else if (model == "IIS") {
      log_lik <- lik_fun(y, x, zeta_ppd, params)
      zeta_sam <- sample_zeta_iis(log_lik, zeta_ppd)
    } else if (model == "AIS") {
      log_lik <- lik_fun(y, x, zeta_ppd, params)
      zeta_sam <- sample_zeta_ais(log_lik, zeta_ppd, 
                                  ais_cache$mu_part, 
                                  ais_cache$Sigma_inv_dif, 
                                  ais_cache$log_sqrt_det_ratio, 
                                  R = R)
    }
    # Note: If model == "plugin", zeta_sam remains fixed across iterations.
    
    # Step C: Save Draws
    if (j > burn_in) {
      idx <- j - burn_in
      for (p in param_names) {
        if (!is.null(storage[[p]])) {
          val <- params[[p]]
          if (is.matrix(val)) {
            storage[[p]][idx, , ] <- val
          } else {
            storage[[p]][idx, ] <- val
          }
        }
      }
      storage[["zeta"]][idx, ] <- zeta_sam
    }
  }
  
  return(storage)
}


# Gaussian: Parameter Initializer
init_params_normal = function(y, x = NULL, zeta_ppd) {
  if (is.null(x)) {
    beta_val <- numeric(0)
  } else {
    beta_val <- rep(0, ncol(as.matrix(x)))
  }
  
  list(
    beta = beta_val,
    theta = 0,
    sigma_y_sq = 1,
    tau = 1
  )
}

# Gaussian: Stage-2 Parameter Gibbs Sampler
sample_params_normal = function(y, x = NULL, zeta_sam, params, prior_prec = 0.001, alpha_0 = 0.01, beta_0 = 0.01) {
  n <- length(y)
  
  # User supplies x (with/without intercept column) or x is NULL
  X <- if (is.null(x)) as.matrix(zeta_sam) else cbind(as.matrix(x), zeta_sam)
  beta_params <- ncol(X)
  
  # Posterior beta & theta update
  Xar <- crossprod(X)
  Q_B <- params$tau * Xar + diag(prior_prec, beta_params)
  l_B <- params$tau * crossprod(X, y)
  ch_Q <- chol(Q_B)
  
  beta_full <- backsolve(ch_Q, forwardsolve(t(ch_Q), l_B) + rnorm(beta_params))
  
  # Residual variance update
  y_hat <- X %*% beta_full
  sse <- sum((y - y_hat)^2)
  tau_new <- rgamma(1, alpha_0 + n/2, beta_0 + 0.5 * sse)
  
  # Separate beta and theta depending on whether x was provided
  if (is.null(x)) {
    beta_out <- numeric(0)
    theta_out <- beta_full[1]
  } else {
    beta_out <- beta_full[1:(beta_params - 1)]
    theta_out <- beta_full[beta_params]
  }
  
  list(
    beta = beta_out,
    theta = theta_out,
    sigma_y_sq = 1 / tau_new,
    tau = tau_new
  )
}

# Gaussian: Log-Likelihood Evaluator
lik_fun_normal = function(y, x = NULL, zeta_ppd, params) {
  S <- nrow(zeta_ppd)
  n <- ncol(zeta_ppd)
  
  if (is.null(x) || length(params$beta) == 0) {
    y_hat_fixed <- 0
  } else {
    y_hat_fixed <- as.vector(as.matrix(x) %*% params$beta)
  }
  
  means_matrix <- sweep(zeta_ppd * params$theta, 2, y_hat_fixed, "+")
  Y_mat <- matrix(y, nrow = S, ncol = n, byrow = TRUE)
  
  dnorm(means_matrix, mean = Y_mat, sd = sqrt(params$sigma_y_sq), log = TRUE)
}

# IIS Sampling Module

sample_zeta_iis = function(log_lik, zeta_ppd) {
  S <- nrow(zeta_ppd)
  n <- ncol(zeta_ppd)
  
  col_maxes <- matrixStats::colMaxs(log_lik)
  prob_zeta <- exp(sweep(log_lik, 2, col_maxes, "-"))
  cum_probs <- matrixStats::colCumsums(prob_zeta)
  
  zeta_indices <- integer(n)
  for (h in 1:n) {
    u_h <- runif(1) * cum_probs[S, h]
    # pmin guards against rare floating-point boundary overflow 
    zeta_indices[h] <- findInterval(u_h, cum_probs[, h]) + 1 
  }
  
  return(zeta_ppd[cbind(zeta_indices, 1:n)])
}


# AIS Sampling Module
sample_zeta_ais = function(log_lik, zeta_ppd, mu_part, Sigma_inv_dif, log_sqrt_det_ratio, R = 300) {
  S <- nrow(zeta_ppd)
  n <- ncol(zeta_ppd)
  
  col_maxes <- matrixStats::colMaxs(log_lik)
  prob_zeta <- exp(sweep(log_lik, 2, col_maxes, "-"))
  cum_probs <- matrixStats::colCumsums(prob_zeta)
  
  zeta_iis_indices <- matrix(NA_integer_, nrow = R, ncol = n)
  for (h in 1:n) {
    u <- runif(R) * cum_probs[S, h]
    zeta_iis_indices[, h] <- findInterval(u, cum_probs[, h]) + 1
  }
  
  zeta_iis_all <- matrix(zeta_ppd[cbind(as.vector(zeta_iis_indices), 
                                        rep(1:n, each = R))], nrow = R, ncol = n)
  
  v <- sweep(zeta_iis_all, 2, mu_part, "-")
  quad_forms <- rowSums((v %*% Sigma_inv_dif) * v)
  log_ratio_AIS_vec <- log_sqrt_det_ratio - 0.5 * quad_forms
  
  log_ratio_mod <- log_ratio_AIS_vec - max(log_ratio_AIS_vec)
  adj_weight_idx <- sample.int(R, 1, prob = exp(log_ratio_mod))
  
  return(zeta_iis_all[adj_weight_idx, ])
}


# Helper function: Transform candidate values to Z-space via ECDFs
transform_to_z = function(zeta_mat, sorted_zeta_ppd) {
  N <- nrow(sorted_zeta_ppd)
  n <- ncol(sorted_zeta_ppd)
  R <- nrow(zeta_mat)
  
  # Standard empirical probability grid points: (i - 0.5) / N
  u_grid <- (seq_len(N) - 0.5) / N
  
  U_mat <- matrix(NA_real_, R, n)
  
  for (h in 1:n) {
    # Continuous piecewise-linear empirical CDF mapping
    ecdf_smooth <- stats::approxfun(
      x = sorted_zeta_ppd[, h], 
      y = u_grid, 
      method = "linear", 
      rule = 2,  # Clamps out-of-bounds candidates safely to grid limits
      ties = "ordered" #don't check for ties, should be better for ecdf # ties = mean is another option
    )
    
    u_raw <- ecdf_smooth(zeta_mat[, h])
    
    # Boundary guard to prevent infinite values in qnorm()
    U_mat[, h] <- pmin(pmax(u_raw, 1 / (2 * N)), 1 - 1 / (2 * N))
  }
  
  # Transform continuous uniform probabilities to standard normal Z-scores
  return(stats::qnorm(U_mat))
}

# Helper AIS_copula: Vecchia AIS Sampler Module
sample_zeta_ais_copula_vecchia = function(log_lik, zeta_ppd, sorted_zeta_ppd, locs, 
                                          covparms, NNarray, R = 300, m = 20) {
  S <- nrow(zeta_ppd)
  n <- ncol(zeta_ppd)
  
  # --- Layer 1: Discrete Proposals ---
  col_maxes <- matrixStats::colMaxs(log_lik)
  prob_zeta <- exp(sweep(log_lik, 2, col_maxes, "-"))
  cum_probs <- matrixStats::colCumsums(prob_zeta)
  
  zeta_iis_indices <- matrix(NA_integer_, nrow = R, ncol = n)
  for (h in 1:n) {
    u <- runif(R) * cum_probs[S, h]
    zeta_iis_indices[, h] <- pmin(findInterval(u, cum_probs[, h]) + 1, S)
  }
  
  zeta_iis_all <- matrix(zeta_ppd[cbind(as.vector(zeta_iis_indices), 
                                        rep(1:n, each = R))], nrow = R, ncol = n)
  
  # --- Layer 2: Copula + Vecchia Re-weighting ---
  Z_cand <- transform_to_z(zeta_iis_all, sorted_zeta_ppd)
  log_ratio_AIS_vec <- numeric(R)
  
  for (r in 1:R) {
    z_r <- Z_cand[r, ]
    
    # Extract $loglik from the returned GpGp list object
    val_vecchia <- GpGp::vecchia_meanzero_loglik(
      covparms = covparms, 
      covfun_name = "exponential_isotropic", 
      y = z_r, 
      locs = locs,
      NNarray = NNarray
    )
    
    # Robust extraction for both list and numeric returns
    log_joint_vecchia <- if (is.list(val_vecchia)) val_vecchia$loglik else as.numeric(val_vecchia)
    
    # Independent N(0,1) log-density denominator
    log_indep_z <- sum(dnorm(z_r, mean = 0, sd = 1, log = TRUE))
    
    log_ratio_AIS_vec[r] <- log_joint_vecchia - log_indep_z
  }
  
  # Sample proportional to exponentiated adjusted weights
  log_ratio_mod <- log_ratio_AIS_vec - max(log_ratio_AIS_vec)
  adj_weight_idx <- sample.int(R, 1, prob = exp(log_ratio_mod))
  
  return(zeta_iis_all[adj_weight_idx, ])
}


## Function to set up copula cache for AIS_copula

setup_copula_cache <- function(zeta_ppd, locs, m_vecchia = 20, nugget_prop = 0.10) {
  n <- ncol(zeta_ppd)
  
  # --- 1. Pre-sort Stage 1 Draws (for Fast Interpolated PIT) ---
  sorted_zeta_ppd <- apply(zeta_ppd, 2, sort)
  
  # --- 2. Pre-compute Vecchia Nearest-Neighbor Structure ---
  NNarray <- GpGp::find_ordered_nn(locs, m = m_vecchia)
  
  # --- 3. Fit Spatial Covariance Parameters on Stage 1 Means ---
  zeta_mean <- colMeans(zeta_ppd)
  fit_gp <- GpGp::fit_model(
    y = zeta_mean, 
    locs = locs, 
    X = matrix(1, nrow = n, ncol = 1),
    covfun_name = "exponential_isotropic",
    silent = TRUE
  )
  
  # Extract range parameter (covparms vector in GpGp is c(variance, range, nugget))
  fitted_range <- fit_gp$covparms[2]
  
  # --- 4. Standardize Variance & Inject Nugget Floor ---
  # In Z-space, total marginal variance must equal 1.0.
  # Partition variance into: Partial Sill (1 - nugget_prop) + Nugget (nugget_prop)
  covparms_adjusted <- c(
    variance = 1.0 - nugget_prop,  # Partial sill (by default, 0.90)
    range    = fitted_range,       # Fitted spatial scale
    nugget   = nugget_prop         # Nugget floor (by default, 0.10)
  )
  
  list(
    sorted_zeta_ppd = sorted_zeta_ppd,
    NNarray = NNarray,
    covparms = covparms_adjusted
  )
}

