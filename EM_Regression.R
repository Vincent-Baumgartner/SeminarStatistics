#-------------------------------------------------------------------------------
# 
#                         Linear Regression by EM
#                       
#          [with general MISSINGNESS PATTERNS in predictor variables]
#         
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#                            Helperfunctions
#-------------------------------------------------------------------------------

new_params <- function(mu, alpha, beta, Sigmax, sigmay2) {
  return(list(mu      = mu,
              alpha   = alpha,
              beta    = beta,
              Sigmax  = Sigmax,
              sigmay2 = sigmay2))
}

max_delta <- function(old, new) {
  return(max(abs(c(new$mu- old$mu,
                   new$alpha - old$alpha,
                   new$beta - old$beta,
                   as.vector(new$Sigmax - old$Sigmax),
                   new$sigmay2 - old$sigmay2))))
}

# initialises starting parameter estimate θ_(0)
em_init <- function(y, X) {
  
  p   <- ncol(X)
  mu0 <- colMeans(X, na.rm = TRUE)
  X_imp <- X
  for (j in 1:p) X_imp[is.na(X_imp[, j]), j] <- mu0[j]
  
  fit <- lm(y ~ X_imp)
  alpha0 <- coef(fit)[1]
  beta0  <- coef(fit)[-1]
  sigmay20 <- mean(residuals(fit)^2)
  Sigmax0 <- cov(X_imp)
  
  new_params(mu0, alpha0, beta0, Sigmax0, sigmay20)
}

#-------------------------------------------------------------------------------
#                       calculational functions
#-------------------------------------------------------------------------------

# posterior per observation (y_i, x_i)
posterior_missing <- function(y_i, x_i, obs_idx, params) {
  
  p <- length(x_i)
  M_idx <- setdiff(1:p, obs_idx)
  
  if (!length(M_idx)) # complete obs.
    return(list(x_hat = x_i, C = NULL, mis_idx = 0))

  # current joint normal (x_i, yi)
  mu_full <- c(params$mu, params$alpha + t(params$beta) %*% params$mu)
  Sigma_full <- rbind(
    cbind(params$Sigmax, 
          params$Sigmax %*% params$beta),
    
    cbind(t(params$beta) %*% params$Sigmax,
          t(params$beta) %*% params$Sigmax %*% params$beta + params$sigmay2)
  )
  
  # Indices for obs. (xio, yi)
  W_idx <- c(obs_idx, p + 1)
  
  # split everything into missing and observed parts 
  Sigma_MW <- Sigma_full[M_idx, W_idx, drop = FALSE]
  Sigma_WW <- Sigma_full[W_idx, W_idx, drop = FALSE]
  Sigma_WW_inv <- solve(Sigma_WW)
  Sigma_MM <- Sigma_full[M_idx, M_idx, drop = FALSE]
  
  w_i <- c(x_i[obs_idx], y_i)
  mu_M <- params$mu[M_idx]
  mu_W <- mu_full[W_idx]
  
  ## conditional estimates
  mu_cond  <- mu_M + Sigma_MW %*% Sigma_WW_inv %*% (w_i - mu_W)
  Sigma_cond <- Sigma_MM - Sigma_MW %*% Sigma_WW_inv %*% t(Sigma_MW)
  
  x_hat <- x_i
  x_hat[M_idx] <- mu_cond
  
  return(list(x_hat = x_hat, C = Sigma_cond, mis_idx = M_idx))
}

# E-Step
e_step <- function(y, X, params) {
  
  n <- length(y); p <- ncol(X)
  Sx  <- numeric(p)
  Sxx <- matrix(0, p, p)
  Sxy <- numeric(p)
  
  for (i in 1:n) {
    
    obs_idx <- which(!is.na(X[i, ]))
    post    <- posterior_missing(y[i], X[i, ], obs_idx, params)
    
    Sx  <- Sx  + post$x_hat
    E_xx <- tcrossprod(post$x_hat)
    
    if (length(post$mis_idx)){
      E_xx[post$mis_idx, post$mis_idx] <-
      E_xx[post$mis_idx, post$mis_idx] + post$C}
    
    Sxx <- Sxx + E_xx
    Sxy <- Sxy + y[i] * post$x_hat
    
  }
  
  return(list(Sx = Sx, Sxx = Sxx, Sxy = Sxy))
}

# M-Step
m_step <- function(stats, const) {
  
  n   <- const$n
  Sx  <- stats$Sx
  Sxx <- stats$Sxx
  Sxy <- stats$Sxy
  Sy  <- const$Sy
  Syy <- const$Syy
  
  mu_new     <- Sx / n
  A_mat      <- Sxx - (1 / n) * tcrossprod(Sx)   
  Sigma_new  <- A_mat / n
  
  beta_new   <- solve(A_mat, Sxy - (Sy / n) * Sx)
  alpha_new  <- (Sy - crossprod(beta_new, Sx)) / n
  
  quad <- as.numeric(t(beta_new) %*% Sxx %*% beta_new)
  sigma_y2_new <- (Syy - 2 * alpha_new * Sy - 2 * crossprod(beta_new, Sxy) +
                     n * alpha_new^2 + 2 * alpha_new * crossprod(beta_new, Sx) +
                     quad) / n
  
  return(new_params(mu_new, alpha_new, beta_new, Sigma_new, sigma_y2_new))
}

#-------------------------------------------------------------------------------
#                             Wrapperfunction
#-------------------------------------------------------------------------------

lm_em <- function(y, X, max_iter = 200, tol = 1e-6, verbose = FALSE) {
  if (!is.numeric(y) || !is.matrix(X) || !is.numeric(X))
    stop("y must be numeric, X must be numeric matrix")
  if (length(y) != nrow(X))
    stop("y and X must have same row length")
  
  n <- length(y)
  const <- list(n = n, Sy = sum(y), Syy = sum(y^2))
  
  params <- em_init(y, X)
  
  for (it in 1:max_iter) {
    stats       <- e_step(y, X, params)
    params_new  <- m_step(stats, const)
    
    if (verbose)
      cat(sprintf("Iter%3d  Delta = %.3e\n",
                  it, max_delta(params, params_new)))
    if (max_delta(params, params_new) < tol) {
      params <- params_new; break
    }
    params <- params_new
  }
  params$iterations <- it
  return(params)
}
