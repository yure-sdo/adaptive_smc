sample_prior <- function(d) {
  log_lengthscales <- rnorm(d)
  log_sigmaf <- rnorm(1, 0, 1)
  log_sigma  <- rnorm(1, 0, 1)
  return(c(log_lengthscales, log_sigmaf, log_sigma))
}

init_particles <- function(N, d) {
  return(t(replicate(N, sample_prior(d))))
}

covariance_matrix <- function(particles, weights) {
  mu <- colSums(particles * weights)
  diffs <- sweep(particles, 2, mu)
  return(t(diffs) %*% (diffs * weights))
}

# GP log marginal likelihood
gp_log_lik <- function(X, y, theta) {
  n <- nrow(X)
  d <- ncol(X)
  
  log_lengthscales <- theta[1:d]
  log_signal_sd <- theta[d+1]
  log_noise_sd <- theta[d+2]
  
  lengthscales <- exp(log_lengthscales)
  signal_var <- exp(2 * log_signal_sd)
  noise_var <- exp(2 * log_noise_sd)
  
  X_scaled <- sweep(X, 2, lengthscales, "/")
  dists <- as.matrix(dist(X_scaled))^2
  K <- signal_var * exp(-0.5 * dists)
  Ky <- K + noise_var * diag(n)
  
  L <- chol(Ky)
  alpha <- backsolve(t(L), forwardsolve(L, y))
  
  log_det <- 2 * sum(log(diag(L)))
  
  lml <- -0.5 * sum(y * alpha) - 0.5 * log_det - 0.5 * n * log(2 * pi)
  
  return(lml)
}
