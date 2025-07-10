annealed_importance_sampling_gp <- function(X, y, N = 500, T_steps = 20, threshold_ess = 0.5) {
  d <- ncol(X)
  particles <- init_particles(N, d)
  
  weights <- rep(1/N, N)
  phi_seq <- seq(0, 1, length.out = T_steps)
  
  logZ <- 0
  logliks <- apply(particles, 1, function(theta) gp_log_lik(X, y, theta))
  
  start <- Sys.time()
  
  for (t in 2:T_steps) {
    # random walk proposal
    cov_prop <- covariance_matrix(particles, weights)
    res <- rwmh_move(particles, logliks, weights, X, y, phi_seq[t], cov_prop, n_steps = 1)
    particles <- res$particles
    logliks <- res$logliks
    
    # update weights
    inc_weights <- exp((phi_seq[t] - phi_seq[t-1]) * logliks)
    weights <- weights * inc_weights
    weights <- weights / sum(weights)
    
    logZ <- logZ + log(mean(inc_weights))
  }
  
  time <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  
  ess <- 1 / sum(weights^2)
  
  return(list(logZ = logZ, ess = ess, time = time, ess_per_sec = ess / time))
}