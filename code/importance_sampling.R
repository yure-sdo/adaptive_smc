importance_sampling_gp <- function(X, y, N = 1000) {
  d <- ncol(X)
  
  particles <- init_particles(N, d)
  
  start <- Sys.time()
  
  logliks <- apply(particles, 1, function(theta) gp_log_lik(X, y, theta))
  
  max_loglik <- max(logliks)
  
  logZ <- max_loglik + log(mean(exp(logliks - max_loglik)))
  
  time <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  
  weights <- exp(logliks - max_loglik)
  weights <- weights / sum(weights)
  
  ess <- 1 / sum(weights^2)
  
  return(list(logZ = logZ, ess = ess, time = time, ess_per_sec = ess / time))
}