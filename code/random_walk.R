rwmh_move <- function(particles, logliks, weights, X, y, phi_t, cov_proposal, n_steps=1) {
  N <- nrow(particles)
  
  for (step in 1:n_steps) {
    for (i in 1:N) {
      theta_curr <- particles[i, ]
      theta_prop <- mvrnorm(1, mu = theta_curr, Sigma = cov_proposal)
      
      log_target_curr <- sum(dnorm(theta_curr, 0, 1, log = TRUE)) + phi_t * logliks[i]
      loglik_prop <- gp_log_lik(X, y, theta_prop)
      log_target_prop <- sum(dnorm(theta_prop, 0, 1, log = TRUE)) + phi_t * loglik_prop
      
      log_alpha <- log_target_prop - log_target_curr
      
      if (log(runif(1)) < log_alpha) {
        particles[i, ] <- theta_prop
        logliks[i] <- loglik_prop
      }
    }
  }
  return(list(particles = particles, logliks = logliks))
}