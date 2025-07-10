generate_gp_data <- function(n, d, signal_sd = 1, noise_sd = 1) {
  X <- matrix(rnorm(n * d), nrow = n, ncol = d)
  
  true_lengthscales = runif(d, 100, 200)
  X_scaled <- sweep(X, 2, true_lengthscales, "/")
  
  dists <- as.matrix(dist(X_scaled))^2
  
  K <- signal_sd^2 * exp(-0.5 * dists)
  
  Ky <- K + noise_sd^2 * diag(n)
  
  y <- as.numeric(mvrnorm(1, mu = rep(0, n), Sigma = Ky))
  
  return(list(X = X, y = y, 
              true_params = list(log_lengthscales = log(true_lengthscales),
                                 log_signal_sd = log(signal_sd),
                                 log_noise_sd = log(noise_sd))))
}