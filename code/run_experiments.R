library(MASS)
library(Matrix)
library(dplyr)
library(tidyr)
library(purrr)
library(knitr)

source("data_generation.R")
source("random_walk.R")
source("adaptive_smc.R")
source("importance_sampling.R")
source("annealed_is.R")
source("utils.R")

sim <- generate_gp_data(n = 300, d = 50, signal_sd = 1, noise_sd = 1)
X <- sim$X
y <- sim$y

n_runs <- 100

results <- list(SMC = list(), IS = list(), AIS = list())

for (i in 1:n_runs) {
  smc  <- adaptive_smc_gp(X, y, N = 500, T_steps = 5)
  is   <- importance_sampling_gp(X, y, N = 500)
  ais  <- annealed_importance_sampling_gp(X, y, N = 500, T_steps = 5)
  
  results$SMC[[i]] <- smc
  results$IS[[i]]  <- is
  results$AIS[[i]] <- ais
  
  cat(sprintf("Completed %d runs\n", i))
}

df_results <- map_dfr(results, bind_rows, .id = "Method")

summary_stats <- df_results %>%
  group_by(Method) %>%
  summarise(
    across(
      .cols = c(logZ, ess, time, ess_per_sec),
      .fns = list(
        mean = ~mean(.x, na.rm = TRUE),
        q05 = ~quantile(.x, 0.05, na.rm = TRUE),
        q95 = ~quantile(.x, 0.95, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>%
  mutate(Method = factor(Method, levels = c("IS", "AIS", "SMC"))) %>%
  arrange(Method)


table1_data <- summary_stats %>%
  transmute(
    Method,
    `Mean(logZ)` = sprintf("%.2f", logZ_mean),
    `90% Sample Interval` = sprintf("[%.2f, %.2f]", logZ_q05, logZ_q95)
  )

print(kable(table1_data, format = "pipe", align = 'lrr'))

table2_data <- summary_stats %>%
  transmute(
    Method,
    `Time (s)` = sprintf("%.2f [%.2f, %.2f]", time_mean, time_q05, time_q95),
    `ESS` = sprintf("%.1f [%.1f, %.1f]", ess_mean, ess_q05, ess_q95),
    `ESS / sec` = sprintf("%.2f [%.2f, %.2f]", ess_per_sec_mean, ess_per_sec_q05, ess_per_sec_q95)
  )

for (i in 1:nrow(table2_data)) {
  cat(paste(table2_data[i,], collapse = " | "), "|\n")
}