# Adaptive Sequential Monte Carlo for High-Dimensional Bayesian Inference

This repository contains the code for the assignment of Computational Statistics, a PhD-level course offered at the School of Applied Mathematics, Getulio Vargas Foundation.

## About the project

This project explores the use of Adaptive Sequential Monte Carlo (SMC) methods for robust posterior inference and model selection in high-dimensional Bayesian models. Standard sampling techniques like Importance Sampling (IS) and even Annealed Importance Sampling (AIS) often struggle in these settings, leading to unreliable estimates.

This work implements and compares three methods on a high-dimensional Gaussian Process regression problem:

- Standard Importance Sampling (IS)

- Annealed Importance Sampling (AIS)

- Adaptive Sequential Monte Carlo (SMC) with adaptive resampling and proposal kernels.
