
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gsaot

<!-- badges: start -->

[![R-CMD-check](https://github.com/pietrocipolla/gsaot/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pietrocipolla/gsaot/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The package `gsaot` provides a set of tools to compute and plot Optimal
Transport (OT) based sensitivity indices. The core functions of the
package are:

- `ot_indices()`: compute OT indices for multivariate outputs using
  different solvers for OT (network simplex, Sinkhorn, and so on).

- `ot_indices_wb()`: compute OT indices for univariate or multivariate
  outputs using the Wasserstein-Bures semi-metric.

- `ot_indices_1d()`: compute OT indices for univariate outputs using OT
  solution in one dimension. The package provides also functions to plot
  the resulting indices and the inner statistics.

## Installation

You can install the development version of gsaot from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("pietrocipolla/gsaot")
```

## Example

A basic application of the functions implemented in the package:

``` r
library(gsaot)

N <- 1000

# Define the inputs distribution
mx <- c(1, 1, 1)
Sigmax <- matrix(data = c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), nrow = 3)

# Sample the inputs
x1 <- rnorm(N)
x2 <- rnorm(N)
x3 <- rnorm(N)

x <- cbind(x1, x2, x3)
x <- mx + x %*% chol(Sigmax)

# Define the (linear)  model and the output
A <- matrix(data = c(4, -2, 1, 2, 5, -1), nrow = 2, byrow = TRUE)
y <- t(A %*% t(x))

x <- data.frame(x)

M <- 25

# Compute the sensitivity indices using Sinkhorn's solver and default parameters
sensitivity_indices <- ot_indices(x, y, M)
#> Using default values for solver sinkhorn
sensitivity_indices
#> Method: sinkhorn 
#> 
#> Indices:
#>        X1        X2        X3 
#> 0.5642179 0.6255448 0.2736724 
#> 
#> Upper bound: 89.13198

# Compute the sensitivity indices using the Network Simplex solver and default parameters
sensitivity_indices <- ot_indices(x, y, M, solver = "transport")
#> Using default values for solver transport
sensitivity_indices
#> Method: transport 
#> 
#> Indices:
#>        X1        X2        X3 
#> 0.4771722 0.5267831 0.1709329 
#> 
#> Upper bound: 89.13198

# Compute the Wasserstein-Bures indices
sensitivity_indices <- ot_indices_wb(x, y, M)
sensitivity_indices
#> Method: wasserstein-bures 
#> 
#> Indices:
#>        X1        X2        X3 
#> 0.4561290 0.5017250 0.1238633 
#> 
#> Advective component:
#>        X1        X2        X3 
#> 0.2777030 0.3160325 0.1072954 
#> 
#> Diffusive component:
#>        X1        X2        X3 
#> 0.1784259 0.1856925 0.0165679 
#> 
#> Upper bound: 89.13198

# Compute the sensitivity map using 1-dimensional solver
sensitivity_indices <- ot_indices_smap(x, y, M)
sensitivity_indices
#>             X1         X2        X3
#> [1,] 0.5710719 0.03378987 0.1580119
#> [2,] 0.2833240 0.70813481 0.1332598
```
