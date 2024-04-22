
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gsaot

<!-- badges: start -->

[![R-CMD-check](https://github.com/pietrocipolla/gsaot/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pietrocipolla/gsaot/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The package `gsaot` provides a set of tools to compute and plot Optimal
Transport (OT) based sensitivity indices. The core functions of the
package are: \* `ot_indices()`: compute OT indices for multivariate
outputs using different solvers for OT (network simplex, Sinkhorn, and
so on). \* `ot_indices_wb()`: compute OT indices for univariate or
multivariate outputs using the Wasserstein-Bures semi-metric. \*
`ot_indices_1d()`: compute OT indices for univariate outputs using OT
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

This is a basic example which shows you how to solve a common problem:

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
#> 0.6104298 0.6432362 0.3063337 
#> 
#> Upper bound: 96.54529
```
