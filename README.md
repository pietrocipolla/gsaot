
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gsaot

<!-- badges: start -->

[![R-CMD-check](https://github.com/pietrocipolla/gsaot/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pietrocipolla/gsaot/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The package `gsaot` provides a set of tools to compute and plot Optimal
Transport (OT) based sensitivity indices. The core functions of the
package are: - `ot_indices()`: compute OT indices for multivariate
outputs using different solvers for OT (network simplex, Sinkhorn, and
so on). - `ot_indices_wb()`: compute OT indices for univariate or
multivariate outputs using the Wasserstein-Bures semi-metric. -
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
#> 0.5904896 0.6259631 0.2804490 
#> 
#> Upper bound: 95.31012

# Compute the sensitivity indices using the Network Simplex solver and default parameters
sensitivity_indices <- ot_indices(x, y, M, solver = "wasserstein")
#> Using default values for solver wasserstein
sensitivity_indices
#> Method: wasserstein 
#> 
#> Indices:
#>        X1        X2        X3 
#> 0.5068366 0.5323364 0.1828405 
#> 
#> Upper bound: 95.31012

# Compute the Wasserstein-Bures indices
sensitivity_indices <- ot_indices_wb(x, y, M)
sensitivity_indices
#> Method: wasserstein-bures 
#> 
#> Indices:
#>        X1        X2        X3 
#> 0.4873015 0.5075676 0.1395621 
#> 
#> Advective component:
#>        X1        X2        X3 
#> 0.2967554 0.3209685 0.1213170 
#> 
#> Diffusive component:
#>         X1         X2         X3 
#> 0.19054611 0.18659905 0.01824503 
#> 
#> Upper bound: 95.31012

# Compute the sensitivity map using 1-dimensional solver
sensitivity_indices <- ot_indices_smap(x, y, M)
sensitivity_indices
#>             X1         X2        X3
#> [1,] 0.5891068 0.04332126 0.1881632
#> [2,] 0.3145604 0.72062949 0.1395558
```
