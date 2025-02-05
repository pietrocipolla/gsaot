
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gsaot <a href="https://pietrocipolla.github.io/gsaot/"><img src="man/figures/logo.png" align="right" height="139" alt="gsaot website" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/pietrocipolla/gsaot/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pietrocipolla/gsaot/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/pietrocipolla/gsaot/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/pietrocipolla/gsaot/actions/workflows/test-coverage.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/gsaot)](https://CRAN.R-project.org/package=gsaot)
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

``` r
install.packages("gsaot")
```

You can install the development version of gsaot from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("pietrocipolla/gsaot")
```

### :exclamation: :exclamation: Installation note

The `sinkhorn` and `sinkhorn_log` solvers in `gsaot` greatly benefit
from optimization in compilation. To add this option (before package
installation), edit your `.R/Makevars` file with the desired flags. Even
though different compilers have different options, a common flag to
enable a safe level of optimization is

    CXXFLAGS+=-O2

More detailed information on how to customize the R packages compilation
can be found in the [R
guide](https://cran.r-project.org/doc/manuals/R-admin.html#Customizing-package-compilation).

## Example

We can use a gaussian toy model with three outputs as an example:

``` r
library(gsaot)

N <- 1000

mx <- c(1, 1, 1)
Sigmax <- matrix(data = c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), nrow = 3)

x1 <- rnorm(N)
x2 <- rnorm(N)
x3 <- rnorm(N)

x <- cbind(x1, x2, x3)
x <- mx + x %*% chol(Sigmax)

A <- matrix(data = c(4, -2, 1, 2, 5, -1), nrow = 2, byrow = TRUE)
y <- t(A %*% t(x))

x <- data.frame(x)
```

After having defined the number of partitions, we compute the
sensitivity indices using different solvers. First, Sinkhorn solver and
default parameters:

``` r
M <- 25

sensitivity_indices <- ot_indices(x, y, M)
sensitivity_indices
#> Method: sinkhorn 
#> 
#> Solver Options:
#> $numIterations
#> [1] 1000
#> 
#> $epsilon
#> [1] 0.01
#> 
#> $maxErr
#> [1] 1e-09
#> 
#> 
#> Indices:
#>        X1        X2        X3 
#> 0.5921956 0.6361059 0.2896426 
#> 
#> Upper bound: 93.68423
```

Second, Network Simplex solver:

``` r
sensitivity_indices <- ot_indices(x, y, M, solver = "transport")
sensitivity_indices
#> Method: transport 
#> 
#> Solver Options:
#> $fullreturn
#> [1] TRUE
#> 
#> 
#> Indices:
#>        X1        X2        X3 
#> 0.5070913 0.5402425 0.1888146 
#> 
#> Upper bound: 93.68423
```

Third, Wasserstein-Bures solver, with bootstrap:

``` r
sensitivity_indices <- ot_indices_wb(x, y, M, boot = TRUE, R = 100)
sensitivity_indices
#> Method: wasserstein-bures 
#> 
#> Indices:
#>        X1        X2        X3 
#> 0.4762357 0.5060195 0.1260812 
#> 
#> Advective component:
#>        X1        X2        X3 
#> 0.2886873 0.3198006 0.1149821 
#> 
#> Diffusive component:
#>         X1         X2         X3 
#> 0.18754839 0.18621891 0.01109912 
#> 
#> Type of confidence interval: norm 
#> Number of replicates: 100 
#> Confidence level: 0.95 
#> Indices confidence intervals:
#>   Inputs     Index      low.ci    high.ci
#> 1     X1        WB 0.458354687 0.49411672
#> 2     X2        WB 0.489574375 0.52246464
#> 3     X3        WB 0.104332704 0.14782972
#> 4     X1 Advective 0.276044671 0.30132996
#> 5     X2 Advective 0.309620215 0.32998098
#> 6     X3 Advective 0.095967873 0.13399630
#> 7     X1 Diffusive 0.179785566 0.19531121
#> 8     X2 Diffusive 0.178874908 0.19356291
#> 9     X3 Diffusive 0.006317555 0.01588069
#> 
#> Upper bound: 93.6336
```

Fourth, we can use the package to compute the sensitivity map on the
output:

``` r
sensitivity_indices <- ot_indices_smap(x, y, M)
sensitivity_indices
#>             X1         X2        X3
#> [1,] 0.5769613 0.05140922 0.2142022
#> [2,] 0.3119494 0.73023174 0.1353518
```
