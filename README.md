
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
  solution in one dimension.

The package also provides functions to plot the resulting indices and
the separation measures.

## Installation

``` r
install.packages("gsaot")
```

You can install the development version of gsaot from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("pietrocipolla/gsaot")
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
#> 0.5856179 0.6321027 0.2833714 
#> 
#> Upper bound: 93.27558
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
#> 0.4867229 0.5197051 0.1618929 
#> 
#> Upper bound: 93.27558
```

Third, Wasserstein-Bures solver, with bootstrap:

``` r
sensitivity_indices <- ot_indices_wb(x, y, M, boot = TRUE, R = 100)
sensitivity_indices
#> Method: wasserstein-bures 
#> 
#> Indices:
#>         X1         X2         X3 
#> 0.45682255 0.48218997 0.09602267 
#> 
#> Advective component:
#>         X1         X2         X3 
#> 0.27698647 0.30531330 0.09049933 
#> 
#> Diffusive component:
#>          X1          X2          X3 
#> 0.179836078 0.176876664 0.005523346 
#> 
#> Type of confidence interval: norm 
#> Number of replicates: 100 
#> Confidence level: 0.95 
#> Bootstrap statistics:
#>   input component   original        bias       low.ci    high.ci
#> 1    X1        WB 0.46632977 0.009507222 0.4386714616 0.47497364
#> 2    X2        WB 0.49360786 0.011417898 0.4647262667 0.49965367
#> 3    X3        WB 0.11565245 0.019629776 0.0772420360 0.11480331
#> 4    X1 Advective 0.28174446 0.004757991 0.2651150430 0.28885790
#> 5    X2 Advective 0.31059858 0.005285280 0.2941242206 0.31650239
#> 6    X3 Advective 0.10105252 0.010553192 0.0744003614 0.10659829
#> 7    X1 Diffusive 0.18458531 0.004749231 0.1717065389 0.18796562
#> 8    X2 Diffusive 0.18300928 0.006132618 0.1690000331 0.18475330
#> 9    X3 Diffusive 0.01459993 0.009076584 0.0009580535 0.01008864
#> 
#> Upper bound: 93.16529
```

Fourth, we can use the package to compute the sensitivity map on the
output:

``` r
sensitivity_indices <- ot_indices_smap(x, y, M)
sensitivity_indices
#>             X1        X2        X3
#> [1,] 0.5897603 0.0426077 0.1493735
#> [2,] 0.2822771 0.7039011 0.1233232
```
