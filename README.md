
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
#> 0.6117995 0.6195840 0.2792416 
#> 
#> Upper bound: 90.61294
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
#> 0.5205944 0.5171832 0.1701659 
#> 
#> Upper bound: 90.61294
```

Third, Wasserstein-Bures solver, with bootstrap:

``` r
sensitivity_indices <- ot_indices_wb(x, y, M, boot = TRUE, R = 100)
sensitivity_indices
#> Method: wasserstein-bures 
#> 
#> Indices:
#>        X1        X2        X3 
#> 0.4907486 0.4820302 0.1067025 
#> 
#> Advective component:
#>        X1        X2        X3 
#> 0.2973732 0.3077923 0.1006815 
#> 
#> Diffusive component:
#>          X1          X2          X3 
#> 0.193375338 0.174237863 0.006020998 
#> 
#> Type of confidence interval: norm 
#> Number of replicates: 100 
#> Confidence level: 0.95 
#> Bootstrap statistics:
#>   input component   original        bias     low.ci    high.ci
#> 1    X1        WB 0.50130745 0.010558874 0.47205915 0.50943799
#> 2    X2        WB 0.49158608 0.009555888 0.46464462 0.49941578
#> 3    X3        WB 0.12584905 0.019146555 0.08836225 0.12504274
#> 4    X1 Advective 0.30313340 0.005760165 0.28592711 0.30881936
#> 5    X2 Advective 0.31226320 0.004470863 0.29666440 0.31892026
#> 6    X3 Advective 0.11089625 0.010214758 0.08537492 0.11598808
#> 7    X1 Diffusive 0.19817405 0.004798709 0.18461781 0.20213287
#> 8    X2 Diffusive 0.17932289 0.005085025 0.16640655 0.18206918
#> 9    X3 Diffusive 0.01495279 0.008931797 0.00109468 0.01094731
#> 
#> Upper bound: 90.31479
```

Fourth, we can use the package to compute the sensitivity map on the
output:

``` r
sensitivity_indices <- ot_indices_smap(x, y, M)
sensitivity_indices
#>             X1         X2        X3
#> [1,] 0.5975251 0.04299889 0.1842869
#> [2,] 0.3201639 0.71046059 0.1222535
```
