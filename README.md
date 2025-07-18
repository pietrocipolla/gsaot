
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
#> Indices:
#>        X1        X2        X3 
#> 0.5925141 0.6584782 0.2949564
```

Second, Network Simplex solver:

``` r
sensitivity_indices <- ot_indices(x, y, M, solver = "transport")
sensitivity_indices
#> Method: transport 
#> 
#> Indices:
#>        X1        X2        X3 
#> 0.4999788 0.5534392 0.1829258
```

Third, Wasserstein-Bures solver, with bootstrap:

``` r
sensitivity_indices <- ot_indices_wb(x, y, M, boot = TRUE, R = 100)
sensitivity_indices
#> Method: wass-bures 
#> 
#> Indices:
#>        X1        X2        X3 
#> 0.4697990 0.5198658 0.1207588 
#> 
#> Advective component:
#>        X1        X2        X3 
#> 0.2902096 0.3296279 0.1091971 
#> 
#> Diffusive component:
#>        X1        X2        X3 
#> 0.1795894 0.1902379 0.0115617 
#> 
#> Type of confidence interval: norm 
#> Number of replicates: 100 
#> Confidence level: 0.95 
#> Bootstrap statistics:
#>   input  component   original        bias      low.ci    high.ci
#> 1    X1 wass-bures 0.47978725 0.009988218 0.450824897 0.48877316
#> 2    X2 wass-bures 0.52917218 0.009306359 0.504610876 0.53512076
#> 3    X3 wass-bures 0.13698299 0.016224178 0.100790712 0.14072692
#> 4    X1  advective 0.29535740 0.005147760 0.276749775 0.30366951
#> 5    X2  advective 0.33352122 0.003893343 0.320432301 0.33882346
#> 6    X3  advective 0.11759222 0.008395100 0.092129728 0.12626451
#> 7    X1  diffusive 0.18442984 0.004840458 0.172016597 0.18716217
#> 8    X2  diffusive 0.19565095 0.005413016 0.182888068 0.19758781
#> 9    X3  diffusive 0.01939077 0.007829078 0.006402971 0.01672042
```

Fourth, we can use the package to compute the sensitivity map on the
output:

``` r
sensitivity_indices <- ot_indices_smap(x, y, M)
sensitivity_indices
#>             X1         X2        X3
#> [1,] 0.5726212 0.04864452 0.1703333
#> [2,] 0.3249486 0.72732618 0.1464068
```
