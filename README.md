
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
#> 0.5774381 0.6170742 0.2532476 
#> 
#> Upper bound: 97.84373
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
#> 0.5003723 0.5318571 0.1652095 
#> 
#> Upper bound: 97.84373
```

Third, Wasserstein-Bures solver, with bootstrap:

``` r
sensitivity_indices <- ot_indices_wb(x, y, M, boot = TRUE, R = 100)
sensitivity_indices
#> Method: wasserstein-bures 
#> 
#> Indices:
#>        X1        X2        X3 
#> 0.4699628 0.4978146 0.1038411 
#> 
#> Advective component:
#>         X1         X2         X3 
#> 0.28890167 0.31629505 0.09761414 
#> 
#> Diffusive component:
#>          X1          X2          X3 
#> 0.181061103 0.181519501 0.006226975 
#> 
#> Type of confidence interval: norm 
#> Number of replicates: 100 
#> Confidence level: 0.95 
#> Bootstrap statistics:
#>   input component   original        bias      low.ci    high.ci
#> 1    X1        WB 0.48149959 0.011536815 0.452453163 0.48747238
#> 2    X2        WB 0.50688549 0.009070932 0.479675009 0.51595410
#> 3    X3        WB 0.12256998 0.018728865 0.085576022 0.12210620
#> 4    X1 Advective 0.29474371 0.005842038 0.277287964 0.30051537
#> 5    X2 Advective 0.32033923 0.004044179 0.305323253 0.32726685
#> 6    X3 Advective 0.10785107 0.010236933 0.082351205 0.11287707
#> 7    X1 Diffusive 0.18675588 0.005694777 0.172723336 0.18939887
#> 8    X2 Diffusive 0.18654625 0.005026753 0.173173777 0.18986522
#> 9    X3 Diffusive 0.01471891 0.008491932 0.001555473 0.01089848
#> 
#> Upper bound: 97.71638
```

Fourth, we can use the package to compute the sensitivity map on the
output:

``` r
sensitivity_indices <- ot_indices_smap(x, y, M)
sensitivity_indices
#>             X1         X2        X3
#> [1,] 0.5800530 0.04109744 0.1571959
#> [2,] 0.3145364 0.71390309 0.1282103
```
