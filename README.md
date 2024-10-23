
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gsaot <a href="https://pietrocipolla.github.io/gsaot/"><img src="man/figures/logo.png" align="right" height="139" alt="gsaot website" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/pietrocipolla/gsaot/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pietrocipolla/gsaot/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/pietrocipolla/gsaot/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/pietrocipolla/gsaot/actions/workflows/test-coverage.yaml)
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
#> Using default values for solver sinkhorn
sensitivity_indices
#> Method: sinkhorn 
#> 
#> Indices:
#>        X1        X2        X3 
#> 0.5968109 0.6523353 0.2959787 
#> 
#> Upper bound: 97.05654
```

Second, Network Simplex solver:

``` r
sensitivity_indices <- ot_indices(x, y, M, solver = "transport")
#> Using default values for solver transport
sensitivity_indices
#> Method: transport 
#> 
#> Indices:
#>        X1        X2        X3 
#> 0.5029419 0.5455398 0.1817795 
#> 
#> Upper bound: 97.05654
```

Third, Wasserstein-Bures solver, with bootstrap:

``` r
sensitivity_indices <- ot_indices_wb(x, y, M, boot = TRUE, R = 100)
sensitivity_indices
#> Method: wasserstein-bures 
#> 
#> Indices:
#>        X1        X2        X3 
#> 0.4715790 0.5106960 0.1217466 
#> 
#> Advective component:
#>        X1        X2        X3 
#> 0.2898702 0.3241905 0.1129709 
#> 
#> Diffusive component:
#>          X1          X2          X3 
#> 0.181708726 0.186505506 0.008775646 
#> 
#> Type of confidence interval: norm 
#> Number of replicates: 100 
#> Confidence level: 0.95 
#> Indices confidence intervals:
#>   Inputs     Index      low.ci    high.ci
#> 1     X1        WB 0.451141101 0.49201684
#> 2     X2        WB 0.496408255 0.52498375
#> 3     X3        WB 0.100694207 0.14279898
#> 4     X1 Advective 0.276749928 0.30299056
#> 5     X2 Advective 0.314497738 0.33388326
#> 6     X3 Advective 0.094646496 0.13129540
#> 7     X1 Diffusive 0.171963032 0.19145442
#> 8     X2 Diffusive 0.180299983 0.19271103
#> 9     X3 Diffusive 0.004281849 0.01326944
#> 
#> Upper bound: 97.17103
```

Fourth, we can use the package to compute the sensitivity map on the
output:

``` r
sensitivity_indices <- ot_indices_smap(x, y, M)
sensitivity_indices
#>             X1         X2        X3
#> [1,] 0.5749624 0.04188444 0.1852877
#> [2,] 0.3261860 0.72632849 0.1404875
```
