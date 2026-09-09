# gsaot

The package `gsaot` provides a set of tools to compute and plot Optimal
Transport (OT) based sensitivity indices. The core functions of the
package are:

- [`ot_indices()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices.md):
  compute OT indices for multivariate outputs using different solvers
  for OT (network simplex, Sinkhorn, and so on).

- [`ot_indices_wb()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices_wb.md):
  compute OT indices for univariate or multivariate outputs using the
  Wasserstein-Bures semi-metric.

- [`ot_indices_1d()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices_1d.md):
  compute OT indices for univariate outputs using OT solution in one
  dimension.

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

### ❗ ❗ Installation note

The `sinkhorn` and `sinkhorn_stable` solvers in `gsaot` greatly benefit
from optimization in compilation. To add this option (before package
installation), edit your `.R/Makevars` file with the desired flags. Even
though different compilers have different options, a common flag to
enable a safe level of optimization is

``` R
CXXFLAGS+=-O2
```

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
#> 0.7230669 0.7725176 0.5109855
```

Second, Network Simplex solver:

``` r

sensitivity_indices <- ot_indices(x, y, M, solver = "transport")
sensitivity_indices
#> Method: transport 
#> 
#> Indices:
#>        X1        X2        X3 
#> 0.4932602 0.5307910 0.1787914
```

Third, Wasserstein-Bures solver, with bootstrap:

``` r

sensitivity_indices <- ot_indices_wb(x, y, M, boot = TRUE, R = 100)
sensitivity_indices
#> Method: wass-bures 
#> 
#> Indices:
#>        X1        X2        X3 
#> 0.4586257 0.4905809 0.1132465 
#> 
#> Advective component:
#>        X1        X2        X3 
#> 0.2810481 0.3110967 0.1048114 
#> 
#> Diffusive component:
#>          X1          X2          X3 
#> 0.177577600 0.179484207 0.008435054 
#> 
#> Type of confidence interval: norm 
#> Number of replicates: 100 
#> Confidence level: 0.95 
#> Bootstrap statistics:
#>   input  component   original        bias      low.ci    high.ci
#> 1    X1 wass-bures 0.47094416 0.012318475 0.439759468 0.47749190
#> 2    X2 wass-bures 0.50234531 0.011764409 0.473281120 0.50788069
#> 3    X3 wass-bures 0.13081694 0.017570472 0.095125939 0.13136699
#> 4    X1  advective 0.28762307 0.006574989 0.268271901 0.29382427
#> 5    X2  advective 0.31613211 0.005035409 0.300270370 0.32192303
#> 6    X3  advective 0.11369616 0.008884746 0.089916886 0.11970594
#> 7    X1  diffusive 0.18332109 0.005743487 0.169803189 0.18535201
#> 8    X2  diffusive 0.18621321 0.006729000 0.171303894 0.18766452
#> 9    X3  diffusive 0.01712078 0.008685726 0.003249537 0.01362057
```

Fourth, we can use the package to compute the sensitivity map on the
output:

``` r

sensitivity_indices <- ot_indices_smap(x, y, M)
sensitivity_indices
#>             X1         X2        X3
#> [1,] 0.5814812 0.04629205 0.1866779
#> [2,] 0.2987823 0.71591879 0.1291407
```
