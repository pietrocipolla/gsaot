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

The `sinkhorn` and `sinkhorn_log` solvers in `gsaot` greatly benefit
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
#> 0.6061460 0.6337609 0.2922763
```

Second, Network Simplex solver:

``` r
sensitivity_indices <- ot_indices(x, y, M, solver = "transport")
sensitivity_indices
#> Method: transport 
#> 
#> Indices:
#>        X1        X2        X3 
#> 0.5049639 0.5188542 0.1682235
```

Third, Wasserstein-Bures solver, with bootstrap:

``` r
sensitivity_indices <- ot_indices_wb(x, y, M, boot = TRUE, R = 100)
sensitivity_indices
#> Method: wass-bures 
#> 
#> Indices:
#>        X1        X2        X3 
#> 0.4751984 0.4826512 0.1044734 
#> 
#> Advective component:
#>         X1         X2         X3 
#> 0.28905226 0.30849658 0.09739885 
#> 
#> Diffusive component:
#>         X1         X2         X3 
#> 0.18614614 0.17415462 0.00707455 
#> 
#> Type of confidence interval: norm 
#> Number of replicates: 100 
#> Confidence level: 0.95 
#> Bootstrap statistics:
#>   input  component   original        bias      low.ci    high.ci
#> 1    X1 wass-bures 0.48387515 0.008676749 0.458657687 0.49173911
#> 2    X2 wass-bures 0.49374040 0.011089198 0.466735426 0.49856697
#> 3    X3 wass-bures 0.12212017 0.017646775 0.083725558 0.12522124
#> 4    X1  advective 0.29324822 0.004195965 0.278010956 0.30009356
#> 5    X2  advective 0.31343443 0.004937852 0.298375896 0.31861727
#> 6    X3  advective 0.10624660 0.008847758 0.080308465 0.11448923
#> 7    X1  diffusive 0.19062692 0.004480784 0.178795145 0.19349713
#> 8    X2  diffusive 0.18030596 0.006151346 0.166726877 0.18158236
#> 9    X3  diffusive 0.01587357 0.008799018 0.001765732 0.01238337
```

Fourth, we can use the package to compute the sensitivity map on the
output:

``` r
sensitivity_indices <- ot_indices_smap(x, y, M)
sensitivity_indices
#>             X1         X2        X3
#> [1,] 0.5898173 0.03863799 0.1744461
#> [2,] 0.3094364 0.70845462 0.1216729
```
