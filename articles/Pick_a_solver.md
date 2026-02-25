# Pick an OT solver

When first navigating the world of Optimal Transport, a reasonable
question may be: how should I solve this optimization problem? This
vignette provides a practical guide to selecting the right solver and
function in `gsaot`, based on output dimensionality, data size, and the
kind of information you want from the indices. The guide is intended for
users who are new to OT-based sensitivity analysis and want a clear
starting point for their analyses.

``` r
library(gsaot)
```

## A quick decision guide

Use the following three dimensions to guide your choice:

1.  **Output dimensionality**
    - If the output is one-dimensional, use
      [`ot_indices_1d()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices_1d.md)
      or
      [`ot_indices_wb()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices_wb.md).
    - If the output is multi-dimensional, use
      [`ot_indices()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices.md)
      or
      [`ot_indices_wb()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices_wb.md).
2.  **Data size and bootstrapping**
    - For small to medium data sets (roughly up to a few thousand
      points),
      [`ot_indices()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices.md)
      with `solver = "transport"` is a safe and accurate choice.
    - For large data sets, or when bootstrapping is needed, prefer
      [`ot_indices()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices.md)
      with `solver = "sinkhorn"` to reduce computational cost.
3.  **Type of information desired**
    - If you want the Wasserstein-Bures decomposition (advective and
      diffusive components), use
      [`ot_indices_wb()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices_wb.md).

The next sections expand on each dimension and provide minimal examples.

## Dimension 1: output is 1D vs multi-D

If the output is one-dimensional,
[`ot_indices_1d()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices_1d.md)
is specialized and faster. It avoids building large cost matrices and is
the recommended default for scalar outputs.

``` r
set.seed(1)
N <- 800

x <- data.frame(
  x1 = rnorm(N),
  x2 = rnorm(N)
)
y1 <- 2 * x$x1 - x$x2 + rnorm(N)

res_1d <- ot_indices_1d(x, y1, M = 20)
res_1d
#> Method: 1-dimensional (p=2) 
#> 
#> Indices:
#>        x1        x2 
#> 0.4230686 0.1135296
```

For multi-dimensional outputs, use
[`ot_indices()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices.md)
for the standard OT indices or
[`ot_indices_wb()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices_wb.md)
if you need the Wasserstein-Bures information.

``` r
Y <- cbind(y1, 0.5 * x$x1 + rnorm(N))

res_md <- ot_indices(x, Y, M = 20)
res_md
#> Method: sinkhorn 
#> 
#> Indices:
#>        x1        x2 
#> 0.4992386 0.2322164
```

## Dimension 2: solver choice for `ot_indices()`

The function
[`ot_indices()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices.md)
supports three solvers:

- `solver = "transport"`: solves the non-regularized OT problem via the
  `transport` package. This is a solid default for small to medium data
  sets when you want a close match to the exact OT solution.
- `solver = "sinkhorn"`: solves the entropic-regularized OT problem.
  This is faster and more scalable, especially when you want
  bootstrapping.
- `solver = "sinkhorn_log"`: numerically more stable for very small
  regularization, at a higher cost.

### Small to medium data: `solver = "transport"`

``` r
res_transport <- ot_indices(
  x, Y, M = 20,
  solver = "transport",
  solver_optns = list(method = "networkflow")
)
res_transport
#> Method: transport 
#> 
#> Indices:
#>        x1        x2 
#> 0.3868914 0.1212555
```

### Large data or bootstrapping: `solver = "sinkhorn"`

For big experiments or bootstrap confidence intervals,
`solver = "sinkhorn"` is usually the best trade-off.

``` r
res_sinkhorn <- ot_indices(
  x, Y, M = 20,
  solver = "sinkhorn",
  solver_optns = list(epsilon = 0.05, numIterations = 1000),
  boot = TRUE, R = 200
)
res_sinkhorn
#> Method: sinkhorn 
#> 
#> Indices:
#>        x1        x2 
#> 0.7305059 0.4970830 
#> 
#> Type of confidence interval: norm 
#> Number of replicates: 200 
#> Confidence level: 0.95 
#> Bootstrap statistics:
#>   input  original        bias    low.ci   high.ci
#> 1    x1 0.7388703 0.008364458 0.7138344 0.7471773
#> 2    x2 0.5098147 0.012731695 0.4696061 0.5245600
```

### Regularization and the speed-accuracy trade-off

The main knob for Sinkhorn is the regularization parameter `epsilon`:

- Larger `epsilon` yields faster convergence but more smoothing, so
  indices may be less sensitive to fine structure.
- Smaller `epsilon` yields solutions closer to non-regularized OT but is
  slower and may require `sinkhorn_log` for stability.

You can explore this trade-off by running the same problem with
different `epsilon` values.

``` r
res_eps_fast <- ot_indices(
  x, Y, M = 20,
  solver = "sinkhorn",
  solver_optns = list(epsilon = 0.1)
)

res_eps_precise <- ot_indices(
  x, Y, M = 20,
  solver = "sinkhorn_log",
  solver_optns = list(epsilon = 0.005, numIterations = 2000)
)

res_eps_fast
#> Method: sinkhorn 
#> 
#> Indices:
#>        x1        x2 
#> 0.8460829 0.6780586
res_eps_precise
#> Method: sinkhorn_log 
#> 
#> Indices:
#>        x1        x2 
#> 0.4432757 0.1757531
```

## Dimension 3: Wasserstein-Bures indices

When you want more structure in the indices,
[`ot_indices_wb()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices_wb.md)
provides a decomposition into advective and diffusive components. This
is useful when the goal is not only ranking inputs but also
understanding how the output distribution changes. The function is
defined for any output dimensionality.

``` r
res_wb <- ot_indices_wb(x, Y, M = 20)
res_wb
#> Method: wass-bures 
#> 
#> Indices:
#>         x1         x2 
#> 0.35592501 0.08260068 
#> 
#> Advective component:
#>         x1         x2 
#> 0.28010831 0.07353109 
#> 
#> Diffusive component:
#>          x1          x2 
#> 0.075816697 0.009069586
```

## Final checklist

- Output is 1D:
  [`ot_indices_1d()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices_1d.md)
  or
  [`ot_indices_wb()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices_wb.md)
- Output is multi-D and you want standard indices:
  [`ot_indices()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices.md)
- Output is multi-D and you want Wasserstein-Bures components:
  [`ot_indices_wb()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices_wb.md)
- Small to medium data: `ot_indices(..., solver = "transport")`
- Large data or bootstrap: `ot_indices(..., solver = "sinkhorn")`

This workflow matches the intended use of the solvers in `gsaot` and
should cover most practical scenarios.
