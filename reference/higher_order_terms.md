# Higher order terms for optimal transport sensitivity indices

Compute the higher order terms as the difference between the output of
[`ot_indices()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices.md)
and the output of
[`ot_indices_wb()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices_wb.md)
computed on the same sample.

## Usage

``` r
higher_order_terms(ot_result, wb_result)
```

## Arguments

- ot_result:

  An object returned by
  [`ot_indices()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices.md).

- wb_result:

  An object returned by
  [`ot_indices_wb()`](https://pietrocipolla.github.io/gsaot/reference/ot_indices_wb.md).

## Value

An object of class `gsaot_indices` containing the higher order terms of
the Wasserstein-Bures decomposition.

## Details

The helper only computes the point estimate difference between two
already computed results. The function does not check that the ground
cost used for the `ot_result` object is the squared Euclidean one
(default). The user should therefore pay attention to this aspect when
using the function.

## Examples

``` r
dat <- gaussian_fun(1000)
ot_result <- ot_indices(dat$x, dat$y, 10)
wb_result <- ot_indices_wb(dat$x, dat$y, 10)
higher_order_terms(ot_result, wb_result)
#> Method: higher order terms 
#> 
#> Indices:
#>        X1        X2        X3 
#> 0.2218041 0.2330969 0.3075903 
```
