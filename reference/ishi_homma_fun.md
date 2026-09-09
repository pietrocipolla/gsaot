# Ishigami-Homma function evaluation

Evaluates the Ishigami-Homma function. Input samples are drawn from a
uniform distribution over \\\[-\pi, \pi\]^3\\

## Usage

``` r
ishi_homma_fun(N, A = 2, B = 1)
```

## Arguments

- N:

  Number of input samples to generate.

- A:

  (default: `2`) Numeric, amplitude of the second sine component .

- B:

  (default: `1`) Numeric, coefficient of the interaction term.

## Value

A list with two elements:

- `x`: a numeric matrix of size `N x 8` containing the input samples.

- `y`: a numeric vector of length `N` with the corresponding function
  outputs.

## Details

The Ishigami-Homma function is defined as: \$\$Y = \sin(X_1) + A \cdot
\sin^2(X_2) + B \cdot X_3^4 \cdot \sin(X_1)\$\$ where \\X_i \sim
\mathcal{U}(-\pi, \pi)\\.

## See also

[`sobol_fun`](https://pietrocipolla.github.io/gsaot/reference/sobol_fun.md),
[`gaussian_fun`](https://pietrocipolla.github.io/gsaot/reference/gaussian_fun.md)

## Examples

``` r
result <- ishi_homma_fun(1000)
head(result$x)
#>              X1         X2         X3
#> [1,]  2.5694307  2.0266196 -3.0235712
#> [2,] -3.1331797 -2.5348989 -2.4485520
#> [3,]  0.7957511 -1.9834278  2.1374600
#> [4,] -1.6037587 -2.3952594  0.5626499
#> [5,] -2.5238342  3.0391805 -0.5978263
#> [6,] -0.1244027  0.6199593  1.0293164
head(result$y)
#> [1] 47.4061162  0.3393381 17.3044482 -0.1776725 -0.6322897  0.4117596
```
