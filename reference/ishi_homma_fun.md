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
#>               X1         X2          X3
#> [1,] -2.28423759 -0.9399756  0.33191817
#> [2,] -0.09561026 -0.3103324  2.08030330
#> [3,] -2.28376501  0.7043859  2.89825661
#> [4,] -2.22049303  2.2053860 -1.09533096
#> [5,]  2.22482608  2.8740950 -0.02519212
#> [6,] -2.77463554  2.8473635  0.09935267
head(result$y)
#> [1]   0.5389623  -1.6968814 -53.2895980  -0.6453498   0.9333679  -0.1906092
```
