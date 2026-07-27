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
#>               X1           X2         X3
#> [1,] -0.05443348 -1.061058398  0.3978556
#> [2,] -3.06827400 -1.312432673 -2.8215873
#> [3,]  1.33332317 -1.192814747 -2.7640902
#> [4,]  0.78278887 -1.297454234 -0.1415366
#> [5,] -2.19194536  0.006835429  1.5755733
#> [6,] -0.10317137  2.491888130  2.4603557
head(result$y)
#> [1]  1.468043 -2.846818 59.433922  2.559795 -5.824497 -3.144867
```
