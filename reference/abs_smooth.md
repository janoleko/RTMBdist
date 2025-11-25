# Smooth approximation to the absolute value function

Smooth approximation to the absolute value function

## Usage

``` r
abs_smooth(x, epsilon = 1e-06)
```

## Arguments

- x:

  vector of evaluation points

- epsilon:

  smoothing constant

## Value

Smooth absolute value of `x`.

## Details

We approximate the absolute value here as \$\$\vert x \vert \approx
\sqrt{x^2 + \epsilon}\$\$

## Examples

``` r
abs(0)
#> [1] 0
abs_smooth(0, 1e-4)
#> [1] 0.01
```
