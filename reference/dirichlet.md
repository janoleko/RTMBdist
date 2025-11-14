# Dirichlet distribution

Density and and random generation for the Dirichlet distribution.

## Usage

``` r
ddirichlet(x, alpha, log = FALSE)

rdirichlet(n, alpha)
```

## Arguments

- x:

  vector or matrix of quantiles. If `x` is a vector, it needs to sum to
  one. If `x` is a matrix, each row should sum to one.

- alpha:

  vector or matrix of positive shape parameters

- log:

  logical; if `TRUE`, densities \\p\\ are returned as \\\log(p)\\.

- n:

  number of random values to return.

## Value

`ddirichlet` gives the density, `rdirichlet` generates random deviates.

## Details

This implementation of `ddirichlet` allows for automatic differentiation
with `RTMB`.

## Examples

``` r
# single alpha
alpha <- c(1,2,3)
x <- rdirichlet(1, alpha)
d <- ddirichlet(x, alpha)
# vectorised over alpha
alpha <- rbind(alpha, 2*alpha)
x <- rdirichlet(2, alpha)
```
