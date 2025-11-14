# Reparameterised von Mises-Fisher distribution

Density, distribution function, and random generation for the von
Mises-Fisher distribution.

## Usage

``` r
dvmf2(x, theta, log = FALSE)

rvmf2(n, theta)
```

## Arguments

- x:

  unit vector or matrix (with each row being a unit vector) of
  evaluation points

- theta:

  direction and concentration vector. The direction of `theta`
  determines the mean direction on the sphere. The norm of `theta` is
  the concentration parameter of the distribution.

- log:

  logical; if `TRUE`, densities are returned on the log scale.

- n:

  number of random values to return.

## Value

`dvmf` gives the density and `rvm` generates random deviates.

## Details

In this parameterisation, \\\theta = \kappa \mu\\, where \\\mu\\ is a
unit vector and \\\kappa\\ is the concentration parameter.

`dvmf2` allows for automatic differentiation with `RTMB`. `rvmf2` is
imported from
[`movMF::rmovMF`](https://rdrr.io/pkg/movMF/man/movMF_distribution.html).

## Examples

``` r
set.seed(123)
# single parameter set
theta <- c(1,2,3)
x <- rvmf2(1, theta)
d <- dvmf2(x, theta)

# vectorised over parameters
theta <- matrix(theta, nrow = 1)
theta <- theta[rep(1,10), ]
x <- rvmf2(10, theta)
d <- dvmf2(x, theta)
```
