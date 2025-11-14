# Wishart distribution

Density and random generation for the wishart distribution

## Usage

``` r
dwishart(x, nu, Sigma, log = FALSE)

rwishart(n, nu, Sigma)
```

## Arguments

- x:

  positive definite p x p matrix of evaluation points

- nu:

  degrees of freedom, needs to be greater than `p - 1`

- Sigma:

  scale matrix, needs to be positive definite and match the dimension of
  `x`.

- log:

  logical; if `TRUE`, densities \\p\\ are returned as \\\log(p)\\.

- n:

  number of random deviates to return

## Value

`dwishart` gives the density, `rwishart` generates random deviates
(matrix for `n = 1`, array for `n > 1`)

## Examples

``` r
x <- rwishart(1, nu = 5, Sigma = diag(3))
d <- dwishart(x, nu = 5, Sigma = diag(3))
```
