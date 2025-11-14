# Reparameterised zero-inflated beta distribution

Density, distribution function, and random generation for the
zero-inflated beta distribution reparameterised in terms of mean and
concentration.

## Usage

``` r
dzibeta2(x, mu, phi, zeroprob = 0, log = FALSE)

pzibeta2(q, mu, phi, zeroprob = 0, lower.tail = TRUE, log.p = FALSE)

rzibeta2(n, mu, phi, zeroprob = 0)
```

## Arguments

- x, q:

  vector of quantiles

- mu:

  mean parameter, must be in the interval from 0 to 1.

- phi:

  concentration parameter, must be positive.

- zeroprob:

  zero-inflation probability between 0 and 1.

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P\[X \le x\]\\,
  otherwise \\P\[X \> x\]\\.

- n:

  number of random values to return.

- p:

  vector of probabilities

## Value

`dzibeta2` gives the density, `pzibeta2` gives the distribution
function, and `rzibeta2` generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

## Examples

``` r
set.seed(123)
x <- rzibeta2(1, 0.5, 1, 0.5)
d <- dzibeta2(x, 0.5, 1, 0.5)
p <- pzibeta2(x, 0.5, 1, 0.5)
```
