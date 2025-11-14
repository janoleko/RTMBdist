# Zero-inflated log normal distribution

Density, distribution function, and random generation for the
zero-inflated log normal distribution.

## Usage

``` r
dzilnorm(x, meanlog = 0, sdlog = 1, zeroprob = 0, log = FALSE)

pzilnorm(q, meanlog = 0, sdlog = 1, zeroprob = 0,
         lower.tail = TRUE, log.p = FALSE)

rzilnorm(n, meanlog = 0, sdlog = 1, zeroprob = 0)

plnorm(q, meanlog = 0, sdlog = 1, lower.tail = TRUE, log.p = FALSE)
```

## Arguments

- x, q:

  vector of quantiles

- meanlog, sdlog:

  mean and standard deviation of the distribution on the log scale with
  default values of 0 and 1 respectively.

- zeroprob:

  zero-inflation probability between 0 and 1.

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \le x\]\\, otherwise,
  \\P\[X \> x\]\\.

- n:

  number of random values to return

## Value

`dzilnorm` gives the density, `pzilnorm` gives the distribution
function, and `rzilnorm` generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

## Examples

``` r
x <- rzilnorm(1, 1, 1, 0.5)
d <- dzilnorm(x, 1, 1, 0.5)
p <- pzilnorm(x, 1, 1, 0.5)
```
