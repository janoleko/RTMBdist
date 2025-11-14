# Box–Cox t distribution (BCT)

Density, distribution function, quantile function, and random generation
for the Box–Cox t distribution.

## Usage

``` r
dbct(x, mu = 5, sigma = 0.1, nu = 1, tau = 2, log = FALSE)

pbct(q, mu = 5, sigma = 0.1, nu = 1, tau = 2, lower.tail = TRUE, log.p = FALSE)

qbct(p, mu = 5, sigma = 0.1, nu = 1, tau = 2, lower.tail = TRUE, log.p = FALSE)

rbct(n, mu = 5, sigma = 0.1, nu = 1, tau = 2)
```

## Arguments

- x, q:

  vector of quantiles

- mu:

  location parameter, must be positive.

- sigma:

  scale parameter, must be positive.

- nu:

  skewness parameter (real).

- tau:

  degrees of freedom, must be positive.

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P\[X \le x\]\\,
  otherwise \\P\[X \> x\]\\.

- p:

  vector of probabilities

- n:

  number of random values to return

## Value

`dbct` gives the density, `pbct` gives the distribution function, `qbct`
gives the quantile function, and `rbct` generates random deviates.

## Details

This implementation of `dbct` and `pbct` allows for automatic
differentiation with `RTMB` while the other functions are imported from
`gamlss.dist` package. See
`gamlss.dist::`[`BCT`](https://rdrr.io/pkg/gamlss.dist/man/BCt.html) for
more details.

## References

Rigby, R. A., Stasinopoulos, D. M., Heller, G. Z., and De Bastiani, F.
(2019) Distributions for modeling location, scale, and shape: Using
GAMLSS in R, Chapman and Hall/CRC, doi:10.1201/9780429298547. An older
version can be found in https://www.gamlss.com/.

## Examples

``` r
x <- rbct(1, mu = 10, sigma = 0.2, nu = 0.5, tau = 4)
d <- dbct(x, mu = 10, sigma = 0.2, nu = 0.5, tau = 4)
p <- pbct(x, mu = 10, sigma = 0.2, nu = 0.5, tau = 4)
q <- qbct(p, mu = 10, sigma = 0.2, nu = 0.5, tau = 4)
```
