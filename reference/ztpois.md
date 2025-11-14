# Zero-truncated Poisson distribution

Probability mass function, distribution function, and random generation
for the zero-truncated Poisson distribution.

## Usage

``` r
dztpois(x, lambda, log = FALSE)

pztpois(q, lambda, lower.tail = TRUE, log.p = FALSE)

rztpois(n, lambda)
```

## Arguments

- x, q:

  integer vector of counts

- lambda:

  vector of (non-negative) means

- log, log.p:

  logical; return log-density if TRUE

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \le x\]\\, otherwise,
  \\P\[X \> x\]\\.

- n:

  number of random values to return.

## Value

`dztpois` gives the probability mass function, `pztpois` gives the
distribution function, and `rztpois` generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

By definition, this distribution only has support on the positive
integers (1, 2, ...). Any zero-truncated distribution is defined as
\$\$P(X=x \| X\>0) = P(X=x) / (1 - P(X=0)),\$\$ where \\P(X=x)\\ is the
probability mass function of the corresponding untruncated distribution.

## Examples

``` r
set.seed(123)
x <- rztpois(1, 0.5)
d <- dztpois(x, 0.5)
p <- pztpois(x, 0.5)
```
