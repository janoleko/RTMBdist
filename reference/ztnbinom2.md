# Reparameterised zero-truncated negative binomial distribution

Probability mass function, distribution function, quantile function, and
random generation for the zero-truncated negative binomial distribution
reparameterised in terms of mean and size.

## Usage

``` r
dztnbinom2(x, mu, size, log = FALSE)

pztnbinom2(q, mu, size, lower.tail = TRUE, log.p = FALSE)

rztnbinom2(n, mu, size)
```

## Arguments

- x, q:

  integer vector of counts

- mu:

  mean parameter, must be positive

- size:

  size/dispersion parameter, must be positive

- log, log.p:

  logical; return log-density if TRUE

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \le x\]\\, otherwise,
  \\P\[X \> x\]\\.

- n:

  number of random values to return.

## Value

`dztnbinom2` gives the probability mass function, `pztnbinom2` gives the
distribution function, and `rztnbinom2` generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

By definition, this distribution only has support on the positive
integers (1, 2, ...). Any zero-truncated distribution is defined as
\$\$P(X=x \| X\>0) = P(X=x) / (1 - P(X=0)),\$\$ where \\P(X=x)\\ is the
probability mass function of the corresponding untruncated distribution.

## Examples

``` r
set.seed(123)
x <- rztnbinom2(1, mu = 2, size = 1)
d <- dztnbinom2(x, mu = 2, size = 1)
p <- pztnbinom2(x, mu = 2, size = 1)
```
