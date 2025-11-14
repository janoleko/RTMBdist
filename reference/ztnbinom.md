# Zero-truncated Negative Binomial distribution

Probability mass function, distribution function, and random generation
for the zero-truncated Negative Binomial distribution.

## Usage

``` r
dztnbinom(x, size, prob, log = FALSE)

pztnbinom(q, size, prob, lower.tail = TRUE, log.p = FALSE)

rztnbinom(n, size, prob)
```

## Arguments

- x, q:

  integer vector of counts

- size:

  target for number of successful trials, or dispersion parameter (the
  shape parameter of the gamma mixing distribution). Must be strictly
  positive, need not be integer.

- prob:

  probability of success in each trial. 0 \< prob \<= 1.

- log, log.p:

  logical; return log-density if TRUE

- lower.tail:

  logical; if `TRUE`, probabilities are \\P\[X \le x\]\\, otherwise,
  \\P\[X \> x\]\\.

- n:

  number of random values to return.

## Value

`dztnbinom` gives the probability mass function, `pztnbinom` gives the
distribution function, and `rztnbinom` generates random deviates.

## Details

This implementation allows for automatic differentiation with `RTMB`.

By definition, this distribution only has support on the positive
integers (1, 2, ...). Any zero-truncated distribution is defined as
\$\$P(X=x \| X\>0) = P(X=x) / (1 - P(X=0)),\$\$ where \\P(X=x)\\ is the
probability mass function of the corresponding untruncated distribution.

## Examples

``` r
set.seed(123)
x <- rztnbinom(1, size = 2, prob = 0.5)
d <- dztnbinom(x, size = 2, prob = 0.5)
p <- pztnbinom(x, size = 2, prob = 0.5)
```
