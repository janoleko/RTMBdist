
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RTMBdist

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/janoleko/RTMBdist/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/janoleko/RTMBdist/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The `RTMBdist` package extends the functionality of the `RTMB` framework
by providing a collection of non-standard probability distributions that
are compatible with automatic differentiation (AD). While `RTMB` enables
flexible and efficient modelling - including random effects - its
built-in support is limited to standard distributions. This package
fills that gap by offering additional, AD-compatible distributions,
broadening the range of models that can be implemented and estimated
using `RTMB`.

## Installation

You can install the development version of RTMBdist from
[GitHub](https://github.com/) with:

``` r
devtools::install_github("janoleko/RTMBdist")
```

## Example

Let’s pretend we want to do numerical maximum likelihood estimation
(MLE) for a gamma distribution that is parametrised in terms of mean and
standard deviation.

``` r
library(RTMBdist)
#> Loading required package: RTMB

# draw data
x <- rgamma2(100, mean = 5, sd = 2)

# negative log-likelihood function
nll <- function(par) {
  mu <- exp(par[1]); ADREPORT(mu)
  sigma <- exp(par[2]); ADREPORT(sigma)
  -sum(dgamma2(x, mu, sigma, log = TRUE))
}

# automatically differentiable objective function object
obj <- MakeADFun(nll, c(log(5), log(2)), silent = TRUE)

# model fitting
opt <- nlminb(obj$par, obj$fn, obj$gr)

summary(sdreport(obj))
#>        Estimate Std. Error
#> par   1.6471688 0.03694199
#> par   0.6513479 0.07840871
#> mu    5.1922584 0.19181235
#> sigma 1.9181246 0.15039767
```
