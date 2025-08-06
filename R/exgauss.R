#' Exponentially modified Gaussian distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the exponentially modified Gaussian distribution.
#'
#' @details
#' This implementation of \code{dexgauss} and \code{pexgauss} allows for automatic differentiation with \code{RTMB}.
#' \code{qexgauss} is only a wrapper for \code{gamlss.dist::qexGAUS}.
#'
#' If \eqn{X \sim N(\mu, \sigma^2)} and \eqn{Y \sim \text{Exp}(\lambda)}, then
#' \eqn{Z = X + Y} follows the exponentially modified Gaussian distribution with parameters \eqn{\mu}, \eqn{\sigma}, and \eqn{\lambda}.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param mu mean parameter of the Gaussian part
#' @param sigma standard deviation parameter of the Gaussian part, must be positive.
#' @param lambda rate parameter of the exponential part, must be positive.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dexgauss} gives the density, \code{pexgauss} gives the distribution function, \code{qexgauss} gives the quantile function, and \code{rexgauss} generates random deviates.
#'
#' @examples
#' x <- rexgauss(1, 1, 2, 2)
#' d <- dexgauss(x, 1, 2, 2)
#' p <- pexgauss(x, 1, 2, 2)
#' q <- qexgauss(p, 1, 2, 2)
#' @name exgauss
NULL

#' @rdname exgauss
#' @export
dexgauss <- function(x, mu = 0, sigma = 1, lambda = 1, log = FALSE) {
  if (!ad_context()) {
    # ensure sigma > 0, lambda > 0
    if (any(sigma <= 0)) stop("sigma must be > 0")
    if (any(lambda <= 0)) stop("lambda must be > 0")
  }

  # potentially escape to RNG or CDF
  if (inherits(x, "simref")) {
    return(dGenericSim("dexgauss", x = x, mu = mu, sigma = sigma, lambda = lambda, log = log))
  }
  if (inherits(x, "osa")) {
    return(dGenericOSA("dexgauss", x = x, mu = mu, sigma = sigma, lambda = lambda, log = log))
  }

  # logdens <- - nu + (mu - x) / nu + sigma^2 / (2 * nu^2) +
  #   log(RTMB::pnorm((x - mu)/sigma - sigma / nu))
  nu <- 1 / lambda

  # one of the rare cases for which computation of non-log-scale seems to work better
  dens <- 1 / nu * exp((mu - x) / nu + sigma^2 / (2 * nu^2)) *
    RTMB::pnorm((x - mu) / sigma - sigma / nu)

  if(log) return(log(dens))
  return(dens)
}
#' @rdname exgauss
#' @export
#' @importFrom RTMB pnorm
pexgauss <- function(q, mu = 0, sigma = 1, lambda = 1, lower.tail = TRUE, log.p = FALSE) {
  if (!ad_context()) {
    # ensure sigma > 0, lambda > 0
    if (any(sigma <= 0)) stop("sigma must be > 0")
    if (any(lambda <= 0)) stop("lambda must be > 0")
  }

  p <- RTMB::pnorm(q, mu, sigma) -
    0.5 * exp(lambda/2 * (2*mu + lambda * sigma^2 - 2*q)) *
    erfc((mu + lambda * sigma^2 - q) / (sqrt(2) * sigma))

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  return(p)
}
#' @rdname exgauss
#' @export
#' @importFrom gamlss.dist qexGAUS
qexgauss <- function(p, mu = 0, sigma = 1, lambda = 1, lower.tail = TRUE, log.p = FALSE) {
  if (!ad_context()) {
    # ensure sigma > 0, lambda > 0
    if (any(sigma <= 0)) stop("sigma must be > 0")
    if (any(lambda <= 0)) stop("lambda must be > 0")
  }

  nu <- 1 / lambda
  gamlss.dist::qexGAUS(p, mu = mu, sigma = sigma, nu = nu,
                       lower.tail = lower.tail, log.p = log.p)
}
#' @rdname exgauss
#' @export
#' @importFrom stats rnorm rexp
rexgauss <- function(n, mu = 0, sigma = 1, lambda = 1) {
  # ensure sigma > 0, lambda > 0
  if (any(sigma <= 0)) stop("sigma must be > 0")
  if (any(lambda <= 0)) stop("lambda must be > 0")

  # Generate n random values from the exponentially modified Gaussian distribution
  rnorm(n, mu, sigma) + rexp(n, rate = lambda)
}



