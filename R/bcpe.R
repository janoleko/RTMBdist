#' Box-Cox Power Exponential distribution (BCPE)
#'
#' Density, distribution function, quantile function, and random generation for
#' the Box-Cox Power Exponential distribution.
#'
#' @details
#' This implementation of \code{dbcpe} allows for automatic differentiation with \code{RTMB} while the other functions are imported from \code{gamlss.dist} package.
#' See \code{gamlss.dist::\link[gamlss.dist]{BCPE}} for more details.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param mu location parameter, must be positive.
#' @param sigma scale parameter, must be positive.
#' @param nu vector of \code{nu} parameter values.
#' @param tau vector of \code{tau} parameter values, must be positive.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dbcpe} gives the density, \code{pbcpe} gives the distribution function, \code{qbcpe} gives the quantile function, and \code{rbcpe} generates random deviates.
#'
#' @examples
#' x <- rbcpe(1, mu = 5, sigma = 0.1, nu = 1, tau = 1)
#' d <- dbcpe(x, mu = 5, sigma = 0.1, nu = 1, tau = 1)
#' p <- pbcpe(x, mu = 5, sigma = 0.1, nu = 1, tau = 1)
#' q <- qbcpe(p, mu = 5, sigma = 0.1, nu = 1, tau = 1)
#' @name bcpe
NULL

#' @rdname bcpe
#' @export
#' @import RTMB
dbcpe <- function(x, mu = 5, sigma = 0.1, nu = 1, tau = 2, log = FALSE) {

  if(!ad_context()) {
    if (any(x <= 0)) stop("BCPE is only defined for x > 0")
    if (mu <= 0 || sigma <= 0 || tau <= 0) stop("mu, sigma, tau must be > 0")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dbcpe", x=x, mu=mu, sigma=sigma, nu=nu, tau=tau, log=log))
  }
  if(inherits(x, "osa")) {
    stop("Currently, GAMLSS distributions don't support OSA residuals.")
  }

  # constant for scaling the PE part
  logc <- 0.5 * ((-2 / tau) * log(2) + lgamma(1 / tau) - lgamma(3 / tau))

  # standardized z
  iszero_nu <- iszero(nu)
  z <- iszero_nu * log(x / mu) / sigma +
    (1 - iszero_nu) * ((x / mu)^nu - 1) / (nu + .Machine$double.xmin * sigma)

  # log-density
  logdens <- (nu - 1) * log(x) - nu * log(mu) - log(sigma) -
    logc - (1 + 1 / tau) * log(2) - lgamma(1 / tau) +
    log(tau) - 0.5 * abs(z / exp(logc))^tau

  if(log) return(logdens)
  return(exp(logdens))
}
#' @rdname bcpe
#' @export
#' @usage pbcpe(q, mu = 5, sigma = 0.1, nu = 1, tau = 2, lower.tail = TRUE, log.p = FALSE)
#' @importFrom gamlss.dist pBCPE
pbcpe <- function(q, mu = 5, sigma = 0.1, nu = 1, tau = 2, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    if (any(x <= 0)) stop("BCPE is only defined for x > 0")
    if (mu <= 0 || sigma <= 0 || tau <= 0) stop("mu, sigma, tau must be > 0")
  }

  gamlss.dist::pBCPE(q, mu = mu, sigma = sigma, nu = nu, tau = tau,
                     lower.tail = lower.tail, log.p = log.p)
}
#' @rdname bcpe
#' @export
#' @usage qbcpe(p, mu = 5, sigma = 0.1, nu = 1, tau = 2, lower.tail = TRUE, log.p = FALSE)
#' @importFrom gamlss.dist qBCPE
qbcpe <- function(p, mu = 5, sigma = 0.1, nu = 1, tau = 2, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    if (any(x <= 0)) stop("BCPE is only defined for x > 0")
    if (mu <= 0 || sigma <= 0 || tau <= 0) stop("mu, sigma, tau must be > 0")
  }

  gamlss.dist::qBCPE(p, mu = mu, sigma = sigma, nu = nu, tau = tau,
                     lower.tail = lower.tail, log.p = log.p)
}
#' @rdname bcpe
#' @export
#' @importFrom gamlss.dist pBCPE
rbcpe <- function(n, mu = 5, sigma = 0.1, nu = 1, tau = 2) {

  if(!ad_context()) {
    if (any(x <= 0)) stop("BCPE is only defined for x > 0")
    if (mu <= 0 || sigma <= 0 || tau <= 0) stop("mu, sigma, tau must be > 0")
  }

  gamlss.dist::rBCPE(n, mu = mu, sigma = sigma, nu = nu, tau = tau)
}
