#' Box–Cox Cole and Green distribution (BCCG)
#'
#' Density, distribution function, quantile function, and random generation for
#' the Box–Cox Cole and Green distribution.
#'
#' @details
#' This implementation of \code{dbccg} allows for automatic differentiation with \code{RTMB} while the other functions are imported from \code{gamlss.dist} package.
#' See \code{gamlss.dist::\link[gamlss.dist]{BCCG}} for more details.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param mu location parameter, must be positive.
#' @param sigma scale parameter, must be positive.
#' @param nu skewness parameter (real).
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{dbccg} gives the density, \code{pbccg} gives the distribution function, \code{qbccg} gives the quantile function, and \code{rbccg} generates random deviates.

#'
#' @examples
#' x <- rbccg(5, mu = 10, sigma = 0.2, nu = 0.5)
#' d <- dbccg(x, mu = 10, sigma = 0.2, nu = 0.5)
#' p <- pbccg(x, mu = 10, sigma = 0.2, nu = 0.5)
#' q <- qbccg(p, mu = 10, sigma = 0.2, nu = 0.5)
#' @name bccg
NULL

#' @rdname bccg
#' @export
#' @import RTMB
dbccg <- function(x, mu = 1, sigma = 0.1, nu = 1, log = FALSE) {

  if (!ad_context()) {
    if (any(x < 0)) stop("BCCG is only defined for x => 0")
    if (mu <= 0 || sigma <= 0) stop("mu and sigma must be > 0")
  }

  # allow simulated references / OSA conventions consistent with your BCPE
  if (inherits(x, "simref")) {
    return(dGenericSim("dbccg", x = x, mu = mu, sigma = sigma, nu = nu, log = log))
  }
  if (inherits(x, "osa")) {
    stop("Currently, GAMLSS distributions don't support OSA residuals.")
  }

  # Standardized z without branching:
  # z = ((x/mu)^nu - 1)/(nu*sigma)  if nu != 0
  #   =  log(x/mu)/sigma            if nu == 0
  iz <- iszero(nu)
  z <- iz * (log(x / mu) / sigma) +
    (1 - iz) * (((x / mu)^nu - 1) / ((nu + .Machine$double.xmin) * sigma))

  # Log-density: (nu-1)log x - nu log mu - log sigma - 0.5 z^2 - 0.5 log(2*pi)
  logdens <- (nu - 1) * log(x) - nu * log(mu) - log(sigma) -
    0.5 * z * z - 0.5 * log(2 * pi)

  if (log) return(logdens)
  exp(logdens)
}

#' @rdname bccg
#' @export
#' @usage pbccg(q, mu = 1, sigma = 0.1, nu = 1, lower.tail = TRUE, log.p = FALSE)
#' @importFrom gamlss.dist pBCCG
pbccg <- function(q, mu = 1, sigma = 0.1, nu = 1, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(q <= 0)) stop("BCCG is only defined for q > 0")
    if (mu <= 0 || sigma <= 0) stop("mu and sigma must be > 0")
  }

  gamlss.dist::pBCCG(q, mu = mu, sigma = sigma, nu = nu,
                     lower.tail = lower.tail, log.p = log.p)
}

#' @rdname bccg
#' @export
#' @usage qbccg(p, mu = 1, sigma = 0.1, nu = 1, lower.tail = TRUE, log.p = FALSE)
#' @importFrom gamlss.dist qBCCG
qbccg <- function(p, mu = 1, sigma = 0.1, nu = 1, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(p < 0 | p > 1)) stop("p must be in [0, 1]")
    if (mu <= 0 || sigma <= 0) stop("mu and sigma must be > 0")
  }

  gamlss.dist::qBCCG(p, mu = mu, sigma = sigma, nu = nu,
                     lower.tail = lower.tail, log.p = log.p)
}

#' @rdname bccg
#' @export
#' @importFrom gamlss.dist rBCCG
rbccg <- function(n, mu = 1, sigma = 0.1, nu = 1) {

  if (!ad_context()) {
    if (length(n) != 1 || !is.finite(n) || n < 0) stop("n must be a non-negative scalar")
    if (mu <= 0 || sigma <= 0) stop("mu and sigma must be > 0")
  }

  gamlss.dist::rBCCG(n, mu = mu, sigma = sigma, nu = nu)
}
