#' Power Exponential distribution (PE and PE2)
#'
#' Density, distribution function, quantile function, and random generation for
#' the Power Exponential distribution (two versions).
#'
#' @details
#' This implementation of \code{dpe} and \code{dpe2} allows for automatic differentiation with \code{RTMB} while the other functions are imported from \code{gamlss.dist} package.
#'
#' For \code{PE}, \code{mu} is the mean and \code{sigma} is the standard deviation while this does not hold for \code{PE2}.
#'
#' See \code{gamlss.dist::\link[gamlss.dist]{PE}} for more details.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param mu location parameter
#' @param sigma scale parameter, must be positive
#' @param nu shape parameter (real)
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}
#'
#' @return
#' \code{dpe} gives the density, \code{ppe} gives the distribution function, \code{qpe} gives the quantile function, and \code{rpe} generates random deviates.
#'
#' @examples
#' # PE
#' x <- rpe(5, mu = 0, sigma = 1, nu = 2)
#' d <- dpe(x, mu = 0, sigma = 1, nu = 2)
#' p <- ppe(x, mu = 0, sigma = 1, nu = 2)
#' q <- qpe(p, mu = 0, sigma = 1, nu = 2)
#'
#' # PE2
#' x <- rpe2(5, mu = 0, sigma = 1, nu = 2)
#' d <- dpe2(x, mu = 0, sigma = 1, nu = 2)
#' p <- ppe2(x, mu = 0, sigma = 1, nu = 2)
#' q <- qpe2(p, mu = 0, sigma = 1, nu = 2)
#' @name pe
NULL

#' @rdname pe
#' @export
#' @import RTMB
dpe <- function(x, mu = 0, sigma = 1, nu = 2, log = FALSE) {

  if (!ad_context()) {
    if (sigma <= 0) stop("sigma must be > 0")
  }

  # constant c
  logc <- (lgamma(1 / nu) - lgamma(3 / nu)) / 2

  # standardized variable
  z <- (x - mu) / (exp(logc) * sigma)

  # log-density
  log_const <- log(nu) - log(2) - log(sigma) - logc - lgamma(1 / nu)
  logdens <- log_const - abs(z)^nu

  if (log) return(logdens)
  exp(logdens)
}

#' @rdname pe
#' @export
#' @usage ppe(q, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, log.p = FALSE)
#' @importFrom gamlss.dist pPE
ppe <- function(q, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (sigma <= 0) stop("sigma must be > 0")
  }

  gamlss.dist::pPE(q, mu = mu, sigma = sigma, nu = nu,
                  lower.tail = lower.tail, log.p = log.p)
}

#' @rdname pe
#' @export
#' @usage qpe(p, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, log.p = FALSE)
#' @importFrom gamlss.dist qPE
qpe <- function(p, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(p < 0 | p > 1)) stop("p must be in [0, 1]")
    if (sigma <= 0) stop("sigma must be > 0")
  }

  gamlss.dist::qPE(p, mu = mu, sigma = sigma, nu = nu,
                   lower.tail = lower.tail, log.p = log.p)
}

#' @rdname pe
#' @export
#' @importFrom gamlss.dist rPE
rpe <- function(n, mu = 0, sigma = 1, nu = 2) {

  if (!ad_context()) {
    if (length(n) != 1 || !is.finite(n) || n < 0) stop("n must be a non-negative scalar")
    if (sigma <= 0) stop("sigma must be > 0")
  }

  gamlss.dist::rPE(n, mu = mu, sigma = sigma, nu = nu)
}

#' @rdname pe
#' @export
#' @import RTMB
dpe2 <- function(x, mu = 0, sigma = 1, nu = 2, log = FALSE) {

  if (!ad_context()) {
    if (sigma <= 0) stop("sigma must be > 0")
  }

  # standardized variable
  z <- (x - mu) / sigma

  # log-density
  log_const <- log(nu) - log(2) - log(sigma) - lgamma(1 / nu)
  logdens <- log_const - abs(z)^nu

  if (log) return(logdens)
  exp(logdens)
}

#' @rdname pe
#' @export
#' @usage ppe2(q, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, log.p = FALSE)
#' @importFrom gamlss.dist pPE2
ppe2 <- function(q, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (sigma <= 0) stop("sigma must be > 0")
  }

  gamlss.dist::pPE2(q, mu = mu, sigma = sigma, nu = nu,
                    lower.tail = lower.tail, log.p = log.p)
}

#' @rdname pe
#' @export
#' @usage qpe2(p, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, log.p = FALSE)
#' @importFrom gamlss.dist qPE2
qpe2 <- function(p, mu = 0, sigma = 1, nu = 2, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) {
    if (any(p < 0 | p > 1)) stop("p must be in [0, 1]")
    if (sigma <= 0) stop("sigma must be > 0")
  }

  gamlss.dist::qPE2(p, mu = mu, sigma = sigma, nu = nu,
                    lower.tail = lower.tail, log.p = log.p)
}

#' @rdname pe
#' @export
#' @importFrom gamlss.dist rPE2
rpe2 <- function(n, mu = 0, sigma = 1, nu = 2) {

  if (!ad_context()) {
    if (length(n) != 1 || !is.finite(n) || n < 0) stop("n must be a non-negative scalar")
    if (sigma <= 0) stop("sigma must be > 0")
  }

  gamlss.dist::rPE2(n, mu = mu, sigma = sigma, nu = nu)
}
