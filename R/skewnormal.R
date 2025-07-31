#' Skew normal distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the skew normal distribution.
#'
#' @details
#' This implementation of \code{dskewnorm} allows for automatic differentiation with \code{RTMB} while the other functions are imported from the \code{sn} package.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param xi location parameter
#' @param omega scale parameter, must be positive.
#' @param alpha skewness parameter, +/- \code{Inf} is allowed.
#' @param log logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param ... additional parameters to be passed to the \code{sn} package functions for \code{pskewnorm} and \code{qskewnorm}.
#'
#' @return
#' \code{dskewnorm} gives the density, \code{pskewnorm} gives the distribution function, \code{qskewnorm} gives the quantile function, and \code{rskewnorm} generates random deviates.
#'
#' @examples
#' # alpha is skew parameter
#' x = rskewnorm(1, alpha = 1)
#' d = dskewnorm(x, alpha = 1)
#' p = pskewnorm(x, alpha = 1)
#' q = qskewnorm(p, alpha = 1)
#' @name skewnorm
NULL

#' @rdname skewnorm
#' @export
#' @importFrom RTMB dnorm
#' @importFrom RTMB pnorm
dskewnorm <- function(x, xi = 0, omega = 1, alpha = 0, log = FALSE) {
  z = (x - xi) / omega # standardised observation

  log_normal_density <- RTMB::dnorm(z, log = TRUE)
  log_skew_component <- log(2) - log(omega) + log(RTMB::pnorm(alpha * z))

  log_density = log_normal_density + log_skew_component

  if(log) {
    return(log_density)
  } else{
    return(exp(log_density))
  }
}

#' @rdname skewnorm
#' @export
#' @importFrom sn psn
pskewnorm <- function(q, xi = 0, omega = 1, alpha = 0, ...) {
  sn::psn(x = q, xi = xi, omega = omega, alpha = alpha, ...)
}

#' @rdname skewnorm
#' @export
#' @importFrom sn qsn
qskewnorm <- function(p, xi = 0, omega = 1, alpha = 0, ...) {
  qsn(p = p, xi = xi, omega = omega, alpha = alpha, ...)
}

#' @rdname skewnorm
#' @export
#' @importFrom sn rsn
rskewnorm <- function(n, xi = 0, omega = 1, alpha = 0) {
  rsn(n = n, xi = xi, omega = omega, alpha = alpha)
}
