#' Reparameterised skew normal distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the skew normal distribution reparameterised in terms of mean, standard deviation and skew magnitude
#'
#' @details
#' This implementation of \code{dskewnorm2} allows for automatic differentiation with \code{RTMB} while the other functions are imported from the \code{sn} package.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param mean mean parameter
#' @param sd standard deviation, must be positive.
#' @param alpha skewness parameter, +/- \code{Inf} is allowed.
#' @param log logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param ... additional parameters to be passed to the \code{sn} package functions for \code{pskewnorm} and \code{qskewnorm}.
#'
#' @return
#' \code{dskewnorm2} gives the density, \code{pskewnorm2} gives the distribution function, \code{qskewnorm2} gives the quantile function, and \code{rskewnorm2} generates random deviates.
#'
#' @examples
#' # alpha is skew parameter
#' x <- rskewnorm2(1, alpha = 1)
#' d <- dskewnorm2(x, alpha = 1)
#' p <- pskewnorm2(x, alpha = 1)
#' q <- qskewnorm2(p, alpha = 1)
#' @name skewnorm2
NULL

#' @rdname skewnorm2
#' @export
dskewnorm2 <- function(x, mean = 0, sd = 1, alpha = 0, log = FALSE) {

  # parameter transformation
  delta <- alpha / sqrt(1 + alpha^2)
  omega <- sd / sqrt(1 - (2 * delta^2) / pi)
  xi <- mean - omega * delta * sqrt(2 / pi)

  dskewnorm(x = x, xi = xi, omega = omega, alpha = alpha, log = log)
}

#' @rdname skewnorm2
#' @export
pskewnorm2 <- function(q, mean = 0, sd = 1, alpha = 0, ...) {

  # parameter transformation
  delta <- alpha / sqrt(1 + alpha^2)
  omega <- sd / sqrt(1 - (2 * delta^2) / pi)
  xi <- mean - omega * delta * sqrt(2 / pi)

  pskewnorm(q = q, xi = xi, omega = omega, alpha = alpha, ...)
}

#' @rdname skewnorm2
#' @export
qskewnorm2 <- function(p, mean = 0, sd = 1, alpha = 0, ...) {
  # parameter transformation
  delta <- alpha / sqrt(1 + alpha^2)
  omega <- sd / sqrt(1 - (2 * delta^2) / pi)
  xi <- mean - omega * delta * sqrt(2 / pi)

  qskewnorm(p = p, xi = xi, omega = omega, alpha = alpha, ...)
}

#' @rdname skewnorm2
#' @export
rskewnorm2 <- function(n, mean = 0, sd = 1, alpha = 0) {
  # parameter transformation
  delta <- alpha / sqrt(1 + alpha^2)
  omega <- sd / sqrt(1 - (2 * delta^2) / pi)
  xi <- mean - omega * delta * sqrt(2 / pi)

  rskewnorm(n = n, xi = xi, omega = omega, alpha = alpha)
}
