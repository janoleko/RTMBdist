#' Reparameterised beta distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the beta distribution reparameterised in terms of mean and concentration.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mu mean parameter, must be in the interval from 0 to 1.
#' @param phi concentration parameter, must be positive.
#' @param log logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#'
#' @return
#' \code{dbeta2} gives the density, \code{pbeta2} gives the distribution function, \code{qbeta2} gives the quantile function, and \code{rbeta2} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x = rbeta2(1, 0.5, 1)
#' d = dbeta2(x, 0.5, 1)
#' p = pbeta2(x, 0.5, 1)
#' q = qbeta2(p, 0.5, 1)
#' @name beta2
NULL
#' @rdname beta2
#' @export
#' @importFrom RTMB dbeta
dbeta2 <- function(x, mu, phi, log = FALSE) {
  # ensure mu in [0,1]
  if (any(mu < 0 | mu > 1)) stop("mu must be in the interval [0, 1].")
  # ensure phi > 0
  if (any(phi <= 0)) stop("phi must be strictly positive.")

  shape1 <- mu * phi
  shape2 <- (1 - mu) * phi
  RTMB::dbeta(x, shape1 = shape1, shape2 = shape2, log = log)
}
#' @rdname beta2
#' @export
#' @importFrom RTMB pbeta
pbeta2 <- function(q, mu, phi) {
  # ensure mu in [0,1]
  if (any(mu < 0 | mu > 1)) stop("mu must be in the interval [0, 1].")
  # ensure phi > 0
  if (any(phi <= 0)) stop("phi must be strictly positive.")

  shape1 <- mu * phi
  shape2 <- (1 - mu) * phi
  RTMB::pbeta(q, shape1 = shape1, shape2 = shape2)
}
#' @rdname beta2
#' @export
#' @importFrom RTMB qbeta
qbeta2 <- function(p, mu, phi) {
  # ensure mu in [0,1]
  if (any(mu < 0 | mu > 1)) stop("mu must be in the interval [0, 1].")
  # ensure phi > 0
  if (any(phi <= 0)) stop("phi must be strictly positive.")

  shape1 <- mu * phi
  shape2 <- (1 - mu) * phi
  RTMB::qbeta(p, shape1 = shape1, shape2 = shape2)
}
#' @rdname beta2
#' @export
#' @importFrom stats rbeta
rbeta2 <- function(n, mu, phi) {
  # ensure mu in [0,1]
  if (any(mu < 0 | mu > 1)) stop("mu must be in the interval [0, 1].")
  # ensure phi > 0
  if (any(phi <= 0)) stop("phi must be strictly positive.")

  shape1 <- mu * phi
  shape2 <- (1 - mu) * phi
  stats::rbeta(n, shape1 = shape1, shape2 = shape2)
}


