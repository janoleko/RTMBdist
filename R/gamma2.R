#' Reparametrised gamma distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the gamma distribution reparametrised in terms of mean and standard deviation.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mean mean parameter, must be positive.
#' @param sd standard deviation parameter, must be positive.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dgamma2} gives the density, \code{pgamma2} gives the distribution function, \code{qgamma2} gives the quantile function, and \code{rgamma2} generates random deviates.
#'
#' @examples
#' x = rgamma2(1)
#' d = dgamma2(x)
#' p = pgamma2(x)
#' q = qgamma2(p)
#' @name gamma2
NULL

#' @rdname gamma2
#' @export
#' @importFrom RTMB dgamma
dgamma2 = function(x, mean = 1, sd = 1, log = FALSE) {
  # ensure mean, sd > 0
  if (any(mean <= 0)) stop("mean must be strictly positive.")
  if (any(sd <= 0)) stop("sd must be strictly positive.")

  shape = mean^2 / sd^2
  scale = sd^2 / mean
  dgamma(x = x, shape = shape, scale = scale, log = log)
}

#' @rdname gamma2
#' @export
#' @importFrom RTMB pgamma
pgamma2 = function(q, mean = 1, sd = 1, lower.tail = TRUE, log.p = FALSE) {
  # ensure mean, sd > 0
  if (any(mean <= 0)) stop("mean must be strictly positive.")
  if (any(sd <= 0)) stop("sd must be strictly positive.")

  shape = mean^2 / sd^2
  scale = sd^2 / mean
  pgamma(q = q, shape = shape, scale = scale, lower.tail = lower.tail, log.p = log.p)
}

#' @rdname gamma2
#' @export
#' @importFrom RTMB qgamma
qgamma2 = function(p, mean = 1, sd = 1, lower.tail = TRUE, log.p = FALSE) {
  # ensure mean, sd > 0
  if (any(mean <= 0)) stop("mean must be strictly positive.")
  if (any(sd <= 0)) stop("sd must be strictly positive.")

  shape = mean^2 / sd^2
  scale = sd^2 / mean
  qgamma(p = p, shape = shape, scale = scale, lower.tail = lower.tail, log.p = log.p)
}

#' @rdname gamma2
#' @export
#' @importFrom stats rgamma
rgamma2 = function(n, mean = 1, sd = 1) {
  # ensure mean, sd > 0
  if (any(mean <= 0)) stop("mean must be strictly positive.")
  if (any(sd <= 0)) stop("sd must be strictly positive.")

  shape = mean^2 / sd^2
  scale = sd^2 / mean
  rgamma(n = n, shape = shape, scale = scale)
}
