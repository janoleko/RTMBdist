#' Reparameterised zero-truncated negative binomial distribution
#'
#' Probability mass function, distribution function, quantile function, and random generation for
#' the zero-truncated negative binomial distribution reparameterised in terms of mean and size.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' By definition, this distribution only has support on the positive integers (1, 2, ...).
#' Any zero-truncated distribution is defined as
#' \deqn{P(X=x | X>0) = P(X=x) / (1 - P(X=0)),}
#' where \eqn{P(X=x)} is the probability mass function of the corresponding untruncated distribution.
#'
#' @param x,q integer vector of counts
#' @param n number of random values to return.
#' @param mu mean parameter, must be positive
#' @param size size/dispersion parameter, must be positive
#' @param log,log.p logical; return log-density if TRUE
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dztnbinom2} gives the probability mass function, \code{pztnbinom2} gives the distribution function, and \code{rztnbinom2} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rztnbinom2(1, mu = 2, size = 1)
#' d <- dztnbinom2(x, mu = 2, size = 1)
#' p <- pztnbinom2(x, mu = 2, size = 1)
#' @name ztnbinom2
NULL

#' @rdname ztnbinom2
#' @export
dztnbinom2 <- function(x, mu, size, log = FALSE) {
  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args)
    if (any(mu <= 0)) stop("mu must be > 0")
    if (any(size <= 0)) stop("size must be > 0")
  }

  log_1m_zprob <- log1p(-dnbinom2(0, mu = mu, size = size))  # log(1-P(X=0))
  logdens <- dnbinom2(x, mu = mu, size = size, log = TRUE)
  logdens <- logdens - log_1m_zprob + log(ispos_strict(x))

  if (log) return(logdens)
  exp(logdens)
}

#' @rdname ztnbinom2
#' @export
pztnbinom2 <- function(q, mu, size, lower.tail = TRUE, log.p = FALSE) {
  if(!ad_context()) {
    q <- floor(q)
    if (any(mu <= 0)) stop("mu must be > 0")
    if (any(size <= 0)) stop("size must be > 0")
  }

  cdf <- pnbinom2(q, mu = mu, size = size)
  p0 <- dnbinom2(0, mu = mu, size = size)
  p <- pmax.ad(cdf - p0, 0) / (1 - p0)

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  p
}

#' @rdname ztnbinom2
#' @export
#' @importFrom stats runif qnbinom
rztnbinom2 <- function(n, mu, size) {
  if(!ad_context()) {
    if (any(mu <= 0)) stop("mu must be > 0")
    if (any(size <= 0)) stop("size must be > 0")
  }

  u <- runif(n)
  p0 <- dnbinom2(0, mu = mu, size = size)
  qnbinom2(p0 + (1 - p0) * u, mu = mu, size = size)
}
