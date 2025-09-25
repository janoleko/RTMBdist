#' Zero-truncated Poisson distribution
#'
#' Probability mass function, distribution function, and random generation for
#' the zero-truncated Poisson distribution.
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
#' @param lambda vector of (non-negative) means
#' @param log,log.p logical; return log-density if TRUE
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dztpois} gives the probability mass function, \code{pztpois} gives the distribution function, and \code{rztpois} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rztpois(1, 0.5)
#' d <- dztpois(x, 0.5)
#' p <- pztpois(x, 0.5)
#' @name ztpois
NULL
#' @rdname ztpois
#' @export
#' @importFrom RTMB dpois
dztpois <- function(x, lambda, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure lambda >= 0
    if (any(lambda < 0)) stop("lambda must be >= 0")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")){
    return(dGenericSim("dztpois", x = x, lambda = lambda, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dztpois", x = x, lambda = lambda, log=log))
  }

  log_1m_zprob <- log1p(-exp(-lambda))  # log(1 - P(X=0))
  logdens <- dpois(x, lambda = lambda, log = TRUE)

  logdens <- logdens - log_1m_zprob + log(ispos_strict(x))

  if (log) return(logdens)
  return(exp(logdens))
}
#' @rdname ztpois
#' @importFrom RTMB ppois dpois
#' @export
pztpois <- function(q, lambda, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    # ensure lambda >= 0
    if (any(lambda < 0)) stop("lambda must be >= 0")
    q <- floor(q)  # make sure it's integer-valued
  }

  cdf <- ppois(q, lambda)
  p0 <- exp(-lambda)
  p <- pmax.ad(cdf - p0, 0) / (1 - p0)

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  return(p)
}
#' @rdname ztpois
#' @importFrom stats runif qpois
#' @export
rztpois <- function(n, lambda) {
  if (any(lambda < 0)) stop("lambda must be >= 0")

  u <- runif(n)
  p0 <- exp(-lambda)
  x <- qpois(p0 + (1 - p0) * u, lambda)
  x
}
