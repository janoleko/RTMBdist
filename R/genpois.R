#' Generalised Poisson distribution
#'
#' Probability mass function, distribution function, and random generation for
#' the generalised Poisson distribution.
#'
#' @details
#' This implementation of \code{dgenpois} allows for automatic differentiation with \code{RTMB}.
#' The other functions are imported from \code{gamlss.dist::GPO}.
#'
#' The distribution has mean \eqn{\lambda} and variance \eqn{\lambda(1 + \phi \lambda)^2}.
#' For \eqn{\phi = 0} it reduces to the Poisson distribution, however \eqn{\phi} must be strictly positive here.
#'
#' @param x,q integer vector of counts
#' @param p vector of probabilities
#' @param n number of random values to return.
#' @param lambda vector of positive means
#' @param phi vector of non-negative dispersion parameters
#' @param log,log.p logical; return log-density if TRUE
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#' @param max.value a constant, set to the default value of 10000 for how far the algorithm should look for \code{q}.
#'
#' @return
#' \code{dgenpois} gives the probability mass function, \code{pgenpois} gives the distribution function, \code{qgenpois} gives the quantile function, and \code{rgenpois} generates random deviates.
#'
#' @examples
#' set.seed(123)
#' x <- rgenpois(1, 2, 3)
#' d <- dgenpois(x, 2, 3)
#' p <- pgenpois(x, 2, 3)
#' q <- qgenpois(p, 2, 3)
#' @name genpois
NULL
#' @rdname genpois
#' @export
#' @import RTMB
dgenpois <- function(x, lambda = 1, phi = 1, log = FALSE) {

  if(!ad_context()) {
    # ensure lambda, phi > 0
    if (any(lambda <= 0)) stop("lambda must be > 0")
    if (any(phi <= 0)) stop("phi must be > 0")
    # Check x is integer >= 0
    if(any(x < 0 | x != floor(x))) {
      stop("x must be a non-negative integer")
    }
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")){
    return(dGenericSim("dgenpois", x = x, lambda = lambda, phi = phi, log=log))
  }
  if(inherits(x, "osa")) {
    stop("Currently, generalised Poisson does not support OSA residuals.")
  }

  phi_lambda_p1 <- 1 + phi * lambda
  phi_x_p1 <- 1 + phi * x

  logdens <- x * (log(lambda) - log(phi_lambda_p1)) +
    (x-1) * log(phi_x_p1) - lgamma(x + 1) -
    (lambda * (phi_x_p1)) / (phi_lambda_p1)

  if(log) return(logdens)
  return(exp(logdens))
}
#' @rdname genpois
#' @export
#' @importFrom gamlss.dist pGPO
pgenpois <- function(q, lambda = 1, phi = 1, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    # ensure lambda, phi > 0
    if (any(lambda <= 0)) stop("lambda must be > 0")
    if (any(phi <= 0)) stop("phi must be > 0")
    # Check q is integer >= 0
    if(any(q < 0 | q != floor(q))) {
      stop("q must be a non-negative integer")
    }
  }

  p <- gamlss.dist::pGPO(q, mu = lambda, sigma = phi)

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)

  return(p)
}
#' @rdname genpois
#' @export
#' @usage qgenpois(lambda = 1, phi = 1,
#'          lower.tail = TRUE, log.p = FALSE, max.value = 1e4)
#' @importFrom gamlss.dist qGPO
qgenpois <- function(p, lambda = 1, phi = 1, lower.tail = TRUE, log.p = FALSE, max.value = 1e4) {

  if(!ad_context()) {
    # ensure lambda, phi > 0
    if (any(lambda <= 0)) stop("lambda must be > 0")
    if (any(phi <= 0)) stop("phi must be > 0")
    # Check p is in [0,1]
    if(any(p < 0 | p > 1)) {
      stop("p must be in [0,1]")
    }
  }

  gamlss.dist::qGPO(p, mu = lambda, sigma = phi,
                    lower.tail = lower.tail, log.p = log.p, max.value = max.value)
}
#' @rdname genpois
#' @export
#' @importFrom gamlss.dist rGPO
rgenpois <- function(n, lambda = 1, phi = 1, max.value = 1e4) {

  if(!ad_context()) {
    # ensure lambda, phi > 0
    if (any(lambda <= 0)) stop("lambda must be > 0")
    if (any(phi <= 0)) stop("phi must be > 0")
  }

  gamlss.dist::rGPO(n, mu = lambda, sigma = phi)
}

