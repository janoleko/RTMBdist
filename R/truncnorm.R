#' Truncated normal distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the truncated normal distribution.
#'
#' @details
#' This implementation of \code{dtruncnorm} allows for automatic differentiation with \code{RTMB}.
#'
#' \strong{Caution:} \code{x} should not be parameter dependent as this introduces a non-differentiability.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return.
#' @param mean mean parameter, must be positive.
#' @param sd standard deviation parameter, must be positive.
#' @param min,max truncation bounds.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dtruncnorm} gives the density, \code{ptruncnorm} gives the distribution function, \code{qtruncnorm} gives the quantile function, and \code{rtruncnorm} generates random deviates.
#'
#' @examples
#' x <- rtruncnorm(1, mean = 2, sd = 2, min = -1, max = 5)
#' d <- dtruncnorm(x, mean = 2, sd = 2, min = -1, max = 5)
#' p <- ptruncnorm(x, mean = 2, sd = 2, min = -1, max = 5)
#' q <- qtruncnorm(p, mean = 2, sd = 2, min = -1, max = 5)
#' @name truncnorm
NULL

#' @rdname truncnorm
#' @export
#' @importFrom RTMB dnorm pnorm
dtruncnorm <- function(x, mean, sd, min = -Inf, max = Inf, log = FALSE) {

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dtruncnorm", x=x, mean=mean, sd=sd, min=min, max=max, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dtruncnorm", x=x, mean=mean, sd=sd, min=min, max=max, log=log))
  }

  # normalisation constant: probability of being within [min, max]
  denom <- RTMB::pnorm(max, mean, sd) - RTMB::pnorm(min, mean, sd)

  # initialise log-density vector
  logdens <- rep(NaN, length(x))

  # logical vector for values inside the truncation bounds
  inside <- (x >= min) & (x <= max)

  # compute density only for values inside the bounds
  logdens[inside] <- RTMB::dnorm(x[inside], mean, sd, log = TRUE) - log(denom)

  # return log-density if requested
  if(log) {
    logdens[!inside] <- -Inf
    return(logdens)
  } else {
    dens <- exp(logdens)
    dens[!inside] <- 0
    return(dens)
  }
}

#' @rdname truncnorm
#' @export
#' @usage
#' ptruncnorm(q, mean = 0, sd = 1, min = -Inf, max = Inf,
#'            lower.tail = TRUE, log.p = FALSE)
#' @importFrom RTMB pnorm
ptruncnorm <- function(q, mean = 0, sd = 1, min = -Inf, max = Inf, lower.tail = TRUE, log.p = FALSE) {
  # if (sd <= 0) stop("Standard deviation 'sd' must be positive.")

  # CDF for values below the lower bound is 0 (or -Inf if log.p)
  below <- q < min
  above <- q > max
  inside <- !below & !above

  # normalisation constant: probability of being within [min, max]
  denom <- RTMB::pnorm(max, mean, sd) - RTMB::pnorm(min, mean, sd)

  # Compute standardized CDF
  p <- numeric(length(q))
  p[below] <- 0
  p[above] <- 1
  p[inside] <- (RTMB::pnorm(q[inside], mean, sd) - RTMB::pnorm(min, mean, sd)) / denom

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)

  return(p)
}

#' @rdname truncnorm
#' @export
#' @usage
#' qtruncnorm(p, mean = 0, sd = 1, min = -Inf, max = Inf,
#'            lower.tail = TRUE, log.p = FALSE)
#' @importFrom RTMB pnorm qnorm
qtruncnorm <- function(p, mean = 0, sd = 1, min = -Inf, max = Inf, lower.tail = TRUE, log.p = FALSE) {
  if (sd <= 0) stop("Standard deviation 'sd' must be positive.")

  # Handle log.p and lower.tail
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  # Check that probabilities are in [0, 1]
  if (any(p < 0 | p > 1)) stop("Probabilities must be between 0 and 1.")

  # Normalization constant
  denom <- RTMB::pnorm(max, mean, sd) - RTMB::pnorm(min, mean, sd)

  # Transform p into quantiles of untruncated normal
  p_untrunc <- p * denom + RTMB::pnorm(min, mean, sd)

  # Invert the untruncated normal CDF
  q <- stats::qnorm(p_untrunc, mean, sd)

  return(q)
}

#' @rdname truncnorm
#' @export
#' @importFrom RTMB pnorm qnorm
#' @importFrom stats rnorm runif
rtruncnorm <- function(n, mean = 0, sd = 1, min = -Inf, max = Inf) {
  u <- runif(n)
  left <- pnorm((min - mean) / sd)
  right <- pnorm((max - mean) / sd) - pnorm((min - mean) / sd)
  x <- qnorm(left + u * right) * sd + mean
  return(x)
}


