#' Zero-inflated Poisson distribution
#'
#' @param x integer vector of counts
#' @param n number of random values to return.
#' @param lambda vector of (non-negative) means
#' @param zeroprob zero-inflation probability between 0 and 1
#' @param log logical; return log-density if TRUE
#'
#' @examples
#' set.seed(123)
#' x = rzipois(1, 0.5, 1)
#' d = dzipois(x, 0.5, 1)
#' @name zipois
NULL
#' @rdname zipois
#' @export
#' @importFrom RTMB logspace_add dpois
dzipois <- function(x, lambda, zeroprob, log = FALSE) {

  # if (any(lambda < 0)) stop("lambda must be >= 0")
  # if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")

  logdens <- numeric(length(x))
  zero_idx <- (x == 0)

  # Zero inflation part
  logdens[zero_idx] <- logspace_add(log(zeroprob), log(1-zeroprob) - lambda)
  logdens[!zero_idx] <- log(1 - zeroprob) + dpois(x[!zero_idx], lambda, log = TRUE)

  if (log) return(logdens)
  return(exp(logdens))
}
#' @rdname zipois
#' @importFrom stats runif rpois
#' @export
rzipois <- function(n, lambda, zeroprob) {
  if (any(lambda < 0)) stop("lambda must be >= 0")
  if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")
  u <- runif(n)
  res <- ifelse(u < zeroprob, 0, rpois(n, lambda))
  return(res)
}
