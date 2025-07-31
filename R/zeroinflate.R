#' Zero-inflated density constructer
#'
#' @description
#' Constructs a zero-inflated density function from a given probability density function
#'
#' @details
#' The definition of zero-inflation is different for discrete and continuous distributions.
#' For discrete distributions with p.m.f. \eqn{f} and zero-inflation probability \eqn{p}, we have
#' \deqn{\Pr(X = 0) = p + (1 - p) \cdot f(0),} and
#' \deqn{\Pr(X = x) = (1 - p) \cdot f(x), \quad x > 0.}
#'
#' For continuous distributions with p.d.f. \eqn{f}, we have
#' \deqn{f_{\text{zinfl}}(x) = p \cdot \delta_0(x) + (1 - p) \cdot f(x),}
#' where \eqn{\delta_0} is the Dirac delta function at zero.
#'
#' @param dist either a probability density function or a probability mass function
#' @param discrete logical; if \code{TRUE}, the density for \code{x = 0} will be \code{zeroprob + (1-zeroprob) * dist(0, ...)}. Otherwise it will just be \code{zeroprob}.
#' In standard cases, this will be determined automatically. For non-standard cases, set this to \code{TRUE} or \code{FALSE} depending on the type of \code{dist}. See details.
#'
#' @returns zero-inflated density function with first argument \code{x}, second argument \code{zeroprob}, and additional arguments \code{...} that will be passed to \code{dist}.
#' @importFrom RTMB ADoverload
#' @export
#'
#' @examples
#' # Zero-inflated normal distribution
#' dzinorm <- zero_inflate(dnorm)
#' dzinorm(c(NA, 0, 2), 0.5, mean = 1, sd = 1)
#'
#' # Zero-inflated Poisson distribution
#' zipois <- zero_inflate(dpois)
#' zipois(c(NA, 0, 1), 0.5, 1)
#'
#' # Non-standard case: Zero-inflated reparametrised beta distribution
#' dzibeta2 <- zero_inflate(dbeta2, discrete = FALSE)
zero_inflate <- function(dist, discrete = NULL) {
  dist_name <- deparse(substitute(dist))

  # Auto-detect whether the function is discrete
  discrete_dists <- c("dpois", "dbinom", "dnbinom", "dgeom", "dmultinom", "dhyper")
  if (is.null(discrete)) {
    discrete <- dist_name %in% discrete_dists
  }

  # Extract function argument names, excluding `x`
  dist_args <- setdiff(names(formals(dist)), "x")

  # Define core function logic for both discrete and continuous cases
  core_function <- function(x, zeroprob, dist, discrete, dist_args, ...) {
    "[<-" <- ADoverload("[<-")  # Ensure correct assignment behavior

    args <- list(...)

    # Initialize output vector
    out <- numeric(length(x))

    # Handle NA values properly
    is_na <- is.na(x)
    is_zero <- x == 0 & !is_na
    is_nonzero <- x != 0 & !is_na  # Exclude NA and 0

    # Common logic: assign probability at zero
    out[is_zero] <- zeroprob

    if (discrete) {
      # Logic for discrete distributions: Handle x = 0 case with dist(0)
      out[is_zero] <- out[is_zero] + (1 - zeroprob) * do.call(dist, c(list(0), args))
    }

    # Logic for non-zero values
    out[is_nonzero] <- (1 - zeroprob) * do.call(dist, c(list(x[is_nonzero]), args))

    # Ensure NAs remain as NA
    out[is_na] <- NA

    return(out)
  }

  # Return function with core logic encapsulated
  function(x, zeroprob, ...) {
    core_function(x, zeroprob, dist, discrete, dist_args, ...)
  }
}
