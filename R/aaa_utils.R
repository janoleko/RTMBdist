# getting one step ahead resiudals and simulation functions from RTMB (not exported)
dGenericOSA <- get("dGenericOSA", envir = asNamespace("RTMB"), inherits = FALSE)
dGenericSim <- get("dGenericSim", envir = asNamespace("RTMB"), inherits = FALSE)

# getting ad_context from RTMB (not exported)
ad_context <- get("ad_context", envir = asNamespace("RTMB"), inherits = FALSE)

#' AD-compatible error function and complementary error function
#'
#' @param x vector of evaluation points
#'
#' @returns \code{erf(x)} returns the error function and \code{erfc(x)} returns the complementary error function.
#'
#' @examples
#' erf(1)
#' erfc(1)
#' @name erf
NULL
#' @rdname erf
#' @export
#' @importFrom RTMB pnorm
erf <- function(x) {
  2 * RTMB::pnorm(x * sqrt(2)) - 1 # + eps
}
#' @rdname erf
#' @export
erfc <- function(x) {
  1 - erf(x) # + eps
}

## AD pmin/pmax helpers that work for both ad and numeric:
pmin.ad <- function(x, y) apply(cbind(x,y), 1, min)
pmax.ad <- function(x, y) apply(cbind(x,y), 1, max)

## AD-indicator constructors
# 1 if x == 0, 0 otherwise
iszero <- function(x) {
  1 - sign(x)^2
}
# 1 if x > 0, 0 otherwise
ispos <- function(x) {
  s <- sign(x)
  0.5 * (s + abs(s))
}
# 1 if x < val, 0 otherwise
smaller <- function(x, val) {
  s <- sign(x - val)
  0.5 * (abs(s) - s)
}
# 1 if x > val, 0 otherwise
greater <- function(x, val) {
  1 - smaller(x, val)
}

