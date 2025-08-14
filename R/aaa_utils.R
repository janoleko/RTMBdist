# getting OSA residual and simulation functions from RTMB (not exported)
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
iszero <- RTMB::ADjoint(f = function(x) as.numeric(x==0),
                        df = function(x,y,dy) rep(0, length(x)),
                        name = "iszero")
# 1 if x != 0, 0 otherwise
isnonzero <- function(x) {
  1 - iszero(x)
}
# 1 if x => 0, 0 otherwise
ispos <- function(x) {
  s <- sign(x)
  0.5 * (s + abs(s))
}
# 1 if x < 0, 0 otherwise
isneg <- function(x) {
  s <- sign(x)
  - 0.5 * (s - abs(s))
}
# 1 if x > 0, 0 otherwise
ispos_strict <- function(x) {
  ispos(x) * isnonzero(x)
}
# 1 if x < val, 0 otherwise
smaller <- function(x, val) {
  s <- sign(x - val)
  0.5 * (abs(s) - s)
}
# 1 if x > val, 0 otherwise
greater <- function(x, val) {
  s <- sign(val - x)
  0.5 * (abs(s) - s)
}
# turns +Inf into largest finite value
as.finite <- function(x) {
  pmin.ad(x, .Machine$double.xmax)
}
as.finite.neg <- function(x) {
  pmax.ad(x, -.Machine$double.xmax)
}


## Logarithm of zero-inflated density/ pmf
# x == 0: p0
# x > 0: (1-p0) * pdf(x)
log_zi <- function(x, logdens, zeroprob) {
  logdens <- as.finite(logdens) # turn + Inf into finite
  logdens <- RTMB::logspace_add(
    log(iszero(x)) + log(zeroprob),
    # log(ispos_strict(x)) + log1p(-zeroprob) + logdens
    log(isnonzero(x)) + log1p(-zeroprob) + logdens
  )
}
# x == 0: p0 + pmf(0)
# x > 0: (1-p0) * pmf(x)
log_zi_discrete <- function(x, logdens, zeroprob) {
  RTMB::logspace_add(
    log(iszero(x)) + log(zeroprob),
    log1p(-zeroprob) + logdens
    )
}
