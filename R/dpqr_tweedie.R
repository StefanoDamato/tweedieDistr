#' Tweedie Distribution Functions
#'
#' @description
#' Density, distribution function, quantile function and random generation for
#' the Tweedie distribution with mean equal to `mean`, dispersion equal to
#' `dispersion`, and power equal to `power`.
#'
#' @details
#' If `mean`, `dispersion`, or `power` are not specified they assume the
#' default values of `1`, `1`, and `1.5`, respectively.
#'
#' The Tweedie distribution used here follows the compound Poisson-Gamma
#' parameterisation with power parameter in \eqn{(1, 2)}. It has
#' \eqn{\mathbb{E}[X] = \mu} and
#' \eqn{\mathrm{Var}(X) = \phi\mu^p}, where \eqn{\mu} is `mean`,
#' \eqn{\phi} is `dispersion`, and \eqn{p} is `power`.
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If `length(n) > 1`, the length is taken
#'   to be the number required.
#' @param mean vector of means.
#' @param dispersion vector of dispersion parameters.
#' @param power vector of power parameters.
#' @param log,log.p logical; if `TRUE`, probabilities `p` are given as `log(p)`.
#' @param lower.tail logical; if `TRUE` (default), probabilities are
#'   \eqn{P[X \le x]}; otherwise, \eqn{P[X > x]}.
#'
#' @return
#' `dtweedie` gives the density, `ptweedie` gives the distribution
#'  function, `qtweedie` gives the quantile function, and `rtweedie`
#'  generates random samples.
#'
#' The length of the result is determined by `n` for `rtweedie`, and is the
#' maximum of the lengths of the numerical arguments for the other functions.
#'
#' The numerical arguments other than `n` are recycled to the length of the
#' result. Only the first elements of the logical arguments are used.
#'
#' @references
#' Dunn, P. K., & Smyth, G. K. (2005). Series evaluation of Tweedie
#' exponential dispersion model densities. *Statistics and Computing*,
#' 15(4), 267--280. \doi{10.1007/s11222-005-4070-y}.
#'
#' @name tweedie
#' @rdname tweedie
#' @aliases dtweedie ptweedie qtweedie rtweedie
#' @export
rtweedie <- function(n, mean = 1, dispersion = 1, power = 1.5) {
  lambda <- (mean^(2 - power)) / (dispersion * (2 - power))
  alpha <- (2 - power) / (power - 1)
  beta <- 1 / (dispersion * (power - 1) * (mean^(power - 1)))

  m <- rpois(n, lambda)
  rgamma(n, m * alpha, beta)
}

#' @rdname tweedie
#' @export
dtweedie <- function(x, mean = 1, dispersion = 1, power = 1.5, log = FALSE) {
  out <- as.vector(tweedieDensity(x, mean, dispersion, power, log))
  # tweedieDensity returns NaN for x = +Inf; the true limit is 0
  inf_pos <- is.infinite(x) & x > 0
  out[inf_pos] <- if (log) -Inf else 0
  out
}

#' @rdname tweedie
#' @export
ptweedie <- function(q, mean = 1, dispersion = 1, power = 1.5, lower.tail = TRUE, log.p = FALSE) {
  cdf <- as.vector(tweedieCDF(q, mean, dispersion, power))
  if (!lower.tail) {
    cdf <- 1 - cdf
  }
  if (log.p) {
    cdf <- log(cdf)
  }
  cdf
}

#' @rdname tweedie
#' @export
qtweedie <- function(p, mean = 1, dispersion = 1, power = 1.5, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) {
    p <- exp(p)
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  out <- as.vector(tweedieInvCDF(p, mean, dispersion, power))
  # tweedieInvCDF cannot converge to +Inf; handle the p = 1 boundary explicitly
  out[p >= 1] <- Inf
  out
}
