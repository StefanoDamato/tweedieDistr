#' Tweedie Distribution
#'
#' Construct a Tweedie distribution object using the compound Poisson--Gamma
#' parameterisation with power parameter in \eqn{(1, 2)}. The Tweedie family
#' is a subclass of exponential dispersion models that naturally produces exact
#' zeros (via the Poisson count component) mixed with continuous positive
#' values (via the Gamma severity component), making it well suited to
#' intermittent demand data.
#'
#' The density is evaluated using the series expansion of Dunn & Smyth (2005),
#' implemented in C++ for performance.
#'
#' @param mean Mean parameter \eqn{\mu > 0}.
#' @param dispersion Dispersion parameter \eqn{\phi > 0}.
#' @param power Power parameter \eqn{p \in (1, 2)}.
#'
#' @return A `distributional` distribution object of class `dist_tweedie`.
#'
#' @references
#' Dunn, P. K., & Smyth, G. K. (2005). Series evaluation of Tweedie
#' exponential dispersion model densities. *Statistics and Computing*,
#' 15(4), 267--280. \doi{10.1007/s11222-005-4070-y}.
#'
#' @export
#'
#' @importFrom rlang abort
#' @importFrom distributional new_dist covariance
#' @importFrom stats rpois rgamma
#'
#' @examples
#' d <- dist_tweedie(mean = 2, dispersion = 0.8, power = 1.5)
#' d |> mean()
#' d |> quantile(c(0.5, 0.9))
#' d |> density(c(0, 1.5, 3))
#' d |> distributional::variance()
#' d |> distributional::generate(10)
dist_tweedie <- function(mean = 1, dispersion = 1, power = 1.5) {
  mean <- as.double(mean)
  dispersion <- as.double(dispersion)
  power <- as.double(power)

  if (any(mean <= 0, na.rm = TRUE)) {
    abort("The mean parameter of a Tweedie distribution must be strictly positive.")
  }
  if (any(dispersion <= 0, na.rm = TRUE)) {
    abort("The dispersion parameter of a Tweedie distribution must be strictly positive.")
  }
  if (any(power <= 1 | power >= 2, na.rm = TRUE)) {
    abort("The power parameter of a Tweedie distribution must be in (1, 2).")
  }

  new_dist(mu = mean, phi = dispersion, p = power, class = "dist_tweedie")
}

#' @noRd
#' @export
format.dist_tweedie <- function(x, digits = 2, ...) {
  sprintf(
    "Tweedie(%s, %s, %s)",
    format(x[["mu"]], digits = digits, ...),
    format(x[["phi"]], digits = digits, ...),
    format(x[["p"]], digits = digits, ...)
  )
}

#' @importFrom stats density
#' @exportS3Method distributional::density
#' @export
#' @noRd
density.dist_tweedie <- function(x, at, ...) {
  dtweedie(at,
    mean = x[["mu"]],
    dispersion = x[["phi"]],
    power = x[["p"]],
    log = FALSE
  )
}

#' @importFrom distributional generate
#' @exportS3Method distributional::generate
#' @noRd
generate.dist_tweedie <- function(x, times, ...) {
  rtweedie(times,
    mean = x[["mu"]],
    dispersion = x[["phi"]],
    power = x[["p"]]
  )
}

#' @exportS3Method distributional::cdf
#' @noRd
cdf.dist_tweedie <- function(x, q, lower.tail = TRUE, log.p = FALSE, ...) {
  ptweedie(q,
    mean = x[["mu"]],
    dispersion = x[["phi"]],
    power = x[["p"]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @exportS3Method distributional::quantile
#' @noRd
quantile.dist_tweedie <- function(x, p, lower.tail = TRUE, log.p = FALSE, ...) {
  qtweedie(p,
    mean = x[["mu"]],
    dispersion = x[["phi"]],
    power = x[["p"]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @export
#' @noRd 
mean.dist_tweedie <- function(x, ...) {
  x[["mu"]]
}

#' @export
#' @noRd
covariance.dist_tweedie <- function(x, ...) {
  x[["phi"]] * x[["mu"]]^x[["p"]]
}

#' @importFrom distributional skewness
#' @exportS3Method distributional::skewness
#' @noRd
skewness.dist_tweedie <- function(x, ...) {
  # Skewness of a Tweedie EDM: third cumulant / variance^(3/2)
  # k3 = p * phi^2 * mu^(2p-1); variance = phi * mu^p
  # Result: p * sqrt(phi) * mu^(p/2 - 1)
  x[["p"]] * sqrt(x[["phi"]]) * x[["mu"]]^(x[["p"]] / 2 - 1)
}

#' @importFrom distributional kurtosis
#' @exportS3Method distributional::kurtosis
#' @noRd
kurtosis.dist_tweedie <- function(x, ...) {
  # Excess kurtosis: fourth cumulant / variance^2
  # k4 = p * (2p-1) * phi^3 * mu^(3p-2); variance = phi * mu^p
  # Result: p * (2p - 1) * phi * mu^(p - 2)
  x[["p"]] * (2 * x[["p"]] - 1) * x[["phi"]] * x[["mu"]]^(x[["p"]] - 2)
}

