# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

test_that("dist_tweedie creates an object of the right class", {
  d <- dist_tweedie(mean = 2, dispersion = 0.5, power = 1.5)
  expect_s3_class(d, "distribution")
  expect_true(inherits(unclass(d)[[1L]], "dist_tweedie"))
})

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

test_that("dist_tweedie rejects non-positive mean", {
  expect_error(dist_tweedie(mean = 0),  "mean")
  expect_error(dist_tweedie(mean = -1), "mean")
})

test_that("dist_tweedie rejects non-positive dispersion", {
  expect_error(dist_tweedie(dispersion = 0),  "dispersion")
  expect_error(dist_tweedie(dispersion = -1), "dispersion")
})

test_that("dist_tweedie rejects power outside (1, 2)", {
  expect_error(dist_tweedie(power = 1),   "power")
  expect_error(dist_tweedie(power = 2),   "power")
  expect_error(dist_tweedie(power = 0.5), "power")
  expect_error(dist_tweedie(power = 2.5), "power")
})

# ---------------------------------------------------------------------------
# format
# ---------------------------------------------------------------------------

test_that("format.dist_tweedie produces a readable string", {
  d   <- dist_tweedie(mean = 2, dispersion = 0.5, power = 1.5)
  fmt <- format(d)
  expect_true(grepl("Tweedie", fmt))
  expect_true(grepl("2", fmt))
})

# ---------------------------------------------------------------------------
# density / cdf / quantile / generate
# ---------------------------------------------------------------------------

test_that("density() method matches dtweedie()", {
  mu  <- 2
  phi <- 0.5
  p   <- 1.5
  d   <- dist_tweedie(mean = mu, dispersion = phi, power = p)
  at  <- c(0, 0.5, 1, 3)
  expect_equal(
    density(d, at)[[1]],
    dtweedie(at, mean = mu, dispersion = phi, power = p),
    tolerance = 1e-8
  )
})

test_that("cdf() method matches ptweedie()", {
  mu  <- 2
  phi <- 0.5
  p   <- 1.5
  d   <- dist_tweedie(mean = mu, dispersion = phi, power = p)
  at  <- c(0, 0.5, 1, 3)
  expect_equal(
    distributional::cdf(d, at)[[1]],
    ptweedie(at, mean = mu, dispersion = phi, power = p),
    tolerance = 1e-8
  )
})

test_that("quantile() method matches qtweedie()", {
  mu   <- 2
  phi  <- 0.5
  p    <- 1.5
  d    <- dist_tweedie(mean = mu, dispersion = phi, power = p)
  prob <- c(0.25, 0.5, 0.75, 0.9)
  expect_equal(
    quantile(d, p = prob)[[1]],
    qtweedie(prob, mean = mu, dispersion = phi, power = p),
    tolerance = 1e-8
  )
})

test_that("generate() method returns non-negative values of the right length", {
  d   <- dist_tweedie(mean = 2, dispersion = 0.5, power = 1.5)
  out <- distributional::generate(d, 50)[[1]]
  expect_length(out, 50)
  expect_true(all(out >= 0))
})

# ---------------------------------------------------------------------------
# mean / variance / covariance / median
# ---------------------------------------------------------------------------

test_that("mean() returns mu", {
  d <- dist_tweedie(mean = 3, dispersion = 0.8, power = 1.6)
  expect_equal(mean(d), 3)
})

test_that("variance() returns phi * mu^p", {
  mu  <- 2
  phi <- 0.5
  p   <- 1.5
  d   <- dist_tweedie(mean = mu, dispersion = phi, power = p)
  expect_equal(
    as.numeric(distributional::variance(d)),
    phi * mu^p,
    tolerance = 1e-8
  )
})

test_that("covariance() returns phi * mu^p (scalar, univariate)", {
  mu  <- 2
  phi <- 0.5
  p   <- 1.5
  d   <- dist_tweedie(mean = mu, dispersion = phi, power = p)
  expect_equal(
    as.numeric(distributional::covariance(d)),
    phi * mu^p,
    tolerance = 1e-8
  )
})

test_that("median() is consistent with qtweedie(0.5)", {
  mu  <- 2
  phi <- 0.5
  p   <- 1.5
  d   <- dist_tweedie(mean = mu, dispersion = phi, power = p)
  expect_equal(
    as.numeric(median(d)),
    qtweedie(0.5, mean = mu, dispersion = phi, power = p),
    tolerance = 1e-8
  )
})

# ---------------------------------------------------------------------------
# skewness / kurtosis
# ---------------------------------------------------------------------------

test_that("skewness() equals p * sqrt(phi) * mu^(p/2 - 1)", {
  mu  <- 2
  phi <- 0.5
  p   <- 1.5
  d   <- dist_tweedie(mean = mu, dispersion = phi, power = p)
  expect_equal(
    as.numeric(distributional::skewness(d)),
    p * sqrt(phi) * mu^(p / 2 - 1),
    tolerance = 1e-8
  )
})

test_that("skewness() is positive and decreases as mu increases", {
  phi <- 0.5
  p   <- 1.5
  s1  <- as.numeric(
    distributional::skewness(dist_tweedie(mean = 1, dispersion = phi, power = p))
  )
  s2  <- as.numeric(
    distributional::skewness(dist_tweedie(mean = 3, dispersion = phi, power = p))
  )
  expect_gt(s1, 0)
  expect_gt(s2, 0)
  expect_gt(s1, s2)
})

test_that("kurtosis() equals p * (2p - 1) * phi * mu^(p - 2)", {
  mu  <- 2
  phi <- 0.5
  p   <- 1.5
  d   <- dist_tweedie(mean = mu, dispersion = phi, power = p)
  expect_equal(
    as.numeric(distributional::kurtosis(d)),
    p * (2 * p - 1) * phi * mu^(p - 2),
    tolerance = 1e-8
  )
})

test_that("kurtosis() is positive", {
  d <- dist_tweedie(mean = 2, dispersion = 0.5, power = 1.5)
  expect_gt(as.numeric(distributional::kurtosis(d)), 0)
})

# ---------------------------------------------------------------------------
# support
# ---------------------------------------------------------------------------

test_that("support() is [0, Inf): closed at 0 (point mass), open at Inf", {
  d <- dist_tweedie(mean = 2, dispersion = 0.5, power = 1.5)
  s <- distributional::support(d)
  expect_s3_class(s, "support_region")
  expect_equal(format(s)[[1]], "[0,Inf)")
})
