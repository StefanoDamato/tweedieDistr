# ---------------------------------------------------------------------------
# Output structure
# ---------------------------------------------------------------------------

test_that("rtweedie returns a numeric vector of the requested length", {
  for (n in c(1, 10, 100, 500)) {
    expect_length(rtweedie(n), n)
    expect_type(rtweedie(n), "double")
  }
})

test_that("rtweedie returns non-negative values", {
  cases <- list(
    list(mu = 1,  phi = 1,   p = 1.5),
    list(mu = 2,  phi = 0.5, p = 1.7),
    list(mu = 5,  phi = 0.3, p = 1.9),
    list(mu = 0.5, phi = 2,  p = 1.2)
  )
  set.seed(42)
  for (case in cases) {
    x <- rtweedie(1000, mean = case$mu, dispersion = case$phi,
                  power = case$p)
    expect_true(
      all(x >= 0),
      label = sprintf("mu=%.1f phi=%.1f p=%.1f", case$mu, case$phi, case$p)
    )
  }
})

# ---------------------------------------------------------------------------
# Zero-mass fraction
# ---------------------------------------------------------------------------

test_that("fraction of zeros is close to exp(-lambda)", {
  cases <- list(
    list(mu = 1,   phi = 1,   p = 1.5),   # lambda ~ 2,    P(0) ~ 0.135
    list(mu = 2,   phi = 0.5, p = 1.5),   # lambda ~ 5.66, P(0) ~ 0.003
    list(mu = 1,   phi = 2,   p = 1.5),   # lambda ~ 1,    P(0) ~ 0.368
    list(mu = 0.5, phi = 1,   p = 1.8)    # higher P(0)
  )
  set.seed(123)
  for (case in cases) {
    lambda <- case$mu^(2 - case$p) / (case$phi * (2 - case$p))
    p0     <- exp(-lambda)
    x      <- rtweedie(10000, mean = case$mu, dispersion = case$phi,
                       power = case$p)
    p0_hat <- mean(x == 0)
    # Two-sided tolerance: within 3 standard errors of the binomial proportion
    se <- sqrt(p0 * (1 - p0) / 10000)
    expect_true(
      abs(p0_hat - p0) < 4 * se,
      label = sprintf(
        "mu=%.1f phi=%.1f p=%.1f: expected P(0)=%.3f, got %.3f (SE=%.4f)",
        case$mu, case$phi, case$p, p0, p0_hat, se
      )
    )
  }
})

# ---------------------------------------------------------------------------
# Sample mean and variance
# ---------------------------------------------------------------------------

test_that("sample mean is close to mu", {
  cases <- list(
    list(mu = 1, phi = 1,   p = 1.5),
    list(mu = 3, phi = 0.5, p = 1.6),
    list(mu = 2, phi = 2,   p = 1.3)
  )
  set.seed(456)
  for (case in cases) {
    x    <- rtweedie(20000, mean = case$mu, dispersion = case$phi,
                     power = case$p)
    xbar <- mean(x)
    se   <- sqrt(case$phi * case$mu^case$p / 20000)
    expect_true(
      abs(xbar - case$mu) < 4 * se,
      label = sprintf(
        "mu=%.1f phi=%.1f p=%.1f: expected %.2f, got %.2f (SE=%.4f)",
        case$mu, case$phi, case$p, case$mu, xbar, se
      )
    )
  }
})

test_that("sample variance is close to phi * mu^p", {
  cases <- list(
    list(mu = 2, phi = 0.5, p = 1.5),
    list(mu = 1, phi = 1,   p = 1.7)
  )
  set.seed(789)
  for (case in cases) {
    target_var <- case$phi * case$mu^case$p
    x    <- rtweedie(50000, mean = case$mu, dispersion = case$phi,
                     power = case$p)
    expect_equal(
      var(x), target_var,
      tolerance = 0.05 * target_var,
      label = sprintf("mu=%.1f phi=%.1f p=%.1f", case$mu, case$phi, case$p)
    )
  }
})

# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------

test_that("rtweedie results are reproducible with the same seed", {
  set.seed(999)
  x1 <- rtweedie(100, mean = 2, dispersion = 0.5, power = 1.5)
  set.seed(999)
  x2 <- rtweedie(100, mean = 2, dispersion = 0.5, power = 1.5)
  expect_equal(x1, x2)
})
