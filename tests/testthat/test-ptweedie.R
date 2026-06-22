# ---------------------------------------------------------------------------
# Boundary values
# ---------------------------------------------------------------------------

test_that("ptweedie(0) equals exp(-lambda) for several parameterisations", {
  cases <- list(
    list(mu = 1,   phi = 1,   p = 1.5),
    list(mu = 2,   phi = 0.5, p = 1.7),
    list(mu = 0.5, phi = 2,   p = 1.2),
    list(mu = 5,   phi = 0.3, p = 1.9),
    list(mu = 3,   phi = 1,   p = 1.1)
  )
  for (case in cases) {
    lambda <- case$mu^(2 - case$p) / (case$phi * (2 - case$p))
    expect_equal(
      ptweedie(0, mean = case$mu, dispersion = case$phi, power = case$p),
      exp(-lambda),
      tolerance = 1e-8,
      label = sprintf("mu=%.1f phi=%.1f p=%.1f", case$mu, case$phi, case$p)
    )
  }
})

test_that("ptweedie returns 0 for q < 0", {
  x_neg <- c(-100, -10, -2, -0.5, -0.001)
  for (x in x_neg) {
    expect_equal(
      ptweedie(x, mean = 2, dispersion = 0.5, power = 1.5),
      0,
      label = paste0("q = ", x)
    )
  }
})

test_that("ptweedie(Inf) = 1", {
  cases <- list(
    list(mu = 1, phi = 1,   p = 1.5),
    list(mu = 3, phi = 0.5, p = 1.7)
  )
  for (case in cases) {
    expect_equal(
      ptweedie(Inf, mean = case$mu, dispersion = case$phi, power = case$p),
      1,
      label = sprintf("mu=%.1f phi=%.1f p=%.1f", case$mu, case$phi, case$p)
    )
  }
})

# ---------------------------------------------------------------------------
# Monotonicity and range
# ---------------------------------------------------------------------------

test_that("ptweedie is non-decreasing on a fine grid", {
  q <- c(0, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10, 20, 50)
  cases <- list(
    list(mu = 1.5, phi = 0.8, p = 1.6),
    list(mu = 0.5, phi = 2,   p = 1.3),
    list(mu = 5,   phi = 0.2, p = 1.9)
  )
  for (case in cases) {
    cdf <- ptweedie(q, mean = case$mu, dispersion = case$phi, power = case$p)
    expect_true(
      all(diff(cdf) >= 0),
      label = sprintf("mu=%.1f phi=%.1f p=%.1f", case$mu, case$phi, case$p)
    )
  }
})

test_that("ptweedie values are in [0, 1]", {
  q <- c(seq(-2, 0, by = 0.5), seq(0, 20, by = 0.5))
  cases <- list(
    list(mu = 2, phi = 0.5, p = 1.5),
    list(mu = 1, phi = 2,   p = 1.1)
  )
  for (case in cases) {
    cdf <- ptweedie(q, mean = case$mu, dispersion = case$phi, power = case$p)
    expect_true(all(cdf >= 0))
    expect_true(all(cdf <= 1))
  }
})

# ---------------------------------------------------------------------------
# lower.tail and log.p flags
# ---------------------------------------------------------------------------

test_that("ptweedie lower.tail = FALSE returns 1 - CDF", {
  q  <- c(0, 0.1, 0.5, 1, 2, 5, 10)
  cases <- list(
    list(mu = 2, phi = 0.5, p = 1.5),
    list(mu = 1, phi = 1,   p = 1.7)
  )
  for (case in cases) {
    lt <- ptweedie(q, mean = case$mu, dispersion = case$phi,
                   power = case$p, lower.tail = TRUE)
    ut <- ptweedie(q, mean = case$mu, dispersion = case$phi,
                   power = case$p, lower.tail = FALSE)
    expect_equal(
      lt + ut, rep(1, length(q)),
      tolerance = 1e-8,
      label = sprintf("mu=%.1f phi=%.1f p=%.1f", case$mu, case$phi, case$p)
    )
  }
})

test_that("ptweedie log.p = TRUE returns log(CDF)", {
  q <- c(0, 0.1, 0.5, 1, 2, 5)
  cases <- list(
    list(mu = 2, phi = 0.5, p = 1.5),
    list(mu = 3, phi = 1,   p = 1.3)
  )
  for (case in cases) {
    cdf     <- ptweedie(q, mean = case$mu, dispersion = case$phi,
                        power = case$p, log.p = FALSE)
    log_cdf <- ptweedie(q, mean = case$mu, dispersion = case$phi,
                        power = case$p, log.p = TRUE)
    expect_equal(
      log_cdf, log(cdf),
      tolerance = 1e-8,
      label = sprintf("mu=%.1f phi=%.1f p=%.1f", case$mu, case$phi, case$p)
    )
  }
})

# ---------------------------------------------------------------------------
# Consistency with dtweedie (numerical integration)
# ---------------------------------------------------------------------------

test_that("ptweedie(b) - ptweedie(a) equals integral of dtweedie over (a,b)", {
  cases <- list(
    list(mu = 1, phi = 1,   p = 1.5, a = 0.5, b = 3),
    list(mu = 2, phi = 0.5, p = 1.7, a = 1,   b = 5),
    list(mu = 3, phi = 2,   p = 1.3, a = 0.2, b = 4)
  )
  for (case in cases) {
    integ <- integrate(
      dtweedie, lower = case$a, upper = case$b,
      mean = case$mu, dispersion = case$phi, power = case$p,
      subdivisions = 500L, rel.tol = 1e-10
    )
    dp <- ptweedie(case$b, mean = case$mu, dispersion = case$phi,
                   power = case$p) -
          ptweedie(case$a, mean = case$mu, dispersion = case$phi,
                   power = case$p)
    expect_equal(
      integ$value, dp,
      tolerance = 1e-8,
      label = sprintf(
        "mu=%.1f phi=%.1f p=%.1f a=%.1f b=%.1f",
        case$mu, case$phi, case$p, case$a, case$b
      )
    )
  }
})

# ---------------------------------------------------------------------------
# Vectorisation
# ---------------------------------------------------------------------------

test_that("ptweedie is consistent when vectorised over q vs scalar calls", {
  q      <- c(0, 0.1, 0.5, 1, 2, 5, 10)
  p_vec  <- ptweedie(q, mean = 2, dispersion = 0.5, power = 1.5)
  p_loop <- vapply(q, ptweedie, numeric(1),
                   mean = 2, dispersion = 0.5, power = 1.5)
  expect_equal(p_vec, p_loop, tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# Cross-validation against the CRAN tweedie package
# ---------------------------------------------------------------------------

local({
  params <- list(
    list(mu = 1,    phi = 1,   xi = 1.5),
    list(mu = 2,    phi = 0.5, xi = 1.5),
    list(mu = 3,    phi = 2,   xi = 1.2),
    list(mu = 0.5,  phi = 0.8, xi = 1.8),
    list(mu = 5,    phi = 0.3, xi = 1.3),
    list(mu = 0.2,  phi = 1.5, xi = 1.7),
    list(mu = 10,   phi = 0.1, xi = 1.9),
    list(mu = 1,    phi = 2,   xi = 1.1)
  )
  x_grid <- c(0, 0.1, 0.5, 1, 2, 3, 5, 10, 20)

  for (pm in params) {
    test_that(
      sprintf(
        "ptweedie matches tweedie::ptweedie (mu=%.1f phi=%.1f xi=%.1f)",
        pm$mu, pm$phi, pm$xi
      ), {
        skip_if_not_installed("tweedie")
        expected <- tweedie::ptweedie(
          x_grid, xi = pm$xi, mu = pm$mu, phi = pm$phi
        )
        got <- ptweedie(
          x_grid, mean = pm$mu, dispersion = pm$phi, power = pm$xi
        )
        expect_equal(got, expected, tolerance = 1e-8)
      }
    )
  }

  test_that(
    "ptweedie lower.tail=FALSE matches 1 - tweedie::ptweedie for all params",
    {
      skip_if_not_installed("tweedie")
      q_test <- c(0.1, 0.5, 1, 2, 5, 10)
      for (pm in params) {
        ref <- tweedie::ptweedie(
          q_test, xi = pm$xi, mu = pm$mu, phi = pm$phi
        )
        got <- ptweedie(
          q_test, mean = pm$mu, dispersion = pm$phi, power = pm$xi,
          lower.tail = FALSE
        )
        expect_equal(
          got, 1 - ref,
          tolerance = 1e-8,
          label = sprintf(
            "mu=%.1f phi=%.1f xi=%.1f", pm$mu, pm$phi, pm$xi
          )
        )
      }
    }
  )
})
