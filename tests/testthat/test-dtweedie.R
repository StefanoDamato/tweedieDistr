# ---------------------------------------------------------------------------
# Point mass at 0
# ---------------------------------------------------------------------------

test_that("dtweedie(0) equals exp(-lambda) for several parameterisations", {
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
      dtweedie(0, mean = case$mu, dispersion = case$phi, power = case$p),
      exp(-lambda),
      tolerance = 1e-8,
      label = sprintf("mu=%.1f phi=%.1f p=%.1f", case$mu, case$phi, case$p)
    )
  }
})

# ---------------------------------------------------------------------------
# Support: zero and negative x
# ---------------------------------------------------------------------------

test_that("dtweedie is 0 for x < 0", {
  x_neg <- c(-10, -2, -1, -0.5, -0.001)
  for (x in x_neg) {
    expect_equal(
      dtweedie(x, mean = 2, dispersion = 0.5, power = 1.5),
      0,
      label = paste0("x = ", x)
    )
  }
})

test_that("dtweedie is 0 at +Inf and returns 0 (not NaN)", {
  expect_equal(dtweedie(Inf,  mean = 2, dispersion = 0.5, power = 1.5), 0)
  expect_equal(dtweedie(-Inf, mean = 1, dispersion = 1,   power = 1.5), 0)
})

# ---------------------------------------------------------------------------
# Non-negativity on the positive axis
# ---------------------------------------------------------------------------

test_that("dtweedie >= 0 on a fine grid for several parameterisations", {
  x <- c(0, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10, 20, 50)
  cases <- list(
    list(mu = 1,   phi = 1,   p = 1.5),
    list(mu = 2,   phi = 0.5, p = 1.3),
    list(mu = 0.5, phi = 2,   p = 1.8),
    list(mu = 5,   phi = 0.2, p = 1.9),
    list(mu = 1,   phi = 0.5, p = 1.1)
  )
  for (case in cases) {
    d <- dtweedie(x, mean = case$mu, dispersion = case$phi, power = case$p)
    expect_true(
      all(d >= 0, na.rm = TRUE),
      label = sprintf("mu=%.1f phi=%.1f p=%.1f", case$mu, case$phi, case$p)
    )
  }
})

# ---------------------------------------------------------------------------
# Continuous part integrates to 1 - P(X = 0)
# ---------------------------------------------------------------------------

test_that("dtweedie integrates to 1 - exp(-lambda) over (0, Inf)", {
  cases <- list(
    list(mu = 1, phi = 1,   p = 1.5),
    list(mu = 2, phi = 0.5, p = 1.7),
    list(mu = 3, phi = 2,   p = 1.3)
  )
  for (case in cases) {
    lambda <- case$mu^(2 - case$p) / (case$phi * (2 - case$p))
    integ <- integrate(
      dtweedie, lower = 0, upper = Inf,
      mean = case$mu, dispersion = case$phi, power = case$p,
      subdivisions = 500L, rel.tol = 1e-10
    )
    expect_equal(
      integ$value,
      1 - exp(-lambda),
      tolerance = 1e-8,
      label = sprintf("mu=%.1f phi=%.1f p=%.1f", case$mu, case$phi, case$p)
    )
  }
})

# ---------------------------------------------------------------------------
# log = TRUE
# ---------------------------------------------------------------------------

test_that("dtweedie log = TRUE equals log(dtweedie) for several cases", {
  cases <- list(
    list(mu = 1,   phi = 1,   p = 1.5),
    list(mu = 2,   phi = 0.5, p = 1.7),
    list(mu = 0.5, phi = 2,   p = 1.2)
  )
  x <- c(0.01, 0.1, 0.5, 1, 2, 5, 10)
  for (case in cases) {
    d     <- dtweedie(x, mean = case$mu, dispersion = case$phi,
                      power = case$p, log = FALSE)
    d_log <- dtweedie(x, mean = case$mu, dispersion = case$phi,
                      power = case$p, log = TRUE)
    expect_equal(
      d_log, log(d),
      tolerance = 1e-8,
      label = sprintf("mu=%.1f phi=%.1f p=%.1f", case$mu, case$phi, case$p)
    )
  }
})

test_that("dtweedie log = TRUE returns -Inf at x = +Inf", {
  expect_equal(
    dtweedie(Inf, mean = 2, dispersion = 0.5, power = 1.5, log = TRUE),
    -Inf
  )
})

# ---------------------------------------------------------------------------
# Vectorisation
# ---------------------------------------------------------------------------

test_that("dtweedie is consistent when vectorised over x vs scalar calls", {
  x <- c(0, 0.1, 0.5, 1, 2, 5, 10)
  d_vec    <- dtweedie(x, mean = 2, dispersion = 0.5, power = 1.5)
  d_scalar <- vapply(x, dtweedie, numeric(1),
                     mean = 2, dispersion = 0.5, power = 1.5)
  expect_equal(d_vec, d_scalar, tolerance = 1e-8)
})

test_that("dtweedie is consistent when vectorised over mean", {
  means <- c(0.5, 1, 2, 3, 5)
  d_vec <- dtweedie(1, mean = means, dispersion = 0.5, power = 1.5)
  d_loop <- vapply(means, function(m) {
    dtweedie(1, mean = m, dispersion = 0.5, power = 1.5)
  }, numeric(1))
  expect_equal(d_vec, d_loop, tolerance = 1e-8)
})

test_that("dtweedie is consistent when vectorised over dispersion", {
  phis <- c(0.2, 0.5, 1, 2, 5)
  d_vec <- dtweedie(1, mean = 2, dispersion = phis, power = 1.5)
  d_loop <- vapply(phis, function(phi) {
    dtweedie(1, mean = 2, dispersion = phi, power = 1.5)
  }, numeric(1))
  expect_equal(d_vec, d_loop, tolerance = 1e-8)
})

test_that("dtweedie is consistent when vectorised over power", {
  powers <- c(1.1, 1.3, 1.5, 1.7, 1.9)
  d_vec  <- dtweedie(1, mean = 2, dispersion = 0.5, power = powers)
  d_loop <- vapply(powers, function(p) {
    dtweedie(1, mean = 2, dispersion = 0.5, power = p)
  }, numeric(1))
  expect_equal(d_vec, d_loop, tolerance = 1e-8)
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
  x_grid <- c(0, 0.01, 0.05, 0.1, 0.5, 1, 2, 3, 5, 10, 20)

  for (pm in params) {
    test_that(
      sprintf(
        "dtweedie matches tweedie::dtweedie (mu=%.1f phi=%.1f xi=%.1f)",
        pm$mu, pm$phi, pm$xi
      ), {
        skip_if_not_installed("tweedie")
        expected <- tweedie::dtweedie(
          x_grid, xi = pm$xi, mu = pm$mu, phi = pm$phi
        )
        got <- dtweedie(
          x_grid, mean = pm$mu, dispersion = pm$phi, power = pm$xi
        )
        expect_equal(got, expected, tolerance = 1e-8)
      }
    )
  }

  test_that(
    "dtweedie log=TRUE matches log(tweedie::dtweedie) on positive x", {
      skip_if_not_installed("tweedie")
      x_pos <- c(0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10)
      for (pm in params) {
        ref <- tweedie::dtweedie(
          x_pos, xi = pm$xi, mu = pm$mu, phi = pm$phi
        )
        got_log <- dtweedie(
          x_pos, mean = pm$mu, dispersion = pm$phi, power = pm$xi,
          log = TRUE
        )
        expect_equal(
          got_log, log(ref),
          tolerance = 1e-8,
          label = sprintf(
            "mu=%.1f phi=%.1f xi=%.1f", pm$mu, pm$phi, pm$xi
          )
        )
      }
    }
  )
})
