# ---------------------------------------------------------------------------
# Boundary values
# ---------------------------------------------------------------------------

test_that("qtweedie(0) = 0 and qtweedie(1) = Inf", {
  cases <- list(
    list(mu = 1,  phi = 1,   p = 1.5),
    list(mu = 2,  phi = 0.5, p = 1.7),
    list(mu = 5,  phi = 0.3, p = 1.9),
    list(mu = 0.5, phi = 2,  p = 1.1)
  )
  for (case in cases) {
    expect_equal(
      qtweedie(0, mean = case$mu, dispersion = case$phi, power = case$p),
      0,
      label = sprintf("q=0, mu=%.1f phi=%.1f p=%.1f", case$mu, case$phi, case$p)
    )
    expect_equal(
      qtweedie(1, mean = case$mu, dispersion = case$phi, power = case$p),
      Inf,
      label = sprintf("q=1, mu=%.1f phi=%.1f p=%.1f", case$mu, case$phi, case$p)
    )
  }
})

test_that("qtweedie returns 0 for p at or below the point mass P(X = 0)", {
  cases <- list(
    list(mu = 1, phi = 1,   p = 1.5),
    list(mu = 2, phi = 0.5, p = 1.7),
    list(mu = 3, phi = 2,   p = 1.2)
  )
  for (case in cases) {
    lambda <- case$mu^(2 - case$p) / (case$phi * (2 - case$p))
    p0 <- exp(-lambda)
    expect_equal(
      qtweedie(p0, mean = case$mu, dispersion = case$phi, power = case$p),
      0,
      label = sprintf("at p0, mu=%.1f phi=%.1f p=%.1f",
                      case$mu, case$phi, case$p)
    )
    expect_equal(
      qtweedie(p0 / 2, mean = case$mu, dispersion = case$phi,
               power = case$p),
      0,
      label = sprintf("below p0, mu=%.1f phi=%.1f p=%.1f",
                      case$mu, case$phi, case$p)
    )
  }
})

# ---------------------------------------------------------------------------
# Round-trip: qtweedie(ptweedie(x)) = x on the positive axis
# ---------------------------------------------------------------------------

test_that("ptweedie(qtweedie(p)) = p on the interior of the distribution", {
  cases <- list(
    list(mu = 2,   phi = 0.5, p = 1.5),
    list(mu = 1,   phi = 1,   p = 1.7),
    list(mu = 3,   phi = 2,   p = 1.3),
    list(mu = 0.5, phi = 0.8, p = 1.6),
    list(mu = 5,   phi = 0.3, p = 1.9)
  )
  for (case in cases) {
    lambda <- case$mu^(2 - case$p) / (case$phi * (2 - case$p))
    p0     <- exp(-lambda)
    probs  <- p0 + (0.99 - p0) * c(0.1, 0.25, 0.5, 0.75, 0.9)
    q_vals <- qtweedie(probs, mean = case$mu, dispersion = case$phi,
                       power = case$p)
    expect_equal(
      ptweedie(q_vals, mean = case$mu, dispersion = case$phi, power = case$p),
      probs,
      tolerance = 1e-8,
      label = sprintf("mu=%.1f phi=%.1f p=%.1f", case$mu, case$phi, case$p)
    )
  }
})

# ---------------------------------------------------------------------------
# Monotonicity
# ---------------------------------------------------------------------------

test_that("qtweedie is non-decreasing in p", {
  p_seq <- seq(0.01, 0.99, by = 0.01)
  cases <- list(
    list(mu = 2,  phi = 0.5, p = 1.5),
    list(mu = 1,  phi = 2,   p = 1.2),
    list(mu = 10, phi = 0.1, p = 1.8)
  )
  for (case in cases) {
    q <- qtweedie(p_seq, mean = case$mu, dispersion = case$phi,
                  power = case$p)
    expect_true(
      all(diff(q) >= 0),
      label = sprintf("mu=%.1f phi=%.1f p=%.1f", case$mu, case$phi, case$p)
    )
  }
})

# ---------------------------------------------------------------------------
# lower.tail and log.p flags
# ---------------------------------------------------------------------------

test_that("qtweedie lower.tail=FALSE is consistent with lower.tail=TRUE", {
  cases <- list(
    list(mu = 2, phi = 0.5, p = 1.5),
    list(mu = 1, phi = 1,   p = 1.3)
  )
  x_vals <- c(0.5, 1, 2, 5)
  for (case in cases) {
    p_lt <- ptweedie(x_vals, mean = case$mu, dispersion = case$phi,
                     power = case$p)
    p_ut <- 1 - p_lt
    q_lt <- qtweedie(p_lt, mean = case$mu, dispersion = case$phi,
                     power = case$p, lower.tail = TRUE)
    q_ut <- qtweedie(p_ut, mean = case$mu, dispersion = case$phi,
                     power = case$p, lower.tail = FALSE)
    expect_equal(
      q_lt, q_ut,
      tolerance = 1e-8,
      label = sprintf("mu=%.1f phi=%.1f p=%.1f", case$mu, case$phi, case$p)
    )
  }
})

test_that("qtweedie log.p=TRUE accepts log-probabilities", {
  cases <- list(
    list(mu = 2, phi = 0.5, p = 1.5),
    list(mu = 1, phi = 1,   p = 1.8),
    list(mu = 3, phi = 2,   p = 1.2)
  )
  x_vals <- c(0.5, 1, 2, 4)
  for (case in cases) {
    prob  <- ptweedie(x_vals, mean = case$mu, dispersion = case$phi,
                      power = case$p)
    q_lin <- qtweedie(prob, mean = case$mu, dispersion = case$phi,
                      power = case$p, log.p = FALSE)
    q_log <- qtweedie(log(prob), mean = case$mu, dispersion = case$phi,
                      power = case$p, log.p = TRUE)
    expect_equal(
      q_lin, q_log,
      tolerance = 1e-8,
      label = sprintf("mu=%.1f phi=%.1f p=%.1f", case$mu, case$phi, case$p)
    )
  }
})

# ---------------------------------------------------------------------------
# Vectorisation
# ---------------------------------------------------------------------------

test_that("qtweedie is consistent when vectorised over p vs scalar calls", {
  prob   <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  q_vec  <- qtweedie(prob, mean = 2, dispersion = 0.5, power = 1.5)
  q_loop <- vapply(prob, qtweedie, numeric(1),
                   mean = 2, dispersion = 0.5, power = 1.5)
  expect_equal(q_vec, q_loop, tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# Cross-validation against the CRAN tweedie package
# ---------------------------------------------------------------------------
#
# tweedie::qtweedie requires scalar xi; mu and phi must be scalar or the same
# length as p.  We call it element-by-element via vapply to be safe.

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
  p_grid <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)

  for (pm in params) {
    test_that(
      sprintf(
        "qtweedie matches tweedie::qtweedie (mu=%.1f phi=%.1f xi=%.1f)",
        pm$mu, pm$phi, pm$xi
      ), {
        skip_if_not_installed("tweedie")
        lambda <- pm$mu^(2 - pm$xi) / (pm$phi * (2 - pm$xi))
        p0     <- exp(-lambda)
        probs  <- p_grid[p_grid > p0]
        expected <- vapply(
          probs, tweedie::qtweedie, numeric(1),
          xi = pm$xi, mu = pm$mu, phi = pm$phi
        )
        got <- qtweedie(
          probs, mean = pm$mu, dispersion = pm$phi, power = pm$xi
        )
        expect_equal(got, expected, tolerance = 1e-8)
      }
    )
  }

  test_that(
    "ptweedie(tweedie::qtweedie(p)) = p for all param sets", {
      skip_if_not_installed("tweedie")
      p_grid <- c(0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
      for (pm in params) {
        lambda <- pm$mu^(2 - pm$xi) / (pm$phi * (2 - pm$xi))
        p0     <- exp(-lambda)
        probs  <- p_grid[p_grid > p0]
        if (length(probs) == 0) next
        ref_q <- vapply(
          probs, tweedie::qtweedie, numeric(1),
          xi = pm$xi, mu = pm$mu, phi = pm$phi
        )
        got <- ptweedie(
          ref_q, mean = pm$mu, dispersion = pm$phi, power = pm$xi
        )
        expect_equal(
          got, probs,
          tolerance = 1e-7,
          label = sprintf(
            "mu=%.1f phi=%.1f xi=%.1f", pm$mu, pm$phi, pm$xi
          )
        )
      }
    }
  )

  test_that(
    "tweedie::ptweedie(qtweedie(p)) = p for all param sets", {
      skip_if_not_installed("tweedie")
      p_grid <- c(0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
      for (pm in params) {
        lambda <- pm$mu^(2 - pm$xi) / (pm$phi * (2 - pm$xi))
        p0     <- exp(-lambda)
        probs  <- p_grid[p_grid > p0]
        if (length(probs) == 0) next
        our_q <- qtweedie(
          probs, mean = pm$mu, dispersion = pm$phi, power = pm$xi
        )
        got <- tweedie::ptweedie(
          our_q, xi = pm$xi, mu = pm$mu, phi = pm$phi
        )
        expect_equal(
          got, probs,
          tolerance = 1e-7,
          label = sprintf(
            "mu=%.1f phi=%.1f xi=%.1f", pm$mu, pm$phi, pm$xi
          )
        )
      }
    }
  )
})
