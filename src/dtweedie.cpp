#include "utils.h"
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// Terms beyond LOG_TOL log-units below the peak contribute less than machine
// epsilon (exp(-37) ~ 8.5e-17) and are dropped from the series.
static const double LOG_TOL = 37.0;


// Compute the series expansion in log(A)
static arma::vec log_A(const arma::vec& y, const arma::vec& phi,
                        const arma::vec& rho) {
  const arma::uword n = y.n_elem;
  arma::vec out(n);

  // Make a vector of unique rho values for each item
  arma::vec urho = arma::unique(rho);
  for (arma::uword g = 0; g < urho.n_elem; ++g) {
    const double r     = urho(g);
    const double alpha = (2.0 - r) / (r - 1.0);
    const arma::uvec idx = arma::find(rho == r);

    // For each item also compute log(z) and j_max
    const arma::vec yg   = y(idx);
    const arma::vec phig = phi(idx);
    const arma::vec log_z = alpha * arma::log(yg) - alpha * std::log(r - 1.0) -
                            std::log(2.0 - r) - (1.0 + alpha) * arma::log(phig);
    const arma::vec j_max = arma::pow(yg, 2.0 - r) / (phig * (2.0 - r));

    // Precompute the lgamma tables over the union of log-gamma windows
    const double jmax_lo = j_max.min(), jmax_hi = j_max.max();
    const double hw_lo = std::sqrt(2.0 * LOG_TOL * jmax_lo / (1.0 + alpha));
    const double hw_hi = std::sqrt(2.0 * LOG_TOL * jmax_hi / (1.0 + alpha));
    const arma::uword j_lo = std::max<double>(1.0, std::floor(jmax_lo - hw_lo) - 5.0);
    const arma::uword j_hi = static_cast<arma::uword>(std::ceil(jmax_hi + hw_hi)) + 5;

    // Fill the table with the values
    const arma::uword tbl_n = j_hi - j_lo + 1;
    arma::vec lg1(tbl_n), lga(tbl_n);
    for (arma::uword t = 0; t < tbl_n; ++t) {
      const double j = static_cast<double>(j_lo + t);
      lg1(t) = std::lgamma(j + 1.0);
      lga(t) = std::lgamma(j * alpha);
    }

    // Take w_j from the table for each observation
    auto term = [&](arma::uword j, double lz) -> double {
      double lg1j, lgaj;
      if (j >= j_lo && j <= j_hi) {
        lg1j = lg1(j - j_lo);
        lgaj = lga(j - j_lo);
      } else {

        // Use a fallback when the value is not stored
        lg1j = std::lgamma(static_cast<double>(j) + 1.0);
        lgaj = std::lgamma(static_cast<double>(j) * alpha);
      }
      return static_cast<double>(j) * lz - lg1j - lgaj;
    };

    for (arma::uword k = 0; k < idx.n_elem; ++k) {
      const double lz = log_z(k);
      arma::uword jstar = static_cast<arma::uword>(
          std::max<double>(1.0, std::llround(j_max(k))));

      // Compute the sum in log-scale subcrating its maximum
      double M = -arma::datum::inf, S = 0.0;
      auto accumulate = [&](double w) {
        if (w > M) { S = S * std::exp(M - w) + 1.0; M = w; }
        else       { S += std::exp(w - M); }
      };

      // Iterate the computation of terms of the sum on the right
      for (arma::uword j = jstar; ; ++j) {
        const double w = term(j, lz);
        if (j > jstar && w < M - LOG_TOL) break;
        accumulate(w);
      }

      // Iterate the computation of terms of the sum on the left
      for (arma::uword j = jstar; j-- > 1; ) {
        const double w = term(j, lz);
        if (w < M - LOG_TOL) break;
        accumulate(w);
      }

      // Compute the sum
      out(idx(k)) = M + std::log(S) - std::log(yg(k));
    }
  }
  return out;
}


// Compute the density function
// [[Rcpp::export]]
arma::vec tweedieDensity(arma::vec x, arma::vec mean, arma::vec dispersion,
                         arma::vec power, bool log) {

  // Reshape all vectors to the same size
  int l = std::max({static_cast<int>(x.n_elem), static_cast<int>(mean.n_elem),
                    static_cast<int>(dispersion.n_elem),
                    static_cast<int>(power.n_elem)});
  x          = recycle_to_length(x,          l, "x");
  mean       = recycle_to_length(mean,       l, "mean");
  dispersion = recycle_to_length(dispersion, l, "dispersion");
  power      = recycle_to_length(power,      l, "power");

  // Initialise to NaN so that non-finite values propagate
  arma::vec log_p(l);
  log_p.fill(arma::datum::nan);

  // For negative or infinite values the density is 0
  arma::uvec neg_idx = arma::find(x < 0);
  if (!neg_idx.is_empty()) log_p(neg_idx).fill(-arma::datum::inf);
  arma::uvec posinf_idx = arma::find(x == arma::datum::inf);
  if (!posinf_idx.is_empty()) log_p(posinf_idx).fill(-arma::datum::inf);

  // In 0 the density is -lambda
  arma::uvec zero_idx = arma::find(x == 0);
  if (!zero_idx.is_empty()) {
    log_p(zero_idx) =
        -(arma::pow(mean(zero_idx), 2 - power(zero_idx))) /
        (dispersion(zero_idx) % (2 - power(zero_idx)));
  }

  // For positive values use the series expansion
  arma::uvec pos_idx = arma::find((x > 0) % (x < arma::datum::inf));
  if (!pos_idx.is_empty()) {
    arma::vec x_p    = x(pos_idx);
    arma::vec mu_p   = mean(pos_idx);
    arma::vec phi_p  = dispersion(pos_idx);
    arma::vec rho_p  = power(pos_idx);
    log_p(pos_idx) =
        log_A(x_p, phi_p, rho_p) +
        ((x_p % (arma::pow(mu_p, 1 - rho_p) / (1 - rho_p))) -
         (arma::pow(mu_p, 2 - rho_p)) / (2 - rho_p)) / phi_p;
  }

  // Resturn the results
  return log ? log_p : arma::exp(log_p);
}
