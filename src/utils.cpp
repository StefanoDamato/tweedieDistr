#include "utils.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// Change the shape of a vector repeating its items
arma::vec recycle_to_length(const arma::vec& x, int l, const std::string& name) {
  if (x.n_elem == static_cast<arma::uword>(l)) return x;
  if (x.n_elem == 1) {
    arma::vec out(l);
    out.fill(x[0]);
    return out;
  }
  stop("`%s` must have length 1 or length %d.", name.c_str(), l);
}

// Compute the lambda of the rate distirbution
arma::vec cpg_lambda(const arma::vec& mean, const arma::vec& dispersion,
                     const arma::vec& power) {
  return arma::pow(mean, 2 - power) / (dispersion % (2 - power));
}


// Compute the shape of the Gamma distribution
arma::vec cpg_alpha(const arma::vec& power) {
  return (2 - power) / (power - 1);
}

// Compute the size of the Gamma distirbution
arma::vec cpg_beta(const arma::vec& mean, const arma::vec& dispersion,
                   const arma::vec& power) {
  return 1 / (dispersion % (power - 1) % arma::pow(mean, power - 1));
}


// Compute the CDF for a batch of points via a single pass per point, walking
// outward from the Poisson mode and accumulating on the fly
arma::vec eval_cdf_direct(const arma::vec& x, const arma::vec& lambda,
                          const arma::vec& alpha, const arma::vec& beta) {
  arma::uword n = x.n_elem;
  arma::vec cdf(n);

  for (arma::uword i = 0; i < n; ++i) {
    double lam   = lambda(i);
    double a     = alpha(i);
    double scale = 1.0 / beta(i);
    double xi    = x(i);

    double l_mode = std::floor(lam);
    if (l_mode < 1) l_mode = 1;
    double log_p_mode = R::dpois(l_mode, lam, 1);

    double acc = std::exp(-lam) +
      std::exp(log_p_mode) * R::pgamma(xi, a * l_mode, scale, 1, 0);

    // Keep adding terms downward until the Poisson weight is too small to matter
    for (double k = l_mode - 1; k > 0; k -= 1) {
      double logp = R::dpois(k, lam, 1);
      if (log_p_mode - logp >= LOG_TOL) break;
      acc += std::exp(logp) * R::pgamma(xi, a * k, scale, 1, 0);
    }

    // Keep adding terms upward until the Poisson weight is too small to matter
    for (double k = l_mode + 1; ; k += 1) {
      double logp = R::dpois(k, lam, 1);
      if (log_p_mode - logp >= LOG_TOL) break;
      acc += std::exp(logp) * R::pgamma(xi, a * k, scale, 1, 0);
    }

    cdf(i) = std::min(std::max(acc, 0.0), 1.0);
  }
  return cdf;
}
