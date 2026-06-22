#include "utils.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::vec tweedieCDF(arma::vec x, arma::vec mean, arma::vec dispersion,
                     arma::vec power) {
  int l = std::max({static_cast<int>(x.n_elem), static_cast<int>(mean.n_elem),
                    static_cast<int>(dispersion.n_elem),
                    static_cast<int>(power.n_elem)});
  x          = recycle_to_length(x,          l, "x");
  mean       = recycle_to_length(mean,       l, "mean");
  dispersion = recycle_to_length(dispersion, l, "dispersion");
  power      = recycle_to_length(power,      l, "power");

  arma::vec lambda = arma::pow(mean, 2 - power) / (dispersion % (2 - power));
  arma::vec alpha  = (2 - power) / (power - 1);
  arma::vec beta   = 1 / (dispersion % (power - 1) % arma::pow(mean, power - 1));
  arma::vec cdf(l, arma::fill::zeros);

  // x < 0: CDF is 0 (already initialised to 0)
  // x = 0: CDF equals the point mass exp(-lambda)
  arma::uvec zero_idx = arma::find(x == 0);
  if (!zero_idx.is_empty()) cdf(zero_idx) = arma::exp(-lambda(zero_idx));

  // x > 0: compound Poisson-Gamma summation
  arma::uvec pos_idx = arma::find(x > 0);
  if (!pos_idx.is_empty()) {
    cdf(pos_idx) = compound_Poisson_Gamma(
        x(pos_idx), lambda(pos_idx), alpha(pos_idx), beta(pos_idx));
  }

  return cdf;
}
