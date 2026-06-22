#include "utils.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// Stirling-type weight W(j) on the log scale.
static arma::vec get_log_W(const arma::vec& alpha, int j,
                            const arma::vec& constant_log_W) {
  return j * (constant_log_W - (1 + alpha) * log(j)) - log(2 * M_PI) -
         0.5 * arma::log(alpha) - log(j);
}


// Series A(y) from Dunn & Smyth (2005), evaluated on the log scale.
static arma::vec log_A(const arma::vec& y, const arma::vec& phi,
                        const arma::vec& rho) {
  arma::vec alpha = (2 - rho) / (rho - 1);
  arma::vec log_z = alpha % arma::log(y) - alpha % arma::log(rho - 1) -
                    arma::log(2 - rho) - (1 + alpha) % arma::log(phi);

  arma::vec j_max        = arma::pow(y, 2 - rho) / (phi % (2 - rho));
  arma::vec log_W_max    = j_max % (1 + alpha) - log(2 * M_PI) -
                           0.5 * arma::log(alpha) - arma::log(j_max);
  arma::vec constant_log_W = log_z + (1 + alpha) - alpha % arma::log(alpha);

  // Expand upper summation bound until the tail is negligible
  int j_U = std::max(1., ceil(arma::max(j_max)));
  arma::vec log_W_U = get_log_W(alpha, j_U, constant_log_W);
  while (any(log_W_max - log_W_U < 37)) {
    j_U     = j_U + 1;
    log_W_U = get_log_W(alpha, j_U, constant_log_W);
  }

  // Expand lower summation bound until the tail is negligible
  int j_L = std::max(1., floor(arma::min(j_max)));
  arma::vec log_W_L = get_log_W(alpha, j_L, constant_log_W);
  while (any(log_W_max - log_W_L < 37) & (j_L > 1)) {
    j_L     = j_L - 1;
    log_W_L = get_log_W(alpha, j_L, constant_log_W);
  }

  arma::vec j = arma::linspace(j_L, j_U, j_U - j_L + 1);
  arma::mat log_W = j * log_z.t();
  log_W.each_col() -= arma::lgamma(j + 1);
  log_W           -= arma::lgamma(j * alpha.t());

  // Log-sum-exp for numerical stability
  arma::rowvec max_log_W   = arma::max(log_W, 0);
  arma::mat exp_stabilized = arma::exp(log_W.each_row() - max_log_W);
  arma::rowvec log_sum_w   = max_log_W + arma::log(arma::sum(exp_stabilized, 0));

  return log_sum_w.t() - arma::log(y);
}


// [[Rcpp::export]]
arma::vec tweedieDensity(arma::vec x, arma::vec mean, arma::vec dispersion,
                         arma::vec power, bool log) {
  int l = std::max({static_cast<int>(x.n_elem), static_cast<int>(mean.n_elem),
                    static_cast<int>(dispersion.n_elem),
                    static_cast<int>(power.n_elem)});
  x          = recycle_to_length(x,          l, "x");
  mean       = recycle_to_length(mean,       l, "mean");
  dispersion = recycle_to_length(dispersion, l, "dispersion");
  power      = recycle_to_length(power,      l, "power");
  arma::vec log_p(l, arma::fill::none);

  // x < 0: density is 0
  arma::uvec neg_idx = arma::find(x < 0);
  if (!neg_idx.is_empty()) log_p(neg_idx).fill(-arma::datum::inf);

  // x = 0: density is exp(-lambda)
  arma::uvec zero_idx = arma::find(x == 0);
  if (!zero_idx.is_empty()) {
    log_p(zero_idx) =
        -(arma::pow(mean(zero_idx), 2 - power(zero_idx))) /
        (dispersion(zero_idx) % (2 - power(zero_idx)));
  }

  // x > 0: series expansion (Dunn & Smyth 2005)
  arma::uvec pos_idx = arma::find(x > 0);
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

  return log ? log_p : arma::exp(log_p);
}
