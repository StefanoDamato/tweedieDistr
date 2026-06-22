#include "utils.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


arma::vec recycle_to_length(const arma::vec& x, int l, const std::string& name) {
  if (x.n_elem == static_cast<arma::uword>(l)) return x;
  if (x.n_elem == 1) {
    arma::vec out(l);
    out.fill(x[0]);
    return out;
  }
  stop("`%s` must have length 1 or length %d.", name.c_str(), l);
}


arma::vec dpois_vec(const arma::vec& k, const arma::vec& lambda, bool log_p) {
  arma::vec out(k.n_elem, arma::fill::none);
  for (arma::uword i = 0; i < k.n_elem; ++i)
    out[i] = R::dpois(k[i], lambda[i], log_p);
  return out;
}


arma::vec pgamma_vec(const arma::vec& x, const arma::vec& shape,
                     const arma::vec& beta) {
  arma::vec out(x.n_elem, arma::fill::none);
  for (arma::uword i = 0; i < x.n_elem; ++i)
    out[i] = R::pgamma(x[i], shape[i], 1.0 / beta[i], true, false);
  return out;
}


arma::vec compound_Poisson_Gamma(const arma::vec& x,
                                  const arma::vec& lambda,
                                  const arma::vec& alpha,
                                  const arma::vec& beta) {
  arma::uword n = x.n_elem;
  arma::vec l_mode = arma::floor(lambda);
  l_mode.elem(arma::find(l_mode < 1)).fill(1);
  arma::vec log_p_mode = dpois_vec(l_mode, lambda, true);

  // Seed with the point mass at 0 and the mode term
  arma::vec cdf = arma::exp(-lambda);
  cdf += arma::exp(log_p_mode) % pgamma_vec(x, alpha % l_mode, beta);

  // Iterate downward from the mode
  arma::vec l_low = l_mode - 1;
  arma::uvec idx_low = arma::find(l_low > 0);
  while (!idx_low.is_empty()) {
    arma::vec l_low_active = l_low.elem(idx_low);
    arma::vec log_p_low    = dpois_vec(l_low_active, lambda.elem(idx_low), true);
    arma::uvec keep = arma::find(log_p_mode.elem(idx_low) - log_p_low < 37.0);
    if (keep.is_empty()) break;

    arma::uvec active   = idx_low.elem(keep);
    arma::vec l_active  = l_low.elem(active);
    cdf.elem(active) += arma::exp(log_p_low.elem(keep)) %
      pgamma_vec(x.elem(active), alpha.elem(active) % l_active, beta.elem(active));
    l_low.elem(active) -= 1;
    idx_low = arma::find(l_low > 0);
  }

  // Iterate upward from the mode
  arma::vec l_high = l_mode + 1;
  arma::uvec idx_high = arma::regspace<arma::uvec>(0, n - 1);
  while (!idx_high.is_empty()) {
    arma::vec l_high_active = l_high.elem(idx_high);
    arma::vec log_p_high    = dpois_vec(l_high_active, lambda.elem(idx_high), true);
    arma::uvec keep = arma::find(log_p_mode.elem(idx_high) - log_p_high < 37.0);
    if (keep.is_empty()) break;

    arma::uvec active  = idx_high.elem(keep);
    arma::vec l_active = l_high.elem(active);
    cdf.elem(active) += arma::exp(log_p_high.elem(keep)) %
      pgamma_vec(x.elem(active), alpha.elem(active) % l_active, beta.elem(active));
    l_high.elem(active) += 1;
    idx_high = active;
  }

  return arma::clamp(cdf, 0.0, 1.0);
}
