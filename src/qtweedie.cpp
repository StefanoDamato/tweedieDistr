#include "utils.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// Build the Poisson terms shared by the CDF and PDF, once per call
static void build_poisson_terms(const arma::vec& lambda,
                                arma::vec& k, arma::vec& logw,
                                arma::uvec& term_start, arma::vec& base_cdf) {
  arma::uword n = lambda.n_elem;

  // Find the lower and upper bounds for the sufficiently large Poisson weights
  arma::vec l_lo(n), l_hi(n);
  for (arma::uword i = 0; i < n; ++i) {
    double lam = lambda(i);
    double l_mode = std::floor(lam);
    if (l_mode < 1) l_mode = 1;
    double log_p_mode = R::dpois(l_mode, lam, 1);

    l_lo(i) = l_mode;
    while (l_lo(i) > 1 && log_p_mode - R::dpois(l_lo(i) - 1, lam, 1) < LOG_TOL) l_lo(i) -= 1;
    l_hi(i) = l_mode;
    while (log_p_mode - R::dpois(l_hi(i) + 1, lam, 1) < LOG_TOL) l_hi(i) += 1;
  }

  // Compute the starting index of each point in the concatenated vector
  term_start.set_size(n + 1);
  term_start(0) = 0;
  for (arma::uword i = 0; i < n; ++i)
    term_start(i + 1) = term_start(i) + static_cast<arma::uword>(l_hi(i) - l_lo(i) + 1);
  k.set_size(term_start(n));
  logw.set_size(term_start(n));
  base_cdf = arma::exp(-lambda);

  // Fill the concatenated vector with the Poisson weights
  for (arma::uword i = 0; i < n; ++i) {
    double lam = lambda(i);
    double kk  = l_lo(i);
    for (arma::uword t = term_start(i); t < term_start(i + 1); ++t) {
      k(t)    = kk;
      logw(t) = R::dpois(kk, lam, 1);
      kk += 1;
    }
  }
}


// Compute the CDF for a batch of points using the saved Poisson weights
static arma::vec eval_cdf_vec(const arma::vec& k, const arma::vec& logw,
                              const arma::uvec& term_start, const arma::vec& base_cdf,
                              const arma::uvec& idx, const arma::vec& x,
                              const arma::vec& alpha, const arma::vec& beta) {
  arma::vec out(idx.n_elem);
  for (arma::uword j = 0; j < idx.n_elem; ++j) {
    arma::uword i = idx(j);
    double scale  = 1.0 / beta(j);
    double cdf    = base_cdf(i);

    // Keep adding terms until the Poisson weight is too small to matter
    for (arma::uword t = term_start(i); t < term_start(i + 1); ++t)
      cdf += std::exp(logw(t)) * R::pgamma(x(j), alpha(j) * k(t), scale, 1, 0);
    out(j) = std::min(std::max(cdf, 0.0), 1.0);
  }
  return out;
}


// Evaluate the density and CDF of the compound Poisson-Gamma given the Poisson weights
static arma::vec eval_cdf_pdf(const arma::vec& k, const arma::vec& logw,
                              const arma::uvec& term_start, const arma::vec& base_cdf,
                              arma::uword i, double x, double alpha, double beta) {
  double scale = 1.0 / beta;
  double cdf   = base_cdf(i);
  double pdf   = 0.0;

  // Keep adding terms until the Poisson weight is too small to matter
  for (arma::uword t = term_start(i); t < term_start(i + 1); ++t) {
    double w     = std::exp(logw(t));
    double shape = alpha * k(t);
    cdf += w * R::pgamma(x, shape, scale, 1, 0);
    pdf += w * R::dgamma(x, shape, scale, 0);
  }
  return arma::vec({std::min(std::max(cdf, 0.0), 1.0), pdf});
}


// Find the root by bisection when Newton-Raphson does not converge
static arma::vec solver_bisection(const arma::vec& k, const arma::vec& logw,
                                   const arma::uvec& term_start, const arma::vec& base_cdf,
                                   const arma::uvec& idx,
                                   const arma::vec& q, const arma::vec& alpha,
                                   const arma::vec& beta, const arma::vec& x0) {
  int n = q.n_elem;
  const double tol      = 1e-8;
  const int    max_iter = 100;

  arma::vec f_L = eval_cdf_vec(k, logw, term_start, base_cdf, idx, x0, alpha, beta) - q;
  arma::vec f_U = f_L;
  arma::vec x_L = x0;
  arma::vec x_U = x0;

  // Expand the upper bound until it exceeds q everywhere
  arma::uvec neg_idx = arma::find(f_U < 0);
  while (!neg_idx.is_empty()) {
    x_L.elem(neg_idx) = x_U.elem(neg_idx);
    f_L.elem(neg_idx) = f_U.elem(neg_idx);
    x_U.elem(neg_idx) *= 2;
    f_U.elem(neg_idx) = eval_cdf_vec(
        k, logw, term_start, base_cdf, idx.elem(neg_idx), x_U.elem(neg_idx),
        alpha.elem(neg_idx), beta.elem(neg_idx)) - q.elem(neg_idx);
    neg_idx = arma::find(f_U < 0);
  }

  // Expand the lower bound until it is below q everywhere
  arma::uvec pos_idx = arma::find(f_L > 0);
  while (!pos_idx.is_empty()) {
    x_U.elem(pos_idx) = x_L.elem(pos_idx);
    f_U.elem(pos_idx) = f_L.elem(pos_idx);
    x_L.elem(pos_idx) /= 2;
    f_L.elem(pos_idx) = eval_cdf_vec(
        k, logw, term_start, base_cdf, idx.elem(pos_idx), x_L.elem(pos_idx),
        alpha.elem(pos_idx), beta.elem(pos_idx)) - q.elem(pos_idx);
    pos_idx = arma::find(f_L > 0);
  }

  // Bisect the bracket until convergence
  arma::uvec active = arma::regspace<arma::uvec>(0, n - 1);
  arma::vec x_mid(n, arma::fill::zeros);
  for (int iter = 0; iter < max_iter; ++iter) {
    arma::vec x_mid_a = 0.5 * (x_L.elem(active) + x_U.elem(active));
    x_mid.elem(active) = x_mid_a;
    arma::vec f_mid = eval_cdf_vec(
        k, logw, term_start, base_cdf, idx.elem(active), x_mid_a,
        alpha.elem(active), beta.elem(active)) - q.elem(active);

    arma::vec bracket = x_U.elem(active) - x_L.elem(active);
    arma::uword na = active.n_elem;
    arma::uvec conv_mask(na);
    for (arma::uword i = 0; i < na; ++i) {
      conv_mask(i) = (std::abs(f_mid(i)) < tol) ||
                    (bracket(i) < tol * (1.0 + x_mid_a(i)));
    }

    arma::uvec mid_neg = arma::find(f_mid < 0);
    if (!mid_neg.is_empty()) x_L.elem(active.elem(mid_neg)) = x_mid_a.elem(mid_neg);
    arma::uvec mid_pos = arma::find(f_mid > 0);
    if (!mid_pos.is_empty()) x_U.elem(active.elem(mid_pos)) = x_mid_a.elem(mid_pos);

    arma::uvec converged = arma::find(conv_mask);
    if (!converged.is_empty()) {
      active = active.elem(arma::find(conv_mask == 0));
      if (active.is_empty()) break;
    }
  }

  // Mark the points that did not converge as missing
  arma::vec result = x_mid;
  if (!active.is_empty()) result.elem(active).fill(R_NaReal);
  return result;
}


// Find the root by Newton-Raphson, falling back to bisection on failure
static arma::vec solver_Newton_Raphson(const arma::vec& q,
                                        const arma::vec& lambda,
                                        const arma::vec& alpha,
                                        const arma::vec& beta,
                                        const arma::vec& mean,
                                        const arma::vec& dispersion,
                                        const arma::vec& power) {
  int n = q.n_elem;
  const double tol      = 1e-8;
  const int    max_iter = 20;

  // Compute an initial guess for x
  arma::vec x0 = mean + 2 * arma::sqrt(dispersion % arma::pow(mean, power)) % (2 * q - 1);
  arma::uvec neg0 = arma::find(x0 <= 0);
  if (!neg0.is_empty()) x0.elem(neg0) = mean.elem(neg0);
  x0 = arma::clamp(x0, tol, 1 / tol);
  arma::vec x = x0;

  arma::uvec active = arma::regspace<arma::uvec>(0, n - 1);
  arma::vec converged(n, arma::fill::zeros);

  // Build the cached Poisson terms
  arma::vec  k, logw, base_cdf;
  arma::uvec term_start;
  build_poisson_terms(lambda, k, logw, term_start, base_cdf);

  // Iterate the Newton-Raphson update
  for (int iter = 0; iter < max_iter && !active.is_empty(); ++iter) {
    arma::vec f(active.n_elem);
    arma::vec f_prime(active.n_elem);

    for (arma::uword j = 0; j < active.n_elem; ++j) {
      arma::uword i = active(j);
      arma::vec cdf_pdf = eval_cdf_pdf(k, logw, term_start, base_cdf, i, x(i), alpha(i), beta(i));
      f(j)       = cdf_pdf(0) - q(i);
      f_prime(j) = cdf_pdf(1) == 0 ? tol : cdf_pdf(1);
    }

    x.elem(active) -= f / f_prime;
    x.elem(active)  = arma::clamp(x.elem(active), tol, 1 / std::sqrt(tol));

    arma::uvec newly_conv = arma::find(arma::abs(f) < tol);
    if (!newly_conv.is_empty()) {
      converged.elem(active.elem(newly_conv)).ones();
      arma::uvec still_active = arma::find(converged.elem(active) == 0);
      active = active.elem(still_active);
    }
  }

  // Fall back to bisection for points that did not converge
  if (!active.is_empty()) {
    x.elem(active) = solver_bisection(
        k, logw, term_start, base_cdf, active, q.elem(active), alpha.elem(active),
        beta.elem(active), x0.elem(active));
  }
  return x;
}


// Compute the quantile function
// [[Rcpp::export]]
arma::vec tweedieInvCDF(arma::vec q, arma::vec mean, arma::vec dispersion,
                         arma::vec power) {

  // Reshape all input vectors to the same length
  int l = std::max({static_cast<int>(q.n_elem), static_cast<int>(mean.n_elem),
                    static_cast<int>(dispersion.n_elem),
                    static_cast<int>(power.n_elem)});
  q          = recycle_to_length(q,          l, "q");
  mean       = recycle_to_length(mean,       l, "mean");
  dispersion = recycle_to_length(dispersion, l, "dispersion");
  power      = recycle_to_length(power,      l, "power");

  // Change to the compound Gamma-Poisson parametrization
  arma::vec lambda = cpg_lambda(mean, dispersion, power);
  arma::vec alpha  = cpg_alpha(power);
  arma::vec beta   = cpg_beta(mean, dispersion, power);
  arma::vec invcdf(l, arma::fill::zeros);

  // Solve for x where q exceeds the point mass at 0
  arma::uvec pos_idx = arma::find(q > arma::exp(-lambda));
  if (!pos_idx.is_empty()) {
    invcdf.elem(pos_idx) = solver_Newton_Raphson(
        q(pos_idx), lambda(pos_idx), alpha(pos_idx), beta(pos_idx),
        mean(pos_idx), dispersion(pos_idx), power(pos_idx));
  }

  // Return the values
  return invcdf;
}
