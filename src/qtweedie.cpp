#include "utils.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Defined in tweedie_density.cpp — needed by solver_Newton_Raphson.
arma::vec tweedieDensity(arma::vec x, arma::vec mean, arma::vec dispersion,
                         arma::vec power, bool log);


// Find the root by bisection when Newton-Raphson does not converge
static arma::vec solver_bisection(const arma::vec& q, const arma::vec& lambda,
                                   const arma::vec& alpha, const arma::vec& beta,
                                   const arma::vec& x0) {
  int n = q.n_elem;
  const double tol      = 1e-8;
  const int    max_iter = 100;
  arma::vec f_L = compound_Poisson_Gamma(x0, lambda, alpha, beta) - q;
  arma::vec f_U = f_L;
  arma::vec x_L = x0;
  arma::vec x_U = x0;

  // Expand the upper bound until it exceeds q everywhere
  arma::uvec neg_idx = arma::find(f_U < 0);
  while (!neg_idx.is_empty()) {
    x_L.elem(neg_idx) = x_U.elem(neg_idx);
    f_L.elem(neg_idx) = f_U.elem(neg_idx);
    x_U.elem(neg_idx) *= 2;
    f_U.elem(neg_idx) = compound_Poisson_Gamma(
        x_U.elem(neg_idx), lambda.elem(neg_idx),
        alpha.elem(neg_idx), beta.elem(neg_idx)) - q.elem(neg_idx);
    neg_idx = arma::find(f_U < 0);
  }

  // Expand the lower bound until it is below q everywhere
  arma::uvec pos_idx = arma::find(f_L > 0);
  while (!pos_idx.is_empty()) {
    x_U.elem(pos_idx) = x_L.elem(pos_idx);
    f_U.elem(pos_idx) = f_L.elem(pos_idx);
    x_L.elem(pos_idx) /= 2;
    f_L.elem(pos_idx) = compound_Poisson_Gamma(
        x_L.elem(pos_idx), lambda.elem(pos_idx),
        alpha.elem(pos_idx), beta.elem(pos_idx)) - q.elem(pos_idx);
    pos_idx = arma::find(f_L > 0);
  }

  // Bisect the bracket until convergence
  arma::uvec active = arma::regspace<arma::uvec>(0, n - 1);
  arma::vec x_mid(n, arma::fill::zeros);
  for (int iter = 0; iter < max_iter; ++iter) {
    arma::vec x_mid_a = 0.5 * (x_L.elem(active) + x_U.elem(active));
    x_mid.elem(active) = x_mid_a;
    arma::vec f_mid = compound_Poisson_Gamma(
        x_mid_a, lambda.elem(active),
        alpha.elem(active), beta.elem(active)) - q.elem(active);

    // Check for convergence
    arma::vec bracket = x_U.elem(active) - x_L.elem(active);
    arma::uword na = active.n_elem;
    arma::uvec conv_mask(na);
    for (arma::uword i = 0; i < na; ++i) {
      conv_mask(i) = (std::abs(f_mid(i)) < tol) ||
                    (bracket(i) < tol * (1.0 + x_mid_a(i)));
    }

    // Update the bracket for the next iteration
    arma::uvec mid_neg = arma::find(f_mid < 0);
    if (!mid_neg.is_empty()) x_L.elem(active.elem(mid_neg)) = x_mid_a.elem(mid_neg);
    arma::uvec mid_pos = arma::find(f_mid > 0);
    if (!mid_pos.is_empty()) x_U.elem(active.elem(mid_pos)) = x_mid_a.elem(mid_pos);

    // Update the active set too
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

  // Compute an initial guess
  arma::vec x0 = mean + 2 * arma::sqrt(dispersion % arma::pow(mean, power)) % (2 * q - 1);
  arma::uvec neg0 = arma::find(x0 <= 0);
  if (!neg0.is_empty()) x0.elem(neg0) = mean.elem(neg0);
  x0 = arma::clamp(x0, tol, 1 / tol);
  arma::vec x = x0;
  arma::uvec active = arma::regspace<arma::uvec>(0, n - 1);
  arma::vec converged(n, arma::fill::zeros);

  // Iterate the Newton-Raphson update
  for (int iter = 0; iter < max_iter && !active.is_empty(); ++iter) {
    arma::vec x_a = x.elem(active);
    arma::vec cdf_x = compound_Poisson_Gamma(
        x_a, lambda.elem(active), alpha.elem(active), beta.elem(active));
    arma::vec pdf_x = tweedieDensity(
        x_a, mean.elem(active), dispersion.elem(active), power.elem(active), false);
    arma::vec f       = cdf_x - q.elem(active);
    arma::vec f_prime = pdf_x;
    f_prime.elem(arma::find(f_prime == 0)).fill(tol);
    x.elem(active) -= f / f_prime;
    x.elem(active)  = arma::clamp(x.elem(active), tol, 1 / std::sqrt(tol));

    // Check for convergence and udpate the active set
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
        q.elem(active), lambda.elem(active), alpha.elem(active),
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
  arma::vec lambda = arma::pow(mean, 2 - power) / (dispersion % (2 - power));
  arma::vec alpha  = (2 - power) / (power - 1);
  arma::vec beta   = 1 / (dispersion % (power - 1) % arma::pow(mean, power - 1));
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
