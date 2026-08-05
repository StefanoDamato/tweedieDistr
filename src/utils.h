#pragma once
#include <RcppArmadillo.h>

// Set the tolerance for the truncation of infinite sums
static const double LOG_TOL = 37.0;

// Recycle a vector to match a given length
arma::vec recycle_to_length(const arma::vec& x, int l, const std::string& name);

// The compound Poisson-Gamma parameters
arma::vec cpg_lambda(const arma::vec& mean, const arma::vec& dispersion, const arma::vec& power);
arma::vec cpg_alpha(const arma::vec& power);
arma::vec cpg_beta(const arma::vec& mean, const arma::vec& dispersion, const arma::vec& power);

// Compute the CDF for a batch of points via a single pass per point; shared
// by tweedieCDF, which evaluates each point once and has no cache to reuse
arma::vec eval_cdf_direct(const arma::vec& x, const arma::vec& lambda,
                          const arma::vec& alpha, const arma::vec& beta);
