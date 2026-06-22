#pragma once
#include <RcppArmadillo.h>

// Recycle x to length l; errors if x is neither length 1 nor length l.
arma::vec recycle_to_length(const arma::vec& x, int l, const std::string& name);

// Vectorised wrappers around R's scalar dpois / pgamma.
arma::vec dpois_vec(const arma::vec& k, const arma::vec& lambda, bool log_p = false);
arma::vec pgamma_vec(const arma::vec& x, const arma::vec& shape, const arma::vec& beta);

// Compound Poisson-Gamma CDF; shared by tweedieCDF and tweedieInvCDF.
arma::vec compound_Poisson_Gamma(const arma::vec& x,
                                  const arma::vec& lambda,
                                  const arma::vec& alpha,
                                  const arma::vec& beta);
