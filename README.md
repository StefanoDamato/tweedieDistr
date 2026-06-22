
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tweedieDistr <a href="https://github.com/StefanoDamato/tweedieDistr/"><img src="man/figures/logo.png" align="right" height="150" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/StefanoDamato/tweedieDistr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/StefanoDamato/tweedieDistr/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/tweedieDistr)](https://CRAN.R-project.org/package=tweedieDistr)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

`tweedieDistr` provides density, distribution function, quantile
function, and random generation for the **Tweedie distribution** under
the compound Poisson–Gamma parameterisation with power parameter
$p \in (1, 2)$. A
[`distributional`](https://github.com/mitchelloharawild/distributional)-compatible
constructor is also provided for use in tidy modelling workflows.

The Tweedie family naturally combines a point mass at zero with a
continuous positive component, making it well suited to **intermittent
demand** data and any setting where exact zeros occur alongside strictly
positive observations.

## Exported functions

| Function | Description |
|----|----|
| `dtweedie()` | Probability density function via the series expansion of Dunn & Smyth (2005), implemented in C++ for performance. |
| `ptweedie()` | Cumulative distribution function using a truncated compound Poisson–Gamma summation. |
| `qtweedie()` | Quantile function via a hybrid Newton–Raphson / bisection root-finding algorithm. |
| `rtweedie()` | Random generation via the exact compound Poisson–Gamma representation. |
| `dist_tweedie()` | [`distributional`](https://github.com/mitchelloharawild/distributional) object with full S3 method support. |

The `dist_tweedie()` object supports the standard `distributional`
interface: `density()`, `CDF()`, `quantile()`, `generate()`, `mean()`,
`variance()`, `distributional::skewness()`,
`distributional::kurtosis()`, and `distributional::support()`.

## Installation

You can install the **stable** version from CRAN:

``` r
install.packages("tweedieDistr")
```

You can install the **development** version from
[GitHub](https://github.com/StefanoDamato/tweedieDistr):

``` r
# install.packages("devtools")
devtools::install_github("StefanoDamato/tweedieDistr")
```

## Usage

### Standard `stats`-style interface

The four `dtweedie` / `ptweedie` / `qtweedie` / `rtweedie` functions
mirror the conventions of base R distribution functions and support
vectorised arguments.

``` r
library(tweedieDistr)

# density at a few points
dtweedie(c(0, 1, 2, 3), mean = 2, dispersion = 0.5, power = 1.5)
#> [1] 0.003493489 0.330271685 0.324015308 0.180836595

# cumulative probabilities
ptweedie(c(0, 1, 2, 3), mean = 2, dispersion = 0.5, power = 1.5)
#> [1] 0.003493489 0.210260258 0.559994415 0.813739540

# quantiles
qtweedie(c(0.25, 0.5, 0.75, 0.9), mean = 2, dispersion = 0.5, power = 1.5)
#> [1] 1.117522 1.820328 2.686368 3.604913

# random samples
set.seed(42)
rtweedie(8, mean = 2, dispersion = 0.5, power = 1.5)
#> [1] 3.4363596 2.8968086 0.8407508 4.9660779 1.8928915 1.1482404 1.2166902
#> [8] 1.7718016
```

### `distributional` interface

`dist_tweedie()` creates a fully-featured `distributional` object,
compatible with tidy modelling frameworks such as
[`fable`](https://fable.tidyverts.org/).

``` r
library(tweedieDistr)

d <- dist_tweedie(mean = 2, dispersion = 0.5, power = 1.5)
d
#> <distribution[1]>
#> [1] Tweedie(2, 0.5, 1.5)

# moments
mean(d)
#> [1] 2
distributional::variance(d)
#> [1] 1.414214
distributional::skewness(d)
#> [1] 0.8919053
distributional::kurtosis(d)
#> [1] 1.06066

# density and CDF
density(d, at = c(0, 1, 2, 3))
#> [[1]]
#> [1] 0.003493489 0.330271685 0.324015308 0.180836595
distributional::cdf(d, q = c(1, 2, 3))
#> [[1]]
#> [1] 0.2102603 0.5599944 0.8137395

# quantiles
quantile(d, p = c(0.5, 0.9))
#> [[1]]
#> [1] 1.820328 3.604913

# random generation
distributional::generate(d, times = 6)
#> [[1]]
#> [1] 3.6691746 1.6225007 2.3630598 3.3738106 4.1138676 0.2777309
```

### Vectorised distributions

`dist_tweedie()` supports vectorised arguments, so a whole column of
distribution objects can be constructed at once — useful when each row
of a forecast table has its own parameters.

``` r
library(tweedieDistr)

d_vec <- dist_tweedie(
  mean       = c(1, 2, 5),
  dispersion = c(0.5, 1, 0.8),
  power      = c(1.3, 1.5, 1.8)
)
d_vec
#> <distribution[3]>
#> [1] Tweedie(1, 0.5, 1.3) Tweedie(2, 1, 1.5)   Tweedie(5, 0.8, 1.8)
mean(d_vec)
#> [1] 1 2 5
quantile(d_vec, p = 0.9)
#> [1]  1.959805  4.290072 10.124295
```

## Mathematical background

The Tweedie distribution with power $p \in (1, 2)$ is a compound
Poisson–Gamma variable: $X = \sum_{i=1}^{N} G_i$, where

$$N \sim \mathrm{Poisson}(\lambda), \qquad
G_i \sim \mathrm{Gamma}(\alpha, \beta),$$

with

$$\lambda = \frac{\mu^{2-p}}{\phi(2-p)}, \quad
\alpha = \frac{2-p}{p-1}, \quad
\beta = \frac{1}{\phi(p-1)\mu^{p-1}}.$$

The distribution has mean $\mu$ and variance $\phi\mu^p$. When $N = 0$
the sum is conventionally defined as zero, giving a point mass
$P(X = 0) =
e^{-\lambda}$. The density for $x > 0$ is evaluated via the series
expansion of Dunn & Smyth (2005), implemented in C++ through
[Rcpp](https://www.rcpp.org/) and
[RcppArmadillo](https://github.com/RcppCore/RcppArmadillo).

## Contributors

<!-- prettier-ignore-start -->

<!-- markdownlint-disable -->

<table>

<tbody>

<tr>

<td align="center" valign="top" width="20%">

<a href="#">
<img src="https://github.com/StefanoDamato.png" width="100px;" alt="Stefano Damato" style="border-radius:50%;border:1px solid #646464;"/><br />
<sub><b>Stefano Damato</b></sub></a><br /> <sub>(Maintainer)</sub><br />
<a href="mailto:stefano.damato@supsi.ch?subject=[tweedieDistr package]">Email</a>
</td>

</tr>

</tbody>

</table>

<!-- markdownlint-restore -->

<!-- prettier-ignore-end -->

## Getting help

If you encounter a bug, please file a minimal reproducible example on
[GitHub](https://github.com/StefanoDamato/tweedieDistr/issues).

## References

Dunn, P. K., & Smyth, G. K. (2005). Series evaluation of Tweedie
exponential dispersion model densities. *Statistics and Computing*,
15(4), 267–280. <https://doi.org/10.1007/s11222-005-4070-y>.
