#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib tweedieDistr, .registration = TRUE
## usethis namespace: end
NULL

.onLoad <- function(libname, pkgname) {
  register <- function(generic, class, fun, pkg = "distributional") {
    if (requireNamespace(pkg, quietly = TRUE)) {
      registerS3method(generic, class, fun, envir = asNamespace(pkg))
    }
  }
  register("skewness", "dist_tweedie", skewness.dist_tweedie)
  register("kurtosis", "dist_tweedie", kurtosis.dist_tweedie)
}
