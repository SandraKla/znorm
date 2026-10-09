#' Box-Cox Standardization of Laboratory Measurements
#'
#' Transform positive measurements into standardized scores relative to lower
#' and upper reference limits, or transform scores back to measurements.
#'
#' @param x A numeric vector. For `znorm()`, positive, finite measurements;
#'   for `iznorm()`, finite standardized scores. Missing values are allowed.
#' @param limits A numeric vector of length two giving positive lower and
#'   upper reference limits, or a numeric matrix with two columns and either
#'   one row or `length(x)` rows. Non-missing limits must be finite and each
#'   lower limit must be smaller than its corresponding upper limit. A missing
#'   limit produces a missing result for that measurement.
#' @param lambda A numeric vector of length one or `length(x)` containing
#'   finite Box-Cox power parameters. Missing parameters produce missing
#'   results. Zero selects a logarithmic transformation and one a linear
#'   transformation. The default is zero.
#' @param probs Two distinct, finite probabilities strictly between zero and
#'   one. They are sorted before use. The default `c(0.025, 0.975)` describes
#'   a central 95 percent reference interval.
#'
#' @details
#' The Box-Cox transformation is
#' \deqn{g(x, \lambda) = (x^\lambda - 1) / \lambda}
#' for nonzero \eqn{\lambda}, and \eqn{g(x, 0) = \log(x)} otherwise.
#' For lower and upper reference limits \eqn{L} and \eqn{U}, define
#' \eqn{a = g(L, \lambda)}, \eqn{b = g(U, \lambda)}, and
#' \eqn{d = \Phi^{-1}(p_2) - \Phi^{-1}(p_1)}, where \eqn{p_1 < p_2}
#' are the sorted `probs` and \eqn{\Phi^{-1}} is the standard normal
#' quantile function. The score is
#' \deqn{z = \frac{g(x, \lambda) - (a + b) / 2}{b - a} d.}
#' `iznorm()` applies the inverse Box-Cox transformation to
#' \eqn{z (b - a) / d + (a + b) / 2}.
#'
#' The transformed midpoint always maps to zero. The reference limits map
#' to \eqn{-d/2} and \eqn{d/2}. For complementary probabilities, such as
#' the defaults, these equal the requested normal quantiles. For asymmetric
#' probabilities the midpoint convention is retained; the limits do not
#' generally map to the individual normal quantiles.
#'
#' For complementary probabilities, `lambda = 0` and `lambda = 1` agree
#' with the zlog and z standardization. Positive inputs are required even for
#' `lambda = 1`, so the domain is consistent across power parameters.
#' The parameter is supplied by the user; it is not estimated from `x`.
#' The transformation alone does not guarantee normally distributed data.
#'
#' For nonzero powers, not every score has a positive inverse. If
#' \eqn{1 + \lambda g} is nonpositive, `iznorm()` returns `NaN` for that
#' element and issues a warning. Extreme powers or reference limits that
#' cannot be distinguished on the transformed scale produce an error.
#' As with other numeric transformations, extreme measurements can overflow
#' or underflow in floating-point arithmetic.
#'
#' @return A numeric vector of the same length and with the same names as
#'   `x`. `znorm()` returns standardized scores; `iznorm()` returns
#'   measurements. Missing inputs propagate to the corresponding results.
#'
#' @references
#' Hoffmann G, Klawonn F, Lichtinghagen R, Orth M (2017).
#' "The Zlog-Value as Basis for the Standardization of Laboratory Results."
#' *LaboratoriumsMedizin*, 41(1), 23-32.
#' \doi{10.1515/labmed-2016-0087}.
#'
#' Box GEP, Cox DR (1964). "An Analysis of Transformations."
#' *Journal of the Royal Statistical Society, Series B*, 26(2), 211-252.
#' \doi{10.1111/j.2517-6161.1964.tb00553.x}.
#'
#' @export
#' @examples
#' # Logarithmic, linear, and square-root-like transformations
#' znorm(1:10, limits = c(2, 8), lambda = 0)
#' znorm(1:10, limits = c(2, 8), lambda = 1)
#' albumin <- c(42, 34, 38, 43, 50, 42, 27, 31, 24)
#' scores <- znorm(albumin, limits = c(35, 52), lambda = 0.5)
#' iznorm(scores, limits = c(35, 52), lambda = 0.5)
#'
#' # A different reference interval and power for each measurement
#' x <- c(albumin = 42, bilirubin = 8)
#' limits <- rbind(c(35, 52), c(2, 21))
#' znorm(x, limits = limits, lambda = c(1, 0.33))
#'
#' # The default reference limits map to normal quantiles
#' znorm(c(2, 8), limits = c(2, 8), lambda = 0.5)
#' qnorm(c(0.025, 0.975))
znorm <- function(x, limits, lambda = 0, probs = c(0.025, 0.975)) {
  if (!is.numeric(x))
    stop("'x' has to be numeric.")
  if (any(x <= 0, na.rm = TRUE))
    stop("'x' has to be positive for a Box-Cox transformation.")
  if (!is.numeric(probs) || length(probs) != 2L)
    stop("'probs' has to be a numeric of length 2.")
  if (!is.numeric(lambda) || !(length(lambda) %in% c(1L, length(x))))
    stop("'lambda' has to be numeric of length 1 or length(x).")
  
  if (is.null(dim(limits))) {
    if (length(limits) != 2L)
      stop("'limits' has to be of length 2 or a two-column matrix.")
    limits <- matrix(limits, nrow = 1L, ncol = 2L)
  }
  if (ncol(limits) != 2L)
    stop("'limits' has to have two columns (lower and upper limit).")
  if (!(nrow(limits) %in% c(1L, length(x))))
    stop("'limits' has to have a single row or as many rows as elements in 'x'.")
  
  lower <- limits[, 1L]
  upper <- limits[, 2L]
  if (any(lower <= 0 | upper <= 0, na.rm = TRUE))
    stop("'limits' have to be positive for a Box-Cox transformation.")
  if (any(lower >= upper, na.rm = TRUE))
    stop("lower limit has to be smaller than the upper limit.")
  
  spread <- diff(stats::qnorm(sort(probs)))
  
  gx <- .boxcox(x, lambda)
  gl <- .boxcox(lower, lambda)
  gu <- .boxcox(upper, lambda)
  
  (gx - (gl + gu) / 2) * spread / (gu - gl)
}

#' Inverse Box-Cox Transformation: Calculate Laboratory Values from znorm Values
#'
#' Inverse function of \code{znorm}, analogous to \code{iz}/\code{izlog}.
#'
#' @param x `numeric`, znorm values.
#' @param limits `numeric` vector or matrix, lower and upper reference limits.
#' @param lambda `numeric`, Box-Cox parameter, see \code{znorm}.
#' @param probs `numeric`, see \code{znorm}.
#'
#' @return `numeric`, laboratory values.
#' @export
#'
#' @examples
#' x <- znorm(1:10, limits = c(2, 8), lambda = 0.5)
iznorm <- function(x, limits, lambda = 0, probs = c(0.025, 0.975)) {
  if (is.null(dim(limits))) {
    if (length(limits) != 2L)
      stop("'limits' has to be of length 2 or a two-column matrix.")
    limits <- matrix(limits, nrow = 1L, ncol = 2L)
  }
  lower <- limits[, 1L]
  upper <- limits[, 2L]
  
  spread <- diff(stats::qnorm(sort(probs)))
  
  gl <- .boxcox(lower, lambda)
  gu <- .boxcox(upper, lambda)
  
  g <- x * (gu - gl) / spread + (gl + gu) / 2
  
  ifelse(lambda == 0, exp(g), (g * lambda + 1) ^ (1 / lambda))
}

.numeric_vector <- function(x) {
  is.numeric(x) && !is.complex(x) && is.null(dim(x))
}

#' Box-Cox-Transformation
#'
#' @param x `numeric` positive values
#' @param lambda `numeric`, Box-Cox-Parameter
#' @return `numeric`
#' @noRd
.boxcox <- function(x, lambda) {
  ifelse(lambda == 0, log(x), (x ^ lambda - 1) / lambda)
}