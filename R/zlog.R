#' Calculate standardized logarithmic values
#'
#' Converts positive laboratory measurements into standardized logarithmic
#' values using lower and upper reference limits.
#'
#' @param x A numeric vector containing positive laboratory measurements.
#' @param limits A numeric vector of length two giving positive lower and
#'   upper reference limits. The upper limit must be greater than the lower.
#'
#' @details
#' The scaling factor is calculated as `diff(stats::qnorm(c(0.025, 0.975)))`.
#'
#' @return A numeric vector containing the calculated zlog values, with the
#'   same length and names as `x`.
#'   Invalid measurements are returned as `NA`.
#'
#' @examples
#' zlog(4, limits = c(2, 8))
#' zlog(c(2, 4, 8), limits = c(2, 8))
#'
#' @export
zlog <- function(x, limits) {
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }

  if (!is.numeric(limits) || is.complex(limits) ||
      !is.null(dim(limits)) || length(limits) != 2L) {
    stop("limits must be a numeric vector of length 2.")
  }

  L <- limits[1L]
  U <- limits[2L]
  result <- rep(NA_real_, length(x))
  names(result) <- names(x)

  if (
    !is.finite(L) ||
    !is.finite(U) ||
    L <= 0 ||
    U <= 0 ||
    U <= L
  ) {
    return(result)
  }

  valid <- is.finite(x) & x > 0

  logl <- log(L)
  logu <- log(U)
  mu.log <- (logl + logu) / 2
  sigma.log <- (logu - logl) / diff(stats::qnorm(c(0.025, 0.975)))

  result[valid] <- (log(x[valid]) - mu.log) / sigma.log

  result
}

#' Invert standardized logarithmic values
#'
#' Converts standardized logarithmic scores back to laboratory measurements
#' using the same reference limits as [zlog()].
#'
#' @param x A numeric vector of finite standardized scores. Negative scores
#'   and zero are allowed.
#' @param limits A numeric vector of length two giving positive lower and
#'   upper reference limits. The upper limit must be greater than the lower.
#'
#' @details
#' The inverse is `exp(x * sigma.log + mu.log)`, where
#' `mu.log = (log(L) + log(U)) / 2` and
#' `sigma.log = (log(U) - log(L)) / diff(stats::qnorm(c(0.025, 0.975)))`.
#'
#' @return A numeric vector of recovered measurements with the same length
#'   and names as `x`. Missing or non-finite scores and invalid reference
#'   limits produce `NA`.
#' @export
#' @examples
#' izlog(c(-1.959964, 0, 1.959964), limits = c(2, 8))
#' values <- c(2, 4, 8, 12)
#' izlog(zlog(values, limits = c(2, 8)), limits = c(2, 8))
izlog <- function(x, limits) {
  if (!is.numeric(x) || is.complex(x) || !is.null(dim(x))) {
    stop("x must be a numeric vector.")
  }
  if (!is.numeric(limits) || is.complex(limits) ||
      !is.null(dim(limits)) || length(limits) != 2L) {
    stop("limits must be a numeric vector of length 2.")
  }
  
  result <- rep(NA_real_, length(x))
  names(result) <- names(x)
  L <- limits[1L]
  U <- limits[2L]
  if (!is.finite(L) || !is.finite(U) || L <= 0 || U <= L) {
    return(result)
  }
  
  mu.log <- (log(L) + log(U)) / 2
  sigma.log <- (log(U) - log(L)) / diff(stats::qnorm(c(0.025, 0.975)))
  valid <- is.finite(x)
  result[valid] <- exp(x[valid] * sigma.log + mu.log)
  result
}