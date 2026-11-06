# math_primitives.R -- Numerically stable elementary functions
#
# PURPOSE
# -------
# Floating-point naive implementations of log(1+exp(x)), log(1-exp(-x)), and
# probability clamping produce incorrect results near the boundaries of
# double precision.  These primitives use piece-wise stable formulas that are
# safe across the full range of inputs encountered in survival analysis
# (very small or very large hazard rates, near-zero probabilities).
#
# All three functions are exported so users can call them directly when
# constructing custom likelihoods or diagnostics.

#' Numerically stable log(1 + exp(x))
#'
#' Computes the "softplus" function
#' \deqn{\mathrm{log1pexp}(x) = \log\!\bigl(1 + e^{x}\bigr)}
#' without overflow.  A naive `log(1 + exp(x))` overflows to `Inf` for
#' \eqn{x \gtrsim 710} even though the true value is \eqn{\approx x}.  The
#' identity \eqn{\log(1 + e^{x}) = x + \log(1 + e^{-x})} for \eqn{x > 0} keeps
#' every intermediate finite.
#'
#' @param x Numeric vector.
#' @return Numeric vector with element-wise \eqn{\log(1 + e^{x})}.
#' @references
#' \enc{Mächler}{Machler} M (2012). *Accurately Computing \eqn{\log(1 - e^{-|a|})}.* R package
#' Rmpfr vignette.
#' \url{https://CRAN.R-project.org/package=Rmpfr/vignettes/log1mexp-note.pdf}
#' @seealso [hzr_log1mexp()] for the complementary \eqn{\log(1 - e^{-x})},
#'   [hzr_clamp_prob()] for boundary-safe probabilities.
#' @examples
#' hzr_log1pexp(c(-50, 0, 50))
#' @export
hzr_log1pexp <- function(x) {
  # Anything but a plain numeric vector, matrix or array takes ifelse()
  # itself: a class's own `[` or `>` (a Matrix, a data frame) can make the
  # form below differ from it.
  if (is.object(x) || !is.numeric(x)) {
    return(ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x))))
  }
  # For the rest, ifelse() written out as it runs, but forming each branch
  # only where it is taken: ifelse() formed both over every element, and
  # this is among the most called functions of a multiphase fit (#631). The
  # result is ifelse()'s: the attributes of x > 0, and NA where x is NA.
  ans <- x > 0
  pos <- which(ans)
  neg <- which(!ans)
  if (length(pos)) ans[pos] <- x[pos] + log1p(exp(-x[pos]))
  if (length(neg)) ans[neg] <- log1p(exp(x[neg]))
  ans
}

#' Numerically stable log(1 - exp(-x)) for x > 0
#'
#' Computes
#' \deqn{\mathrm{log1mexp}(x) = \log\!\bigl(1 - e^{-x}\bigr), \qquad x > 0}
#' accurately across the full range.  Following \enc{Mächler}{Machler} (2012), the evaluation
#' switches at \eqn{x = \log 2}: for small \eqn{x} it uses
#' \eqn{\log(-\mathrm{expm1}(-x))} (avoiding cancellation as \eqn{x \to 0}), and
#' for larger \eqn{x} it uses \eqn{\mathrm{log1p}(-e^{-x})} (avoiding loss of
#' precision as \eqn{x \to \infty}).  Non-positive or non-finite inputs return
#' `NA`.
#'
#' @param x Numeric vector with positive values.
#' @return Numeric vector with element-wise \eqn{\log(1 - e^{-x})}.
#' @references
#' \enc{Mächler}{Machler} M (2012). *Accurately Computing \eqn{\log(1 - e^{-|a|})}.* R package
#' Rmpfr vignette.
#' \url{https://CRAN.R-project.org/package=Rmpfr/vignettes/log1mexp-note.pdf}
#' @seealso [hzr_log1pexp()] for the softplus \eqn{\log(1 + e^{x})},
#'   [hzr_clamp_prob()] for boundary-safe probabilities.
#' @examples
#' hzr_log1mexp(c(0.01, 0.5, 5))
#' @export
hzr_log1mexp <- function(x) {
  out <- rep(NA_real_, length(x))
  cutoff <- log(2)

  valid <- is.finite(x) & (x > 0)
  small <- valid & (x <= cutoff)
  large <- valid & (x > cutoff)

  out[small] <- log(-expm1(-x[small]))
  out[large] <- log1p(-exp(-x[large]))
  out
}

#' Clamp probabilities away from 0 and 1
#'
#' Bounds a probability vector to \eqn{[\epsilon,\, 1 - \epsilon]} so that
#' downstream \eqn{\log p} or \eqn{\log(1 - p)} evaluations stay finite.  Used
#' throughout the likelihood and prediction code to guard survival and CDF
#' values against exact 0 or 1.
#'
#' @param p Numeric vector of probabilities.
#' @param eps Small positive tolerance in \eqn{(0, 0.5)}; default `1e-12`.
#' @return Numeric vector bounded to `[eps, 1 - eps]`.
#' @seealso [hzr_log1pexp()] and [hzr_log1mexp()] for the companion
#'   stable-logarithm primitives.
#' @examples
#' hzr_clamp_prob(c(0, 0.5, 1))
#' @export
hzr_clamp_prob <- function(p, eps = 1e-12) {
  if (!is.numeric(eps) || length(eps) != 1 || !is.finite(eps) || eps <= 0 || eps >= 0.5) {
    stop("'eps' must be a finite scalar in (0, 0.5).", call. = FALSE)
  }

  pmin(1 - eps, pmax(eps, p))
}
