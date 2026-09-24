# boundary-detect.R -- is a fitted phase indistinguishable from a step? (#448)
#
# The G1 decomposition already refuses `nu == 0` with `m >= 0` outright: that
# limit "is degenerate (G collapses to a step function), so there is no usable
# form" (decomposition.R).  But a fit can settle at nu = -1.4e-16, take the
# ordinary `m > 0 && nu < 0` branch, and return converged = TRUE.  The guard
# sits at the exact point; the pathology is a neighbourhood of it.
#
# The test below deliberately does NOT threshold `nu`.  It asks a question of
# the DATA: does this phase complete its entire rise, from numerically 0 to
# numerically 1, inside a single gap between adjacent observed times?  If so
# the phase is a step AT THIS DATA'S RESOLUTION -- a fact about this fit and
# these times, not a tuning choice.  The magnitude is reported so the reader
# judges severity; nothing is suppressed on its behalf.

# Numerically-0 / numerically-1 scale.  sqrt(.Machine$double.eps) is the
# conventional "these are not distinguishable" scale, not a value tuned to any
# dataset here.
.hzr_step_resolution <- sqrt(.Machine$double.eps)

#' Is a fitted phase indistinguishable from a step at this data's resolution?
#'
#' The G1 decomposition refuses `nu == 0` with `m >= 0` outright, because that
#' limit is degenerate: \eqn{G} collapses to a step. But a fit can settle at
#' `nu = -1.4e-16`, take the ordinary `m > 0 && nu < 0` branch, and return
#' `converged = TRUE`. The guard sits at the exact point while the pathology is
#' a neighbourhood of it (#448).
#'
#' This thresholds no parameter. It asks whether THE DATA CAN RESOLVE THE RISE:
#' the phase must span the range, below `tol` somewhere and above `1 - tol`
#' somewhere, while at most ONE distinct observed time falls strictly inside
#' the transition. A genuine curve puts many times inside it; a step admits at
#' most the one sitting on it. That one is not an edge case -- an observed time
#' landing exactly on `t_half`, where \eqn{G} is 0.5, is the mechanism #448
#' reports.
#'
#' @param time,time_lower,time_upper Observed times. All three are read: an
#'   interval-censored row is observed as \[lower, upper\] and a
#'   left-truncated fit carries the entry time in `time_lower`, so both are
#'   points the likelihood evaluates at and both can fall below `min(time)`.
#' @param t_half,nu The fitted shape values, reported in the detail string.
#' @param g_fn Function of a time vector returning \eqn{G}. Injected so this
#'   is testable without a fitted object.
#' @param tol The numerically-0 / numerically-1 scale.
#' @return `NULL` when the phase is not a step at this resolution, otherwise a
#'   list with `parameter` and `detail`.
#' @keywords internal
#' @noRd
.hzr_phase_step_detail <- function(time, t_half, nu, g_fn,
                                   time_lower = NULL, time_upper = NULL,
                                   tol = .hzr_step_resolution) {
  if (!is.function(g_fn)) return(NULL)
  # EVERY observed time, not just `time`.  An interval-censored row is observed
  # as [time_lower, time_upper], and on a left-truncated fit time_lower is the
  # entry time, so both can fall below min(time) and both are points the
  # likelihood evaluates at.  The criterion is about the resolution of the
  # times the phase is actually asked about, so it has to see all of them; on
  # `time` alone it would miss the degeneracy on precisely those fits.
  time <- c(time, time_lower, time_upper)
  ok <- is.finite(time) & time > 0
  if (sum(ok) < 2L) return(NULL)
  if (!is.finite(t_half) || !is.finite(nu)) return(NULL)

  ut <- sort(unique(time[ok]))
  if (length(ut) < 2L) return(NULL)

  g <- tryCatch(g_fn(ut), error = function(e) NULL)
  if (is.null(g) || length(g) != length(ut) || !all(is.finite(g))) return(NULL)

  # THE DATA CANNOT RESOLVE THE RISE.  The phase must span the whole range --
  # somewhere below tol and somewhere above 1 - tol -- while at most ONE
  # distinct observed time falls strictly inside the transition.  A genuine
  # curve puts many times inside it; a step admits at most the one sitting on
  # it.
  #
  # The earlier form of this test asked for one ADJACENT PAIR to bracket the
  # entire rise, and it declined on #448's own fit: an observed time lands
  # exactly on t_half, where G is 0.5, so the rise splits across two gaps
  # (-0 -> 0.5, then 0.5 -> 1).  That tie AT t_half is the defect's signature
  # -- the issue's five tied event times -- not an edge case, so a criterion
  # that a tie defeats is the wrong one.
  if (!any(g < tol) || !any(g > 1 - tol)) return(NULL)
  inside <- which(g > tol & g < 1 - tol)
  if (length(inside) > 1L) return(NULL)

  # The bracketing pair, for the message: last time below, first time above.
  j <- max(which(g < tol))
  jj <- min(which(g > 1 - tol))
  n_at <- sum(ok & time >= ut[j] & time <= ut[jj])
  list(
    parameter = "nu",
    detail = sprintf(
      paste0("The early-phase shape is a step at this data's resolution: the ",
             "phase rises from 0 to 1 entirely between two adjacent observed ",
             "times (%.6g and %.6g), with nu = %.3g and t_half = %.6g. %d ",
             "observation(s) fall in that interval, so the log-likelihood is ",
             "discontinuous there and a one-ulp change in nu can move it. ",
             "Treat the shape parameters as unidentified rather than estimated."),
      ut[j], ut[jj], nu, t_half, n_at
    )
  )
}

#' Build a `$boundary` record for one phase, or `NULL`
#'
#' Kept out of `hazard_api.R` so the shared check's body stays small. It is
#' called once per named phase, BEFORE the unbounded-type filter, because a
#' `"cdf"` phase is not an unbounded type and would otherwise never be
#' examined.
#'
#' Only the G1-based types (`"cdf"`, `"hazard"`) can degenerate this way: the
#' step is the \eqn{\nu \to 0} limit of G1, which is what the decomposition
#' refuses at `nu == 0` exactly. `"g3"` has its own parameterisation and
#' `"constant"` has no shape.
#'
#' @param name,type The phase's name and type.
#' @param theta The fitted parameter vector, phase-name prefixed.
#' @param time,time_lower,time_upper Observed times.
#' @return A `$boundary` record, or `NULL`.
#' @keywords internal
#' @noRd
.hzr_phase_step_record <- function(name, type, theta, time,
                                   time_lower = NULL, time_upper = NULL) {
  if (!length(type) || !type %in% c("cdf", "hazard")) return(NULL)
  need <- paste0(name, c(".log_t_half", ".nu", ".m"))
  if (!all(need %in% names(theta))) return(NULL)
  t_half <- exp(unname(theta[[need[[1L]]]]))
  nu <- unname(theta[[need[[2L]]]])
  m <- unname(theta[[need[[3L]]]])
  if (!is.finite(t_half) || !is.finite(nu) || !is.finite(m)) return(NULL)
  g_fn <- function(x) {
    hzr_decompos(x, t_half = t_half, nu = nu, m = m)$G
  }
  d <- .hzr_phase_step_detail(time, t_half = t_half, nu = nu, g_fn = g_fn,
                              time_lower = time_lower, time_upper = time_upper)
  if (is.null(d)) return(NULL)
  list(mechanism = "phase_discontinuity", phase = name,
       parameter = d$parameter, detail = d$detail)
}
