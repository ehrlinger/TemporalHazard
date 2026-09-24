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

# Returns NULL when the phase is not a step at this resolution, otherwise a
# list(parameter, detail) describing it.  `g_fn` takes a time vector and
# returns G; injecting it keeps this testable without a fitted object.
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

  # An adjacent pair that brackets the whole rise: below tol, then above 1-tol.
  lo <- g[-length(g)]
  hi <- g[-1L]
  jump <- which(lo < tol & hi > 1 - tol)
  if (!length(jump)) return(NULL)
  j <- jump[[1L]]

  n_at <- sum(ok & time >= ut[j] & time <= ut[j + 1L])
  list(
    parameter = "nu",
    detail = sprintf(
      paste0("The early-phase shape is a step at this data's resolution: the ",
             "phase rises from 0 to 1 entirely between two adjacent observed ",
             "times (%.6g and %.6g), with nu = %.3g and t_half = %.6g. %d ",
             "observation(s) fall in that interval, so the log-likelihood is ",
             "discontinuous there and a one-ulp change in nu can move it. ",
             "Treat the shape parameters as unidentified rather than estimated."),
      ut[j], ut[j + 1L], nu, t_half, n_at
    )
  )
}
