# Conservation of Events whose supremum lies where no events remain (#261)
#
# Under CoE the conserved phase's log_mu is solved so that predicted events
# equal observed ones: it absorbs jevent, the events the other phases leave.
# Where the other phases already account for every event (jevent <= 0) there
# is nothing to solve, and .hzr_conserve_events() returns theta unchanged,
# so the objective there falls back to that log_mu's STARTING value. The
# objective is discontinuous at jevent = 0, and a fit can stop short of a
# higher likelihood at that boundary -- the conserved phase switched off --
# while reporting an ordinary-looking optimum (#261). SAS/C refuses the whole
# run instead (SETCOE1200, consrv.c:66-69).
#
# The check is a certificate, as for #418: at the returned estimates, send
# the conserved phase's scale to the boundary and evaluate the fit's own
# objective there. It records a finding only when that beats the reported
# objective, so it can miss a case but cannot report one that is not there.

#' Certify a higher likelihood on the CoE boundary
#'
#' Two points are scored at the returned estimates, both with the conserved
#' phase's `log_mu` sent to `-1e4`, where `exp()` is exactly 0 and the phase
#' contributes nothing at any finite time or linear predictor (checked):
#' holding every other parameter, and with the other phases' intercepts
#' shifted by one common amount so predicted events again equal observed
#' ones. Only intercepts move; shapes, covariate coefficients and fixed
#' parameters keep their fitted values.
#'
#' @param theta The fitted, named full theta.
#' @param value The objective at `theta`, the likelihood at the returned point.
#' @param converged Did the optimizer report convergence?
#' @param objective_fn Function of a full theta returning the objective the
#'   fit maximised, without conservation applied.
#' @param fixmu_phase,fixmu_pos The conserved phase and its `log_mu` slot.
#' @param log_mu_positions Named list of each phase's `log_mu` slot.
#' @param phases,covariate_counts,x_list The fit's phase specification.
#' @param time,status,time_lower,weights The rows the likelihood read.
#' @param total_events Weighted observed events.
#' @param delta The gain the certificate must exceed.
#' @return A boundary record, or `NULL`.
#' @keywords internal
#' @noRd
.hzr_coe_boundary_record <- function(theta, value, converged, objective_fn,
                                     fixmu_phase, fixmu_pos, log_mu_positions,
                                     phases, covariate_counts, x_list, time,
                                     status, time_lower, weights,
                                     total_events, delta = 0.01) {
  if (!isTRUE(converged) || length(value) != 1L || !is.finite(value) ||
        !is.finite(total_events) || total_events <= 0) {
    return(NULL)
  }
  w <- if (is.null(weights)) rep(1, length(time)) else weights
  off <- theta
  # -1e4, not a milder sentinel: exp(-700) times a large time or a large
  # linear predictor is still a live phase, and the record says the phase is
  # off. Verified below rather than assumed.
  off[[fixmu_pos]] <- -1e4
  left <- .hzr_multiphase_cumhaz(time, off, phases, covariate_counts, x_list,
                                 per_phase = TRUE)[[fixmu_phase]]
  if (any(!is.finite(left)) || any(left != 0)) return(NULL)
  # The event hazard too: a shaped phase's cumulative hazard can underflow
  # where its hazard does not. Identical to the hazard with no such phase.
  none <- off
  none[[fixmu_pos]] <- -Inf
  h_off <- .hzr_multiphase_hazard(time, off, phases, covariate_counts, x_list)
  h_none <- .hzr_multiphase_hazard(time, none, phases, covariate_counts,
                                   x_list)
  if (!isTRUE(all(h_off == h_none))) return(NULL)
  candidates <- list(off)

  # The other phases rescaled together so the events are conserved again,
  # on the same entry-time scale the conservation step uses.
  others <- setdiff(names(phases), fixmu_phase)
  cum <- .hzr_multiphase_cumhaz(time, off, phases, covariate_counts, x_list)
  entry <- .hzr_multiphase_entry(time, status, time_lower)
  if (!is.null(entry) && any(entry > 0)) {
    cum <- cum - .hzr_multiphase_cumhaz(entry, off, phases, covariate_counts,
                                        x_list)
  }
  predicted <- sum(w * cum)
  if (length(others) && is.finite(predicted) && predicted > 0) {
    shift <- log(total_events / predicted)
    rescaled <- off
    for (nm in others) {
      pos <- log_mu_positions[[nm]]
      rescaled[pos] <- rescaled[pos] + shift
    }
    candidates[[2L]] <- rescaled
  }

  scored <- lapply(candidates, function(p) {
    v <- tryCatch(objective_fn(p), error = function(e) NA_real_)
    list(theta = p, ll = if (length(v) == 1L && is.finite(v)) v else -Inf)
  })
  best <- scored[[which.max(vapply(scored, `[[`, numeric(1), "ll"))]]
  gain <- best$ll - value
  if (!(gain > delta)) return(NULL)

  list(
    mechanism = "coe_no_events_left",
    phase = fixmu_phase,
    parameter = "log_mu",
    gain = gain,
    certificate_theta = best$theta,
    certificate_loglik = best$ll,
    detail = paste0(
      "Conservation of Events solved phase '", fixmu_phase, "''s log_mu to ",
      "absorb the events the other phases leave, but the log-likelihood is ",
      format(gain, digits = 3), " higher at a point with that phase's scale ",
      "sent to zero, so the reported estimates are not the maximum. A CoE ",
      "fit can stop short of such a point: where the other phases account ",
      "for every event, the conserved log_mu has nothing to solve for and ",
      "falls back to its starting value, so the objective is discontinuous ",
      "there. Consider the model without phase '", fixmu_phase, "', or ",
      "conserve = FALSE."
    )
  )
}

#' Record a conserved phase that has left the model
#'
#' With the gradient of the CoE objective (#565) a fit can run towards the
#' boundary the certificate above points at, not stop short of it: the
#' conserved phase's scale heads for zero and the phase contributes nothing
#' the data can see. The identifiability check already warns about such a
#' phase, whichever phase it is. This records, in the same family as the
#' certificate, that it is the CONSERVED one, so a reader of `$boundary`
#' finds it under `coe_no_events_left` either way. It raises no warning of
#' its own: `warned_by` names the warning that has already said so.
#'
#' @param phase_share The frame from `.hzr_check_phase_identifiability()`.
#' @param fixmu_phase The conserved phase.
#' @param tol The threshold that check used.
#' @return A boundary record, or `NULL`.
#' @keywords internal
#' @noRd
.hzr_coe_vanished_record <- function(phase_share, fixmu_phase, tol) {
  share <- phase_share$share[match(fixmu_phase, phase_share$phase)]
  if (length(share) != 1L || !is.finite(share) || !(share < tol)) {
    return(NULL)
  }
  list(
    mechanism = "coe_no_events_left",
    phase = fixmu_phase,
    parameter = "log_mu",
    share = share,
    tol = tol,
    warned_by = "phase_share",
    detail = paste0(
      "Conservation of Events solved phase '", fixmu_phase, "''s log_mu to ",
      "absorb the events the other phases leave, and they leave almost ",
      "none: that phase contributes at most ", signif(share, 3), " of the ",
      "cumulative hazard at any observed time. The estimates are on the ",
      "boundary where the conserved phase is switched off. Consider the ",
      "model without phase '", fixmu_phase, "', or conserve = FALSE."
    )
  )
}
