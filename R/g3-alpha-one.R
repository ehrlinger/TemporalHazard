# g3 phase with alpha FIXED at 1: the SAS/C reparameterisation (#415)
# ============================================================================
#
# At alpha = 1 the g3 form collapses to G3(t) = (t / tau)^(gamma * eta)
# (hzr_decompos_g3()). Two directions are then unidentified: gamma and eta
# enter only through their product, and tau only through mu * tau^-(gamma*eta).
# A fit that leaves them free walks a ridge and reports convergence at an
# arbitrary point on it (#415 measured gamma-hat from 1.6e6 to 1.04e8).
#
# PROC HAZARD does not fit that model. When ALPHA is FIXED at 1, SETG3 calls
# SETG3_ignore_tau() (hazard/src/model/setg3.c:313-315), which "reexpresses
# [the late phase] in terms of a single exponential function of time"
# (setg3.c:300-303):
#
#   - tau is set to 1 and fixed                          (setg3.c:380-381)
#   - under FIXGE2 gamma and eta are both fixed, at
#     gamma = 2, eta = 1 (or eta = 2, gamma = 1 when
#     ETA = 2 was given)                                 (setg3.c:393-403)
#   - otherwise the product is carried by ONE parameter:
#     both free -> eta is fixed                          (setg3.c:406-407)
#     eta fixed -> gamma = gamma * eta, eta = 1          (setg3.c:409-413)
#     else      -> eta = gamma * eta, gamma = 1          (setg3.c:414-418)
#
# The likelihood is unchanged by this: at alpha = 1 the held and the unheld
# models are the same family, and the hold only removes the directions the
# data cannot see. So this mirrors SAS and says so (warn + record in
# $boundary), rather than refusing. It keys on alpha FIXED at exactly 1, as
# SAS does (`Late.alpha==ONE && hzr_parm_is_fixed(HZ_ALPHA)`); a FREE alpha
# started at 1 is not held, by SAS or here.

#' Hold a g3 phase's redundant shapes when alpha is fixed at 1
#'
#' @param theta Full named starting vector (`<phase>.<parameter>`).
#' @param phases Validated phase list.
#' @return List: `theta` and `phases` after the hold, and `records`, a list
#'   of `$boundary` records (mechanism `"g3_alpha_one"`), one per phase whose
#'   values or fixed set the hold changed. Empty when nothing changed.
#' @keywords internal
#' @noRd
.hzr_g3_alpha_one_hold <- function(theta, phases) {
  records <- list()
  for (k in seq_along(phases)) {
    ph <- phases[[k]]
    nm <- names(phases)[[k]]
    if (!identical(ph$type, "g3")) next
    fixed <- if (is.null(ph$fixed)) character(0) else ph$fixed
    key <- function(p) paste0(nm, ".", p)
    if (!"alpha" %in% fixed || !key("alpha") %in% names(theta)) next
    if (!identical(unname(theta[[key("alpha")]]), 1)) next

    before_theta <- theta
    before_fixed <- fixed
    gamma <- unname(theta[[key("gamma")]])
    eta   <- unname(theta[[key("eta")]])

    theta[[key("log_tau")]] <- 0
    fixed <- union(fixed, "tau")
    if (identical(.hzr_phase_constraint(ph), "eta_gamma")) {
      # eta is derived as 2 / gamma here, so "ETA = 2 given" is gamma = 1.
      gamma_new <- if (isTRUE(all.equal(gamma, 1))) 1 else 2
      theta[[key("gamma")]] <- gamma_new
      theta[[key("eta")]] <- 2 / gamma_new
      fixed <- union(fixed, "gamma")
    } else {
      if (!"gamma" %in% fixed && !"eta" %in% fixed) fixed <- union(fixed, "eta")
      if ("eta" %in% fixed) {
        theta[[key("gamma")]] <- gamma * eta
        theta[[key("eta")]] <- 1
      } else {
        theta[[key("eta")]] <- gamma * eta
        theta[[key("gamma")]] <- 1
      }
    }
    phases[[k]]$fixed <- fixed

    params <- c("log_tau", "gamma", "eta")
    moved <- params[theta[key(params)] != before_theta[key(params)]]
    newly_fixed <- setdiff(fixed, before_fixed)
    if (!length(moved) && !length(newly_fixed)) next

    show <- function(v) format(v, digits = 6)
    listed <- params[sub("^log_", "", params) %in% c(sub("^log_", "", moved),
                                                     newly_fixed)]
    changes <- vapply(listed, function(p) {
      label <- if (p == "log_tau") "tau" else p
      was <- before_theta[[key(p)]]
      now <- theta[[key(p)]]
      if (p == "log_tau") {
        was <- exp(was)
        now <- exp(now)
      }
      paste0(label, " ", show(was), " -> ", show(now),
             if (label %in% newly_fixed) " (now fixed)" else "")
    }, character(1))
    records[[length(records) + 1L]] <- list(
      mechanism = "g3_alpha_one",
      phase = nm,
      parameter = unique(c(sub("^log_", "", moved), newly_fixed)),
      detail = paste0(
        "phase '", nm, "' is g3 with alpha fixed at 1, where G3 = ",
        "(t/tau)^(gamma*eta): tau is confounded with mu, and gamma with eta. ",
        "As PROC HAZARD does, the product is carried by one parameter and ",
        "tau is held at 1 (", paste(changes, collapse = "; "), "). The ",
        "likelihood is unchanged; mu and the product are what is estimated."
      )
    )
  }
  list(theta = theta, phases = phases, records = records)
}


# g3 phase under FIXGE2 (constraint = "eta_gamma"): SAS's alpha rules (#418)
# ============================================================================
#
# Under FIXGE2, SETG3_verify_ge_2() makes gamma * eta = 2 (setg3.c:870-897),
# and SETG3_alpha_fixup() then requires gamma * eta / alpha > 2, that is
# alpha < 1 (setg3.c:815-850):
#
#   - alpha FIXED above 1: PROC HAZARD stops with SETG31040
#     (setg3.c:843-847). Refused here the same way.
#   - alpha FIXED at exactly 1: diverted earlier to SETG3_ignore_tau()
#     (setg3.c:313-315), the #415 hold above.
#   - alpha FREE and started at 1 or more: the start is rewritten to
#     gamma * eta / 3 = 2/3 (setg3.c:848-850). Mirrored, with a warning and
#     a record: a start at alpha = 1 sits on the ridge where gamma is not
#     identified (#415).
#
# NOT covered, and #418 stays open for it: a fit that starts off the ridge
# and climbs toward the gamma -> Inf corner law, where the likelihood has a
# supremum and no maximum. That needs a post-fit profile comparison.

#' Apply PROC HAZARD's FIXGE2 alpha rules to a g3 phase
#'
#' @inheritParams .hzr_g3_alpha_one_hold
#' @return List: `theta` and `records` (mechanism `"g3_fixge2_alpha_start"`).
#'   Stops for a fixed alpha above 1.
#' @keywords internal
#' @noRd
.hzr_g3_fixge2_alpha <- function(theta, phases) {
  records <- list()
  for (k in seq_along(phases)) {
    ph <- phases[[k]]
    nm <- names(phases)[[k]]
    if (!identical(ph$type, "g3")) next
    if (!identical(.hzr_phase_constraint(ph), "eta_gamma")) next
    a_key <- paste0(nm, ".alpha")
    if (!a_key %in% names(theta)) next
    alpha <- unname(theta[[a_key]])
    if (!is.finite(alpha) || alpha < 1) next
    if ("alpha" %in% ph$fixed) {
      if (alpha == 1) next
      stop("phase '", nm, "': alpha is fixed at ", format(alpha, digits = 6),
           " under constraint = \"eta_gamma\" (SAS FIXGE2). With gamma * ",
           "eta = 2 the g3 form needs gamma * eta / alpha > 2, so alpha < 1; ",
           "PROC HAZARD stops on this setup (SETG31040). Fix alpha below 1, ",
           "at 1, or leave it free.", call. = FALSE)
    }
    theta[[a_key]] <- 2 / 3
    records[[length(records) + 1L]] <- list(
      mechanism = "g3_fixge2_alpha_start",
      phase = nm,
      parameter = "alpha",
      detail = paste0(
        "phase '", nm, "' is g3 under constraint = \"eta_gamma\" with a ",
        "free alpha started at ", format(alpha, digits = 6), ". With ",
        "gamma * eta = 2 the form needs alpha < 1, and at alpha = 1 gamma is ",
        "not identified, so the start was moved to 2/3, as PROC HAZARD ",
        "does. The fitted alpha is estimated from there."
      )
    )
  }
  list(theta = theta, records = records)
}
