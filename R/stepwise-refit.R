# stepwise-refit.R -- Scope-mutating refit wrapper for stepwise selection.
#
# Step 8.4 of STEPWISE-DESIGN.md (shared with the forward-step driver).
# Given a current fit plus a single scope mutation (add/drop a variable,
# optionally scoped to a phase), rebuild the hazard() call with the
# updated formula / phases list and return the refitted object.
#
# Design note on warm-start
# -------------------------
# Single-distribution fits: hazard() requires an explicit `theta`
# starting vector (the `fit && !is.null(theta)` guard in hazard_api.R).
# We warm-start by re-using the current theta, appending a zero for the
# newly-added beta or dropping the row for the removed one. This lands
# the optimizer near the MLE and typically converges in a handful of
# BFGS iterations.
#
# Multiphase fits: theta = NULL lets hazard() reassemble starting values
# from the phase specs. Warm-starting the full multiphase vector is
# intricate (per-phase layout with optional fixed shapes) and deferred
# to a future optimisation; the Conservation-of-Events adjustment and
# multi-start loop make re-initialisation cheap enough for v1.
#
# Scope for v1: **main effects only**. A term like a multi-level
# factor or a spline that expands to several coefficients would break
# the one-zero-per-add / one-index-per-drop theta layout. The upstream
# guard in .hzr_candidate_coef_name() errors on such terms before the
# refit fires, so callers via hzr_stepwise() never reach this function
# with an expanding term. Anyone calling .hzr_refit_with_scope()
# directly should pre-expand factors and splines into explicit main
# effects.

#' The objective a fit was estimated under
#'
#' `objective` changes the estimand, not just the arithmetic: on the
#' esophagectomy reference the SAS interval density and the likelihood differ
#' by about 22 log-likelihood units, which is larger than most single-variable
#' effects. Anything that refits a stored fit, or evaluates its likelihood
#' again, has to reuse the same one or it silently compares two estimands ---
#' a populated `delta_logLik` / `aic` computed against a model the base fit
#' was never fitted to. Every such caller reads this one accessor so they
#' cannot drift apart.
#'
#' A fit with no `objective` recorded predates the argument, so it was
#' necessarily estimated under the likelihood; the default is a fact about
#' those objects, not a fallback guess.
#'
#' @param fit A fitted `hazard` object.
#' @return `"likelihood"` or `"sas"`.
#'
#' @keywords internal
#' @noRd
.hzr_fit_objective <- function(fit) {
  fit$spec$objective %||% "likelihood"
}

#' Why a fit cannot be refit with a mutated scope
#'
#' Single decision point for "can `.hzr_refit_with_scope()` handle this
#' fit at all". `hzr_stepwise()` asks the same question up front so that a
#' base fit no candidate could ever be refit from is rejected once, with a
#' message that names the remedy, rather than failing candidate-by-candidate
#' into N warnings and an empty screen that looks like an honest null result
#' (#159). Both callers read this one function so the two answers cannot
#' drift apart.
#'
#' The answer differs by distribution, and deliberately so. A
#' single-distribution scope change adds or drops a term in the **global**
#' formula, so without one there is nothing to mutate. A multiphase scope
#' change rewrites the **phase** formula instead, and the global formula only
#' ever carried the response --- so a vector-interface multiphase fit refits
#' perfectly well from its stored response vectors (#160). Blocking it shut
#' out every translated SAS `SELECTION` job, since SAS's censoring statements
#' map onto this package's `-1/0/1/2` coding, which `survival::Surv()` does
#' not share; 25 of the 26 such jobs also use `ICENSOR`. Requiring a formula
#' there would have forced exactly the status-code round-trip `AGENTS.md`
#' records as having shipped a wrong-answer bug.
#'
#' @param fit A fitted `hazard` object.
#' @return `NULL` when the fit can be refit, otherwise a character scalar
#'   naming the obstruction, phrased to follow "... because".
#'
#' @keywords internal
#' @noRd
.hzr_refit_blocker <- function(fit) {
  if (!is.null(fit$call$formula)) {
    return(NULL)
  }

  if (!identical(fit$spec$dist, "multiphase")) {
    return(paste0(
      "it was built via the vector interface (`time =` / `status =`) ",
      "rather than the formula interface, so it stores no model formula ",
      "to mutate. A multiphase fit does not need one --- its scope lives in ",
      "the phase formulas --- but a single-distribution fit adds and drops ",
      "terms in the global formula"
    ))
  }

  # Multiphase and no formula: the refit rebuilds the call from the stored
  # response vectors, so they have to be there. hazard() always records them,
  # so this fires only for a hand-assembled object -- but an absent `time`
  # would otherwise reach hazard() as a missing argument several frames later.
  if (is.null(fit$data$time) || is.null(fit$data$status)) {
    return(paste0(
      "it was built via the vector interface but stores no `time` / ",
      "`status` vectors to refit from"
    ))
  }

  NULL
}


#' Refit a hazard model with a single scope mutation applied
#'
#' @param current A fitted `hazard` object that was built via the
#'   `formula` / `data` interface.
#' @param action Either `"add"` or `"drop"`.
#' @param var Character scalar naming the variable to add or drop.
#' @param phase For multiphase models, character scalar naming the
#'   phase whose scope changes. Ignored (must be NULL) otherwise.
#' @param data Data frame the original fit was built on. Passed to
#'   `hazard()` for the refit.
#' @param ... Additional named args forwarded to `hazard()` (for
#'   example `control = ...`). `time_windows` and `weights` are pulled
#'   from `current$data` automatically; passing them via `...` will
#'   override.
#'
#' @return A new fitted `hazard` object with `$converged` possibly
#'   FALSE if the refit failed to converge.
#'
#' @keywords internal
#' @noRd
.hzr_refit_with_scope <- function(current, action = c("add", "drop"),
                                   var, phase = NULL, data, ...) {
  action <- match.arg(action)
  if (!inherits(current, "hazard")) {
    stop("`current` must be a fitted `hazard` object.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  if (!is.character(var) || length(var) != 1L || !nzchar(var)) {
    stop("`var` must be a non-empty character scalar.", call. = FALSE)
  }

  dist <- current$spec$dist
  user_args <- list(...)

  # Reuse the original fit's weights / time_windows unless user overrides
  default_weights <- current$data$weights
  default_windows <- current$spec$time_windows

  weights <- if ("weights" %in% names(user_args)) {
    user_args$weights
  } else {
    default_weights
  }
  time_windows <- if ("time_windows" %in% names(user_args)) {
    user_args$time_windows
  } else {
    default_windows
  }

  # `objective` is supplied from the base fit below, so a user-supplied one
  # would match the same formal twice and error inside do.call() with a message
  # about argument matching that says nothing about estimands. Unlike `weights`
  # and `time_windows`, it is NOT user-overridable: the whole point of carrying
  # it is that a candidate refit must be comparable to the model it is being
  # compared against. A redundant value is dropped; a conflicting one is
  # refused, because silently discarding it would return a full result that
  # ignored an explicit argument.
  if ("objective" %in% names(user_args)) {
    fit_objective <- .hzr_fit_objective(current)
    if (!identical(user_args$objective, fit_objective)) {
      stop("`objective` cannot be changed in a refit. The base fit was ",
           "estimated under objective = \"", fit_objective, "\", and refitting ",
           "candidates under ", deparse1(user_args$objective), " would make ",
           "`delta_logLik`, `aic` and `delta_aic` differences between two ",
           "estimands rather than between two models -- a populated `$steps` ",
           "table whose comparisons do not mean what they appear to. Refit the ",
           "base model with that objective and run the selection from there.",
           call. = FALSE)
    }
  }

  extra_args <- user_args[!names(user_args) %in%
                            c("weights", "time_windows", "objective")]

  blocker <- .hzr_refit_blocker(current)
  if (!is.null(blocker)) {
    stop("`current` cannot be refit because ", blocker, ".", call. = FALSE)
  }
  # Recover the original formula; see .hzr_stored_formula() for why this
  # cannot be a deparse. A vector-interface fit has none -- which the blocker
  # above has already established is survivable on the multiphase path only.
  has_formula <- !is.null(current$call$formula)
  current_formula <- if (has_formula) {
    .hzr_stored_formula(current, "`current`")
  } else {
    NULL
  }

  if (dist == "multiphase") {
    if (is.null(phase)) {
      stop("`phase` is required when `current` is a multiphase model.",
           call. = FALSE)
    }
    if (!phase %in% names(current$spec$phases)) {
      stop("Unknown phase: ", sQuote(phase), ".  Available: ",
           paste(sQuote(names(current$spec$phases)), collapse = ", "),
           call. = FALSE)
    }

    new_phases <- current$spec$phases
    new_phases[[phase]] <- .hzr_phase_update_formula(
      new_phases[[phase]], action = action, var = var
    )

    # The scope change above rewrote the PHASE formula; the global formula
    # only ever carried the response. So a vector-interface base fit refits
    # fine, provided its stored response vectors are handed back --
    # time_lower / time_upper included. The formula path gets those from
    # Surv(); a vector refit that dropped them would silently refit a
    # left-truncated cohort as at risk from time 0.
    response_args <- if (has_formula) {
      list(formula = current_formula, data = data)
    } else {
      c(Filter(Negate(is.null),
               current$data[c("time", "status", "time_lower", "time_upper")]),
        list(data = data))
    }

    # Reuse the base fit's objective. Dropping it refits every candidate
    # under the likelihood while the base fit's `objective` is the SAS
    # density, so `delta_logLik` and `aic` would be differenced across two
    # estimands -- a full `$steps` table, no warning, wrong numbers.
    do.call(hazard, c(
      response_args,
      list(
        dist         = "multiphase",
        phases       = new_phases,
        weights      = weights,
        time_windows = time_windows,
        objective    = .hzr_fit_objective(current),
        fit          = TRUE
      ),
      extra_args
    ))
  } else {
    # Single-distribution path: mutate the global formula, warm-start
    # theta by inserting / dropping the relevant beta slot. Unlike multiphase
    # this genuinely cannot proceed without a formula, which is why
    # .hzr_refit_blocker() still refuses that combination -- assert it rather
    # than letting a NULL formula travel into .hzr_formula_update().
    if (!has_formula) {
      stop("Internal: a single-distribution refit reached the formula ",
           "mutation with no formula. .hzr_refit_blocker() should have ",
           "refused this fit.", call. = FALSE)
    }
    new_formula <- .hzr_formula_update(current_formula, action, var)

    n_shape <- .hzr_shape_parameter_count(dist, control = current$spec$control)
    theta_old <- current$fit$theta
    if (is.null(theta_old)) {
      stop("`current` has no fitted theta; refit requires a fitted model.",
           call. = FALSE)
    }

    current_vars <- .hzr_scope_current_vars(current)
    if (action == "add") {
      # Warm-start with an extra zero for the new beta (appended last,
      # matching formula-term ordering).
      if (var %in% current_vars) {
        # Already present; just rebuild theta as-is
        theta_start <- theta_old
      } else {
        theta_start <- c(theta_old, 0)
      }
    } else {
      if (!var %in% current_vars) {
        theta_start <- theta_old
      } else {
        # Drop position: position in beta slot = match index within
        # current_vars, shifted by n_shape.
        drop_idx <- match(var, current_vars) + n_shape
        theta_start <- theta_old[-drop_idx]
      }
    }

    do.call(hazard, c(
      list(
        formula      = new_formula,
        data         = data,
        dist         = dist,
        theta        = theta_start,
        weights      = weights,
        time_windows = time_windows,
        fit          = TRUE
      ),
      extra_args
    ))
  }
}
