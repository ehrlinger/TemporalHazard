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

  extra_args <- user_args[!names(user_args) %in%
                            c("weights", "time_windows")]

  # Recover the original formula; see .hzr_stored_formula() for why this
  # cannot be a deparse.
  raw_formula <- current$call$formula
  if (is.null(raw_formula)) {
    stop("`current` was not built via the formula interface; refit ",
         "requires a formula-based base fit.", call. = FALSE)
  }
  current_formula <- .hzr_stored_formula(current, "`current`")

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

    do.call(hazard, c(
      list(
        formula      = current_formula,
        data         = data,
        dist         = "multiphase",
        phases       = new_phases,
        weights      = weights,
        time_windows = time_windows,
        fit          = TRUE
      ),
      extra_args
    ))
  } else {
    # Single-distribution path: mutate the global formula, warm-start
    # theta by inserting / dropping the relevant beta slot.
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
