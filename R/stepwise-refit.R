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
# Multiphase fits: fitted from two starts, the base's estimates (by
# parameter name, .hzr_multiphase_warm_start()) and the phase specs' default
# values, and the better fit kept (.hzr_refit_best_start()). From the default
# start alone a refit ended below the base it contains, reporting
# converged = TRUE (#551); from the warm start alone it sometimes stopped at
# a lower optimum than the default start reaches.
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
#' again, has to reuse the same one or it silently compares two estimands:
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
#' ever carried the response, so a vector-interface multiphase fit refits
#' perfectly well from its stored response vectors (#160). Blocking it shut
#' out every translated SAS `SELECTION` job, since SAS's censoring statements
#' map onto this package's `-1/0/1/2` coding, which `survival::Surv()` does
#' not share; 25 of the 26 such jobs also use `ICENSOR`. Requiring a formula
#' there would have forced exactly the status-code round-trip `AGENTS.md`
#' records as having shipped a wrong-answer bug.
#'
#' @param fit A fitted `hazard` object.
#' @param stepped Names of the multiphase phases a step can change, or
#'   `NULL` for all of them; passed to `.hzr_inherit_blocker()`.
#' @return `NULL` when the fit can be refit, otherwise a character scalar
#'   naming the obstruction, phrased to follow "... because".
#'
#' @keywords internal
#' @noRd
.hzr_refit_blocker <- function(fit, stepped = NULL) {
  inherit_blocker <- .hzr_inherit_blocker(fit, stepped = stepped)
  if (!is.null(inherit_blocker)) {
    return(inherit_blocker)
  }
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


#' A phase formula the fit ignored
#'
#' A multiphase fit saved before #299 and built without `data` can carry a
#' phase formula it never used. Any refit, a stepwise step or a bootstrap
#' replicate, would start using it, so both refuse such a fit. This is the
#' first check of `.hzr_inherit_blocker()`, kept apart so `hzr_bootstrap()`
#' can apply it without the blocker's stepping-only checks.
#'
#' @param fit A fitted `hazard` object.
#' @return `NULL`, or a character scalar naming the phase and its formula,
#'   carrying its own remedy.
#' @keywords internal
#' @noRd
.hzr_ignored_phase_formula <- function(fit) {
  if (!identical(fit$spec$dist, "multiphase")) {
    return(NULL)
  }
  # A fit saved before #299 can carry a phase formula it never used: built
  # without `data`, the phase took the global `x` or no columns. A refit is
  # given `data`, so it would build that formula's columns into a model the
  # base fit never had. The fit's record decides; a fit saved before the
  # record is judged by its columns, and only when it was built without
  # `data`. `~ 1` counts when the fit had a global `x`: it took `x`, and a
  # refit given `data` would give the phase no columns.
  record <- attr(fit$fit$x_list, "from_formula")
  has_x <- !is.null(fit$data$x) && NCOL(fit$data$x) > 0L
  # "Built without `data`" is read from the evaluated frame, which every fit
  # since 1.1.0 stores (NULL when there was no `data`). The call records
  # syntax, not values: `data = d` with `d` NULL still names `d`, so a test
  # on the call let such a fit through (#310). The call is read only for an
  # object saved before the frame was stored.
  pre_frame <- !"frame" %in% names(fit$data)
  no_data <- if (pre_frame) is.null(fit$call$data) else is.null(fit$data$frame)
  # An object saved before 1.1.0 has neither the frame nor the record, so a
  # vector-interface call that names `data` cannot be judged without
  # evaluating it: `data = dd` reads the same whether `dd` was a data frame or
  # NULL, and by the time a saved object is reloaded `dd` may be gone or
  # rebound. Such a fit is refused whenever its phase looks inherited, which
  # also refuses a fit that genuinely used a phase formula whose columns are
  # the inherited names: a refit request, not a wrong number (#324). A fit
  # with a record, or with no `data` in its call, is decided by the check
  # before it.
  # The vector interface is recognised by `time =` in the call, not by a
  # missing formula: a wrapper writes `formula = fml`, a symbol even when
  # `fml` was NULL, and without `call_env` (none before 1.2.2) a later
  # binding of `fml` would decide the refit instead (#324).
  vector_call <- is.null(fit$call$formula) || "time" %in% names(fit$call)
  undecidable <- pre_frame && vector_call
  for (nm in names(fit$spec$phases)) {
    pf <- fit$spec$phases[[nm]]$formula
    has_terms <- !is.null(pf) && .hzr_phase_formula_has_terms(pf)
    if (!is.null(pf) && (has_x || has_terms) &&
          (!is.null(record) || no_data) &&
          .hzr_phase_inherits_global(fit, nm)) {
      consequence <- if (has_terms) {
        "would add that formula's columns"
      } else {
        "would drop the global `x` columns the phase was fitted with"
      }
      return(paste0(
        "phase '", nm, "' has a formula, `",
        paste(deparse(pf), collapse = " "), "`, that the fit ignored: it ",
        "was built without `data`, and a refit given `data` ", consequence,
        ". Refit the base model with `data =` and retry"
      ))
    }
    if (!is.null(pf) && (has_x || has_terms) && undecidable &&
          .hzr_phase_inherits_global(fit, nm)) {
      return(paste0(
        "phase '", nm, "' has a formula, `",
        paste(deparse(pf), collapse = " "), "`, and the fit stores neither ",
        "its data frame nor a record of whether a phase formula was used, as ",
        "a fit saved before 1.1.0 does (or one whose `data$frame` was ",
        "removed). Its columns are the ones the phase would inherit, so the ",
        "fit could be either model, and a refit given `data` could be the ",
        "other one. Refit the base model with the current version, passing ",
        "`data =`, and retry"
      ))
    }
  }
  NULL
}


#' Why a multiphase fit's inherited design cannot be stepped
#'
#' A phase with no formula of its own inherits the global design, and a
#' stepwise step on it rebuilds that design as a phase formula from the
#' global formula's terms (#284). The rebuild is exact only for plain
#' one-column terms, so three cases are refused here, before any fitting:
#' a design matrix passed directly as `x`, which has no terms at all;
#' time-varying coefficients (`time_windows`), which a phase formula would
#' fit as one constant effect; and a term that expands to more than one
#' column, such as a factor, which a step cannot add or drop as one
#' coefficient. Each silently changed the phase's design before. A fourth,
#' checked first and for every phase: a phase formula the fit ignored because
#' it was built without `data` (#299), which any refit would start using.
#'
#' @param fit A fitted `hazard` object.
#' @param stepped Names of the phases a step can change, or `NULL` for all
#'   of them. The `time_windows` and multi-column checks apply only to these;
#'   a direct `x` is refused for any inheriting phase.
#' @return `NULL`, or a character scalar phrased to follow "... because",
#'   carrying its own remedy.
#' @keywords internal
#' @noRd
.hzr_inherit_blocker <- function(fit, stepped = NULL) {
  if (!identical(fit$spec$dist, "multiphase")) {
    return(NULL)
  }
  ignored <- .hzr_ignored_phase_formula(fit)
  if (!is.null(ignored)) {
    return(ignored)
  }
  inherits <- vapply(fit$spec$phases, function(ph) is.null(ph$formula),
                     logical(1))
  if (!any(inherits)) {
    return(NULL)
  }
  name_phases <- function(p) {
    paste0(if (length(p) > 1L) "phases " else "phase ",
           paste0("'", p, "'", collapse = ", "),
           if (length(p) > 1L) " inherit" else " inherits")
  }
  inheriting <- names(fit$spec$phases)[inherits]
  who <- name_phases(inheriting)
  fix <- paste0(". Give each phase its covariates with ",
                "hzr_phase(formula = ~ ...), with the columns in `data`")

  if (is.null(fit$call$formula)) {
    # Ask the design that was built, not the call: an `x` argument bound to
    # NULL names `x` but built no columns, so an inheriting phase loses none.
    if (!is.null(fit$data$x) && ncol(fit$data$x) > 0L) {
      return(paste0(who, " a design matrix passed directly as `x`, and a ",
                    "refit has no terms to rebuild it from", fix))
    }
    return(NULL)
  }

  # A global `.` is left to hzr_stepwise()'s own refusal (#279).
  tt <- tryCatch(stats::terms(.hzr_stored_formula(fit)),
                 error = function(e) NULL)
  labels <- if (is.null(tt)) character() else attr(tt, "term.labels")
  if (length(labels) == 0L) {
    return(NULL)
  }
  # A refit hands an inheriting phase that is not stepped the same global
  # design back, so these two checks apply only to a phase being stepped.
  # (A direct `x`, above, is lost by every inheriting phase, stepped or not.)
  if (!is.null(stepped)) {
    inheriting <- intersect(stepped, inheriting)
    if (length(inheriting) == 0L) {
      return(NULL)
    }
    who <- name_phases(inheriting)
  }
  reason <- if (!is.null(fit$spec$time_windows)) {
    paste0("time-varying coefficients (`time_windows`), which a phase ",
           "formula would fit as one constant effect")
  } else if (!is.null(fit$data$x) && ncol(fit$data$x) != length(labels)) {
    paste0("a term that expands to more than one column, such as a factor ",
           "with more than two levels or `poly()`, which a step cannot add ",
           "or drop as one coefficient")
  }
  if (is.null(reason)) {
    return(NULL)
  }
  paste0(who, " a global design with ", reason, fix)
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

  # The refit supplies its own start, so a `theta` forwarded from the caller
  # met the same formal twice and failed with an argument-matching error;
  # before #551 a multiphase refit instead used it as every candidate's
  # start, for a model of a different length.
  if ("theta" %in% names(user_args)) {
    stop("`theta` cannot be passed to a stepwise refit: each candidate is ",
         "started from the current model's estimates. Set the starting ",
         "values on the base fit instead.", call. = FALSE)
  }

  extra_args <- user_args[!names(user_args) %in%
                            c("weights", "time_windows", "objective")]

  blocker <- .hzr_refit_blocker(current, stepped = phase)
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
    # A phase with no formula inherits the global terms, and the step starts
    # from them (#284). Read only for such a phase.
    inherited <- if (is.null(new_phases[[phase]]$formula)) {
      .hzr_inherited_rhs(current)
    }
    new_phases[[phase]] <- .hzr_phase_update_formula(
      new_phases[[phase]], action = action, var = var, inherited = inherited
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
    refit_args <- c(
      response_args,
      list(
        dist         = "multiphase",
        phases       = new_phases,
        weights      = weights,
        time_windows = time_windows,
        objective    = .hzr_fit_objective(current)
      ),
      extra_args
    )
    # Two starts, and the better fit kept (#551). The base's estimates with a
    # new coefficient at 0 ARE the base model, so a refit started there cannot
    # end below the base; from the phase specs' default start alone it did,
    # by up to 25 log-likelihood units, reporting converged = TRUE. But the
    # likelihood is multimodal, and the default start sometimes reaches a
    # higher optimum than the warm one, so both are fitted. The default start
    # holds fixed shapes at the base's values, or it fits another model. The
    # unfitted object gives the candidate's parameter names by the fit's own
    # naming; its warnings repeat the fits' and are dropped.
    proto <- suppressWarnings(do.call(hazard,
                                      c(refit_args, list(fit = FALSE))))
    theta_start <- .hzr_multiphase_warm_start(current, proto)
    .hzr_refit_best_start(
      warm    = .hzr_refit_capture(do.call(hazard, c(
        refit_args, list(theta = theta_start, fit = TRUE)
      ))),
      default = .hzr_refit_capture(do.call(hazard, c(
        refit_args,
        list(theta = .hzr_multiphase_default_start(current, proto),
             fit = TRUE)
      )))
    )
  } else {
    # Single-distribution path: mutate the global formula, warm-start
    # theta from the base fit by design column. Unlike multiphase
    # this genuinely cannot proceed without a formula, which is why
    # .hzr_refit_blocker() still refuses that combination -- assert it rather
    # than letting a NULL formula travel into .hzr_formula_update().
    if (!has_formula) {
      stop("Internal: a single-distribution refit reached the formula ",
           "mutation with no formula. .hzr_refit_blocker() should have ",
           "refused this fit.", call. = FALSE)
    }
    new_formula <- .hzr_formula_update(current_formula, action, var)

    # Counted as hazard() counts it when it checks theta (.hzr_check_theta()),
    # without control$shape_param_count: the fit was accepted on that count,
    # and a misleading control would refuse its refit.
    n_shape <- .hzr_shape_parameter_count(dist)
    theta_old <- current$fit$theta
    if (is.null(theta_old)) {
      stop("`current` has no fitted theta; refit requires a fitted model.",
           call. = FALSE)
    }

    # The warm start is laid out against the design the refit will build,
    # read off an unfitted hazard() call with the same arguments. Appending a
    # zero for an added term put it after an interaction that terms() orders
    # last, and hazard() keeps a named theta's names, so the new coefficient
    # was reported under the interaction's name (#489).
    new_design <- suppressWarnings(.hzr_muffle_intercept_warning(
      do.call(hazard, c(
        list(
          formula      = new_formula,
          data         = data,
          dist         = dist,
          weights      = weights,
          time_windows = time_windows,
          fit          = FALSE
        ),
        extra_args
      ))
    ))
    theta_start <- .hzr_refit_warm_start(
      theta_old, n_shape,
      old_cols = colnames(current$data$x),
      new_cols = colnames(new_design$data$x),
      old_windows = current$spec$time_windows,
      new_windows = new_design$spec$time_windows
    )

    .hzr_muffle_intercept_warning(do.call(hazard, c(
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
    )))
  }
}


#' Warm start for a single-distribution refit, matched by design column
#'
#' Each coefficient of the base fit moves to the column of the new design
#' with the same name; a column the base fit did not have starts at 0, and a
#' dropped column's coefficient is discarded. With `time_windows`, the design
#' is one block of columns per window, and the match is made within each
#' block. When the refit uses different windows from the base fit, no block
#' corresponds to another and every coefficient starts at 0.
#'
#' A named theta keeps its names for the coefficients it carries over, and a
#' new column's coefficient is named after the column (`<column>_w<k>` under
#' windows, as the expanded design names it). An unnamed theta stays unnamed.
#' When the refit's rebuilt formula names a column differently, as when
#' `mal:age` becomes `age:mal` or an interaction is reordered after a main
#' effect is dropped, that coefficient is matched by the new column name: it
#' takes the column's name rather than the user's and warm-starts at 0. Its
#' fitted value is unaffected.
#'
#' @param theta_old The base fit's theta: shape parameters, then one
#'   coefficient per design column (per window).
#' @param n_shape Number of leading shape parameters.
#' @param old_cols,new_cols Column names of the base and new unexpanded
#'   designs; `NULL` for a model with no covariates.
#' @param old_windows,new_windows The two fits' `time_windows`, or `NULL`.
#' @return The warm-start theta for the new design.
#' @noRd
.hzr_refit_warm_start <- function(theta_old, n_shape, old_cols, new_cols,
                                  old_windows, new_windows) {
  n_win_old <- length(old_windows) + 1L
  n_win_new <- length(new_windows) + 1L
  p_old <- length(old_cols)
  p_new <- length(new_cols)
  if (length(theta_old) != n_shape + p_old * n_win_old) {
    stop("Internal: the base fit's theta has ", length(theta_old),
         " entries, but its design implies ", n_shape + p_old * n_win_old,
         ", so its coefficients cannot be matched to the refit's columns.",
         call. = FALSE)
  }

  beta_old <- theta_old[-seq_len(n_shape)]
  beta_new <- numeric(p_new * n_win_new)
  names_new <- if (n_win_new > 1L) {
    paste0(rep(new_cols, n_win_new), "_w",
           rep(seq_len(n_win_new), each = p_new))
  } else {
    as.character(new_cols)
  }
  if (identical(old_windows, new_windows)) {
    from <- match(new_cols, old_cols)
    for (k in seq_len(n_win_new)) {
      at_new <- (k - 1L) * p_new + which(!is.na(from))
      at_old <- (k - 1L) * p_old + from[!is.na(from)]
      beta_new[at_new] <- beta_old[at_old]
      if (!is.null(names(beta_old))) {
        names_new[at_new] <- names(beta_old)[at_old]
      }
    }
  }

  theta_start <- c(theta_old[seq_len(n_shape)], beta_new)
  if (is.null(names(theta_old))) {
    return(unname(theta_start))
  }
  names(theta_start) <- c(names(theta_old)[seq_len(n_shape)], names_new)
  theta_start
}

#' Warm-start theta for a multiphase refit
#'
#' Each parameter of the refit model that the base model also has starts at
#' the base's estimate, matched by the fit's own name (`phase.column`, or
#' `phase.log_mu` and the shapes); a coefficient the base does not have, the
#' one a step adds, starts at 0. A dropped coefficient is simply absent from
#' the refit model. Fixed shapes carry their fixed values, which the base
#' holds unchanged. So an entry starts at a point where the refit model
#' reproduces the base log-likelihood, and cannot end below it (#551).
#'
#' @param current The base fit.
#' @param proto The refit model, unfitted (`fit = FALSE`).
#' @return A named numeric theta for `proto`, on the fit's internal scale.
#' @keywords internal
#' @noRd
.hzr_multiphase_warm_start <- function(current, proto) {
  old_names <- .hzr_evaluate_prepare(current)$names
  new_names <- .hzr_evaluate_prepare(proto)$names
  theta_old <- current$fit$theta
  if (length(theta_old) != length(old_names) || anyDuplicated(old_names) ||
        anyDuplicated(new_names)) {
    stop("Internal: the base fit's parameters could not be matched to the ",
         "refit's by name (", length(theta_old), " estimates, ",
         length(old_names), " names).", call. = FALSE)
  }
  theta_start <- stats::setNames(numeric(length(new_names)), new_names)
  shared <- intersect(new_names, old_names)
  theta_start[shared] <- unname(theta_old[match(shared, old_names)])
  theta_start
}

#' Default-start theta for a multiphase refit
#'
#' The phase specs' default starting values -- what the optimizer assembles
#' when handed no theta -- except that every entry the refit does not search
#' over (a fixed or derived shape) carries the base's value. A fixed entry is
#' held at its START value, and a base fitted with a user `theta` holds it at
#' that value rather than the spec's, so a start with the spec's value fits a
#' different model: it won on a changed fixed shape, not on the candidate.
#'
#' @inheritParams .hzr_multiphase_warm_start
#' @return A named numeric theta for `proto`, on the fit's internal scale.
#' @keywords internal
#' @noRd
.hzr_multiphase_default_start <- function(current, proto) {
  prep <- .hzr_evaluate_prepare(proto)
  start <- unlist(lapply(names(prep$phases), function(nm) {
    .hzr_phase_start(prep$phases[[nm]],
                     n_covariates = prep$covariate_counts[[nm]])
  }), use.names = FALSE)
  names(start) <- prep$names
  held <- prep$names[!.hzr_phase_free_mask(prep$phases,
                                            prep$covariate_counts)]
  old_names <- .hzr_evaluate_prepare(current)$names
  held <- intersect(held, old_names)
  start[held] <- unname(current$fit$theta[match(held, old_names)])
  start
}

#' Run one refit, holding its warnings back
#'
#' A multiphase refit is fitted from two starts and only one fit is kept, so
#' the warnings of the discarded fit must not reach the user as if they
#' described the result. The intercept warning the base fit already gave is
#' dropped outright, as `.hzr_muffle_intercept_warning()` drops it.
#'
#' @param expr The refit call, evaluated here.
#' @return A list: `value`, the fit or the error condition, and `warnings`,
#'   the warning conditions it raised.
#' @keywords internal
#' @noRd
.hzr_refit_capture <- function(expr) {
  caught <- list()
  value <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        if (!inherits(w, "hzr_intercept_removed")) {
          caught[[length(caught) + 1L]] <<- w
        }
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  list(value = value, warnings = caught)
}

#' Keep the better of a multiphase refit's two fits
#'
#' A fit is usable when it is a `hazard` object that did not report
#' non-convergence and has a finite objective. Of two usable fits the higher
#' objective wins, and a tie goes to the warm start, which is the one that
#' cannot end below the base. One usable fit wins alone. When neither is
#' usable the default start's outcome is returned exactly as the refit
#' returned it before #551: its error re-raised, or its non-converged fit.
#' The kept fit's warnings are then raised, and `fit$fit$refit_start` records
#' which start won, beside `fit$fit$refit_objectives`, both starts' objectives
#' (`NA` for one that failed).
#'
#' @param warm,default `.hzr_refit_capture()` results.
#' @return The kept fit.
#' @keywords internal
#' @noRd
.hzr_refit_best_start <- function(warm, default) {
  objective_of <- function(r) {
    v <- r$value
    if (inherits(v, "hazard") && !isFALSE(v$fit$converged) &&
          isTRUE(is.finite(v$fit$objective))) v$fit$objective else NA_real_
  }
  obj <- c(warm = objective_of(warm), default = objective_of(default))
  winner <- if (all(is.na(obj))) {
    "default"
  } else if (is.na(obj[["default"]]) ||
               (!is.na(obj[["warm"]]) && obj[["warm"]] >= obj[["default"]])) {
    "warm"
  } else {
    "default"
  }
  kept <- if (winner == "warm") warm else default
  for (w in kept$warnings) warning(w)
  if (inherits(kept$value, "condition")) stop(kept$value)
  fit <- kept$value
  fit$fit$refit_start <- winner
  fit$fit$refit_objectives <- obj
  fit
}


# A refit re-parses the base fit's formula, and the base fit already warned
# that its intercept removal is ignored (#337); do not repeat it per step.
.hzr_muffle_intercept_warning <- function(expr) {
  withCallingHandlers(
    expr,
    hzr_intercept_removed = function(w) invokeRestart("muffleWarning")
  )
}
