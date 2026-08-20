# stepwise.R -- user-facing driver for Phase 4b stepwise covariate
# selection.  Step 8.6 of STEPWISE-DESIGN.md.
#
# The driver combines .hzr_stepwise_forward_step() and
# .hzr_stepwise_backward_step() into a two-way loop with:
#
#   * SAS-style MOVE oscillation guard: when a variable hits
#     `max_move` entries/exits, it is frozen (added to both the
#     force_in and force_out sets used internally).
#   * `max_steps` hard cap that warns on hit.
#   * `force_in` / `force_out` user-supplied constraints.
#   * Optional console trace mirroring the sec.5 output spec.
#
# The returned object inherits from both `hzr_stepwise` and `hazard`
# so it can be passed to `predict()`, `summary()`, `coef()`, etc.
# through the same S3 infrastructure as a plain fit.

#' Stepwise covariate selection for a parametric hazard model
#'
#' Run forward, backward, or two-way stepwise selection on an existing
#' `hazard` fit using score (Q) statistics, Wald p-values, or AIC deltas as
#' the entry / retention criterion.  Phase-specific entry is supported for
#' multiphase models: a covariate can enter one phase and not another.
#'
#' @section Selection direction and criterion:
#'
#' Two arguments shape the search.  `direction` decides which moves are
#' allowed at each step; `criterion` decides how a candidate move is scored
#' and whether it is accepted.
#'
#' \describe{
#'   \item{`direction = "forward"`}{Start from the base model and only
#'     *add* variables --- the best eligible candidate enters each step
#'     until none clears the entry rule.  Variables never leave once in.}
#'   \item{`direction = "backward"`}{Start from the full candidate model and
#'     only *drop* variables --- the weakest term leaves each step until all
#'     survivors clear the retention rule.}
#'   \item{`direction = "both"` (default)}{Two-way stepwise: after each
#'     entry, already-selected variables are re-tested and may be dropped.
#'     This is the SAS `SELECTION = STEPWISE` strategy.  `max_move` caps how
#'     often a single variable may oscillate before it is frozen.}
#' }
#'
#' \describe{
#'   \item{`criterion = "score"` (default)}{Accept moves on SAS-style
#'     significance thresholds, using the score (Q) statistic of the candidate
#'     coefficient --- this reproduces C/SAS HAZARD's `SELECTION` statistic.
#'     Q is evaluated at the *current* model's MLE with the candidate's
#'     coefficient pinned at zero, so **no candidate refit is needed**: the
#'     reduced-model information is inverted once per step and reused across
#'     every candidate.  Only the winner is refit.  A candidate enters if its
#'     p-value is below `slentry`.
#'
#'   Score is an *entry* criterion; the drop path never refit per candidate
#'   in the first place, so removals are tested on the current model's Wald
#'   p-value against `slstay`, as SAS does.
#'
#'   For single-distribution fits, the score criterion computes the observed
#'   information numerically via the suggested \pkg{numDeriv} package and
#'   errors with a clear message if it is not installed; a multiphase fit
#'   uses the analytic Hessian instead and does not need it.
#'
#'   Following SAS, the variance used during *selection* is approximate:
#'   shaping-parameter covariances are ignored.  This affects selection only
#'   --- final-model standard errors are unchanged and still come from the
#'   full Hessian.  Candidates must be single-column numeric main-effect
#'   terms; a factor is rejected with an error rather than skipped.}
#'   \item{`criterion = "wald"`}{Accept moves on SAS-style
#'     significance thresholds, using the Wald \eqn{\chi^2} of the affected
#'     coefficient(s): a candidate enters if its p-value is below `slentry`, and
#'     a term is dropped if its p-value rises above `slstay`.  Entry candidates
#'     are scored from a refit that adds the candidate (so its new coefficient
#'     can be tested); drop candidates are scored from the *current* model's
#'     Wald p-values without a per-candidate refit, and a single refit is run
#'     only after a drop is chosen.  This was the default before version 1.2.0.
#'     It differs algorithmically from C/SAS HAZARD, so the two criteria can
#'     take different step paths --- and select different variable sets ---
#'     even when they converge to a similar final model.}
#'   \item{`criterion = "aic"`}{Accept any move with
#'     \eqn{\Delta\mathrm{AIC} < 0} (a strictly better penalised fit), ignoring
#'     `slentry` / `slstay`.  Entry candidates use the actual
#'     \eqn{\Delta\mathrm{AIC}} from the candidate refit; drop candidates use a
#'     Wald-to-likelihood-ratio approximation,
#'     \eqn{\Delta\mathrm{AIC} \approx W - 2\,\mathrm{df}}, computed from the
#'     current model without a per-candidate refit (the chosen drop is refit
#'     afterwards).  Use this for a non-significance-based,
#'     information-criterion search.}
#' }
#'
#' @param fit A fitted `hazard` object built via the
#'   `formula = Surv(...) ~ predictors, data = df` interface.
#' @param scope Candidate set.  `NULL` (default) uses every data-frame
#'   column not already in the model for every phase.  For
#'   single-distribution fits, pass a one-sided formula
#'   (`~ age + nyha`) or a character vector of names.  For multiphase
#'   fits, pass a named list of one-sided formulas keyed by phase.
#' @param data Data frame the base fit was built on.  Required for
#'   refits.
#' @param direction Search strategy --- one of `"both"` (default),
#'   `"forward"`, or `"backward"`.  Controls whether variables may only
#'   enter, only leave, or both.  See the **Selection direction and
#'   criterion** section.
#' @param criterion Entry / retention rule --- one of `"score"` (default),
#'   `"wald"`, or `"aic"`.  `"score"` and `"wald"` both apply SAS-style
#'   p-value thresholds (`slentry` / `slstay`) but score entry candidates
#'   differently, and can therefore select different variable sets; `"score"`
#'   reproduces C/SAS HAZARD and needs no per-candidate refit.  `"aic"` adds or
#'   drops whenever it lowers the AIC.  See the **Selection direction and
#'   criterion** section.
#' @param slentry Entry p-value threshold for the score / Wald criteria.
#'   Default `0.30` matches SAS `SLENTRY`.
#' @param slstay Retention p-value threshold for the score / Wald criteria.
#'   Default `0.20` matches SAS `SLSTAY`.
#' @param max_steps Hard cap on total accepted actions.  Emits a
#'   `warning()` if hit.  Default `50`.
#' @param max_move Per-variable oscillation cap.  When a variable has
#'   entered + exited more than `max_move` times it is frozen for the
#'   remainder of the run.  Default `4`.
#' @param force_in Character vector of variables that must remain in
#'   the model.  Such variables are still scored and reported in the
#'   selection trace, but are never dropped.
#' @param force_out Character vector of variables that may never be
#'   considered as candidates.
#' @param trace Logical; print step-by-step progress to the console.
#'   Default `TRUE`.
#' @param ... Passed to the underlying `hazard()` refits (e.g.
#'   `control = list(n_starts = 3)`).
#'
#' @return An object of class `c("hzr_stepwise", "hazard")` -- the
#'   final fit augmented with:
#'   \describe{
#'     \item{\code{steps}}{Data frame with one row per accepted /
#'       frozen action; see Details.}
#'     \item{\code{scope}}{Record of the candidate scope, plus
#'       `force_in`, `force_out`, and the frozen set.}
#'     \item{\code{criteria}}{Named list of the threshold / direction
#'       settings actually applied, plus, under `criterion = "score"`,
#'       `n_uncomputable_scores` (how many candidate scores were `NA`),
#'       `uncomputable_reasons` (a named integer vector of *why*) and
#'       `stopped_uncomputable`. Read `uncomputable_reasons` before treating
#'       an unscored candidate as a bad one: `information_indefinite` marks
#'       candidates whose effect is too large for the score test's
#'       approximation at zero, which are typically the strongest variables
#'       on offer rather than degenerate ones. `criterion = "wald"` tests
#'       them.}
#'     \item{\code{trace_msg}}{Character vector of the trace lines,
#'       captured regardless of the `trace` flag.}
#'     \item{\code{elapsed}}{`difftime` from start to finish.}
#'     \item{\code{final_call}}{The call that produced this result.}
#'   }
#'
#' @details
#' The `steps` data frame has columns:
#'
#' \describe{
#'   \item{\code{step_num}}{Integer sequence starting at 1.}
#'   \item{\code{action}}{`"enter"`, `"drop"`, or `"frozen"`.}
#'   \item{\code{variable}}{Variable affected.}
#'   \item{\code{phase}}{Phase name (multiphase) or `NA_character_`.}
#'   \item{\code{criterion}}{The criterion actually applied to this step ---
#'     `"score"`, `"wald"`, or `"aic"`.  Under `criterion = "score"` the drop
#'     rows read `"wald"`, because score is entry-only.}
#'   \item{\code{score}}{Winning score used for the decision.}
#'   \item{\code{stat}, \code{df}}{Test statistic (score Q or Wald) and
#'     degrees of freedom.}
#'   \item{\code{p_value}, \code{delta_aic}}{Always populated when
#'     computable, regardless of the active criterion.}
#'   \item{\code{logLik}, \code{aic}, \code{n_coef}}{Goodness-of-fit
#'     diagnostics of the model *after* this step.}
#' }
#'
#' @examples
#' data(avc)
#' avc <- na.omit(avc)
#' base <- hazard(survival::Surv(int_dead, dead) ~ age,
#'                data = avc, dist = "weibull", fit = TRUE,
#'                theta = c(mu = 0.01, nu = 0.5, 0))
#' \donttest{
#' sw <- hzr_stepwise(base, scope = ~ age + mal,
#'                    data = avc, direction = "forward",
#'                    control = list(n_starts = 1))
#' print(sw)
#' }
#'
#' @seealso [hazard()] for the base model and [hzr_phase()] for multiphase
#'   scopes; [stepwise_trace()] to retrieve the captured selection log.
#' @export
hzr_stepwise <- function(fit,
                         scope     = NULL,
                         data,
                         direction = c("both", "forward", "backward"),
                         criterion = c("score", "wald", "aic"),
                         slentry   = 0.30,
                         slstay    = 0.20,
                         max_steps = 50L,
                         max_move  = 4L,
                         force_in  = character(),
                         force_out = character(),
                         trace     = TRUE,
                         ...) {
  direction <- match.arg(direction)
  criterion <- match.arg(criterion)

  if (!inherits(fit, "hazard")) {
    stop("`fit` must be a `hazard` object.", call. = FALSE)
  }
  if (missing(data) || !is.data.frame(data)) {
    stop("`data` must be a data frame (typically the frame used for the base fit).",
         call. = FALSE)
  }

  # The score (Q) statistic is evaluated at the base model's fitted MLE, so a
  # base that did not converge (or was never fitted, leaving an empty theta)
  # has nothing to score.  Catch that here with an actionable message rather
  # than letting the internal `.hzr_score_free_idx()` guard fire deep in the
  # candidate loop.  `wald` and `aic` refit each candidate from the call and
  # legitimately tolerate a non-converged base, so they are left alone.
  if (criterion == "score" &&
        (!isTRUE(fit$fit$converged) || length(fit$fit$theta) == 0L)) {
    stop("criterion = 'score' requires a converged base model with fitted ",
         "coefficients; this fit did not converge. Supply theta starting ",
         "values to hazard(), or use criterion = 'wald'.",
         call. = FALSE)
  }

  ts_start <- Sys.time()
  call <- match.call()

  extra_args <- list(...)

  steps     <- list()
  trace_msg <- character()

  emit <- function(msg) {
    trace_msg[[length(trace_msg) + 1L]] <<- msg
    if (isTRUE(trace)) cat(msg, "\n", sep = "")
  }

  # The score criterion is entry-only: it exists to remove the per-candidate
  # refit from the forward step, and the drop path never had one. Following
  # SAS, removal is tested on the current model's Wald p-value.
  drop_criterion <- if (criterion == "score") "wald" else criterion

  # Header line (mirrors design sec.5).
  header <- if (criterion == "aic") {
    sprintf(
      "Stepwise selection (direction = %s, criterion = aic)",
      direction
    )
  } else {
    sprintf(
      "Stepwise selection (direction = %s, criterion = %s, slentry = %.2f, slstay = %.2f)",
      direction, criterion, slentry, slstay
    )
  }
  emit(header)
  emit("")

  # Move counter: per-variable tally of entries + exits.  Use a named
  # list rather than a named integer vector -- `lst[[missing]]` returns
  # NULL, whereas `vec[[missing]]` errors with "subscript out of bounds".
  move_counts <- list()
  frozen      <- character()

  current <- fit
  step_no <- 0L
  # The log-likelihood the next step starts from. A forward step produces a
  # model that CONTAINS the current one, so at the optimum the objective
  # cannot fall. When it does, the refit did not converge -- provable in one
  # comparison, which nothing was doing: the objective was written at every
  # step and read at none.
  prev_objective <- current$fit$objective %||% NA_real_
  n_nonmonotone_entries <- 0L
  stopped_by_max_steps <- FALSE
  # Candidates whose score statistic could not be computed, summed over
  # steps.  Tracked so a screen that stopped because nothing was
  # computable is distinguishable from one that stopped because nothing
  # was good enough -- the two produce identical empty steps otherwise.
  n_uncomputable_scores <- 0L
  uncomputable_reasons  <- stats::setNames(integer(0), character(0))
  stopped_uncomputable  <- FALSE

  # `crit` is the criterion actually applied to THIS step, which is not always
  # the run's `criterion`: score is entry-only, so its drops are decided by
  # Wald and must be labelled as such.
  record_step <- function(action, out, crit = criterion) {
    step_no <<- step_no + 1L
    row <- data.frame(
      step_num  = step_no,
      action    = action,
      variable  = out$variable,
      phase     = out$phase,
      criterion = crit,
      score     = out$score,
      stat      = out$stat,
      df        = out$df,
      p_value   = out$p_value,
      delta_aic = out$delta_aic,
      logLik    = current$fit$objective   %||% NA_real_,
      delta_logLik = (current$fit$objective %||% NA_real_) - prev_objective,
      aic       = .hzr_aic(current),
      n_coef    = length(current$fit$theta),
      stringsAsFactors = FALSE
    )
    prev_objective <<- current$fit$objective %||% NA_real_
    steps[[length(steps) + 1L]] <<- row

    score_fmt <- if (criterion == "aic") {
      sprintf("\u0394AIC = %+.2f", out$delta_aic)
    } else {
      sprintf("p = %.3f", out$p_value)
    }
    phase_txt <- if (is.na(out$phase)) {
      ""
    } else if (action == "enter") {
      paste0("  into  ", out$phase)
    } else {
      paste0("  from  ", out$phase)
    }
    emit(sprintf(
      "Step %d: %-6s %s%s   (%s)",
      step_no, toupper(action), out$variable, phase_txt, score_fmt
    ))
  }

  record_freeze <- function(var, phase_hint = NA_character_) {
    step_no <<- step_no + 1L
    row <- data.frame(
      step_num  = step_no,
      action    = "frozen",
      variable  = var,
      phase     = phase_hint,
      criterion = criterion,
      score     = NA_real_,
      stat      = NA_real_,
      df        = NA_integer_,
      p_value   = NA_real_,
      delta_aic = NA_real_,
      logLik    = current$fit$objective   %||% NA_real_,
      delta_logLik = 0,   # freezing changes no parameter
      aic       = .hzr_aic(current),
      n_coef    = length(current$fit$theta),
      stringsAsFactors = FALSE
    )
    steps[[length(steps) + 1L]] <<- row
    emit(sprintf(
      "Step %d: FROZEN %s   (exceeded max_move = %d; OSCILLATING)",
      step_no, var, max_move
    ))
  }

  bump_move <- function(var) {
    move_counts[[var]] <<- (move_counts[[var]] %||% 0L) + 1L
    if (move_counts[[var]] > max_move && !var %in% frozen) {
      frozen <<- c(frozen, var)
      record_freeze(var)
    }
  }

  # Main loop
  repeat {
    if (step_no >= max_steps) {
      warning("Stepwise selection hit max_steps = ", max_steps,
              "; stopping early.", call. = FALSE)
      stopped_by_max_steps <- TRUE
      break
    }

    add_happened  <- FALSE
    drop_happened <- FALSE

    effective_force_out <- unique(c(force_out, frozen))
    effective_force_in  <- unique(c(force_in,  frozen))

    if (direction %in% c("forward", "both")) {
      fwd <- do.call(.hzr_stepwise_forward_step, c(list(
        current   = current,
        scope     = scope,
        data      = data,
        criterion = criterion,
        slentry   = slentry,
        force_out = effective_force_out
      ), extra_args))

      n_uncomputable_scores <- n_uncomputable_scores +
        (fwd$n_uncomputable %||% 0L)
      uncomputable_reasons <- .hzr_merge_reasons(
        uncomputable_reasons, fwd$uncomputable_reasons
      )
      if (identical(fwd$stop_reason, "scores_uncomputable")) {
        stopped_uncomputable <- TRUE
      }

      if (fwd$accepted) {
        # Nested models: the entered model contains the current one, so at the
        # optimum objective_new >= objective_old. A violation is not a
        # statistical result, it is proof the refit failed -- and every later
        # step is then scored against a model that is not at its own optimum.
        # The tolerance keeps optimizer noise from firing this; the cases that
        # matter are whole log-likelihood units, not 1e-10.
        entered_objective <- fwd$fit$fit$objective %||% NA_real_
        obj_tol <- 1e-8 * max(1, abs(prev_objective))
        if (is.finite(entered_objective) && is.finite(prev_objective) &&
              entered_objective < prev_objective - obj_tol) {
          n_nonmonotone_entries <- n_nonmonotone_entries + 1L
          warning("Stepwise forward step ", step_no + 1L, " entered ",
                  fwd$variable,
                  if (!is.na(fwd$phase)) paste0(" (", fwd$phase, ")") else "",
                  " and the log-likelihood FELL, ", format(prev_objective),
                  " -> ", format(entered_objective), " (",
                  format(entered_objective - prev_objective), "). The entered ",
                  "model contains the current one, so this cannot happen at ",
                  "the optimum: the refit did not converge, and every later ",
                  "step is scored against a model that is not at its optimum. ",
                  "See `$steps$delta_logLik`.", call. = FALSE)
        }
        current <- fwd$fit
        record_step("enter", fwd)
        bump_move(fwd$variable)
        add_happened <- TRUE
      }
    }

    if (direction %in% c("backward", "both")) {
      bwd <- do.call(.hzr_stepwise_backward_step, c(list(
        current   = current,
        data      = data,
        criterion = drop_criterion,
        slstay    = slstay,
        force_in  = effective_force_in
      ), extra_args))

      if (bwd$accepted) {
        current <- bwd$fit
        record_step("drop", bwd, crit = drop_criterion)
        bump_move(bwd$variable)
        drop_happened <- TRUE
      }
    }

    if (!add_happened && !drop_happened) {
      emit(sprintf("(no further action after %d step%s)",
                   step_no, if (step_no == 1L) "" else "s"))
      break
    }
  }

  elapsed <- difftime(Sys.time(), ts_start, units = "secs")

  emit("")
  emit(sprintf("Final model: %d covariate%s, logLik = %.2f, AIC = %.2f",
               max(0L, length(current$fit$theta) -
                     .hzr_stepwise_shape_count(current)),
               if (length(current$fit$theta) -
                     .hzr_stepwise_shape_count(current) == 1L) "" else "s",
               current$fit$objective %||% NA_real_,
               .hzr_aic(current)))

  steps_df <- if (length(steps) == 0L) {
    data.frame(
      step_num = integer(), action = character(),
      variable = character(), phase = character(),
      criterion = character(), score = numeric(),
      stat = numeric(), df = integer(),
      p_value = numeric(), delta_aic = numeric(),
      logLik = numeric(), delta_logLik = numeric(),
      aic = numeric(), n_coef = integer(),
      stringsAsFactors = FALSE
    )
  } else {
    do.call(rbind, steps)
  }

  result <- current
  result$steps      <- steps_df
  result$scope      <- list(
    candidates = scope,
    force_in   = force_in,
    force_out  = force_out,
    frozen     = frozen
  )
  result$criteria   <- list(
    direction = direction,
    criterion = criterion,
    slentry   = slentry,
    slstay    = slstay,
    max_steps = max_steps,
    max_move  = max_move,
    hit_max_steps = stopped_by_max_steps,
    n_uncomputable_scores = n_uncomputable_scores,
    uncomputable_reasons  = uncomputable_reasons,
    stopped_uncomputable  = stopped_uncomputable,
    n_nonmonotone_entries = n_nonmonotone_entries
  )

  n_indefinite <- unname(uncomputable_reasons["information_indefinite"])
  if (is.na(n_indefinite)) n_indefinite <- 0L

  if (stopped_uncomputable) {
    warning("Stepwise selection stopped because the score statistic ",
            "could not be computed for any remaining candidate (",
            n_uncomputable_scores, " candidate score(s) were NA across the ",
            "run). This is not the same as no candidate meeting `slentry`: ",
            "the screen stopped without being able to test them.",
            .hzr_format_reasons(uncomputable_reasons), call. = FALSE)
  } else if (n_indefinite > 0L) {
    # The run finished normally, so the branch above stays quiet -- but a
    # candidate the score test could not evaluate at beta = 0 is usually one
    # with a LARGE effect, and it was passed over in favour of candidates that
    # could be scored.  A completed run is where that is least visible and
    # most misleading, so it warns on its own.
    warning("Stepwise selection completed, but ", n_indefinite,
            " candidate score(s) could not be computed because ",
            .hzr_score_reason_text("information_indefinite"),
            ". Those candidates were passed over rather than tested, so the ",
            "selected set may omit strong variables. Re-run with ",
            "`criterion = \"wald\"` to test them; see ",
            "`$criteria$uncomputable_reasons`.", call. = FALSE)
  }
  result$trace_msg  <- trace_msg
  result$elapsed    <- elapsed
  result$final_call <- call

  class(result) <- unique(c("hzr_stepwise", class(result)))
  result
}


#' @rdname hzr_stepwise
#' @param x An `hzr_stepwise` object.
#' @param ... Unused.
#' @return `print.hzr_stepwise` returns `x` invisibly.
#' @export
print.hzr_stepwise <- function(x, ...) {
  cat(paste(x$trace_msg, collapse = "\n"), "\n", sep = "")
  invisible(x)
}


#' @rdname hzr_stepwise
#' @param object An `hzr_stepwise` object.
#' @return `summary.hzr_stepwise` returns a `summary.hzr_stepwise` object
#'   (extends `summary.hazard`) with `$stepwise_steps` and `$stepwise_trace`
#'   appended.
#' @export
summary.hzr_stepwise <- function(object, ...) {
  # Strip the stepwise class so NextMethod dispatches cleanly to
  # summary.hazard.
  class(object) <- setdiff(class(object), "hzr_stepwise")
  out <- NextMethod()
  out$stepwise_steps <- object$steps
  out$stepwise_trace <- object$trace_msg
  class(out) <- unique(c("summary.hzr_stepwise", class(out)))
  out
}


#' @rdname hzr_stepwise
#' @return `print.summary.hzr_stepwise` returns `x` invisibly.
#' @export
print.summary.hzr_stepwise <- function(x, ...) {
  if (!is.null(x$stepwise_trace)) {
    cat(paste(x$stepwise_trace, collapse = "\n"), "\n\n", sep = "")
  }
  class(x) <- setdiff(class(x), "summary.hzr_stepwise")
  NextMethod()
}


#' @rdname hzr_stepwise
#' @return `as.data.frame.hzr_stepwise` returns the `$steps` data frame.
#' @export
as.data.frame.hzr_stepwise <- function(x, ...) {
  x$steps
}


#' Extract the captured console trace from an `hzr_stepwise` fit
#'
#' Every run of [hzr_stepwise()] records the header, per-step lines,
#' and final summary regardless of the `trace` flag.  This accessor
#' returns the full character vector for display or logging.
#'
#' @param fit An `hzr_stepwise` object.
#' @return Character vector, one element per console line.
#' @seealso [hzr_stepwise()], which produces the object this accessor reads.
#' @examples
#' data(avc)
#' avc <- na.omit(avc)
#' base <- hazard(survival::Surv(int_dead, dead) ~ age,
#'                data = avc, dist = "weibull", fit = TRUE,
#'                theta = c(mu = 0.01, nu = 0.5, 0))
#' \donttest{
#' sw <- hzr_stepwise(base, scope = ~ age + mal,
#'                    data = avc, direction = "forward",
#'                    control = list(n_starts = 1))
#' cat(stepwise_trace(sw), sep = "\n")
#' }
#' @export
stepwise_trace <- function(fit) {
  if (!inherits(fit, "hzr_stepwise")) {
    stop("`fit` must be an `hzr_stepwise` object.", call. = FALSE)
  }
  fit$trace_msg
}


# Local shape-count helper that tolerates multiphase as well as single
# distributions (the shared `.hzr_shape_parameter_count` in
# parity-helpers only handles the latter).
.hzr_stepwise_shape_count <- function(fit) {
  if (fit$spec$dist == "multiphase") {
    # For multiphase, "coefficient count" is sum of betas across phases
    # -- i.e. theta length minus all non-beta slots.  Easiest path:
    # count columns in x_list.
    x_list <- fit$fit$x_list
    if (is.null(x_list)) return(length(fit$fit$theta))
    total_betas <- sum(vapply(x_list, function(m) {
      if (is.null(m)) 0L else ncol(m)
    }, integer(1L)))
    length(fit$fit$theta) - total_betas
  } else {
    .hzr_shape_parameter_count(fit$spec$dist,
                                control = fit$spec$control)
  }
}
