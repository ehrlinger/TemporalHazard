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
#'     *add* variables; the best eligible candidate enters each step
#'     until none clears the entry rule.  Variables never leave once in.}
#'   \item{`direction = "backward"`}{Start from the base model, which must
#'     already hold every candidate, and only *drop* variables; the weakest
#'     term leaves each step until all survivors clear the retention rule.
#'     `scope` is not read, so a non-empty one is an error: protect terms
#'     with `force_in`.}
#'   \item{`direction = "both"` (default)}{Two-way stepwise: on every
#'     iteration, whether or not a variable entered, every term in the model,
#'     the base model's included, is re-tested and may be dropped unless it
#'     is in `force_in` or was frozen by `max_move` before the iteration
#'     began; a variable frozen on entry can still be dropped in the same
#'     iteration (see the **Known limitation (the frozen set)** section).
#'     `scope` limits what may enter, not what may leave.
#'     This is the SAS `SELECTION = STEPWISE` strategy.  `max_move` caps how
#'     often a single variable may oscillate before it is frozen.}
#' }
#'
#' \describe{
#'   \item{`criterion = "score"` (default)}{Accept moves on SAS-style
#'     significance thresholds, using the score (Q) statistic of the candidate
#'     coefficient; this reproduces C/SAS HAZARD's `SELECTION` statistic.
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
#'   shaping-parameter covariances are ignored.  This affects selection only;
#'   final-model standard errors are unchanged and still come from the
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
#'     take different step paths (and select different variable sets)
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
#'   fits, pass a named list of one-sided formulas keyed by phase, naming
#'   each phase once.  `scope` lists what may enter; a drop considers every
#'   term in the model except `force_in` and terms frozen by `max_move`
#'   before the iteration began (see the **Known limitation (the frozen
#'   set)** section).  A two-sided formula is an error,
#'   since its left-hand side would never be a candidate, and so is a
#'   non-empty `scope` under `direction = "backward"`, which does not read
#'   it.  An empty scope (`~ 1`, `character()`, or a list of `NULL`s and
#'   `~ 1`s) is accepted there.
#' @param data Data frame the base fit was built on.  Required for
#'   refits.
#' @param direction Search strategy: one of `"both"` (default),
#'   `"forward"`, or `"backward"`.  Controls whether variables may only
#'   enter, only leave, or both.  See the **Selection direction and
#'   criterion** section.
#' @param criterion Entry / retention rule: one of `"score"` (default),
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
#'   remainder of the run.  Default `4`.  In a two-way screen
#'   (`direction = "both"`), a variable frozen on entry can still be
#'   dropped in the same iteration; see the
#'   **Known limitation (the frozen set)** section.
#' @param force_in Character vector of variables that must remain in
#'   the model.  Such variables are still scored and reported in the
#'   selection trace, but are never dropped.
#' @param force_out Character vector of variables that may never be
#'   considered as candidates.
#' @param trace Logical; print step-by-step progress to the console.
#'   Default `TRUE`.
#' @param ... Passed to every candidate refit. Only `control` (e.g.
#'   `control = list(maxit = 500)`) and an `objective` equal to the base
#'   fit's are accepted. Any other name is an error: a misspelling such as
#'   `slentyr` would be stored by `hazard()` without being read, and every
#'   other `hazard()` argument (the response, data, `weights`,
#'   `time_windows`, `dist`, `theta`, `phases`, `fit`) is taken from the base
#'   model, so that each candidate is compared with the model it extends.
#'   The `print()`, `summary()` and `as.data.frame()` methods ignore `...`.
#'
#' @return An object of class `c("hzr_stepwise", "hazard")`, the
#'   final fit augmented with:
#'   \describe{
#'     \item{\code{steps}}{Data frame with one row per accepted /
#'       frozen action; see Details.}
#'     \item{\code{scope}}{Record of the candidate scope, plus
#'       `force_in`, `force_out`, and the frozen set.  In a two-way
#'       screen, `frozen` can name a variable the final model does not
#'       contain; see the **Known limitation (the frozen set)** section.}
#'     \item{\code{criteria}}{Named list of the threshold / direction
#'       settings actually applied, plus
#'       `n_uncomputable_scores` (how many candidate scores were `NA`,
#'       counted once per step: an entry the score test could not score
#'       under `criterion = "score"`; an entry under `criterion = "wald"`, or
#'       a removal under any criterion, whose Wald statistic could not be
#'       computed for want of a variance, `wald_no_variance`; or an entry
#'       under `criterion = "aic"` whose fit had no finite objective,
#'       `nonfinite`),
#'       `uncomputable_reasons` (a named integer vector of *why*),
#'       `wald_untested_removals` and `wald_untested_entries` (the
#'       `"var"` / `"var@phase"` tokens of variables kept in, or left out,
#'       on a step whose Wald test for them could not be computed; a
#'       variable tested at a later step is not listed.  Entries are listed
#'       under `criterion = "wald"` only: under `"score"` an entry no test
#'       could reach is reported by its reason, such as
#'       `fallback_no_variance`) and
#'       `stopped_uncomputable` (`TRUE` when the last iteration had
#'       candidates for entry or for removal and could test none of them).
#'       Read `uncomputable_reasons` before treating
#'       an unscored candidate as a bad one: `information_indefinite` marks
#'       candidates whose effect is too large for the score test's
#'       approximation at zero, which are typically the strongest variables
#'       on offer rather than degenerate ones. Candidates with that cause,
#'       or with `coefficient_diverging`, are refit and tested by Wald
#'       automatically, counted in `n_wald_fallbacks`. A candidate still
#'       reaches `uncomputable_reasons` when that refit fails, or when its
#'       cause is any other, which no refit can rescue. Read
#'       `uncomputable_reasons` for which one it was in any given run.  For
#'       every criterion it also carries
#'       `refit_failures` (the `"var"` / `"var@phase"` tokens of candidate
#'       moves whose refit errored, failed to converge, or was refused
#'       because the move would not change the model -- a drop that removes
#'       no design column, #320),
#'       `refit_failure_reasons` (why each one failed or was refused: the
#'       refit's error message, that it did not converge, or that the move
#'       changes nothing; named by the same tokens),
#'       `n_refit_failures`,
#'       and `stopped_refit_failed` (`TRUE` when the run ended on an
#'       iteration in which a refit failed or a move was refused.  A refit
#'       failure is a screen that could not test its candidates, rather than
#'       one that tested them and liked none; a refusal is determinate -- the
#'       move was tested and would have left the model as it was.  Read
#'       `refit_failure_reasons` for which it was).
#'       Check it before reading a zero-row `steps` as an honest null result.}
#'     \item{\code{trace_msg}}{Character vector of the trace lines,
#'       captured regardless of the `trace` flag.}
#'     \item{\code{elapsed}}{`difftime` from start to finish.}
#'     \item{\code{final_call}}{The call that produced this result.}
#'   }
#'
#' @section Known limitation (the frozen set):
#'
#' In a two-way screen (`direction = "both"`), `$scope$frozen` can name a
#' variable that the final model does not contain.  Each two-way iteration
#' makes a forward step and then a backward step, and the sets of protected
#' variables are fixed at the start of the iteration.  A variable that the
#' forward step freezes can therefore still be dropped by the backward step
#' that follows it, so it is reported as frozen while the final model
#' excludes it.  Nothing warns when this happens.
#'
#' Forward-only and backward-only screens are not affected: a forward-only
#' screen never makes a backward step to drop the frozen variable, and a
#' backward-only screen never makes a forward step to freeze it on entry.
#' In those, `$scope$frozen` and the final model agree.
#'
#' **When the two disagree, trust the final model and `$steps`.**  The
#' final model is what was selected, and `$steps` records both the
#' `"frozen"` row and the `"drop"` that followed it.  Read `$scope$frozen`
#' only as the list of variables that reached the `max_move` cap, not as a
#' list of variables held in the model.
#'
#' This is a known limitation of this release.  The analysis, including why
#' fixing it changes which variables are selected, is in
#' \url{https://github.com/ehrlinger/TemporalHazard/issues/378} and
#' \url{https://github.com/ehrlinger/TemporalHazard/issues/379}.
#'
#' @details
#' The `steps` data frame has columns:
#'
#' \describe{
#'   \item{\code{step_num}}{Integer sequence starting at 1.}
#'   \item{\code{action}}{`"enter"`, `"drop"`, or `"frozen"`.}
#'   \item{\code{variable}}{Variable affected.}
#'   \item{\code{phase}}{Phase name (multiphase) or `NA_character_`.}
#'   \item{\code{criterion}}{The criterion actually applied to this step:
#'     `"score"`, `"wald"`, or `"aic"`.  Under `criterion = "score"` the drop
#'     rows read `"wald"`, because score is entry-only.}
#'   \item{\code{score}}{Winning score used for the decision.}
#'   \item{\code{stat}, \code{df}}{Test statistic and degrees of freedom.}
#'   \item{\code{stat_type}}{What \code{stat} is on this row, and so which
#'     reference distribution recomputes its p-value: \code{"score_q"}
#'     (chi-square on \code{df}), \code{"wald_z"} (standard normal) or
#'     \code{"wald_chisq"} (chi-square on \code{df}).  \code{df} alone does
#'     not distinguish them; a scalar Wald is reported as a \emph{z}, not
#'     as its square, so it and a score Q are both recorded at \code{df = 1}
#'     while calling for different distributions.  It also identifies the
#'     rows the Wald fallback rescued, but only among \emph{entry} rows:
#'     under \code{criterion = "score"} those are the rows with
#'     \code{action == "enter" & stat_type == "wald_z"}.  Drop rows are
#'     always Wald-tested under that criterion (removal follows SAS and
#'     is tested on the current model's Wald p-value), so they read
#'     \code{"wald_z"} whether or not the fallback ever fired.  See
#'     \code{$criteria$n_wald_fallbacks} for the count.}
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
  extra_args <- .hzr_check_forwarded_dots(list(...), "hzr_stepwise",
                                          own = names(formals(hzr_stepwise)),
                                          fit = fit)
  if (missing(data) || !is.data.frame(data)) {
    stop("`data` must be a data frame (typically the frame used for the base fit).",
         call. = FALSE)
  }
  .hzr_refuse_unhonoured_scope(scope, direction)

  # Every accepted step goes through .hzr_refit_with_scope(), so a base fit
  # it cannot refit makes the entire screen a no-op.  Left to fail
  # candidate-by-candidate that produces N warnings and a zero-step result
  # indistinguishable from an honest "nothing met slentry" (#159).  Ask the
  # refit's own predicate once, up front: one message, naming the remedy,
  # before any fitting happens.  Sharing .hzr_refit_blocker() with the refit
  # is what stops the two answers drifting apart.
  # A forward screen steps only the phases its scope gives candidates; a
  # NULL entry offers none, so its phase is never refit. A backward or
  # two-way screen can drop from any phase, so it steps them all (#284).
  stepped <- if (identical(fit$spec$dist, "multiphase") &&
                   direction == "forward" && is.list(scope) &&
                   !is.null(names(scope))) {
    names(scope)[!vapply(scope, is.null, logical(1))]
  }
  refit_blocker <- .hzr_refit_blocker(fit, stepped = stepped)
  if (!is.null(refit_blocker)) {
    # The remedy is distribution-specific. Telling a multiphase caller to
    # "rebuild with the formula interface" would be wrong twice over: it is
    # not what blocked them, and for a translated SAS job it would force the
    # -1/0/1/2 -> Surv() 0/1/2/3 status round-trip this package has already
    # shipped a wrong answer through.
    remedy <- if (!is.null(.hzr_inherit_blocker(fit, stepped = stepped))) {
      # An inherited-design refusal names its own remedy; the generic one
      # below would ask for `time` and `status` the caller already gave.
      ""
    } else if (identical(fit$spec$dist, "multiphase")) {
      paste0("Refit the base model with hazard(), supplying `time` and ",
             "`status` (or a formula), and retry.")
    } else {
      paste0("Rebuild the base fit with the formula interface -- for a ",
             "forward screen that usually means an intercept-only base, ",
             "hazard(Surv(time, status) ~ 1, data = df, ...) -- and retry.")
    }
    stop("`fit` cannot be used as a stepwise base model because ",
         refit_blocker, ". No candidate could be refit as the step ",
         "describes, so the screen cannot run.",
         if (nzchar(remedy)) paste0(" ", remedy),
         call. = FALSE)
  }

  # Read the base model's terms once, before any output, so a `~ .` base
  # formula stops here with its remedy rather than after the header (#279).
  # A multiphase screen is exempt only when every phase has its own formula,
  # which hazard() has already written out (#277). A phase without one
  # inherits the global design, and each candidate refit would re-expand a
  # global `.` against the screen's data.
  own_formulas <- identical(fit$spec$dist, "multiphase") &&
    all(vapply(fit$spec$phases, function(ph) !is.null(ph$formula),
               logical(1)))
  if (!own_formulas) {
    .hzr_formula_rhs_terms(.hzr_stored_formula(fit))
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
  n_wald_fallbacks      <- 0L
  # Variables whose Wald test for entry or removal could not be computed
  # (#389), by "var" / "var@phase" token: an untested removal stays in the
  # model and an untested entry stays out, each as if it had been tested.
  wald_untested_entries  <- character()
  wald_untested_removals <- character()
  # A variable's latest step decides: one untested at step 1 and tested at
  # step 3 was tested.
  wald_tokens <- function(scores, keep) {
    scores <- scores[keep, , drop = FALSE]
    paste0(scores$variable,
           ifelse(is.na(scores$phase), "", paste0("@", scores$phase)))
  }
  update_untested <- function(set, scores, untested) {
    setdiff(union(set, wald_tokens(scores, untested)),
            wald_tokens(scores, !is.na(scores$score)))
  }
  stopped_uncomputable  <- FALSE
  stopped_untestable    <- character()
  # Candidates whose REFIT failed, summed over steps.  The per-candidate
  # warning already fires inside the step, but nothing recorded it on the
  # result, so a screen that could not fit any candidate returned the same
  # empty object as one that fit them all and liked none (#159).  The step
  # functions have returned `refit_failures` since v1; this is the caller
  # finally reading it.
  refit_failures       <- character()
  refit_failure_reasons <- character()
  stopped_refit_failed <- FALSE

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
      stat_type = out$stat_type %||% NA_character_,
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
      stat_type = NA_character_,
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
    # Per-ITERATION, not per-run: what stopped the screen is what the last
    # iteration did, and a failure three steps back is not why it ended.
    iter_refit_failures <- character()
    iter_refit_reasons  <- character()
    # Which half of this iteration had candidates and could test none:
    # "entry", "removal" or both.  Decided per iteration, so a two-way
    # screen that recovers at a later iteration is not reported as stopped.
    iter_untestable     <- character()

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
      n_wald_fallbacks <- n_wald_fallbacks + (fwd$n_wald_fallbacks %||% 0L)
      uncomputable_reasons <- .hzr_merge_reasons(
        uncomputable_reasons, fwd$uncomputable_reasons
      )
      if (identical(fwd$stop_reason, "scores_uncomputable")) {
        iter_untestable <- c(iter_untestable, "entry")
      }
      if (criterion == "wald" && nrow(fwd$all_scores) > 0L) {
        wald_untested_entries <- setdiff(update_untested(
          wald_untested_entries, fwd$all_scores, is.na(fwd$all_scores$score)
        ), fwd$refit_failures %||% character())
      }
      iter_refit_failures <- c(iter_refit_failures,
                               fwd$refit_failures %||% character())
      iter_refit_reasons <- c(iter_refit_reasons,
                              fwd$refit_failure_reasons %||% character())

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

      # A removal whose Wald p-value is NA was not tested, and it stays in
      # the model as if it met `slstay` (#389).  Counted with the forward
      # step's unscored entries.
      n_uncomputable_scores <- n_uncomputable_scores +
        (bwd$n_uncomputable %||% 0L)
      uncomputable_reasons <- .hzr_merge_reasons(
        uncomputable_reasons, bwd$uncomputable_reasons
      )
      if (nrow(bwd$all_scores) > 0L) {
        wald_untested_removals <- update_untested(
          wald_untested_removals, bwd$all_scores,
          !bwd$all_scores$force_in & is.na(bwd$all_scores$score)
        )
      }
      if (identical(bwd$stop_reason, "scores_uncomputable")) {
        iter_untestable <- c(iter_untestable, "removal")
      }
      iter_refit_failures <- c(iter_refit_failures,
                               bwd$refit_failures %||% character())
      iter_refit_reasons <- c(iter_refit_reasons,
                              bwd$refit_failure_reasons %||% character())

      if (bwd$accepted) {
        current <- bwd$fit
        record_step("drop", bwd, crit = drop_criterion)
        bump_move(bwd$variable)
        drop_happened <- TRUE
      }
    }

    refit_failures <- c(refit_failures, iter_refit_failures)
    refit_failure_reasons <- c(refit_failure_reasons, iter_refit_reasons)

    if (!add_happened && !drop_happened) {
      step_txt <- sprintf("%d step%s", step_no,
                          if (step_no == 1L) "" else "s")
      # "no further action" is a claim that candidates were tested and none
      # was good enough.  Say that only when it is true.
      if (length(iter_untestable) > 0L) {
        stopped_uncomputable <- TRUE
        stopped_untestable   <- iter_untestable
      }
      if (length(iter_refit_failures) > 0L) {
        stopped_refit_failed <- TRUE
        emit(sprintf(
          paste0("(stopped after %s: %d candidate move%s could not be ",
                 "completed -- a refit FAILED, or the move was REFUSED as ",
                 "changing nothing -- %s)"),
          step_txt, length(iter_refit_failures),
          if (length(iter_refit_failures) == 1L) "" else "s",
          paste(iter_refit_failures, collapse = ", ")
        ))
      } else if (length(iter_untestable) > 0L) {
        # Name the half: in a two-way screen the other half may have tested
        # its candidates and rejected them.
        emit(sprintf(
          paste0("(stopped after %s: no candidate score could be COMPUTED ",
                 "for %s -- none was tested)"),
          step_txt, paste(iter_untestable, collapse = " or ")
        ))
      } else {
        emit(sprintf("(no further action after %s)", step_txt))
      }
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
      stat = numeric(), stat_type = character(), df = integer(),
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
    n_wald_fallbacks      = n_wald_fallbacks,
    stopped_uncomputable  = stopped_uncomputable,
    wald_untested_removals = wald_untested_removals,
    wald_untested_entries  = wald_untested_entries,
    n_refit_failures      = length(refit_failures),
    refit_failures        = refit_failures,
    refit_failure_reasons = refit_failure_reasons,
    stopped_refit_failed  = stopped_refit_failed,
    n_nonmonotone_entries = n_nonmonotone_entries
  )

  # Both codes mean "no criterion tested this candidate". Keying the warning
  # below on information_indefinite alone would go silent the moment the
  # rescue's own failure was labelled separately -- quieter, for a case that
  # needs to be louder.
  untested_codes <- c("information_indefinite", "fallback_no_variance")
  n_indefinite <- sum(unname(uncomputable_reasons[untested_codes]),
                      na.rm = TRUE)


  if (stopped_uncomputable) {
    warning("Stepwise selection stopped because no remaining candidate ",
            "could be tested for ", paste(stopped_untestable, collapse = " or "),
            ": its score statistic, or its Wald statistic for want of a ",
            "variance, could not be computed (",
            n_uncomputable_scores, " candidate score(s) were NA across the ",
            "run). This is not the same as no candidate meeting `slentry` ",
            "or `slstay`: the screen stopped without being able to test them.",
            .hzr_format_reasons(uncomputable_reasons), call. = FALSE)
  } else if (n_indefinite > 0L) {
    # The run finished normally, so the branch above stays quiet -- but a
    # candidate the score test could not evaluate at beta = 0 is usually one
    # with a LARGE effect, and it was passed over in favour of candidates that
    # could be scored.  A completed run is where that is least visible and
    # most misleading, so it warns on its own.
    warning("Stepwise selection completed, but ", n_indefinite,
            " candidate(s) were tested by NEITHER criterion. The score ",
            "statistic could not be computed for them, and under ",
            "`criterion = \"score\"` they are then refit and Wald-tested ",
            "automatically -- so reaching this means that rescue did not ",
            "produce a test either: it errored or did not converge ",
            "(`information_indefinite`, listed in ",
            "`$criteria$refit_failures`), or it converged but yielded no ",
            "usable variance to test with (`fallback_no_variance`, which ",
            "leaves `refit_failures` empty). Such candidates are typically ",
            "STRONG -- that is what drives the score's information ",
            "indefinite -- so the selected set may omit them. Re-running ",
            "with `criterion = \"wald\"` runs the same refit and fails the ",
            "same way. See `$criteria$uncomputable_reasons` for which ",
            "mechanism applied.", call. = FALSE)
  }
  # Each of these variables was decided with no test at all, which reads
  # exactly like a test it failed (#389).  A half the stop warning above
  # already reports is left out here.
  warn_removals <- if ("removal" %in% stopped_untestable) character() else
    wald_untested_removals
  warn_entries <- if ("entry" %in% stopped_untestable) character() else
    wald_untested_entries
  if (length(c(warn_removals, warn_entries)) > 0L) {
    warning("Stepwise selection decided ",
            length(c(warn_removals, warn_entries)),
            " variable(s) without a Wald test",
            if (length(warn_removals)) {
              paste0("; kept in the model with its removal untested: ",
                     paste(warn_removals, collapse = ", "))
            },
            if (length(warn_entries)) {
              paste0("; left out with its entry untested: ",
                     paste(warn_entries, collapse = ", "))
            },
            ". Cause: ", .hzr_score_reason_text("wald_no_variance"),
            ". See `$criteria$uncomputable_reasons`.", call. = FALSE)
  }
  # A completed run keeps the tally but says nothing about it, and a collision
  # is a naming mistake the user can fix, not a property of the data. The
  # stopped_uncomputable warning above already spells the reason out.
  n_duplicate <- sum(uncomputable_reasons[names(uncomputable_reasons) ==
                                            "duplicate_column"])
  if (!stopped_uncomputable && n_duplicate > 0L) {
    warning("Stepwise selection declined ", n_duplicate, " candidate ",
            "score(s) without testing them: ",
            .hzr_score_reason_text("duplicate_column"), ". See ",
            "`$criteria$uncomputable_reasons`.", call. = FALSE)
  }
  if (stopped_refit_failed) {
    warning("Stepwise selection stopped after ", nrow(steps_df),
            " accepted step(s), and the iteration it stopped on had ",
            "candidate refit failure(s) (", length(refit_failures),
            " across the run: ",
            paste(unique(refit_failures), collapse = ", "),
            "). A candidate whose refit FAILED was never tested, so stopping ",
            "here is NOT evidence that nothing met the entry or retention ",
            "rule; a move that was REFUSED was tested and would not have ",
            "changed the model. See `$criteria$refit_failure_reasons` for ",
            "which each one was.", call. = FALSE)
  }
  result$trace_msg  <- trace_msg
  result$elapsed    <- elapsed
  result$final_call <- call

  class(result) <- unique(c("hzr_stepwise", class(result)))
  result
}

# Check the `...` of hzr_stepwise() / hzr_bootstrap() against what a
# candidate refit may take from it, and return it with abbreviations spelled
# out (#386). Both forward `...` to every refit, and hazard()'s own `...` is
# legacy pass-through that stores any name unread, so a misspelled `slentyr`
# for `slentry` was silently dropped and the screen ran at the default.
# Only `control`, and an `objective` equal to the base fit's, may pass. Every
# other hazard() argument describes the model the candidates are compared
# with, so the refit takes it from the base fit: a forwarded one was ignored
# (time_lower on a formula refit), collided with the refit's own (dist,
# failing every candidate), or changed the estimand of the candidates alone
# (weights, time_windows), so score and AIC differences compared two models
# fitted to different likelihoods. An abbreviation R would match to `control`
# or `objective` is kept, since hazard() used to apply it by partial
# matching. `extra` names the caller consumes itself (hzr_bootstrap()'s
# `trace`).
.hzr_check_forwarded_dots <- function(dots, caller, own, fit,
                                      extra = character()) {
  if (!length(dots)) return(dots)
  nms <- names(dots) %||% rep("", length(dots))
  if (any(!nzchar(nms))) {
    stop(caller, "(): `...` holds an unnamed argument. Everything in `...` ",
         "is forwarded by name to the hazard() refits, so name it.",
         call. = FALSE)
  }
  forwardable <- c("control", "objective")
  hz <- setdiff(names(formals(hazard)), "...")
  tick <- function(x) paste0("`", x, "`", collapse = ", ")
  full <- nms
  unknown <- from_base <- character()
  for (i in seq_along(nms)) {
    if (nms[i] %in% extra) next
    m <- pmatch(nms[i], hz)
    if (is.na(m)) {
      prefixed <- hz[startsWith(hz, nms[i])]
      if (length(prefixed) > 1L) {
        stop(caller, "(): `", nms[i], "` abbreviates more than one ",
             "hazard() argument (", tick(prefixed), "). Spell it out.",
             call. = FALSE)
      }
      unknown <- c(unknown, nms[i])
    } else if (hz[m] %in% forwardable) {
      full[i] <- hz[m]
    } else {
      from_base <- c(from_base, nms[i])
    }
  }
  if (length(unknown)) {
    known <- unique(c(setdiff(own, "..."), forwardable, extra))
    hint <- vapply(unknown, function(nm) {
      dist <- utils::adist(nm, known)[1L, ]
      if (min(dist) > 2L) return("")
      paste0(" Did you mean `", known[which.min(dist)], "`",
             if (length(unknown) > 1L) paste0(" for `", nm, "`") else "", "?")
    }, character(1))
    stop(caller, "(): ", tick(unknown),
         if (length(unknown) > 1L) " are not arguments" else
           " is not an argument",
         " of ", caller, "() or of hazard(). Everything in `...` is forwarded ",
         "to the hazard() refits, and hazard() stores a name it does not ",
         "declare without reading it, so it would have had no effect.",
         paste(hint, collapse = ""), call. = FALSE)
  }
  if (length(from_base)) {
    stop(caller, "(): ", tick(from_base), " cannot be passed through ",
         "`...`. Every candidate refit takes the response, data, weights, ",
         "time windows, dist, theta, phases and fit from the base model, so ",
         "that each candidate is compared with the model it would extend. ",
         "Of hazard()'s arguments, `...` forwards only `control` and an ",
         "`objective` equal to the base fit's. To change anything else, refit ",
         "the base model and screen from that.", call. = FALSE)
  }
  if (anyDuplicated(full[!full %in% extra])) {
    dup <- unique(full[duplicated(full) & !full %in% extra])
    stop(caller, "(): ", tick(dup), " is given more than once in `...`, ",
         "counting abbreviations.", call. = FALSE)
  }
  names(dots) <- full
  if ("control" %in% full && !is.list(dots$control)) {
    stop(caller, "(): `control` must be a list, as hazard() requires; ",
         "every candidate refit would have failed on it.", call. = FALSE)
  }
  # A refit may not change the estimand; refuse here rather than once per
  # candidate, where hzr_bootstrap() would tally it as a failed replicate.
  # match.arg() first, as hazard() does, so "lik" is "likelihood".
  if ("objective" %in% full) {
    fit_objective <- .hzr_fit_objective(fit)
    objective <- tryCatch(
      match.arg(dots$objective, eval(formals(hazard)$objective)),
      error = function(e) dots$objective
    )
    if (!identical(objective, fit_objective)) {
      stop(caller, "(): `objective = ", deparse1(dots$objective),
           "` differs from the base fit's objective = \"", fit_objective,
           "\", and a candidate refit cannot change the estimand. Refit the ",
           "base model under that objective and screen from there.",
           call. = FALSE)
    }
    dots$objective <- objective
  }
  dots
}


#' @rdname hzr_stepwise
#' @param x An `hzr_stepwise` object.
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
