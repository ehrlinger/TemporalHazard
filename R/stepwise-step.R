# stepwise-step.R -- Per-step drivers for hzr_stepwise()
#
# Step 8.4-8.5 of STEPWISE-DESIGN.md.  One forward step or one backward
# step, driven by the unified `.hzr_candidate_score()` wrapper.  Higher
# level control (two-way loop, MOVE cap, force_in, max_steps) lives in
# the hzr_stepwise() driver to be built in sec.8.6.

# ---------------------------------------------------------------------------
# Candidate enumeration
# ---------------------------------------------------------------------------

#' Normalise a user-supplied scope into a flat list of (var, phase) pairs
#'
#' Drops variables already in the model, variables in `force_out`, and
#' candidates whose phase is not present in a multiphase fit.
#'
#' @param fit The current fit.
#' @param scope One of: NULL (data-driven; infer from `data`), a
#'   one-sided formula (single-dist), or a named list of one-sided
#'   formulas keyed by phase (multiphase).
#' @param data Data frame to derive candidates from when `scope = NULL`.
#' @param force_out Character vector of variable names to exclude.
#'
#' @return A list of `list(var = <chr>, phase = <chr or NULL>)`.
#'
#' @keywords internal
#' @noRd
.hzr_stepwise_candidates <- function(fit, scope = NULL, data,
                                      force_out = character()) {
  dist <- fit$spec$dist

  if (dist == "multiphase") {
    phase_names <- names(fit$spec$phases)
    current_per_phase <- .hzr_scope_current_vars(fit)

    if (is.null(scope)) {
      # Default: all data-frame vars (excluding Surv components and
      # force_out) are candidates for every phase.
      stored <- .hzr_stored_formula(fit)
      lhs_vars <- if (is.null(stored)) character() else all.vars(stored[[2L]])
      data_vars <- setdiff(colnames(data), c(lhs_vars, force_out))
      data_vars <- .hzr_modellable_vars(data, data_vars)
      scope <- setNames(
        lapply(phase_names, function(p) {
          rhs_syms <- data_vars
          # Return as a one-sided formula for symmetry with the
          # user-supplied case.
          if (length(rhs_syms) == 0L) return(NULL)
          stats::as.formula(paste("~", paste(rhs_syms, collapse = " + ")))
        }),
        phase_names
      )
    } else {
      if (!is.list(scope) || is.null(names(scope))) {
        stop("`scope` for multiphase fits must be a named list keyed by phase.",
             call. = FALSE)
      }
      unknown <- setdiff(names(scope), phase_names)
      if (length(unknown) > 0L) {
        stop("Unknown phase(s) in scope: ",
             paste(sQuote(unknown), collapse = ", "),
             call. = FALSE)
      }
    }

    candidates <- list()
    for (p in names(scope)) {
      sc <- scope[[p]]
      if (is.null(sc)) next
      terms_p <- .hzr_formula_rhs_terms(sc)
      eligible <- setdiff(terms_p, c(current_per_phase[[p]], force_out))
      for (v in eligible) {
        candidates[[length(candidates) + 1L]] <-
          list(var = v, phase = p)
      }
    }
    return(candidates)
  }

  # Single-distribution
  if (is.null(scope)) {
    f <- .hzr_stored_formula(fit)
    lhs_vars <- if (is.null(f)) character() else all.vars(f[[2L]])
    data_vars <- setdiff(colnames(data), c(lhs_vars, force_out))
    data_vars <- .hzr_modellable_vars(data, data_vars)
  } else {
    if (inherits(scope, "formula")) {
      data_vars <- .hzr_formula_rhs_terms(scope)
    } else if (is.character(scope)) {
      data_vars <- scope
      # A character scope never passes through terms(), so an offset in it
      # reached a refit and failed there with a message that did not say why.
      scope_f <- tryCatch(stats::reformulate(scope), error = function(e) NULL)
      if (!is.null(scope_f)) .hzr_refuse_offset(scope_f, "`scope`")
    } else {
      stop("`scope` must be NULL, a one-sided formula, or a character vector.",
           call. = FALSE)
    }
    data_vars <- setdiff(data_vars, force_out)
  }

  current_vars <- .hzr_scope_current_vars(fit)
  eligible <- setdiff(data_vars, current_vars)
  lapply(eligible, function(v) list(var = v, phase = NULL))
}


# ---------------------------------------------------------------------------
# Forward step
# ---------------------------------------------------------------------------

#' Execute one forward step of stepwise selection
#'
#' Under `"wald"` / `"aic"`, refits the current model with each candidate
#' (variable, phase) appended and scores each via
#' `.hzr_candidate_score(mode = "entry")`.  Under `"score"` there is no
#' candidate refit at all: each candidate is scored by `.hzr_score_q()` with
#' its coefficient pinned at zero.  The reduced-model nuisance block is
#' inverted ONCE per step by `.hzr_score_nuisance()` and reused for every
#' candidate; that reuse is what removes the optimizer from the loop, and
#' letting `.hzr_score_q()` recompute it per candidate would silently give
#' most of the speedup back.  The best candidate is accepted if its score
#' clears the entry threshold.
#'
#' Divergent candidate refits emit a `warning()` naming the failing
#' `(variable, phase)` pair and are excluded from the selection, per
#' the sec.2 Q5 decision in STEPWISE-DESIGN.md.  The score path has no refit
#' to diverge, so it reports no refit failures.
#'
#' @param current Fitted `hazard` object that is the starting point.
#' @param scope Scope specification (see `.hzr_stepwise_candidates`).
#' @param data Data frame for refits.  Must be row-aligned with the fit
#'   under `criterion = "score"`, which reads candidate columns from it.
#' @param criterion One of `"score"`, `"wald"`, or `"aic"`.
#' @param slentry Entry threshold for the score / Wald criteria (ignored when
#'   `criterion = "aic"`; the entry rule there is dAIC < 0).
#' @param force_out Character vector of variables that may never be
#'   considered as candidates.
#' @param ... Forwarded to `.hzr_refit_with_scope()` (and thence to
#'   `hazard()`), e.g. `control = list(...)`.
#'
#' @return A list with:
#' \describe{
#'   \item{accepted}{Logical; did any candidate enter this step?}
#'   \item{fit}{If accepted, the new fit.  Otherwise `current` echoed.}
#'   \item{variable}{Variable that entered, or `NA_character_`.}
#'   \item{phase}{Phase entered, or `NA_character_` (single-dist).}
#'   \item{score}{Winning score (p or dAIC), or `NA_real_`.}
#'   \item{p_value}{Winning p-value.}
#'   \item{delta_aic}{Winning dAIC.}
#'   \item{stat, stat_type, df}{Test statistic of the winner, what it is
#'     (`"score_q"`, `"wald_z"` or `"wald_chisq"`), and its df.}
#'   \item{all_scores}{Tibble-like data frame of every candidate
#'     considered and its score.}
#'   \item{refit_failures}{Character vector of `"var@phase"` tokens for
#'     candidates whose refit diverged.}
#'   \item{refit_failure_reasons}{Why each of those refits failed: the
#'     refit's error message, or that it did not converge. Named by the same
#'     tokens, in the same order.}
#' }
#'
#' @keywords internal
#' @noRd
.hzr_stepwise_forward_step <- function(current, scope = NULL, data,
                                        criterion = c("score", "wald", "aic"),
                                        slentry   = 0.30,
                                        force_out = character(),
                                        ...) {
  criterion <- match.arg(criterion)
  if (!inherits(current, "hazard")) {
    stop("`current` must be a fitted `hazard` object.", call. = FALSE)
  }

  cands <- .hzr_stepwise_candidates(current, scope = scope, data = data,
                                     force_out = force_out)

  null_result <- function() {
    list(
      accepted  = FALSE,
      fit       = current,
      variable  = NA_character_,
      phase     = NA_character_,
      score     = NA_real_,
      p_value   = NA_real_,
      delta_aic = NA_real_,
      stat      = NA_real_,
      stat_type = NA_character_,
      df        = NA_integer_,
      all_scores = data.frame(
        variable  = character(),
        phase     = character(),
        score     = numeric(),
        p_value   = numeric(),
        delta_aic = numeric(),
        stat      = numeric(),
        stat_type = character(),
        df        = integer(),
        stringsAsFactors = FALSE
      ),
      refit_failures = character(),
      refit_failure_reasons = character()
    )
  }

  if (length(cands) == 0L) {
    return(null_result())
  }

  rows     <- vector("list", length(cands))
  failures <- character()
  failure_reasons <- character()

  if (criterion == "score") {
    return(.hzr_stepwise_forward_step_score(
      current = current, cands = cands, data = data,
      slentry = slentry, null_result = null_result, ...
    ))
  }

  for (i in seq_along(cands)) {
    cand <- cands[[i]]
    candidate_fit <- tryCatch(
      .hzr_refit_with_scope(
        current, action = "add",
        var = cand$var, phase = cand$phase,
        data = data, ...
      ),
      error = function(e) e
    )

    failure_token <- if (is.null(cand$phase)) {
      cand$var
    } else {
      paste0(cand$var, "@", cand$phase)
    }

    if (inherits(candidate_fit, "error") ||
          isFALSE(candidate_fit$fit$converged)) {
      reason <- .hzr_refit_failure_reason(candidate_fit)
      warning("Stepwise forward: candidate refit failed for ",
              failure_token, ": ", reason, call. = FALSE)
      failures <- c(failures, failure_token)
      failure_reasons <- c(failure_reasons,
                           stats::setNames(reason, failure_token))
      rows[[i]] <- data.frame(
        variable  = cand$var,
        phase     = cand$phase %||% NA_character_,
        score     = NA_real_,
        p_value   = NA_real_,
        delta_aic = NA_real_,
        stat      = NA_real_,
        stat_type = NA_character_,
        df        = NA_integer_,
        stringsAsFactors = FALSE
      )
      next
    }

    # Coefficient name of the newly-entered variable in the candidate fit,
    # resolved by the column the refit added rather than by `var`, which
    # another term's column can carry (a factor dummy named `flag` beside a
    # logical `flag`, whose own column is `flagTRUE`).
    coef_name <- .hzr_candidate_coef_name(candidate_fit, cand$var,
                                           cand$phase, current = current)

    s <- .hzr_candidate_score(
      criterion = criterion, mode = "entry",
      current = current, candidate = candidate_fit,
      names = coef_name
    )

    rows[[i]] <- data.frame(
      variable  = cand$var,
      phase     = cand$phase %||% NA_character_,
      score     = s$score,
      p_value   = s$p_value,
      delta_aic = s$delta_aic,
      stat      = s$stat,
      stat_type = s$stat_type,
      df        = s$df,
      stringsAsFactors = FALSE
    )
    # Cache the fit on the row so we can recover it for the winner
    attr(rows[[i]], "fit") <- candidate_fit
  }

  all_scores <- do.call(rbind, rows)
  # Strip the per-row fit attributes from the combined frame but keep
  # them in a parallel list keyed by row for winner lookup.
  candidate_fits <- lapply(rows, function(r) attr(r, "fit"))

  valid <- which(!is.na(all_scores$score))
  if (length(valid) == 0L) {
    out <- null_result()
    out$all_scores <- all_scores
    out$refit_failures <- failures
    out$refit_failure_reasons <- failure_reasons
    return(out)
  }

  best_idx <- valid[which.min(all_scores$score[valid])]
  best     <- all_scores[best_idx, ]

  threshold_met <- if (criterion == "wald") {
    best$score < slentry
  } else {
    best$score < 0
  }

  if (!threshold_met) {
    out <- null_result()
    out$all_scores <- all_scores
    out$refit_failures <- failures
    out$refit_failure_reasons <- failure_reasons
    return(out)
  }

  list(
    accepted  = TRUE,
    fit       = candidate_fits[[best_idx]],
    variable  = best$variable,
    phase     = best$phase,
    score     = best$score,
    p_value   = best$p_value,
    delta_aic = best$delta_aic,
    stat      = best$stat,
    stat_type = best$stat_type,
    df        = best$df,
    all_scores = all_scores,
    refit_failures = failures,
    refit_failure_reasons = failure_reasons
  )
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

#' Score-criterion body of the forward step
#'
#' Split out of `.hzr_stepwise_forward_step()` to keep the refit-based and
#' refit-free paths readable side by side.  Returns the same list shape.
#'
#' The winner is refit ONCE, after the decision, mirroring the backward
#' step, which has always worked this way.
#'
#' @param cands List of `list(var, phase)` pairs to score.
#' @param null_result Zero-arg constructor for the "nothing entered" result,
#'   supplied by the caller so both paths return an identical shape.
#'
#' @keywords internal
#' @noRd
.hzr_stepwise_forward_step_score <- function(current, cands, data, slentry,
                                             null_result, ...) {
  # The entire speedup: the reduced-model information does not depend on any
  # candidate, so invert its nuisance block once and hand it to every
  # candidate.  Letting .hzr_score_q() default `nuisance` to NULL would
  # re-invert it per candidate and give most of the benefit back silently.
  nuisance <- .hzr_score_nuisance(current)

  rows <- vector("list", length(cands))
  for (i in seq_along(cands)) {
    cand <- cands[[i]]
    .hzr_score_check_numeric(data, cand$var, cand$phase)

    s_q <- .hzr_score_q(current, cand$var, phase = cand$phase, data = data,
                        nuisance = nuisance)
    s <- .hzr_candidate_score(
      criterion = "score", mode = "entry",
      current = current, score = s_q,
      names = cand$var
    )

    rows[[i]] <- data.frame(
      variable  = cand$var,
      phase     = cand$phase %||% NA_character_,
      score     = s$score,
      p_value   = s$p_value,
      delta_aic = s$delta_aic,
      stat      = s$stat,
      stat_type = s$stat_type,
      df        = s$df,
      reason    = s_q$reason %||% NA_character_,
      stringsAsFactors = FALSE
    )
  }
  all_scores <- do.call(rbind, rows)

  # --- Wald fallback for candidates the score could not test (#130) --------
  # Q is SAS's exactly (src/vars/q1.c), and so is its blind spot: the observed
  # information at beta = 0 goes indefinite precisely when the candidate's
  # effect is LARGE, so the criterion declines candidates in proportion to how
  # predictive they are. Refit those few and test them by Wald rather than
  # dropping the screen's best variables. Only the reasons a refit can
  # actually rescue qualify -- see .hzr_score_fallback_reasons.
  all_scores$fallback <- FALSE
  fallback_failures <- character()
  fallback_reasons  <- character()
  # Keep each fallback's refit, keyed by row, the way the Wald path keeps its
  # candidate fits. The acceptance step below refits the winner with exactly
  # the same arguments, so without this a rescued candidate that goes on to
  # win is fitted twice -- against a criterion whose whole advantage is that
  # it does not refit per candidate.
  fallback_fits <- vector("list", nrow(all_scores))
  for (i in which(is.na(all_scores$score) &
                    all_scores$reason %in% .hzr_score_fallback_reasons)) {
    cand_phase <- if (is.na(all_scores$phase[i])) NULL else all_scores$phase[i]
    refit <- tryCatch(
      .hzr_refit_with_scope(current, action = "add",
                            var = all_scores$variable[i], phase = cand_phase,
                            data = data, ...),
      error = function(e) e
    )
    # A refit that fails or does not converge leaves the row NA with its
    # original reason, so it still counts as uncomputable below rather than
    # quietly becoming a candidate with no score.  Record and warn as every
    # other refit failure in the package does: without this the row is
    # byte-identical to one that was never refit at all, and the only signal
    # left is an uncomputable_reasons count that now means the opposite.
    if (is.null(refit) || inherits(refit, "error") ||
          isFALSE(refit$fit$converged)) {
      fallback_token <- if (is.null(cand_phase)) {
        all_scores$variable[i]
      } else {
        paste0(all_scores$variable[i], "@", cand_phase)
      }
      reason <- .hzr_refit_failure_reason(refit)
      warning("Stepwise forward: Wald-fallback refit failed for ",
              fallback_token, ": ", reason, call. = FALSE)
      fallback_failures <- c(fallback_failures, fallback_token)
      fallback_reasons  <- c(fallback_reasons,
                             stats::setNames(reason, fallback_token))
      next
    }
    w <- .hzr_candidate_score(
      criterion = "wald", mode = "entry", current = current, candidate = refit,
      names = .hzr_candidate_coef_name(refit, all_scores$variable[i],
                                       cand_phase, current = current)
    )
    if (is.na(w$score)) {
      # The refit CONVERGED -- it returned a point estimate -- but its Hessian
      # was not invertible, so there is no standard error and .hzr_wald_p()
      # cannot compute a test (it returns NA for "no vcov or non-finite SE").
      #
      # A bare `next` here recorded nothing at all. The row kept the SCORE's
      # reason, which describes the first of two independent failures and says
      # nothing about the second, and no counter moved -- so a strong
      # candidate that neither criterion could test vanished, and the screen
      # rendered as an honest "nothing met slentry". That indistinguishability
      # is the defect #159 and #130 were both about.
      all_scores$reason[i] <- "fallback_no_variance"
      next
    }
    fallback_fits[[i]]     <- refit
    all_scores$score[i]    <- w$score
    all_scores$p_value[i]  <- w$p_value
    all_scores$stat[i]     <- w$stat
    # `stat` is now a Wald z where every other row on this screen carries a
    # score Q. Leaving stat_type saying "score_q" would make a reader who
    # recomputes pchisq(stat, df) wrong by dozens of orders of magnitude.
    all_scores$stat_type[i] <- w$stat_type
    all_scores$df[i]       <- w$df
    all_scores$fallback[i] <- TRUE
    all_scores$reason[i]   <- NA_character_
  }
  n_wald_fallbacks <- sum(all_scores$fallback)

  valid <- which(!is.na(all_scores$score))

  # A candidate whose Q could not be computed scores NA and drops out of
  # `valid` silently.  Count them so the caller can tell a screen that finished
  # from one that could not run: both otherwise return the same empty step.
  #
  # Carry the REASONS up as well as the count.  The causes are not one thing:
  # a collinear column should be dropped, while an indefinite information
  # matrix usually marks a candidate whose effect is too large for the score
  # test at 0 -- a strong candidate, not a degenerate one.  A bare count
  # invites exactly the wrong reading of the second case.
  n_uncomputable <- sum(is.na(all_scores$score))
  uncomputable_reasons <- .hzr_tally_reasons(
    all_scores$reason[is.na(all_scores$score)]
  )

  if (length(valid) == 0L) {
    out <- null_result()
    out$all_scores     <- all_scores
    out$n_uncomputable <- n_uncomputable
    out$uncomputable_reasons <- uncomputable_reasons
    out$n_wald_fallbacks <- n_wald_fallbacks
    out$refit_failures <- fallback_failures
    out$refit_failure_reasons <- fallback_reasons
    out$stop_reason    <- if (n_uncomputable > 0L) {
      "scores_uncomputable"
    } else {
      "no_candidates"
    }
    return(out)
  }

  best_idx <- valid[which.min(all_scores$score[valid])]
  best     <- all_scores[best_idx, ]

  if (!(best$score < slentry)) {
    out <- null_result()
    out$all_scores     <- all_scores
    out$n_uncomputable <- n_uncomputable
    out$uncomputable_reasons <- uncomputable_reasons
    out$n_wald_fallbacks <- n_wald_fallbacks
    out$refit_failures <- fallback_failures
    out$refit_failure_reasons <- fallback_reasons
    out$stop_reason    <- "no_candidate_met_slentry"
    return(out)
  }

  best_phase <- if (is.na(best$phase)) NULL else best$phase
  # A rescued winner was already fitted, with identical arguments, in the
  # fallback loop above. Reuse it rather than paying for the same fit twice.
  refitted <- fallback_fits[[best_idx]]
  if (is.null(refitted)) {
    refitted <- tryCatch(
      .hzr_refit_with_scope(
        current, action = "add", var = best$variable, phase = best_phase,
        data = data, ...
      ),
      error = function(e) e
    )
  }

  failure_token <- if (is.na(best$phase)) {
    best$variable
  } else {
    paste0(best$variable, "@", best$phase)
  }

  if (inherits(refitted, "error") || isFALSE(refitted$fit$converged)) {
    # The candidate won on Q but the model that would realise it will not fit.
    # Entering it anyway would put a non-converged fit into the chain.
    reason <- .hzr_refit_failure_reason(refitted)
    warning("Stepwise forward: post-entry refit failed for ",
            failure_token, ": ", reason, call. = FALSE)
    out <- null_result()
    out$all_scores     <- all_scores
    out$refit_failures <- c(fallback_failures, failure_token)
    out$refit_failure_reasons <- c(fallback_reasons,
                                   stats::setNames(reason, failure_token))
    out$n_uncomputable <- n_uncomputable
    out$uncomputable_reasons <- uncomputable_reasons
    out$n_wald_fallbacks <- n_wald_fallbacks
    out$stop_reason    <- "refit_failed"
    return(out)
  }

  list(
    accepted  = TRUE,
    fit       = refitted,
    variable  = best$variable,
    phase     = best$phase,
    score     = best$score,
    p_value   = best$p_value,
    delta_aic = best$delta_aic,
    stat      = best$stat,
    stat_type = best$stat_type,
    df        = best$df,
    all_scores = all_scores,
    refit_failures = fallback_failures,
    refit_failure_reasons = fallback_reasons,
    n_uncomputable = n_uncomputable,
    uncomputable_reasons = uncomputable_reasons,
    n_wald_fallbacks = n_wald_fallbacks,
    stop_reason    = "accepted"
  )
}


#' Why a stepwise refit failed, for its warning and `refit_failure_reasons`
#'
#' The refit sites catch the error so one bad candidate cannot end the screen.
#' Catching it used to discard the message too, so a refit that hazard()
#' refused by name (a duplicated design column, a variable missing from
#' `data`) was reported as "refit failed for <var>" and nothing more.
#'
#' @param refit What the refit's `tryCatch()` returned.
#' @return A one-line character string.
#' @noRd
.hzr_refit_failure_reason <- function(refit) {
  if (inherits(refit, "error")) {
    return(conditionMessage(refit))
  }
  if (is.null(refit)) {
    return("the refit returned no fit")
  }
  "the refit did not converge"
}


#' Reject a candidate the score path cannot expand into a single column
#'
#' `.hzr_score_q()` returns a silent `NA` for a non-numeric or absent
#' candidate, which under `criterion = "score"` would mean the variable
#' simply never enters: no error, no warning, just absent from the
#' selected model.  The Wald path surfaces both cases loudly: a factor
#' candidate fails in `.hzr_candidate_coef_name()` (cannot find its expanded
#' coefficient by name), and a candidate missing from `data` fails inside
#' `.hzr_refit_with_scope()`'s refit, which emits the
#' "candidate refit failed for ..." warning. Matching that here keeps the
#' two criteria in agreement about what is selectable.
#'
#' A column absent from `data` warns (matching the Wald refit-failure
#' warning's severity) and is left to the caller's normal `NA` handling;
#' the run continues past it. A present-but-non-numeric column still
#' errors, unchanged.
#'
#' @keywords internal
#' @noRd
.hzr_score_check_numeric <- function(data, var, phase) {
  xcand <- data[[var]]
  if (is.null(xcand)) {
    where <- if (is.null(phase)) "" else paste0(" in phase ", sQuote(phase))
    warning(
      "Stepwise forward: candidate ", sQuote(var), where,
      " not found in `data`; skipping.",
      call. = FALSE
    )
    return(invisible(NULL))
  }
  if (!is.null(.hzr_candidate_numeric(xcand))) {
    return(invisible(NULL))
  }
  where <- if (is.null(phase)) "" else paste0(" in phase ", sQuote(phase))
  stop(
    "Stepwise forward: candidate ", sQuote(var), where, " is not numeric (",
    class(xcand)[1L], "). The score criterion tests the column as it stands ",
    "and cannot expand it. `criterion = \"wald\"` refits once per candidate, ",
    "so it does handle a term that expands to a single design-matrix column ",
    "-- a two-level factor or character column. For anything wider, expand it ",
    "into numeric main effects (e.g. a contrast matrix) and retry.",
    call. = FALSE
  )
}


# ---------------------------------------------------------------------------
# Backward step
# ---------------------------------------------------------------------------

#' Enumerate the currently-in-model (variable, phase) pairs eligible
#' for a drop test.
#'
#' Force-in variables are *included* in the returned list so their
#' scores appear in the trace (sec.2 Q4 of STEPWISE-DESIGN.md); the
#' forward/backward driver filters them out of the argmax pool via the
#' `force_in` flag on each row.
#'
#' @keywords internal
#' @noRd
.hzr_stepwise_drop_candidates <- function(fit) {
  dist <- fit$spec$dist
  if (dist == "multiphase") {
    per_phase <- .hzr_scope_current_vars(fit)
    out <- list()
    for (p in names(per_phase)) {
      for (v in per_phase[[p]]) {
        out[[length(out) + 1L]] <- list(var = v, phase = p)
      }
    }
    return(out)
  }
  vars <- .hzr_scope_current_vars(fit)
  lapply(vars, function(v) list(var = v, phase = NULL))
}


#' Execute one backward step of stepwise selection
#'
#' For each currently-in-model variable, scores drop via
#' `.hzr_candidate_score(mode = "drop")`.  Force-in variables are
#' scored and reported but excluded from the drop decision.
#'
#' Unlike the forward step there is no per-candidate refit: the Wald
#' statistic and its AIC approximation both operate on the current
#' model's vcov.  A single refit fires only after the drop decision
#' via `.hzr_refit_with_scope()`.
#'
#' @param current Fitted `hazard` object.
#' @param data Data frame used to rebuild the dropped model.
#' @param criterion Either `"wald"` or `"aic"`.
#' @param slstay Retention threshold for the Wald criterion (ignored
#'   when `criterion = "aic"`; the drop rule there is dAIC_drop < 0).
#' @param force_in Character vector of variables that may never be
#'   dropped.
#' @param ... Forwarded to `.hzr_refit_with_scope()` for the post-drop
#'   refit.
#'
#' @return A list with the same top-level shape as
#'   `.hzr_stepwise_forward_step()`, except:
#'   * `all_scores` gains a logical `force_in` column.
#'   * The action this represents is a drop, so `accepted = TRUE` means
#'     the variable was removed from the model.
#'   * `refit_failures` and `refit_failure_reasons` carry one more case than
#'     the forward step's: a drop this step REFUSED because it leaves the
#'     design no smaller -- an interaction whose main effect has gone is
#'     recoded, so the "reduced" model is the model it started from (#320).
#'     A multiphase drop is judged on the converged refit; a
#'     single-distribution drop on its reduced design, before any refit
#'     (#323).  The reason then says the drop removes no
#'     column, and `accepted` is `FALSE` with the current fit returned.
#'
#' @keywords internal
#' @noRd
.hzr_stepwise_backward_step <- function(current, data,
                                         criterion = c("wald", "aic"),
                                         slstay    = 0.20,
                                         force_in  = character(),
                                         ...) {
  criterion <- match.arg(criterion)
  if (!inherits(current, "hazard")) {
    stop("`current` must be a fitted `hazard` object.", call. = FALSE)
  }

  cands <- .hzr_stepwise_drop_candidates(current)

  empty_scores <- data.frame(
    variable  = character(),
    phase     = character(),
    force_in  = logical(),
    score     = numeric(),
    p_value   = numeric(),
    delta_aic = numeric(),
    stat      = numeric(),
    stat_type = character(),
    df        = integer(),
    stringsAsFactors = FALSE
  )

  null_result <- function(all_scores = empty_scores) {
    list(
      accepted  = FALSE,
      fit       = current,
      variable  = NA_character_,
      phase     = NA_character_,
      score     = NA_real_,
      p_value   = NA_real_,
      delta_aic = NA_real_,
      stat      = NA_real_,
      stat_type = NA_character_,
      df        = NA_integer_,
      all_scores     = all_scores,
      refit_failures = character(),
      refit_failure_reasons = character()
    )
  }

  if (length(cands) == 0L) {
    return(null_result())
  }

  rows <- vector("list", length(cands))
  for (i in seq_along(cands)) {
    cand <- cands[[i]]
    # Resolved by the variable's term in the stored design, not by name: a
    # factor `fla` with level `g` owns a column named `flag`, while a logical
    # `flag` owns `flagTRUE` (#315).  The design is rebuilt from the data the
    # fit was built on, so a `data` whose column types differ cannot skew it.
    coef_name <- .hzr_candidate_coef_name(current, cand$var, cand$phase)

    s <- .hzr_candidate_score(
      criterion = criterion, mode = "drop",
      current = current, names = coef_name
    )

    rows[[i]] <- data.frame(
      variable  = cand$var,
      phase     = cand$phase %||% NA_character_,
      force_in  = cand$var %in% force_in,
      score     = s$score,
      p_value   = s$p_value,
      delta_aic = s$delta_aic,
      stat      = s$stat,
      stat_type = s$stat_type,
      df        = s$df,
      stringsAsFactors = FALSE
    )
  }
  all_scores <- do.call(rbind, rows)

  eligible <- which(!all_scores$force_in & !is.na(all_scores$score))
  if (length(eligible) == 0L) {
    return(null_result(all_scores))
  }

  best_idx <- eligible[which.min(all_scores$score[eligible])]
  best     <- all_scores[best_idx, ]

  threshold_met <- if (criterion == "wald") {
    best$score < (1 - slstay)          # i.e. p > slstay
  } else {
    best$score < 0                     # dAIC_drop < 0
  }

  if (!threshold_met) {
    return(null_result(all_scores))
  }

  failure_token <- if (is.na(best$phase)) {
    best$variable
  } else {
    paste0(best$variable, "@", best$phase)
  }

  # A drop must remove a column.  Under treatment contrasts an interaction
  # whose main effect has gone is coded with a full set of dummies, so
  # dropping `z` from `~ z + z:f` turns `z, z:fb` into `z:fa, z:fb`: the same
  # column space and the same likelihood, a "reduced" model that is the model
  # it started from (#320).  Accepting it recorded a drop whose p-value
  # described a variable the fit still carries.  The forward step refuses the
  # mirror of this, a candidate that adds no column
  # (`.hzr_entered_coef_name()`, #306).
  refuse_no_column <- function(old_cols, new_cols) {
    reason <- paste0(
      "removes no column: the refit's design (",
      paste(sQuote(new_cols), collapse = ", "),
      ") has no fewer columns than the current one (",
      paste(sQuote(old_cols), collapse = ", "),
      "), so the model is not reduced"
    )
    warning("Stepwise backward: dropping ", failure_token, " ", reason, ".",
            call. = FALSE)
    out <- null_result(all_scores)
    out$refit_failures <- failure_token
    out$refit_failure_reasons <- stats::setNames(reason, failure_token)
    out
  }

  # A single-distribution refit warm-starts from `theta_old[-drop_idx]`, one
  # element shorter than the design a drop that removes no column leaves, so
  # that refit fails to conform and would report the arithmetic rather than
  # the cause (#323).  Its reduced design is therefore decided here, before
  # refitting, from the same formula parse `hazard()` uses.  Both designs are
  # built on `data`, so a factor level absent from it cannot pass for a
  # removed column.  A design that cannot be built is left to the refit,
  # which reports why.
  if (is.na(best$phase) && !is.null(current$call$formula)) {
    designs <- tryCatch({
      old_formula <- .hzr_stored_formula(current, "`current`")
      new_formula <- .hzr_formula_update(old_formula, "drop", best$variable)
      list(old = colnames(.hzr_parse_formula(old_formula, data)$x),
           new = colnames(.hzr_parse_formula(new_formula, data)$x))
    }, error = function(e) NULL)
    if (!is.null(designs) && length(designs$new) >= length(designs$old)) {
      return(refuse_no_column(designs$old, designs$new))
    }
  }

  refitted <- tryCatch(
    .hzr_refit_with_scope(
      current, action = "drop",
      var = best$variable,
      phase = if (is.na(best$phase)) NULL else best$phase,
      data = data, ...
    ),
    error = function(e) e
  )

  if (inherits(refitted, "error") || isFALSE(refitted$fit$converged)) {
    reason <- .hzr_refit_failure_reason(refitted)
    warning("Stepwise backward: post-drop refit failed for ",
            failure_token, ": ", reason, call. = FALSE)
    out <- null_result(all_scores)
    out$refit_failures <- failure_token
    out$refit_failure_reasons <- stats::setNames(reason, failure_token)
    return(out)
  }

  # The same check on the multiphase path, where the refit conforms and its
  # phase design is read off the result.
  if (!is.na(best$phase)) {
    old_cols <- colnames(current$fit$x_list[[best$phase]])
    new_cols <- colnames(refitted$fit$x_list[[best$phase]])
    if (is.null(old_cols)) {
      stop("Internal: phase ", sQuote(best$phase), " has no design columns, ",
           "so a drop from it cannot be checked.", call. = FALSE)
    }
    if (length(new_cols) >= length(old_cols)) {
      return(refuse_no_column(old_cols, new_cols))
    }
  }

  list(
    accepted  = TRUE,
    fit       = refitted,
    variable  = best$variable,
    phase     = best$phase,
    score     = best$score,
    p_value   = best$p_value,
    delta_aic = best$delta_aic,
    stat      = best$stat,
    stat_type = best$stat_type,
    df        = best$df,
    all_scores     = all_scores,
    refit_failures = character(),
    refit_failure_reasons = character()
  )
}


#' Name under which a newly-entered variable appears in coef(fit)
#'
#' Canonical naming differs between fit kinds:
#'   multiphase: phase-prefixed formula names (e.g. `"early.age"`).
#'   single-dist: positional `"betaN"` from `.hzr_parameter_names()`,
#'     where N is the index in `colnames(fit$data$x)` of the column `var`'s
#'     term builds (or of the column named `var`, for a fit with no stored
#'     design).
#'
#' This matches the naming `summary.hazard()` prints and the canonical
#' name `.hzr_wald_p()` uses for coefficient lookup.
#'
#' Stepwise v1 is scoped to main-effect terms only.  Multi-column
#' expansions (factors with > 2 levels, splines, interactions)
#' trigger an error so the caller can rebuild their formula with the
#' expansion baked in rather than silently scoring just one coefficient.
#'
#' @param current For a forward step, the fit `fit` was refit from.  The
#'   candidate is then resolved by the column the refit added, not by name.
#' @param data Data frame the fit's stored design is rebuilt against to
#'   resolve `var` by its term (`.hzr_term_coef_name()`).  Only a fit with no
#'   stored design falls back to looking `var` up by name.
#'
#' @keywords internal
#' @noRd
.hzr_candidate_coef_name <- function(fit, var, phase, current = NULL,
                                     data = fit$data$frame) {
  if (!is.null(current)) {
    return(.hzr_entered_coef_name(fit, current, var, phase))
  }
  by_term <- .hzr_term_coef_name(fit, var, phase, data)
  if (!is.null(by_term)) {
    return(by_term)
  }
  # No stored design (the `time =` / `x =` interface, or a fit saved before
  # the design was kept): the bare name is all there is to go on.
  if (fit$spec$dist == "multiphase") {
    target <- paste0(phase, ".", var)
    coef_names <- names(stats::coef(fit))
    if (!is.null(coef_names) && !target %in% coef_names) {
      expanded <- grep(
        paste0("^", .hzr_regex_escape(target)),
        coef_names, value = TRUE
      )
      if (length(expanded) > 1L) {
        stop(
          "Variable ", sQuote(var), " in phase ", sQuote(phase),
          " expands to multiple coefficients (",
          paste(sQuote(expanded), collapse = ", "),
          ").  Stepwise v1 supports main-effect terms only; ",
          "rebuild your candidate as pre-expanded main effects and retry.",
          call. = FALSE
        )
      }
      # A single expansion is the term under another name -- a logical becomes
      # `varTRUE`, a two-level factor `varb`. The name model.matrix() gave it is
      # the one the coefficient vector carries, so return that rather than the
      # name we guessed.
      if (length(expanded) == 1L) {
        return(expanded)
      }
    }
    return(target)
  }
  xcols <- colnames(fit$data$x)
  idx <- which(xcols == var)
  if (length(idx) == 0L) {
    # Check for a multi-column expansion masquerading as one variable.
    prefix_hits <- grep(paste0("^", .hzr_regex_escape(var)), xcols)
    if (length(prefix_hits) > 1L) {
      stop(
        "Variable ", sQuote(var),
        " expands to multiple coefficients (",
        paste(sQuote(xcols[prefix_hits]), collapse = ", "),
        ").  Stepwise v1 supports main-effect terms only; ",
        "rebuild your candidate as an already-expanded set of ",
        "main effects (e.g. a numeric contrast matrix) and retry.",
        call. = FALSE
      )
    }
    stop("Variable ", sQuote(var),
         " not found in the design matrix.",
         call. = FALSE)
  }
  if (length(idx) > 1L) {
    stop("Variable ", sQuote(var),
         " matches multiple design-matrix columns -- ",
         "stepwise v1 supports main-effect terms only.",
         call. = FALSE)
  }
  paste0("beta", idx)
}


#' Name of an in-model variable's coefficient, found by its term
#'
#' A drop has no refit to compare columns against, so the variable is found
#' through the design the fit stored: rebuilt against `data`, the model
#' matrix's `assign` attribute maps each column to its term, and `var` is a
#' term label.  Looking `var` up among the column names instead finds another
#' term's column when that column carries the name: a logical `flag` becomes
#' `flagTRUE`, while a factor `fla` with level `g` owns `flag` (#315).
#'
#' @return The coefficient name (`beta<k>` for a single distribution,
#'   `<phase>.<column>` for multiphase), or `NULL` when the fit stored no
#'   design, which leaves only the name lookup.
#'
#' @keywords internal
#' @noRd
.hzr_term_coef_name <- function(fit, var, phase, data) {
  multiphase <- fit$spec$dist == "multiphase"
  if (multiphase) {
    design <- if (.hzr_phase_inherits_global(fit, phase)) {
      fit$data$x_design
    } else {
      fit$fit$x_design[[phase]]
    }
    cols <- colnames(fit$fit$x_list[[phase]])
  } else {
    design <- fit$data$x_design
    cols <- colnames(fit$data$x)
  }
  if (is.null(design) || is.null(data)) {
    return(NULL)
  }
  where <- if (multiphase) paste0(" in phase ", sQuote(phase)) else ""

  mf <- stats::model.frame(design$terms, data = data, xlev = design$xlevels,
                           na.action = stats::na.pass)
  mm <- stats::model.matrix(design$terms, data = mf,
                            contrasts.arg = design$contrasts)
  # Column names are unique within a design (.hzr_refuse_duplicate_columns()),
  # so the stored columns can be found in the rebuild by name.  The stored
  # set lacks the intercept, and under `time_windows` it is window-expanded.
  pos <- match(cols, colnames(mm))
  if (length(cols) == 0L || anyNA(pos)) {
    stop("Cannot map variable ", sQuote(var), where, " to its coefficient: ",
         "the fit's design columns (", paste(sQuote(cols), collapse = ", "),
         ") are not those its stored formula builds (",
         paste(sQuote(colnames(mm)), collapse = ", "), ").",
         call. = FALSE)
  }
  term <- match(var, attr(design$terms, "term.labels"))
  if (is.na(term)) {
    stop("Variable ", sQuote(var), where, " is not a term of the fitted ",
         "model's formula.", call. = FALSE)
  }
  idx <- which(attr(mm, "assign")[pos] == term)
  if (length(idx) > 1L) {
    stop(
      "Variable ", sQuote(var), where,
      " expands to multiple coefficients (",
      paste(sQuote(cols[idx]), collapse = ", "),
      ").  Stepwise v1 supports main-effect terms only; ",
      "rebuild your candidate as pre-expanded main effects and retry.",
      call. = FALSE
    )
  }
  if (length(idx) == 0L) {
    stop("Variable ", sQuote(var), where, " has no column in the fitted ",
         "model's design, so its coefficient cannot be identified.",
         call. = FALSE)
  }
  if (multiphase) paste0(phase, ".", cols[idx]) else paste0("beta", idx)
}


#' Name of the coefficient a forward step's refit added
#'
#' The refit's design differs from the current fit's by the candidate's
#' column, so the difference names it, whatever `model.matrix()` called it.
#' Looking the candidate up by its bare name finds the wrong column when
#' another term's column carries that name: a logical `flag` becomes
#' `flagTRUE`, while a factor `fla` with level `g` owns `flag`.  The score
#' path places its pinned zero by the same difference
#' (`.hzr_score_expand_multiphase()`).
#'
#' @keywords internal
#' @noRd
.hzr_entered_coef_name <- function(fit, current, var, phase) {
  multiphase <- fit$spec$dist == "multiphase"
  design_cols <- function(f) {
    colnames(if (multiphase) f$fit$x_list[[phase]] else f$data$x)
  }
  new_cols <- design_cols(fit)
  old_cols <- design_cols(current)
  added <- which(!new_cols %in% old_cols)
  where <- if (multiphase) paste0(" in phase ", sQuote(phase)) else ""
  if (length(added) > 1L) {
    stop(
      "Variable ", sQuote(var), where,
      " expands to multiple coefficients (",
      paste(sQuote(new_cols[added]), collapse = ", "),
      ").  Stepwise v1 supports main-effect terms only; ",
      "rebuild your candidate as pre-expanded main effects and retry.",
      call. = FALSE
    )
  }
  # One new NAME is not one new column: `z` added to `~ z:f` turns
  # `z:fa, z:fb` into `z, z:fb`, the same column space and likelihood.  The
  # score path requires the count to rise by one for the same reason.
  if (length(new_cols) != length(old_cols) + 1L) {
    stop(
      "Variable ", sQuote(var), where, " does not add a column to the model: ",
      "the refit's design (", paste(sQuote(new_cols), collapse = ", "),
      ") reparameterises the current one (",
      paste(sQuote(old_cols), collapse = ", "),
      "), so there is no coefficient of its own to test.",
      call. = FALSE
    )
  }
  if (length(added) == 0L) {
    stop("Variable ", sQuote(var), where,
         " added no design-matrix column the current fit lacks, so its ",
         "coefficient cannot be identified.",
         call. = FALSE)
  }
  if (multiphase) paste0(phase, ".", new_cols[added]) else paste0("beta", added)
}


#' Escape a character vector for literal use in a regular expression
#'
#' @keywords internal
#' @noRd
.hzr_regex_escape <- function(x) {
  gsub("([.\\\\|*+?^$(){}\\[\\]])", "\\\\\\1", x, perl = TRUE)
}


#' %||%: NULL-coalesce for lazy defaults
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a
