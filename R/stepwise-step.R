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
#' candidate -- that reuse is what removes the optimizer from the loop, and
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
#'   \item{stat, df}{Wald statistic / df of the winner.}
#'   \item{all_scores}{Tibble-like data frame of every candidate
#'     considered and its score.}
#'   \item{refit_failures}{Character vector of `"var@phase"` tokens for
#'     candidates whose refit diverged.}
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
      df        = NA_integer_,
      all_scores = data.frame(
        variable  = character(),
        phase     = character(),
        score     = numeric(),
        p_value   = numeric(),
        delta_aic = numeric(),
        stat      = numeric(),
        df        = integer(),
        stringsAsFactors = FALSE
      ),
      refit_failures = character()
    )
  }

  if (length(cands) == 0L) {
    return(null_result())
  }

  rows     <- vector("list", length(cands))
  failures <- character()

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
      warning("Stepwise forward: candidate refit failed for ",
              failure_token, ".", call. = FALSE)
      failures <- c(failures, failure_token)
      rows[[i]] <- data.frame(
        variable  = cand$var,
        phase     = cand$phase %||% NA_character_,
        score     = NA_real_,
        p_value   = NA_real_,
        delta_aic = NA_real_,
        stat      = NA_real_,
        df        = NA_integer_,
        stringsAsFactors = FALSE
      )
      next
    }

    # Coefficient name(s) of the newly-entered variable in the candidate
    # fit.  For single-dist this is `var` directly; for multiphase the
    # coef is phase-prefixed.
    coef_name <- .hzr_candidate_coef_name(candidate_fit, cand$var,
                                           cand$phase)

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
    df        = best$df,
    all_scores = all_scores,
    refit_failures = failures
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
#' The winner is refit ONCE, after the decision -- mirroring the backward
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
      df        = s$df,
      reason    = s_q$reason %||% NA_character_,
      stringsAsFactors = FALSE
    )
  }
  all_scores <- do.call(rbind, rows)

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
    out$stop_reason    <- "no_candidate_met_slentry"
    return(out)
  }

  best_phase <- if (is.na(best$phase)) NULL else best$phase
  refitted <- tryCatch(
    .hzr_refit_with_scope(
      current, action = "add", var = best$variable, phase = best_phase,
      data = data, ...
    ),
    error = function(e) e
  )

  failure_token <- if (is.na(best$phase)) {
    best$variable
  } else {
    paste0(best$variable, "@", best$phase)
  }

  if (inherits(refitted, "error") || isFALSE(refitted$fit$converged)) {
    # The candidate won on Q but the model that would realise it will not fit.
    # Entering it anyway would put a non-converged fit into the chain.
    warning("Stepwise forward: post-entry refit failed for ",
            failure_token, ".", call. = FALSE)
    out <- null_result()
    out$all_scores     <- all_scores
    out$refit_failures <- failure_token
    out$n_uncomputable <- n_uncomputable
    out$uncomputable_reasons <- uncomputable_reasons
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
    df        = best$df,
    all_scores = all_scores,
    refit_failures = character(),
    n_uncomputable = n_uncomputable,
    uncomputable_reasons = uncomputable_reasons,
    stop_reason    = "accepted"
  )
}


#' Reject a candidate the score path cannot expand into a single column
#'
#' `.hzr_score_q()` returns a silent `NA` for a non-numeric or absent
#' candidate, which under `criterion = "score"` would mean the variable
#' simply never enters -- no error, no warning, just absent from the
#' selected model.  The Wald path surfaces both cases loudly: a factor
#' candidate fails in `.hzr_candidate_coef_name()` (cannot find its expanded
#' coefficient by name), and a candidate missing from `data` fails inside
#' `.hzr_refit_with_scope()`'s refit, which emits the
#' "candidate refit failed for ..." warning. Matching that here keeps the
#' two criteria in agreement about what is selectable.
#'
#' A column absent from `data` warns (matching the Wald refit-failure
#' warning's severity) and is left to the caller's normal `NA` handling --
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
      df        = NA_integer_,
      all_scores     = all_scores,
      refit_failures = character()
    )
  }

  if (length(cands) == 0L) {
    return(null_result())
  }

  rows <- vector("list", length(cands))
  for (i in seq_along(cands)) {
    cand <- cands[[i]]
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

  refitted <- tryCatch(
    .hzr_refit_with_scope(
      current, action = "drop",
      var = best$variable,
      phase = if (is.na(best$phase)) NULL else best$phase,
      data = data, ...
    ),
    error = function(e) e
  )

  failure_token <- if (is.na(best$phase)) {
    best$variable
  } else {
    paste0(best$variable, "@", best$phase)
  }

  if (inherits(refitted, "error") || isFALSE(refitted$fit$converged)) {
    warning("Stepwise backward: post-drop refit failed for ",
            failure_token, ".", call. = FALSE)
    out <- null_result(all_scores)
    out$refit_failures <- failure_token
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
    df        = best$df,
    all_scores     = all_scores,
    refit_failures = character()
  )
}


#' Name under which a newly-entered variable appears in coef(fit)
#'
#' Canonical naming differs between fit kinds:
#'   multiphase  -- phase-prefixed formula names (e.g. `"early.age"`).
#'   single-dist -- positional `"betaN"` from `.hzr_parameter_names()`,
#'     where N is the column index of `var` in `colnames(fit$data$x)`.
#'
#' This matches the naming `summary.hazard()` prints and the canonical
#' name `.hzr_wald_p()` uses for coefficient lookup.
#'
#' Stepwise v1 is scoped to main-effect terms only.  Multi-column
#' expansions (factors with > 2 levels, splines, interactions)
#' trigger an error so the caller can rebuild their formula with the
#' expansion baked in rather than silently scoring just one coefficient.
#'
#' @keywords internal
#' @noRd
.hzr_candidate_coef_name <- function(fit, var, phase) {
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


#' Escape a character vector for literal use in a regular expression
#'
#' @keywords internal
#' @noRd
.hzr_regex_escape <- function(x) {
  gsub("([.\\\\|*+?^$(){}\\[\\]])", "\\\\\\1", x, perl = TRUE)
}


#' %||% -- NULL-coalesce for lazy defaults
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a
