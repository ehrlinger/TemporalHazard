# stepwise-formula.R -- Formula and phase-scope plumbing for stepwise
# covariate selection.
#
# Step 8.3 of STEPWISE-DESIGN.md.  Pure helpers that mutate formulas /
# phase lists / coefficient scopes without touching the optimizer.  The
# actual refit driver that combines these into a new hazard() call
# lands in sec.8.4 alongside the forward/backward step implementations.
#
# Two entry points matter to downstream code:
#
#   .hzr_formula_update(formula, action, var)
#     Add or drop a variable from a one-sided or two-sided formula.
#
#   .hzr_phase_update_formula(phase, action, var)
#     Same mutation, but on the `formula` slot of an `hzr_phase` object;
#     returns an updated `hzr_phase` so the multiphase driver can
#     rebuild the phases list by name.
#
# Both are symmetric under add/drop and idempotent w.r.t. already-present
# or already-absent variables.  Variables currently in scope come from:
#
#   .hzr_scope_current_vars(fit, phase = NULL)
#     Lists the RHS term labels for the global formula (single-dist) or
#     the named phase's formula (multiphase).

# ---------------------------------------------------------------------------
# Formula helpers
# ---------------------------------------------------------------------------

#' Extract RHS term labels from a formula
#'
#' Handles both one-sided (`~ x + y`) and two-sided (`Surv(time, status) ~
#' x + y`) formulas.  Returns an empty character vector for `~ 1`.
#'
#' @keywords internal
#' @noRd
.hzr_formula_rhs_terms <- function(formula) {
  if (is.null(formula)) return(character())
  if (!inherits(formula, "formula")) {
    stop("formula must be a `formula` object.", call. = FALSE)
  }
  # terms() needs a formula with no empty LHS/RHS; ~ 1 has one intercept term.
  tt <- tryCatch(stats::terms(formula),
                 error = function(e) NULL)
  if (is.null(tt)) return(character())
  attr(tt, "term.labels")
}


#' Add or drop a variable from a formula's RHS
#'
#' @param formula Existing formula.  One-sided (`~ x`) or two-sided
#'   (`Surv(time, status) ~ x`) -- the LHS is preserved verbatim.
#' @param action Either `"add"` or `"drop"`.
#' @param var Character scalar naming the variable to add / drop.
#'
#' @return A new formula with the requested change.  Idempotent: adding
#'   a variable that is already present is a no-op; dropping one that
#'   is not present is a no-op.  Dropping the last term leaves `~ 1`.
#'
#' @keywords internal
#' @noRd
.hzr_formula_update <- function(formula, action = c("add", "drop"), var) {
  action <- match.arg(action)
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a `formula` object.", call. = FALSE)
  }
  if (!is.character(var) || length(var) != 1L || !nzchar(var)) {
    stop("`var` must be a non-empty character scalar.", call. = FALSE)
  }

  current_terms <- .hzr_formula_rhs_terms(formula)

  if (action == "add") {
    if (var %in% current_terms) return(formula)
    new_terms <- c(current_terms, var)
  } else {
    if (!var %in% current_terms) return(formula)
    new_terms <- setdiff(current_terms, var)
  }

  rhs <- if (length(new_terms) == 0L) "1" else paste(new_terms, collapse = " + ")

  is_two_sided <- length(formula) == 3L
  lhs <- if (is_two_sided) deparse(formula[[2L]]) else NULL
  env <- environment(formula)

  new_text <- if (is_two_sided) {
    paste(lhs, "~", rhs)
  } else {
    paste("~", rhs)
  }

  stats::as.formula(new_text, env = env)
}


# ---------------------------------------------------------------------------
# Phase-level helpers (multiphase models)
# ---------------------------------------------------------------------------

#' Add or drop a variable from a phase's formula
#'
#' @param phase An `hzr_phase` object.  Its `formula` slot may be NULL
#'   (no phase-specific covariates) -- in the add case a fresh
#'   `~ var` formula is created.
#' @param action Either `"add"` or `"drop"`.
#' @param var Character scalar.
#'
#' @return An updated `hzr_phase` object with the new formula.
#'
#' @keywords internal
#' @noRd
.hzr_phase_update_formula <- function(phase, action = c("add", "drop"), var) {
  action <- match.arg(action)
  if (!inherits(phase, "hzr_phase")) {
    stop("`phase` must be an `hzr_phase` object.", call. = FALSE)
  }
  if (!is.character(var) || length(var) != 1L || !nzchar(var)) {
    stop("`var` must be a non-empty character scalar.", call. = FALSE)
  }

  f <- phase$formula
  if (is.null(f)) {
    if (action == "drop") {
      return(phase)   # nothing to drop
    }
    # Create a fresh one-sided formula in the current calling env so the
    # variable can be resolved later by data-frame column lookup.
    phase$formula <- stats::as.formula(
      paste("~", var),
      env = parent.frame()
    )
    return(phase)
  }

  phase$formula <- .hzr_formula_update(f, action, var)

  # If dropping emptied the RHS (`~ 1`), null the slot out for
  # consistency with hzr_phase(formula = NULL) construction.
  if (action == "drop" &&
        identical(.hzr_formula_rhs_terms(phase$formula), character())) {
    phase$formula <- NULL
  }

  phase
}


# ---------------------------------------------------------------------------
# Scope inspection
# ---------------------------------------------------------------------------

#' List variables currently in a hazard model's scope
#'
#' For single-distribution models returns the RHS term labels of the
#' global formula recorded in `fit$call$formula`.  For multiphase models
#' returns the term labels of the named phase's formula; pass
#' `phase = NULL` to get a named list keyed by phase.
#'
#' @param fit A fitted `hazard` object.
#' @param phase Character scalar phase name (multiphase only), or NULL
#'   to return all phases.
#'
#' @return Character vector of variable names, or a named list of such
#'   vectors.
#'
#' @keywords internal
#' @noRd
.hzr_scope_current_vars <- function(fit, phase = NULL) {
  if (!inherits(fit, "hazard")) {
    stop("`fit` must be a `hazard` object.", call. = FALSE)
  }

  if (fit$spec$dist != "multiphase") {
    if (!is.null(phase)) {
      stop("`phase` is only meaningful for multiphase models.", call. = FALSE)
    }
    f <- .hzr_stored_formula(fit)
    if (is.null(f)) {
      # Non-formula fit (time/x interface); infer from design matrix names.
      xcols <- colnames(fit$data$x)
      return(if (is.null(xcols)) character() else xcols)
    }
    return(.hzr_formula_rhs_terms(f))
  }

  # Multiphase
  per_phase <- lapply(fit$spec$phases, function(ph) {
    if (is.null(ph$formula)) character() else .hzr_formula_rhs_terms(ph$formula)
  })

  if (is.null(phase)) return(per_phase)

  if (!is.character(phase) || length(phase) != 1L) {
    stop("`phase` must be a single character scalar.", call. = FALSE)
  }
  if (!phase %in% names(per_phase)) {
    stop("Unknown phase: ", sQuote(phase), ".  Available: ",
         paste(sQuote(names(per_phase)), collapse = ", "),
         call. = FALSE)
  }
  per_phase[[phase]]
}

# ---------------------------------------------------------------------------
# Phase-scope detection helper
# ---------------------------------------------------------------------------
#' Detect phase-scoped calls in a formula RHS parse tree
#'
#' Walks the call tree of \code{rhs} (the RHS of a formula, i.e.
#' \code{formula[[3L]]}) and returns \code{TRUE} if any call node has a
#' function symbol that exactly matches one of \code{phase_names}.
#'
#' This is stricter than a string-regex approach: a phase named \code{"log"}
#' will NOT produce a false positive when the formula contains \code{log(age)},
#' because \code{log} would also appear in \code{phase_names} only if the user
#' deliberately named a phase \code{"log"}.  Conversely, the function only
#' fires when the call head is an exact match to a known phase name -- standard
#' R functions that happen to share names with phases do not trigger the check
#' unless those names are actually phase names.
#'
#' @param rhs  A language object (the RHS of a formula, typically
#'   \code{formula[[3L]]}).
#' @param phase_names  Character vector of phase names to look for.
#' @return \code{TRUE} if any call in the tree has its function head in
#'   \code{phase_names}; \code{FALSE} otherwise.
#' @noRd
.hzr_formula_has_phase_scope <- function(rhs, phase_names) {
  if (is.null(rhs) || length(phase_names) == 0L) return(FALSE)
  # Walk recursively: check current node, then recurse into sub-expressions.
  .walk <- function(node) {
    if (is.call(node)) {
      # The head of the call (node[[1L]]) is the function symbol or expression.
      head <- node[[1L]]
      if (is.symbol(head) && as.character(head) %in% phase_names) {
        return(TRUE)
      }
      # Recurse into all arguments (node[[2L]], node[[3L]], ...) and the head.
      for (i in seq_along(node)) {
        if (.walk(node[[i]])) return(TRUE)
      }
    }
    FALSE
  }
  .walk(rhs)
}

#' Resolve the model formula a fit stored in its call
#'
#' `match.call()` stores the `formula` argument unevaluated, so it arrives as
#' a `call` when the caller wrote it out at the `hazard()` call and as a
#' *symbol* when they passed a variable holding it. Evaluating in the
#' recorded calling environment covers both; deparsing does not, because
#' `deparse(quote(f))` is `"f"`, which `as.formula()` rejects as "not a call".
#'
#' Every site that needs the stored formula goes through here. Three did it
#' independently once, and two of them were still deparsing after the third
#' was fixed.
#'
#' @param fit A fitted `hazard` object.
#' @param what Character label for the object in error messages.
#' @return The formula, or `NULL` when the fit carries none (the vector
#'   interface). Errors when a stored formula cannot be resolved.
#' @keywords internal
#' @noRd
.hzr_stored_formula <- function(fit, what = "`fit`") {
  raw <- fit$call$formula
  if (is.null(raw)) {
    return(NULL)
  }
  if (inherits(raw, "formula")) {
    return(raw)
  }
  # Evaluating can fail outright when the formula was passed by variable and
  # that binding is gone -- after saveRDS()/readRDS() in a new session, say.
  # Bare, the error reads "object 'f' not found", which names neither the fit
  # at fault nor the remedy.
  out <- tryCatch(
    eval(raw, envir = fit$call_env %||% parent.frame(2L)),
    error = function(e) {
      stop(what, "'s stored model formula (", deparse(raw),
           ") could not be resolved: ", conditionMessage(e),
           ". The formula was passed to `hazard()` by variable and that ",
           "binding is no longer reachable -- for example after saving and ",
           "reloading the fit in a new session. Refit the base model with ",
           "the formula written at the `hazard()` call.", call. = FALSE)
    }
  )
  if (!inherits(out, "formula")) {
    stop(what, "'s stored model formula (", deparse(raw),
         ") did not resolve to a formula. Refit the base model with the ",
         "formula written at the `hazard()` call.", call. = FALSE)
  }
  out
}
