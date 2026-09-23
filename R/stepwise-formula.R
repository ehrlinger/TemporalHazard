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
  # `.` can only be expanded against data, and this reader has none. terms()
  # errors on it, and the tryCatch() below turned that into "no terms", so a
  # `~ .` base model reported a finished screen of zero steps (#279).
  if ("." %in% all.vars(formula)) {
    stop("hzr_stepwise() cannot expand `.` in `",
         paste(deparse(formula), collapse = " "), "`: write the terms out, ",
         "as in `~ age + mal`. A base model written with `.` must be refit ",
         "that way; a `scope` can list its variables.", call. = FALSE)
  }
  # terms() needs a formula with no empty LHS/RHS; ~ 1 has one intercept term.
  tt <- tryCatch(stats::terms(formula),
                 error = function(e) NULL)
  if (is.null(tt)) return(character())
  # term.labels leaves an offset out, so a scope naming one lost it silently.
  .hzr_refuse_offset(formula, paste0(
    "the formula `",
    paste(deparse(formula, width.cutoff = 500L), collapse = " "), "`"
  ))
  attr(tt, "term.labels")
}


#' The term label `terms()` gives a data column
#'
#' The identity a candidate is compared by. `terms()` backquotes a name that
#' is not syntactic, so the column `_X1` is labelled `` `_X1` `` and the
#' column `TRUE` is labelled `` `TRUE` ``. The label is produced by
#' `terms()` itself, from the column's symbol, so no string is parsed and it
#' is exactly the label a model containing the column carries. Distinct
#' columns get distinct labels, and a column's label never equals an
#' expression's: the column `age:mal` is `` `age:mal` ``, the interaction is
#' `age:mal`.
#'
#' @param x Character vector of column names.
#' @return Character vector of labels, the same length. A name no symbol can
#'   carry (`""`) or that `terms()` refuses (`"."`) gets a placeholder no term
#'   label can equal, so it matches only itself.
#' @keywords internal
#' @noRd
.hzr_column_label <- function(x) {
  vapply(x, function(nm) {
    tryCatch(
      attr(stats::terms(stats::as.formula(call("~", as.name(nm)))),
           "term.labels"),
      error = function(e) {
        paste0("<column ", encodeString(nm, quote = "\""), ">")
      }
    )
  }, character(1), USE.NAMES = FALSE)
}

#' Is a column label the placeholder for a name no formula can hold?
#'
#' `.hzr_column_label()` gives `""` and `"."` a placeholder, `<column ".">`,
#' because `terms()` cannot label them (`.` means every other column). Such a
#' column is never a stepwise candidate: pasted into formula text, the
#' placeholder does not parse (#449).
#'
#' @keywords internal
#' @noRd
.hzr_is_label_placeholder <- function(label) {
  startsWith(label, "<column ")
}

#' The term a candidate's refit writes into the formula
#'
#' `.hzr_formula_update()` writes its `var` into the formula TEXT, so the
#' model gains whatever that text parses to. The candidate's spelling is the
#' user's column name, which can parse to a different term: `age ` reads as
#' `age`, the column `age:mal` as the interaction, and `_X1` does not parse
#' (#449, #441). Its identity is the label `terms()` wrote for the resolved
#' column or term, which parses back to exactly that, so that is what the
#' refit is given.
#'
#' @param cand A candidate from `.hzr_stepwise_candidates()`.
#' @return A term label.
#' @keywords internal
#' @noRd
.hzr_candidate_term <- function(cand) {
  cand$id %||% cand$var
}

#' The data column a candidate is, if it is one
#'
#' The score reads the candidate's values from `data`, so it needs the COLUMN
#' the candidate resolved to, which is neither its spelling nor its label in
#' general: a formula `scope` spells `_X1` as its label `` `_X1` ``, which is
#' no column name (#438), and the interaction `age:mal` is spelled like a
#' literal column `age:mal` that it is not (#449). The column is the one
#' whose own `terms()` label is the candidate's identity.
#'
#' @param cand A candidate from `.hzr_stepwise_candidates()`.
#' @param data The screen's data frame.
#' @return The column name, or `NA` when the candidate is not a column.
#' @keywords internal
#' @noRd
.hzr_candidate_column <- function(cand, data) {
  cols <- names(data)
  hit <- match(.hzr_candidate_term(cand), .hzr_column_label(cols))
  if (is.na(hit)) NA_character_ else cols[[hit]]
}

#' Is a string a term label of `data`, exactly as `terms()` writes it?
#'
#' True when the string is the single label `terms()` gives a formula whose
#' right-hand side is that string, and every variable that term reads is a
#' column of `data`. `"age "` and `"age # x"` are not labels (`terms()`
#' writes `age`), nor is `"TRUE"` (a constant, no term), nor `"age*mal"`
#' (three terms), nor a bare name that is not a column.
#'
#' @keywords internal
#' @noRd
.hzr_is_term_label <- function(s, data) {
  isTRUE(tryCatch({
    f <- stats::reformulate(s)
    vars <- all.vars(f)
    identical(attr(stats::terms(f), "term.labels"), s) &&
      length(vars) > 0L && all(vars %in% names(data))
  }, error = function(e) FALSE))
}

#' Resolve user-supplied names to the terms they name
#'
#' `force_in`, `force_out` and a character `scope` are documented as
#' variables, but a candidate is compared by its `terms()` label, and the
#' two spellings differ for a name that is not syntactic (#437). Resolving a
#' string by PARSING it cannot be right: the same text is a raw column name
#' at some sites and a term label at others, and each rule tried for #442
#' merged two things that were different. So a string is resolved by LOOKUP,
#' once, and every later comparison is on the result:
#'
#' 1. exactly a column of `data`: that column, identified by
#'    `.hzr_column_label()`;
#' 2. otherwise exactly a label in `labels`, or a column's label, or, with
#'    `self_label = TRUE`, a string that is itself a term label over columns
#'    of `data` (`.hzr_is_term_label()`): that term;
#' 3. otherwise it names nothing. It is WARNED about, naming it, and dropped.
#'
#' A string that is both a column and a term label names the COLUMN, so with
#' a column literally called `age:mal` the interaction is reachable only
#' through a string that is not a column name.
#'
#' @param x Character vector supplied by the user.
#' @param data The screen's data frame.
#' @param labels Term labels step 2 accepts, besides the columns' own.
#' @param arg The argument's name, for the warning, e.g. `` "`force_in`" ``.
#' @param self_label Accept a string that is itself a term label. Used for a
#'   character `scope`, which introduces its own terms.
#' @return A list: `spelling`, the resolved elements as the user wrote them;
#'   `id`, the label each resolves to, in the same order; and `unresolved`,
#'   the elements that resolved to nothing, which hzr_stepwise() records on
#'   its result so that the warning is not the only trace of them.
#' @keywords internal
#' @noRd
.hzr_resolve_names <- function(x, data, labels = character(), arg,
                               self_label = FALSE) {
  x <- as.character(x)
  id <- rep(NA_character_, length(x))
  cols <- names(data)
  is_col <- x %in% cols
  id[is_col] <- .hzr_column_label(x[is_col])
  known <- c(labels, .hzr_column_label(cols))
  is_lab <- !is_col & x %in% known
  id[is_lab] <- x[is_lab]
  if (self_label) {
    for (i in which(is.na(id))) {
      if (.hzr_is_term_label(x[i], data)) id[i] <- x[i]
    }
  }
  bad <- is.na(id)
  if (any(bad)) {
    warning(arg, " names ",
            paste(encodeString(x[bad], quote = "\""), collapse = ", "),
            ", which is neither a column of `data` nor a term label of the ",
            "model or `scope`; ",
            if (sum(bad) == 1L) "it is" else "they are", " ignored.",
            call. = FALSE)
  }
  list(spelling = x[!bad], id = id[!bad], unresolved = x[bad])
}

#' Add or drop a variable from a formula's RHS
#'
#' @param formula Existing formula.  One-sided (`~ x`) or two-sided
#'   (`Surv(time, status) ~ x`); the LHS is preserved verbatim.
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
#'   (no phase-specific covariates), in which case the phase inherits the
#'   global design: the step starts from `inherited` when that is given, and
#'   otherwise the add case creates a fresh `~ var` formula.
#' @param action Either `"add"` or `"drop"`.
#' @param var Character scalar.
#' @param inherited One-sided formula of the global terms a phase with no
#'   formula inherits (see `.hzr_inherited_rhs()`), or `NULL` when there are
#'   none. Ignored for a phase with a formula of its own.
#'
#' @return An updated `hzr_phase` object with the new formula.
#'
#' @keywords internal
#' @noRd
.hzr_phase_update_formula <- function(phase, action = c("add", "drop"), var,
                                      inherited = NULL) {
  action <- match.arg(action)
  if (!inherits(phase, "hzr_phase")) {
    stop("`phase` must be an `hzr_phase` object.", call. = FALSE)
  }
  if (!is.character(var) || length(var) != 1L || !nzchar(var)) {
    stop("`var` must be a non-empty character scalar.", call. = FALSE)
  }

  f <- phase$formula
  if (is.null(f) && !is.null(inherited)) {
    # A phase with no formula of its own carries the global terms, so a step
    # starts from them: adding `mal` to a phase inheriting `age` gives
    # `~ age + mal`, not `~ mal` (#284). An emptied result stays `~ 1`
    # rather than NULL, which would inherit the global terms again.
    phase$formula <- .hzr_formula_update(inherited, action, var)
    return(phase)
  }
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

  # A drop that empties the RHS leaves `~ 1`, never NULL: a NULL phase
  # formula inherits the global design, so the phase would take back the
  # global covariates the step table says it no longer has (#284).
  phase$formula <- .hzr_formula_update(f, action, var)
  phase
}


#' The global terms a phase with no formula of its own inherits
#'
#' A multiphase phase whose `formula` is `NULL` uses the global design: the
#' right-hand side of the fit's formula. A stepwise step on such a phase has
#' to start from those terms (#284). A vector-interface fit has no formula,
#' so nothing here can say what a directly passed `x` holds;
#' `.hzr_refit_blocker()` refuses the case where a phase inherits one.
#'
#' @param fit A fitted `hazard` object.
#' @return A one-sided formula in the stored formula's environment, or `NULL`
#'   when the global design has no terms.
#' @keywords internal
#' @noRd
.hzr_inherited_rhs <- function(fit) {
  f <- .hzr_stored_formula(fit)
  if (is.null(f) || length(f) < 3L) {
    return(NULL)
  }
  rhs_terms <- .hzr_formula_rhs_terms(f)
  if (length(rhs_terms) == 0L) {
    return(NULL)
  }
  stats::reformulate(rhs_terms, env = environment(f))
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

  # Multiphase. A phase with no formula of its own carries the global terms
  # it inherits (#284); they are read only when some phase does inherit, so
  # a global `.` beside phases that all have formulas is never parsed here.
  inherits_global <- vapply(fit$spec$phases, function(ph) is.null(ph$formula),
                            logical(1))
  inherited_terms <- character()
  if (any(inherits_global)) {
    rhs <- .hzr_inherited_rhs(fit)
    if (!is.null(rhs)) inherited_terms <- .hzr_formula_rhs_terms(rhs)
  }
  per_phase <- lapply(fit$spec$phases, function(ph) {
    if (is.null(ph$formula)) inherited_terms else .hzr_formula_rhs_terms(ph$formula)
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
#' This is stricter than a string-regex approach.  A call such as
#' \code{log(age)} triggers the check only when a phase is actually named
#' \code{"log"}; with phases named \code{"early"} and \code{"constant"} it
#' does not.  When a phase is named \code{"log"}, \code{log(age)} does
#' trigger it and \code{hazard()} refuses the formula (#275); write
#' \code{base::log(age)} to use the function.  The function fires only when a
#' call head exactly matches a known phase name, so a variable such as
#' \code{early_age}, or a bare symbol \code{early} that is not called, does
#' not trigger it.
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

#' Resolve the model formula stored in a fit's call
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
#' @param envir Environment to resolve in when the fit carries no recorded
#'   `call_env` (objects serialised before it was stored). Defaults to the
#'   caller's frame; passing it explicitly is preferable to relying on the
#'   default, which is only correct because it is evaluated lazily.
#' @return The formula, or `NULL` when the fit carries none (the vector
#'   interface). Errors when a stored formula cannot be resolved.
#' @keywords internal
#' @noRd
.hzr_stored_formula <- function(fit, what = "`fit`",
                                envir = parent.frame()) {
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
    eval(raw, envir = fit$call_env %||% envir),
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

#' Columns the package may offer itself as candidates
#'
#' Under `scope = NULL` the *package* enumerates candidates from the data, so a
#' column it cannot model is its own bad choice rather than the caller's, and
#' dropping it beats erroring on it. An explicit `scope` still errors, because
#' there the caller named the column.
#'
#' Logical columns count: they are ordinary 0/1 predictors, and whether a
#' 0/1 field arrives logical or numeric depends on the reader that produced
#' the frame rather than on anything about the variable.
#'
#' @keywords internal
#' @noRd
.hzr_modellable_vars <- function(data, vars) {
  keep <- vapply(
    data[vars],
    function(col) is.numeric(col) || is.logical(col),
    logical(1L)
  )
  vars[keep]
}

#' Coerce a candidate column to the numeric vector the screen models
#'
#' Logical columns are ordinary 0/1 predictors, and `.hzr_modellable_vars()`
#' offers them as candidates under `scope = NULL`. Everything downstream (the
#' score statistic, the design-matrix column) wants a numeric vector, so the
#' translation happens once here rather than teaching each site about logicals.
#'
#' Returns `NULL` for anything that is not modellable as a single numeric
#' column, which is the caller's signal to skip or refuse it.
#'
#' @keywords internal
#' @noRd
.hzr_candidate_numeric <- function(x) {
  if (is.logical(x)) return(as.numeric(x))
  if (is.numeric(x)) return(x)
  NULL
}

#' Refuse a `scope` the screen would not honour
#'
#' A backward screen only drops terms the base model already has, so it
#' never reads `scope`: the variables it leaves out are still dropped, and the
#' ones the base lacks are never tested. A two-sided scope formula is read by
#' its right-hand side only, so its left-hand side is never a candidate, and a
#' phase named twice in a scope list is read at its first entry only. Each ran
#' to a result with no message (#343). `hzr_stepwise()` and `hzr_bootstrap()`
#' call this before any fitting or seeding.
#'
#' @param scope The `scope` argument as given.
#' @param direction The matched `direction`.
#' @param caller `"hzr_stepwise"` or `"hzr_bootstrap"`, which need different
#'   remedies for a backward screen.
#' @return `NULL`, invisibly; otherwise stops.
#' @keywords internal
#' @noRd
.hzr_refuse_unhonoured_scope <- function(scope, direction,
                                         caller = "hzr_stepwise") {
  if (is.null(scope)) {
    return(invisible(NULL))
  }
  # An empty scope offers nothing to enter, which a backward screen honours.
  # An offset is no candidate, but it is refused under "both", so it is not
  # treated as empty here either.
  empty_formula <- function(sc) {
    if (!inherits(sc, "formula") || length(sc) != 2L) return(FALSE)
    tt <- tryCatch(stats::terms(sc), error = function(e) NULL)
    !is.null(tt) && length(attr(tt, "term.labels")) == 0L &&
      is.null(attr(tt, "offset"))
  }
  empty <- if (is.list(scope) && !inherits(scope, "formula")) {
    # A zero-length element offers nothing to enter, as `character()` does for
    # a whole scope; refusing one while accepting the other was arbitrary.
    all(vapply(scope, function(sc) {
      is.null(sc) || length(sc) == 0L || empty_formula(sc)
    }, logical(1)))
  } else {
    length(scope) == 0L || empty_formula(scope)
  }
  one_sided <- function(sc, what) {
    if (inherits(sc, "formula") && length(sc) == 3L) {
      stop(what, " must be one-sided: its left-hand side (`",
           paste(deparse(sc[[2L]]), collapse = " "), "`) would be ignored, ",
           "so it would never be a candidate. Write every candidate on the ",
           "right, as in `~ ", paste(deparse(sc[[3L]]), collapse = " "), "`.",
           call. = FALSE)
    }
  }
  if (is.list(scope) && !inherits(scope, "formula")) {
    dup <- unique(names(scope)[duplicated(names(scope)) & nzchar(names(scope))])
    if (length(dup)) {
      stop("`scope` names ", paste0("`", dup, "`", collapse = ", "),
           " more than once; only the first entry would be read. Give each ",
           "phase one formula listing all its candidates.", call. = FALSE)
    }
    for (i in seq_along(scope)) {
      one_sided(scope[[i]], paste0("`scope$", names(scope)[i], "`"))
    }
  } else {
    one_sided(scope, "`scope`")
  }
  # After the structural checks, so a malformed scope gets the error that
  # names its fault (the left-hand side, a repeated phase) under every
  # direction; the backward refusal is about a well-formed scope.
  if (direction == "backward" && !empty) {
    remedy <- if (caller == "hzr_bootstrap") {
      paste0("For a backward screen on each replicate, pass an empty ",
             "`scope` such as `~ 1`; to screen a candidate set, use ",
             "`direction = \"both\"` or `\"forward\"`. With `scope` unset, ",
             "hzr_bootstrap() does not select at all.")
    } else {
      paste0("Pass the full model as the base fit, protect terms with ",
             "`force_in`, and leave `scope` unset or empty.")
    }
    stop("`scope` has no effect when `direction = \"backward\"`: a backward ",
         "screen only drops terms the base model already has, so a variable ",
         "left out of `scope` is still dropped and one the base lacks is ",
         "never tested. ", remedy, call. = FALSE)
  }
  invisible(NULL)
}
