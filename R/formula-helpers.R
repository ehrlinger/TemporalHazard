#' Parse Surv() formula for hazard modeling
#'
#' Extracts time, status, time_lower, time_upper, and predictors from a formula
#' of the form `Surv(time, status) ~ x1 + x2 + ...`.
#' Supports right-censored, left-censored, interval-censored, and
#' counting-process (start-stop) data.
#'
#' For counting-process (start-stop) data, use `Surv(start, stop, event)`.
#' The start times are returned as `time_lower` and stop times as `time`,
#' enabling the likelihood to compute `H(stop) - H(start)` per epoch.
#'
#' `Surv()` and this package code censoring status differently, so the
#' returned `status` is translated, not passed through:
#'
#' | Meaning   | TemporalHazard | `Surv` "left" | `Surv` "interval" |
#' | --------- | -------------- | ------------- | ----------------- |
#' | left      | `-1`           | `0`           | `2`               |
#' | right     | `0`            | --            | `0`               |
#' | event     | `1`            | `1`           | `1`               |
#' | interval  | `2`            | --            | `3`               |
#'
#' Under `type = "interval"`, `Surv()` reuses the `time2` column to hold the
#' status of any non-interval row, so an upper bound is read only where the
#' row is genuinely interval-censored.
#'
#' @param formula A formula object with Surv() on the LHS.
#' @param data A data frame containing variables referenced in the formula.
#' @return A list with elements: time, status, time_lower, time_upper, x
#' @keywords internal
.hzr_parse_formula <- function(formula, data) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }

  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object.", call. = FALSE)
  }

  # Parse the LHS (should be Surv(...))
  lhs <- formula[[2L]]
  rhs <- formula[[3L]]

  # Extract Surv() call
  if (!is.call(lhs) || !grepl("Surv$", deparse(lhs[[1L]]))) {
    stop("Formula LHS must be a Surv() call.", call. = FALSE)
  }

  # Evaluate the entire Surv() call in the context of data
  # Make sure Surv is available in the evaluation environment
  env_with_surv <- list2env(c(as.list(data), list(Surv = survival::Surv)))
  surv_obj <- eval(lhs, envir = env_with_surv)

  # Surv() returns a Surv object; extract the matrix and attributes
  if (!inherits(surv_obj, "Surv")) {
    stop("Formula LHS must return a Surv object.", call. = FALSE)
  }

  resp <- .hzr_surv_response(surv_obj)

  # Parse RHS (predictors)
  x <- NULL
  if (!is.null(rhs)) {
    # One-sided formula for model.matrix(), with `.` expanded against `data`
    # without the Surv() variables (#273). See .hzr_expand_rhs().
    rhs_formula <- .hzr_expand_rhs(formula, data)
    tryCatch({
      x <- stats::model.matrix(rhs_formula, data = data)
      x_contrasts <- attr(x, "contrasts")
      # Remove intercept column if present
      if (ncol(x) > 0 && colnames(x)[1L] == "(Intercept)") {
        x <- x[, -1L, drop = FALSE]
      }
      if (ncol(x) == 0) {
        x <- NULL
      }
    }, error = function(e) {
      stop("Failed to parse formula RHS: ", e$message, call. = FALSE)
    })
  }

  # What predict(newdata = ) needs to rebuild `x` from new rows: the terms,
  # the factor levels and the contrasts seen at fit time (as predict.lm()
  # keeps them). The terms come from model.frame() so that they carry
  # `predvars`: scale(x) and poly(x, 2) then reuse the fit's centre, scale
  # and basis at new rows instead of recomputing them from those rows.
  x_design <- NULL
  if (!is.null(x)) {
    mf <- stats::model.frame(rhs_formula, data = data)
    x_terms <- attr(mf, "terms")
    x_design <- list(
      terms = x_terms,
      xlevels = stats::.getXlevels(x_terms, mf),
      contrasts = x_contrasts,
      # The formula's variables that were columns of `data`: newdata must
      # supply exactly these. Any other variable (`cutoff` in
      # I(x > cutoff)) comes from the formula's environment.
      # .hzr_mask_symbols(), not all.vars(): the name after `$` is never
      # looked up, so cfg$time must not make a `time` column required.
      data_vars = intersect(
        .hzr_mask_symbols(stats::formula(x_terms)[[2L]]), names(data)
      )
    )
  }

  list(
    time = resp$time,
    status = resp$status,
    time_lower = resp$time_lower,
    time_upper = resp$time_upper,
    x = x,
    x_design = x_design,
    surv_type = resp$surv_type
  )
}


#' Write out `.` in a model formula's right-hand side
#'
#' Returns the right-hand side of a two-sided `Surv(...) ~ ...` formula as a
#' one-sided formula, with `.` expanded to every column of `data` that the
#' left-hand side does not use, as `survival::coxph()` does: `terms()` drops
#' those columns itself. Both the global formula (`.hzr_parse_formula()`,
#' #273) and each `hzr_phase(formula = )` (`hazard()`, #277) go through here,
#' so `.` means the same thing in both.
#'
#' @param formula A two-sided formula with the `Surv()` term on the left.
#' @param data The data frame `.` is expanded against.
#' @return A one-sided formula in `environment(formula)`, with no `.` left.
#' @keywords internal
#' @noRd
.hzr_expand_rhs <- function(formula, data) {
  rhs_formula <- stats::formula(
    stats::delete.response(stats::terms(formula, data = data))
  )
  # When the response uses every column, terms() has nothing to put in
  # place of `.` and leaves it, and model.matrix() would then expand it
  # against all of `data`, response included. Alone, `.` then stands for
  # no column and the model has no covariates. Beside other terms it is
  # refused: replacing the RHS would drop those terms without a word.
  if ("." %in% all.vars(rhs_formula)) {
    if (!identical(all.vars(rhs_formula), ".")) {
      stop("`.` in the formula stands for no column: `data` holds only ",
           "the variables of the Surv() response. Remove `.`, or add the ",
           "covariates to `data`.", call. = FALSE)
    }
    rhs_formula <- stats::reformulate("1", env = environment(rhs_formula))
  }
  rhs_formula
}


#' Read a Surv object into this package's response vectors
#'
#' The one place a `survival::Surv()` object is translated, called by both
#' interfaces: `.hzr_parse_formula()` for the formula's left-hand side, and
#' `hazard()` when a `Surv` is passed as `status`. Keeping a single copy is
#' the point -- the vector path used to take the second column unchanged,
#' which misread left-censored rows as right-censored, and for `"interval"`
#' and `"counting"` is not the status column at all (#226).
#'
#' The translation is driven by `attr(surv_obj, "type")`, never by the codes
#' observed: a right-censored vector and a `"left"` one both hold only 0 and
#' 1, with different meanings. `type = "interval2"` arrives here as
#' `"interval"`, because `Surv()` converts it.
#'
#' @param surv_obj A `Surv` object.
#' @return A list with `time`, `status`, `time_lower`, `time_upper` (either
#'   bound may be `NULL`) and `surv_type`.
#' @keywords internal
#' @noRd
.hzr_surv_response <- function(surv_obj) {
  surv_type <- attr(surv_obj, "type")
  surv_mat <- unclass(surv_obj)

  if (surv_type == "right") {
    # Format: [time, status]
    time <- surv_mat[, 1L]
    status <- surv_mat[, 2L]
    time_lower <- NULL
    time_upper <- NULL
  } else if (surv_type == "left") {
    # Format: [time, status], where Surv codes 1 = event, 0 = left-censored.
    # TemporalHazard codes left-censoring as -1; passing 0 through would
    # silently read these rows as right-censored.
    time <- surv_mat[, 1L]
    status <- ifelse(surv_mat[, 2L] == 0, -1, 1)
    time_lower <- NULL
    time_upper <- surv_mat[, 1L]
  } else if (surv_type == "interval") {
    # Format: [time1, time2, status], where Surv codes
    #   0 = right-censored, 1 = event, 2 = left-censored, 3 = interval.
    # Map onto TemporalHazard's 0 / 1 / -1 / 2.
    surv_status <- surv_mat[, 3L]
    status <- c(0, 1, -1, 2)[surv_status + 1L]
    time <- surv_mat[, 1L]

    # Surv stores the status in `time2` for every non-interval row, so that
    # column is a sentinel except where surv_status == 3.
    is_interval <- surv_status == 3
    time_upper <- ifelse(is_interval, surv_mat[, 2L], time)

    # `time_lower` doubles as the counting-process entry time for status
    # 0 and 1, where the likelihood forms H(stop) - H(start).  Surv
    # "interval" carries no left truncation, so entry time is 0 outside the
    # interval rows; using `time` there would cancel those rows out.
    time_lower <- ifelse(is_interval, surv_mat[, 1L], 0)
  } else if (surv_type == "counting") {
    # Start-stop (counting process) format: Surv(start, stop, event)
    # Used for repeating events / epoch-decomposed longitudinal data.
    # Each epoch contributes H(stop) - H(start) to the likelihood.
    time_lower <- surv_mat[, 1L]  # entry (start) time
    time <- surv_mat[, 2L]        # exit (stop) time
    time_upper <- NULL
    status <- surv_mat[, 3L]
  } else {
    stop("Unsupported Surv() type: ", surv_type, call. = FALSE)
  }

  list(
    time = time,
    status = status,
    time_lower = time_lower,
    time_upper = time_upper,
    surv_type = surv_type
  )
}


#' Refuse `predict(newdata = )` where `newdata`'s `time` would be misread
#'
#' A formula symbol is looked up in `newdata` first, so a model whose
#' formula uses a variable named `time` -- a data column or a constant such
#' as `I(age > time)` -- would read `newdata$time` for it (#270). For the
#' time-based types that column is the prediction time, so any such variable
#' is refused. For the eta-based types (no prediction time) a `time` data
#' column is a genuine covariate and is read as one; only a `time` constant
#' is refused, and only when `newdata` has a `time` column to mask it.
#'
#' @param object A fitted `hazard` object.
#' @param newdata The `newdata` data frame.
#' @param time_based `TRUE` when `newdata$time` is the prediction time.
#' @return `NULL`, invisibly; stops when `newdata$time` would be misread.
#' @keywords internal
#' @noRd
.hzr_check_time_covariate <- function(object, newdata, time_based) {
  design <- object$data$x_design
  x_cols <- colnames(object$data$x)
  from_data <- if (!is.null(design)) design$data_vars else x_cols
  # The symbols a formula looks up; .hzr_mask_symbols() skips the name after
  # `$`, so cfg$time is not a variable `time`.
  rhs_symbols <- function(f) .hzr_mask_symbols(f[[length(f)]])
  used <- if (is.null(design) || .hzr_uses_design_columns(object, newdata)) {
    # Design columns are used as they are and nothing is evaluated, so only
    # a design column literally named `time` collides.
    x_cols
  } else {
    # What model.frame() evaluates: the predvars the fit recorded (where
    # scale()'s centre is already a number), otherwise the formula.
    pv <- attr(design$terms, "predvars")
    if (is.null(pv)) {
      rhs_symbols(stats::formula(design$terms))
    } else {
      .hzr_mask_symbols(pv)
    }
  }
  phases <- object$fit$phases
  if (is.null(phases)) phases <- object$spec$phases
  for (nm in names(phases)) {
    # Only a phase the fit built from its own formula evaluates it; one the
    # fit ignored (a vector-interface fit) cannot misread `time`. The
    # routing helper decides, as it does for predict().
    ph <- phases[[nm]]
    if (!is.null(ph$formula) && !.hzr_phase_inherits_global(object, nm)) {
      used <- c(used, rhs_symbols(ph$formula))
    }
  }
  misread <- if (time_based) {
    "time" %in% used
  } else {
    "time" %in% names(newdata) && "time" %in% setdiff(used, from_data)
  }
  if (misread) {
    stop("The model uses a variable named 'time', but in 'newdata' the ",
         "column 'time' is the prediction time or would stand in for that ",
         "variable, so the model cannot be evaluated there. This stops ",
         "the survival, cumulative hazard and multiphase predictions, and ",
         "hzr_gof() and (for a single-distribution model) hzr_deciles(), ",
         "which build 'newdata' themselves. Rename the variable and refit.",
         call. = FALSE)
  }
  invisible(NULL)
}


#' Covariate design for `predict(newdata = )` on a single-distribution fit
#'
#' Matches `newdata`'s covariates to the fit by name, through
#' `.hzr_global_design()`, so their column order does not matter (#267).
#'
#' @param object A fitted `hazard` object.
#' @param newdata Data frame of new rows, possibly with a `time` column.
#' @return `NULL` when `newdata` has no covariate columns (the time-based
#'   types then evaluate the baseline, every covariate at 0); otherwise a
#'   numeric matrix of the fit's global columns. An object that stored no
#'   `x` leaves position as the only mapping, so its columns are taken as
#'   given.
#' @keywords internal
#' @noRd
.hzr_newdata_design <- function(object, newdata, drop_time = TRUE) {
  covs <- newdata[, names(newdata) != "time", drop = FALSE]
  if (is.null(object$data$x)) {
    return(if (ncol(covs) == 0L) NULL else as.matrix(covs))
  }
  # Where `time` is not the prediction time (the eta-based types without
  # time windows), it may itself be the covariate.
  if (!drop_time) covs <- newdata
  if (ncol(covs) == 0L) {
    return(NULL)
  }
  .hzr_global_design(object, newdata)
}


#' The newdata frame a formula design is rebuilt from
#'
#' Shared by the global and the per-phase rebuild, so the two follow one
#' rule: `newdata` supplies only the model's `data` columns. Every other
#' formula symbol (a constant such as `cutoff` in `I(age > cutoff)`, spline
#' knots, or an object kept outside `data`) resolves from the formula's
#' environment, so a same-named `newdata` column can never stand in for it.
#' A term that took row-level values from outside `data` is then refused by
#' `.hzr_check_equivariant()`.
#'
#' @param newdata Data frame of new rows.
#' @param data_vars The formula's fitting-data variables.
#' @return `newdata` with only the `data_vars` columns.
#' @keywords internal
#' @noRd
.hzr_newdata_frame <- function(newdata, data_vars) {
  newdata[, intersect(names(newdata), data_vars), drop = FALSE]
}


#' Build a newdata design, refusing terms that do not follow newdata's rows
#'
#' A term built from `newdata`'s columns is row-wise: shifting `newdata`'s
#' rows shifts the design's rows the same way. A term that took row-level
#' values from outside `data` (a vector, matrix, list or environment in the
#' formula's environment) does not move with them, so the design is rebuilt
#' with the rows cyclically shifted and any column that does not follow is
#' refused, naming its term. This does not depend on the shape of the
#' outside object. One row cannot be shifted; the row-count backstop
#' (`.hzr_check_design_rows()`) covers it.
#'
#' @param build Function of a data frame of new rows, returning the model
#'   matrix with its `assign` attribute.
#' @param newdata Data frame of new rows.
#' @param labels The terms' labels, indexed by `assign`.
#' @param where Text naming the design, for messages.
#' @return The model matrix built from `newdata`.
#' @keywords internal
#' @noRd
.hzr_check_equivariant <- function(build, newdata, labels, where) {
  mm <- build(newdata)
  .hzr_check_design_rows(mm, newdata, where)
  n <- nrow(newdata)
  if (n < 2L) {
    return(mm)
  }
  s <- c(2:n, 1L)
  ms <- build(newdata[s, , drop = FALSE])
  .hzr_check_design_rows(ms, newdata, where)
  b <- mm[s, , drop = FALSE]
  # Each column's tolerance is scaled to its spread in newdata, so a small
  # variation is still resolved, plus a few ulps of its magnitude for
  # summation-order noise (scale(age) averages reordered rows). The limit:
  # outside-data variation in a column below sqrt(eps) times its spread
  # plus ~64 ulps of its magnitude is not detected (a column wholly from
  # outside `data` is caught down to the ulps; one mixing a data column
  # with an outside part, as in I(age + zt), only down to sqrt(eps) of the
  # spread). The error then left in that column is within the same band.
  tol <- vapply(seq_len(ncol(b)), function(j) {
    f <- b[is.finite(b[, j]), j]
    if (length(f) == 0L) {
      return(0)
    }
    sqrt(.Machine$double.eps) * (max(f) - min(f)) +
      64 * .Machine$double.eps * max(abs(f))
  }, numeric(1))
  same <- abs(ms - b) <= rep(tol, each = nrow(b))
  same[is.na(ms) & is.na(b)] <- TRUE
  same[is.na(same)] <- FALSE
  bad <- which(colSums(!same) > 0L)
  if (length(bad) > 0L) {
    term <- unique(c("(Intercept)", labels)[attr(mm, "assign")[bad] + 1L])
    stop("term ", paste0("'", term, "'", collapse = ", "), " of ", where,
         " uses row-level values taken from outside `data`; ",
         "predict(newdata =) cannot rebuild them for new rows. Move them ",
         "into `data` as columns and refit.", call. = FALSE)
  }
  mm
}


#' Refuse a rebuilt design whose rows are not newdata's
#'
#' A backstop: a formula variable taken from anywhere but `newdata` (a
#' fitting-length vector in the environment, say) would give the design the
#' fitting rows, and every prediction would silently be for the wrong rows.
#'
#' @param mm The rebuilt model matrix.
#' @param newdata Data frame of new rows.
#' @param where Text naming the design, for messages.
#' @return `NULL`, invisibly; stops on a row-count mismatch.
#' @keywords internal
#' @noRd
.hzr_check_design_rows <- function(mm, newdata, where) {
  if (nrow(mm) != nrow(newdata)) {
    stop("The design rebuilt for ", where, " has ", nrow(mm), " rows for ",
         nrow(newdata), " row(s) of 'newdata': a term uses row-level values ",
         "taken from outside `data`, which predict(newdata =) cannot rebuild ",
         "for new rows. Move them into `data` as columns and refit.",
         call. = FALSE)
  }
  invisible(NULL)
}


#' Does a multiphase phase inherit the global design?
#'
#' The fit builds a phase from its own formula only on the formula
#' interface; otherwise the phase takes the global `x` (window-expanded
#' under `time_windows`). `predict(newdata = )` and `hzr_gof()` both need to
#' know which, and both call this one helper so they cannot drift.
#'
#' The fit's record, `attr(fit$x_list, "from_formula")`, decides when it is
#' there. A fit saved before that record (no tag through v1.2.9 has it)
#' falls back to its columns: a phase with no formula inherits, and a phase
#' with one inherits only if its stored columns are exactly the inherited
#' ones -- the window-expanded names under `time_windows`, otherwise
#' `colnames(data$x)`. That keeps a 1.0.3-era fit, which kept neither the
#' record nor its data, on the phase formula it was fitted with. Known
#' limit: a pre-record phase formula whose columns are literally the
#' inherited names; without windows both routes then build the same design.
#'
#' @param object A fitted multiphase `hazard` object.
#' @param nm Phase name.
#' @return A single logical: `TRUE` if the phase inherits the global design.
#' @keywords internal
#' @noRd
.hzr_phase_inherits_global <- function(object, nm) {
  from_formula <- attr(object$fit$x_list, "from_formula")
  if (!is.null(from_formula)) {
    return(!isTRUE(unname(from_formula[nm])))
  }
  phases <- object$fit$phases
  if (is.null(phases)) phases <- object$spec$phases
  if (is.null(phases[[nm]]$formula)) {
    return(TRUE)
  }
  time_windows <- object$spec$time_windows
  inherited <- if (!is.null(time_windows) && !is.null(object$data$x)) {
    colnames(.hzr_expand_time_varying_design(
      x = object$data$x[1, , drop = FALSE], time = 0,
      time_windows = time_windows
    ))
  } else {
    colnames(object$data$x)
  }
  identical(colnames(object$fit$x_list[[nm]]), inherited)
}


#' Rebuild the formula design of a fit saved before it was stored
#'
#' A formula fit from 1.2.10 or earlier kept no `data$x_design`, so
#' `predict(newdata = )` cannot tell a formula variable from an unused
#' column, and refuses any column but the design columns and `time`
#' (`.hzr_uses_design_columns()`). Such a fit usually kept what the design
#' is built from: the call's formula, the bindings that call references
#' (`call_env`) and the data frame it was fitted on (`data$frame`, stored
#' since 1.1.0). The formula is parsed again against that frame, as
#' `hazard()` parsed it, and the result is used only if it reproduces the
#' fitted design `data$x`: the same columns, in the same order, with the
#' same values. A design with a term whose row-level values come from
#' outside `data` reproduces the fitted rows but cannot follow new ones, so
#' it is not used either (`.hzr_check_equivariant()` on the fitting frame).
#' Anything missing, an error, or any difference leaves the fit as it was,
#' so the refusal stands (#301).
#'
#' Only a formula the call itself fixes is used: a formula object stored in
#' the call, a literal `~` expression, or a name bound in `call_env` itself.
#' A name bound only beyond `call_env` (a top-level `f`) means whatever it
#' means now, not what the fit used, and a cutoff with no fitted row between
#' the old and new value reproduces `data$x` exactly. A computed formula
#' (`as.formula(...)`) would run again, side effects and all. What the
#' formula looks up is held to the same rule (`.hzr_formula_trusted()`).
#'
#' A vector-interface fit has no formula and is returned unchanged.
#'
#' @param object A fitted `hazard` object.
#' @return `object`, with `data$x_design` filled in when it was rebuilt.
#' @keywords internal
#' @noRd
.hzr_recover_x_design <- function(object) {
  x_fit <- object$data$x
  frame <- object$data$frame
  env <- object$call_env
  f <- object$call$formula
  if (!is.null(object$data$x_design) || is.null(x_fit) || is.null(f) ||
        !is.environment(env) || !is.data.frame(frame)) {
    return(object)
  }
  literal <- FALSE
  if (is.symbol(f)) {
    if (!exists(as.character(f), envir = env, inherits = FALSE)) {
      return(object)
    }
    f <- get(as.character(f), envir = env, inherits = FALSE)
  } else if (!inherits(f, "formula")) {
    if (!is.call(f) || !identical(f[[1L]], as.name("~"))) {
      return(object)
    }
    # Evaluating a literal `~` looks nothing up: it only makes the formula,
    # with the captured bindings as its environment. (Parsing it below
    # evaluates its terms, as hazard() did.)
    f <- eval(f, env)
    literal <- TRUE
  }
  # The same rule for what the formula itself looks up: only what cannot
  # have changed since the fit. See .hzr_formula_trusted().
  if (!.hzr_formula_trusted(f, frame, if (literal) env)) {
    return(object)
  }
  # Warnings are muffled: whether the rebuild is right is decided by the
  # checks below, not by what re-parsing said on the way.
  parsed <- tryCatch(suppressWarnings(.hzr_parse_formula(f, frame)),
                     error = function(e) NULL)
  design <- parsed$x_design
  x <- parsed$x
  reproduces <- !is.null(design) && is.matrix(x) &&
    identical(dim(x), dim(x_fit)) &&
    identical(colnames(x), colnames(x_fit)) &&
    isTRUE(all(x == x_fit | (is.na(x) & is.na(x_fit))))
  if (reproduces) {
    reproduces <- tryCatch({
      .hzr_check_equivariant(.hzr_design_builder(design), frame,
                             attr(design$terms, "term.labels"), "the model")
      TRUE
    }, error = function(e) FALSE)
  }
  if (reproduces) {
    object$data$x_design <- design
  }
  object
}


#' Does a formula look up only what cannot have changed since the fit?
#'
#' A rebuilt legacy design is checked against the fitted rows, and that
#' check cannot see a lookup that moved between two fitted values: a cutoff
#' `k` in `I(age > k)` moved to a value no fitted age separates from the
#' old one, or a user's `thr()` redefined the same way, reproduces `data$x`
#' exactly and then predicts with the new value at new rows. So every name
#' the right-hand side looks up, other than a data column, must resolve to
#' something fixed: a function or constant of R or a package (base, a
#' namespace, an attached `package:` environment, or `pkg::fn`), or, for a
#' formula written in the call, a value `hazard()` copied into `call_env`
#' when it was fitted. A name found anywhere else -- the workspace, a
#' function frame kept with the formula (saved as it stood when the object
#' was saved, not when it was fitted), a user's function -- is not trusted.
#'
#' @param f A two-sided formula.
#' @param frame The fitting data frame.
#' @param call_env The fit's `call_env` when `f` was written in the call,
#'   otherwise `NULL`.
#' @return A single logical.
#' @keywords internal
#' @noRd
.hzr_formula_trusted <- function(f, frame, call_env = NULL) {
  fenv <- environment(f)
  if (!is.environment(fenv)) {
    return(FALSE)
  }
  where <- function(nm, mode) {
    e <- fenv
    while (!identical(e, emptyenv())) {
      if (exists(nm, envir = e, mode = mode, inherits = FALSE)) {
        return(e)
      }
      e <- parent.env(e)
    }
    NULL
  }
  from_r <- function(e) {
    !is.null(e) && (identical(e, baseenv()) || isNamespace(e) ||
                      startsWith(environmentName(e), "package:"))
  }
  ok_value <- function(nm) {
    e <- where(nm, "any")
    from_r(e) || (!is.null(call_env) && identical(e, call_env) &&
                    !is.function(get(nm, envir = e, inherits = FALSE)))
  }
  ok_function <- function(nm) from_r(where(nm, "function"))
  looks <- .hzr_formula_lookups(f[[length(f)]])
  values <- setdiff(looks$values[nzchar(looks$values)], c(names(frame), "."))
  all(vapply(values, ok_value, logical(1))) &&
    all(vapply(looks$functions, ok_function, logical(1)))
}


#' The names a formula's right-hand side looks up, by kind
#'
#' Like `.hzr_mask_symbols()`, but keeps the names called as functions
#' (`thr` in `thr(age)`), which that helper leaves out, apart from the
#' values. The operand after `$` or `@` is never looked up. A `pkg::fn`
#' call is a package's own function and adds no name.
#'
#' @param e A language object, symbol or constant.
#' @return A list of two character vectors, `values` and `functions`.
#' @keywords internal
#' @noRd
.hzr_formula_lookups <- function(e) {
  if (is.symbol(e)) {
    return(list(values = as.character(e), functions = character(0)))
  }
  if (!is.call(e)) {
    return(list(values = character(0), functions = character(0)))
  }
  head <- e[[1L]]
  parts <- as.list(e)[-1L]
  functions <- character(0)
  if (is.symbol(head)) {
    h <- as.character(head)
    if (h %in% c("::", ":::")) {
      return(list(values = character(0), functions = character(0)))
    }
    functions <- h
    if (h %in% c("$", "@") && length(e) >= 3L) {
      parts <- parts[1L]
    }
  } else {
    parts <- c(list(head), parts)
  }
  sub <- lapply(parts, .hzr_formula_lookups)
  list(
    values = unique(unlist(lapply(sub, `[[`, "values"), use.names = FALSE)),
    functions = unique(c(functions, unlist(lapply(sub, `[[`, "functions"),
                                           use.names = FALSE)))
  )
}


#' The function that builds a formula design at new rows
#'
#' The fit's terms, levels and contrasts, applied to `newdata`'s data
#' columns only (`.hzr_newdata_frame()`). One builder, so a design rebuilt
#' for a legacy fit is checked with the builder `predict()` then uses.
#'
#' @param design A stored `x_design`.
#' @return A function of a data frame of new rows, returning the model
#'   matrix with its `assign` attribute.
#' @keywords internal
#' @noRd
.hzr_design_builder <- function(design) {
  function(x) {
    nd <- .hzr_newdata_frame(x, design$data_vars)
    mf <- stats::model.frame(design$terms, data = nd, xlev = design$xlevels,
                             na.action = stats::na.pass)
    stats::model.matrix(design$terms, data = mf,
                        contrasts.arg = design$contrasts)
  }
}


#' Is the global design taken from `newdata`'s design columns?
#'
#' The design is then used as it is and the formula is not re-evaluated.
#' That needs every fitted design column present by name, and then either
#' the caller declaring the columns design-level (`hzr_deciles()` and
#' `hzr_gof()` set the `hzr_design_columns` attribute: their rows are
#' fitted design rows or their means, which the formula cannot rebuild), or
#' the formula's variables being absent (a fit saved before `x_design`, or
#' newdata given as design columns only). Otherwise the variables win, so a
#' design-named column that contradicts them cannot override them (#272).
#' Some formula variables beside the design columns, with others missing,
#' is an error: neither route could honour what was given. A variable that
#' is itself a design column (numeric `age`) counts as given only when
#' another column is built from it (`I(age^2)`, `age:grpyoung`).
#' Shared by the rebuild and the `time` check so the two cannot disagree
#' about the route.
#'
#' @param object A fitted `hazard` object.
#' @param newdata Data frame of new rows.
#' @return A single logical.
#' @keywords internal
#' @noRd
.hzr_uses_design_columns <- function(object, newdata) {
  cols <- colnames(object$data$x)
  if (is.null(cols) || !all(cols %in% names(newdata))) {
    return(FALSE)
  }
  if (isTRUE(attr(newdata, "hzr_design_columns"))) {
    return(TRUE)
  }
  design <- object$data$x_design
  if (is.null(design)) {
    # No stored design: a vector-interface fit (no formula, so no variable a
    # column could contradict), or a formula fit saved before the design was
    # stored whose design predict() could not rebuild exactly
    # (.hzr_recover_x_design()). For the latter an extra column cannot be
    # told from a formula
    # variable, and taking the design columns beside a contradicting one
    # would be a silent wrong answer, so it is refused -- as it was before,
    # when such a column made the positional match fail.
    extra <- setdiff(names(newdata), c(cols, "time"))
    if (!is.null(object$call$formula) && length(extra) > 0L) {
      stop("This fit was saved by an earlier version of TemporalHazard, ",
           "without a stored formula design, so 'newdata' may hold only its ",
           "design columns (", paste0("'", cols, "'", collapse = ", "),
           ") and 'time'; it also has ",
           paste0("'", extra, "'", collapse = ", "), ". Refit the model ",
           "with the current version, or pass only the design columns.",
           call. = FALSE)
    }
    return(TRUE)
  }
  missing <- setdiff(design$data_vars, names(newdata))
  if (length(missing) == 0L) {
    return(FALSE)
  }
  # The design route would ignore any formula variable given beside the
  # design columns, and without the missing ones the variables cannot be
  # rebuilt, so a mix is refused rather than guessed (#272). A variable
  # that is itself a design column (numeric `age`) counts only if another
  # term is built from it (I(age^2), age:grp): a changed `age` would leave
  # those columns stale.
  feeds_derived <- unlist(lapply(
    attr(design$terms, "term.labels"),
    function(label) {
      v <- all.vars(parse(text = label)[[1L]])
      if (identical(v, label)) character(0) else v
    }
  ))
  counted <- union(setdiff(design$data_vars, cols),
                   intersect(design$data_vars, feeds_derived))
  given <- intersect(counted, names(newdata))
  if (length(given) > 0L) {
    stop("'newdata' gives the formula variable(s) ",
         paste0("'", given, "'", collapse = ", "), " but lacks ",
         paste0("'", missing, "'", collapse = ", "),
         ", while carrying the fitted design columns. Give all of the ",
         "formula's variables, so the design can be rebuilt from them.",
         call. = FALSE)
  }
  TRUE
}


#' Rebuild the global design matrix at new rows
#'
#' Used by `predict()` for a multiphase phase without its own formula, which
#' inherits the global design. Returns the fit's global columns, in the fit's
#' order, whatever other columns `newdata` carries.
#'
#' @param object A fitted `hazard` object.
#' @param newdata Data frame of new rows.
#' @return Numeric matrix with `nrow(newdata)` rows and the columns of
#'   `object$data$x`.
#' @keywords internal
#' @noRd
.hzr_global_design <- function(object, newdata) {
  x_fit <- object$data$x
  cols <- colnames(x_fit)
  design <- object$data$x_design

  # Design columns taken as they are: hzr_deciles() and hzr_gof() pass them
  # (a factor as `grpyoung`, a transform as `log(age)`), and so does a
  # newdata that lacks the formula's variables. When the variables are
  # there, they are rebuilt below instead (#272).
  if (.hzr_uses_design_columns(object, newdata)) {
    return(as.matrix(newdata[, cols, drop = FALSE]))
  }

  # Otherwise newdata gives the formula's variables. A covariate must come
  # from newdata itself: were it looked up in the formula's environment, a
  # stray `mal` in the workspace would silently stand in for a missing
  # column. A fit saved by an earlier version stored no design, so its
  # design columns are all there is to match.
  needed <- if (is.null(design)) cols else design$data_vars
  missing <- setdiff(needed, names(newdata))
  if (length(missing) > 0L) {
    stop("'newdata' lacks the covariate column(s) ",
         paste0("'", missing, "'", collapse = ", "),
         " that the model uses. Columns are matched by name.", call. = FALSE)
  }

  if (!is.null(design)) {
    # Formula interface: the same terms, levels and contrasts as the fit, so
    # a factor given as a single label still codes to the fit's columns.
    # newdata supplies only the data columns, and a term that does not
    # follow its rows is refused; see .hzr_newdata_frame() and
    # .hzr_check_equivariant().
    mm <- .hzr_check_equivariant(.hzr_design_builder(design), newdata,
                                 attr(design$terms, "term.labels"),
                                 "the model")
    return(mm[, cols, drop = FALSE])
  }
  if (!is.null(cols)) {
    # Vector interface with a named `x`: select by name.
    return(as.matrix(newdata[, cols, drop = FALSE]))
  }
  # Vector interface with an unnamed `x`: position is all there is.
  nd <- newdata[, names(newdata) != "time", drop = FALSE]
  if (ncol(nd) != ncol(x_fit)) {
    stop("The fit's global design has ", ncol(x_fit), " unnamed column(s); ",
         "'newdata' has ", ncol(nd), " covariate column(s).", call. = FALSE)
  }
  as.matrix(nd)
}


#' Collect the symbols a masked argument would look up
#'
#' Like [base::all.vars()], but skips the `name` operand of `$` and `@`, which
#' `all.vars()` reports as a variable: `all.vars(quote(df$tt))` is
#' `c("df", "tt")` even though `tt` is never looked up. Counting it makes the
#' ambiguity warning name a column the fit did not use, and makes `data$col`,
#' the remedy that warning prescribes, trigger the warning.
#'
#' @param e A language object, symbol or constant.
#' @return Character vector of symbol names, possibly empty.
#' @keywords internal
#' @noRd
.hzr_mask_symbols <- function(e) {
  if (is.symbol(e)) {
    return(as.character(e))
  }
  if (!is.call(e)) {
    return(character(0))
  }
  head <- e[[1L]]
  if (is.symbol(head) && as.character(head) %in% c("$", "@") &&
        length(e) >= 3L) {
    return(.hzr_mask_symbols(e[[2L]]))
  }
  parts <- as.list(e)[-1L]
  if (!is.symbol(head)) {
    parts <- c(list(head), parts)
  }
  unique(unlist(lapply(parts, .hzr_mask_symbols), use.names = FALSE))
}

#' Is a name bound anywhere between a frame and the global environment?
#'
#' `exists(inherits = FALSE)` sees only the immediate frame, so a wrapper that
#' forwards its own argument (`g <- function(d) hazard(data = d, time = tt)`
#' with `tt` bound one frame out) looks unambiguous when it is not.
#' `inherits = TRUE` goes too far the other way, reaching package namespaces
#' and base, where a column named `c`, `t` or `df` would match on every call.
#' This walks the lexical parents up to and including [globalenv()] and stops
#' before the search path.
#'
#' @param nm Character name to look for.
#' @param env Environment to start from.
#' @return `TRUE` if `nm` is bound in `env` or a lexical parent up to the
#'   global environment, otherwise `FALSE`.
#' @keywords internal
#' @noRd
.hzr_bound_locally <- function(nm, env) {
  while (!identical(env, emptyenv())) {
    is_global <- identical(env, globalenv())
    # Namespaces, attached packages and base are named; a frame or a plain
    # local environment is not. Stop before the search path.
    if (!is_global && (isNamespace(env) || nzchar(environmentName(env)))) {
      return(FALSE)
    }
    if (exists(nm, envir = env, inherits = FALSE)) {
      return(TRUE)
    }
    if (is_global) {
      return(FALSE)
    }
    env <- parent.env(env)
  }
  FALSE
}
