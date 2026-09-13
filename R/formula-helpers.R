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
      data_vars = intersect(all.vars(x_terms), names(data))
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
  for (ph in phases) {
    if (!is.null(ph$formula)) used <- c(used, rhs_symbols(ph$formula))
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


#' Is the global design taken from `newdata`'s design columns?
#'
#' TRUE when `newdata` carries every column of the fitted global design by
#' name (as `hzr_deciles()` and `hzr_gof()` pass it). The design is then
#' used as it is and the formula is not re-evaluated. Shared by the rebuild
#' and the `time` check so the two cannot disagree about the route.
#'
#' @param object A fitted `hazard` object.
#' @param newdata Data frame of new rows.
#' @return A single logical.
#' @keywords internal
#' @noRd
.hzr_uses_design_columns <- function(object, newdata) {
  cols <- colnames(object$data$x)
  !is.null(cols) && all(cols %in% names(newdata))
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

  # newdata that already carries the fit's design columns by name is taken
  # as it is. That is how hzr_deciles() and hzr_gof() pass it, built from
  # object$data$x, so a factor arrives as `grpyoung` and a transform as
  # `log(age)`. For a plain numeric covariate the design column and the
  # variable are the same column, so both routes agree.
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
    mf <- stats::model.frame(design$terms, data = newdata,
                             xlev = design$xlevels, na.action = stats::na.pass)
    mm <- stats::model.matrix(design$terms, data = mf,
                              contrasts.arg = design$contrasts)
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
