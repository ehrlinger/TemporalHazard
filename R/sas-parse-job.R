# sas-parse-job.R -- map SAS censoring statements onto this package's status
# coding. Later job-translation tasks add to this file.
#
# This package codes censoring -1 left, 0 right, 1 event, 2 interval.
# survival::Surv() uses DIFFERENT integers for the same meanings under
# type = "interval" (0/1/2/3 for right/event/left/interval). Carrying one
# coding into the other is a wrong answer with no error, and this package has
# already shipped exactly that bug once (Surv(type = "left") read
# left-censored rows as right-censored). Never mix the two vocabularies.

#' Translate SAS censoring statements to this package's status coding.
#'
#' Builds an unevaluated `status` expression from the `EVENT`/`ICENSOR`/
#' `LCENSOR`/`RCENSOR` operands of a `PROC HAZARD` statements list, ready to
#' drop into a `hazard()` call. Priority, most specific first: `ICENSOR` (only
#' where the event did not occur), then `LCENSOR` (also only where the event
#' did not occur), defaulting to the bare `EVENT` variable, which is already
#' `0`/`1` for right-censored/event and so also covers `RCENSOR` with no
#' extra branch (see below).
#'
#' `EVENT` is optional: the reference `HAZARD` program (`src/hazard/varterm.c`)
#' terminates only when *both* `EVENT` and `ICENSOR` are missing, so a job
#' that specifies `ICENSOR` alone is legitimate. In that case the base status
#' is right-censored (`0`) everywhere, overridden to `2` where `ICENSOR`'s
#' lower bound is non-missing. A job specifying neither is rejected with an
#' error, mirroring HAZARD's own termination.
#'
#' An event outranks a censoring flag: a row that is an event is not left- or
#' interval-censored, so both `LCENSOR` and `ICENSOR` override only where
#' `EVENT == 0`. Treating a malformed row where `EVENT` and a censoring flag
#' both fire as an event is the conservative reading; this is an assumption
#' of this translation, not something confirmed against a SAS reference run.
#'
#' When `LCENSOR` and `ICENSOR` both fire for the same row, `ICENSOR` is
#' applied last and wins, because it is applied after `LCENSOR` below. SAS's
#' own precedence for this combination is undefined; this is this package's
#' documented tie-break.
#'
#' @param statements Named list of SAS statement operands, e.g.
#'   `list(EVENT = "DEAD", TIME = "T", ICENSOR = c("LO", "HI"))`. `ICENSOR`
#'   carries two operands (lower, upper bound variables); `LCENSOR` and
#'   `RCENSOR` carry one. `EVENT` may be absent if `ICENSOR` is present.
#' @return `list(status_expr = <call>, time_lower = <name|NULL>,
#'   time_upper = <name|NULL>, untranslated = <data.frame>)`. `status_expr`
#'   is a `bquote()`-built call, evaluable against an environment/list
#'   holding the named SAS variables.
#' @noRd
.hzr_censor_spec <- function(statements) {
  has_event <- !is.null(statements$EVENT)
  has_icensor <- !is.null(statements$ICENSOR)

  if (!has_event && !has_icensor) {
    stop("The EVENT or ICENSOR variable must be specified.", call. = FALSE)
  }

  ev <- if (has_event) as.name(statements$EVENT)
  expr <- if (has_event) ev else 0

  # RCENSOR needs no branch: right-censoring is code 0, which is exactly what
  # EVENT already carries when the event did not occur. This is deliberate,
  # not an omission -- there is nothing for an RCENSOR flag to change.

  if (!is.null(statements$LCENSOR)) {
    flag <- as.name(statements$LCENSOR)
    # An event outranks a left-censoring flag (see roxygen); LCENSOR only
    # overrides where the event did not occur, symmetric with ICENSOR below.
    expr <- if (has_event) {
      bquote(ifelse(.(ev) == 0 & .(flag) == 1, -1, .(expr)))
    } else {
      bquote(ifelse(.(flag) == 1, -1, .(expr)))
    }
  }

  if (has_icensor) {
    lo <- as.name(statements$ICENSOR[[1L]])
    # Interval censoring only applies where the event did not occur; if it
    # did, EVENT == 1 and the row stays coded as an event, not an interval.
    # Applied last, so ICENSOR wins over LCENSOR where both would fire (see
    # roxygen).
    expr <- if (has_event) {
      bquote(ifelse(.(ev) == 0 & !is.na(.(lo)), 2, .(expr)))
    } else {
      bquote(ifelse(!is.na(.(lo)), 2, .(expr)))
    }
  }

  list(
    status_expr = expr,
    time_lower = if (has_icensor) as.name(statements$ICENSOR[[1L]]),
    time_upper = if (has_icensor) as.name(statements$ICENSOR[[2L]]),
    untranslated = .hzr_untranslated_frame()
  )
}

#' Parse a PROC HAZARD block into a hazard() or hzr_stepwise() call.
#'
#' Every keyword is resolved through .hzr_sas_token() in its lexer context.
#' A keyword that does not resolve is recorded in `untranslated`, never
#' skipped: a dropped option changes the model silently.
#' @noRd
.hzr_parse_hazard <- function(block) {
  st <- strsplit(block$text, ";", fixed = TRUE)[[1L]]
  untr <- .hzr_untranslated_frame()
  seen <- 0L
  mapped <- 0L
  note <- function(kw, reason) {
    untr <<- rbind(untr, .hzr_untranslated_frame(NA_integer_, kw, reason))
  }

  # --- statement 1: the PROC line and its options -------------------------
  toks <- strsplit(trimws(st[[1L]]), " ", fixed = TRUE)[[1L]]
  toks <- toks[nzchar(toks)]
  ctl <- list()
  data_name <- NULL
  outhaz <- NULL

  for (tok in toks) {
    eqp <- .idx(tok, "=")
    key <- if (eqp > 0L) substring(tok, 1L, eqp - 1L) else tok
    val <- if (eqp > 0L) substring(tok, eqp + 1L) else ""
    token <- .hzr_sas_token(key, "HAZARD", "HZRP")
    if (identical(token, "PROC") || identical(token, "HAZARD")) next
    seen <- seen + 1L
    if (is.na(token)) {
      note(key, "unknown PROC HAZARD option")
      next
    }
    mapped <- mapped + 1L
    switch(token,
      DATA        = data_name <- val,
      OUTHAZ      = outhaz <- val,
      MAXITER     = {
        val_num <- suppressWarnings(as.numeric(val))
        if (is.na(val_num)) {
          mapped <- mapped - 1L
          note("MAXITER", "non-numeric value for MAXITER")
        } else {
          ctl$maxit <- val_num
        }
      },
      CONDITION   = {
        val_num <- suppressWarnings(as.numeric(val))
        if (is.na(val_num)) {
          mapped <- mapped - 1L
          note("CONDITION", "non-numeric value for CONDITION")
        } else {
          ctl$condition <- val_num
        }
      },
      CONSERVE    = ctl$conserve <- TRUE,
      NOCONSERVE  = ctl$conserve <- FALSE,
      QUASINEWTON = ctl$method <- "bfgs",
      STEEPEST    = {
        mapped <- mapped - 1L
        note("STEEPEST", "no R equivalent for steepest descent (issue #145)")
      },
      # Print-only options. Mapped, because "no R effect" is the correct
      # translation, not a gap.
      PRINTIT = NULL, NOPRINT = NULL, NOCOR = NULL, NOCOV = NULL,
      NOLOG = NULL, NONOTES = NULL,
      {
        mapped <- mapped - 1L
        note(key, "no hazard() equivalent")
      }
    )
  }

  # --- statements 2..n ----------------------------------------------------
  statements <- list()
  parms_ops <- character(0)
  sel_ops <- NULL
  covars <- list()

  for (i in seq_along(st)[-1L]) {
    stmt_text <- trimws(st[[i]])
    w <- strsplit(stmt_text, " ", fixed = TRUE)[[1L]]
    w <- w[nzchar(w)]
    if (!length(w)) next
    kw <- w[[1L]]
    ops <- w[-1L]
    # EARLY/CONSTANT/LATE operands are a comma-separated VAR=VALUE list, not
    # whitespace-separated bare names -- keep the raw tail text (commas and
    # all) for .hzr_parse_parms()/.hzr_parse_phase_covars() to split.
    ops_text <- trimws(substring(stmt_text, nchar(kw) + 1L))
    token <- .hzr_sas_token(kw, "HAZARD", "STMT")
    seen <- seen + 1L
    if (is.na(token)) {
      note(kw, "unknown HAZARD statement")
      next
    }
    mapped <- mapped + 1L
    switch(token,
      TIME       = statements$TIME <- ops[[1L]],
      EVENT      = statements$EVENT <- ops[[1L]],
      ICENSOR    = statements$ICENSOR <- ops,
      LCENSOR    = statements$LCENSOR <- ops[[1L]],
      RCENSOR    = statements$RCENSOR <- ops[[1L]],
      WEIGHT     = statements$WEIGHT <- ops[[1L]],
      PARAMETERS = parms_ops <- ops,
      STEPWISE   = sel_ops <- ops,
      EARLY      = covars$early <- ops_text,
      CONSTANT   = covars$constant <- ops_text,
      LATE       = covars$late <- ops_text,
      {
        mapped <- mapped - 1L
        note(kw, "no R equivalent")
      }
    )
  }

  # --- assemble -----------------------------------------------------------
  parms <- .hzr_parse_parms(parms_ops, covars = covars)
  untr <- rbind(untr, parms$untranslated)
  cens <- .hzr_censor_spec(statements)
  untr <- rbind(untr, cens$untranslated)

  args <- list()
  if (!is.null(data_name)) args$data <- as.name(data_name)
  args$time <- as.name(statements$TIME)
  args$status <- cens$status_expr
  if (!is.null(cens$time_lower)) args$time_lower <- cens$time_lower
  if (!is.null(cens$time_upper)) args$time_upper <- cens$time_upper
  args$phases <- parms$phases
  args$theta <- parms$theta
  if (!is.null(statements$WEIGHT)) args$weights <- as.name(statements$WEIGHT)

  # Canonical control order, so the emitted call does not depend on the order
  # the options happened to appear in the SAS text.
  ctl <- ctl[intersect(c("maxit", "condition", "conserve", "method"),
                       names(ctl))]
  if (length(ctl)) args$control <- as.call(c(quote(list), ctl))

  head <- quote(hazard)
  if (!is.null(sel_ops)) {
    sel <- .hzr_selection_spec(sel_ops)
    untr <- rbind(untr, sel$untranslated)
    if (!is.null(sel$direction)) {
      head <- quote(hzr_stepwise)
      args$direction <- sel$direction
      if (!is.null(sel$slentry)) args$slentry <- sel$slentry
      if (!is.null(sel$slstay)) args$slstay <- sel$slstay
    }
  }

  list(call = as.call(c(head, args)), outhaz = outhaz,
       untranslated = untr, tokens_seen = seen, tokens_mapped = mapped)
}

#' Translate the DATA step that builds a HAZPRED prediction grid.
#'
#' Returns an unevaluated `data.frame()` call, or `NULL` when the step is not
#' one of the two stereotyped forms this package can read. `NULL` means
#' untranslated, never "no grid": a `predict()` call with no `newdata` is a
#' hollow result -- right shape, empty inside -- so the caller
#' (`.hzr_parse_hazpred()`) must record it in `untranslated`, not treat it as
#' nothing to translate.
#' @noRd
.hzr_parse_grid <- function(txt, name) {
  if (is.null(name) || !nzchar(name)) return(NULL)
  p <- .idx(txt, paste0("DATA ", name, ";"))
  if (p == 0L) return(NULL)

  rest <- substring(txt, p + 1L)
  ends <- c(.idx(rest, "DATA "), .idx(rest, "PROC "), .idx(rest, "%HAZ"))
  ends <- ends[ends > 0L]
  body <- if (length(ends)) substring(rest, 1L, min(ends)) else rest

  time_var <- local({
    m <- regmatches(body, regexec("([A-Z_][A-Z0-9_]*) *= *EXP\\(", body))[[1L]]
    if (length(m) > 1L) m[[2L]] else NA_character_
  })

  # --- log-spaced grid: DO LN_TIME=-5 TO LN_MAX BY INC; t = EXP(LN_TIME) ----
  if (grepl(" DO ", body) && grepl("LOG(", body, fixed = TRUE) &&
      !is.na(time_var)) {
    lo_txt <- local({
      m <- regmatches(body,
             regexec("DO [A-Z_][A-Z0-9_]* *= *(-?[0-9.]+) TO", body))[[1L]]
      if (length(m) > 1L) m[[2L]] else NA_character_
    })
    hi_txt <- local({
      m <- regmatches(body, regexec("MAX *= *([0-9.]+)", body))[[1L]]
      if (length(m) > 1L) m[[2L]] else NA_character_
    })
    if (is.na(lo_txt) || is.na(hi_txt)) return(NULL)
    lo <- suppressWarnings(as.numeric(lo_txt))
    hi <- suppressWarnings(as.numeric(hi_txt))
    if (is.na(lo) || is.na(hi)) return(NULL)
    # SAS writes INC=(5+LN_MAX)/99.9, i.e. 100 points inclusive. Splice the
    # matched text through str2lang(), not the numeric value: a negative
    # bound such as -5 parses (like source code) to a unary-minus call, not
    # a bare negative double, and only str2lang() reproduces that so the
    # emitted call matches what quote()ing the equivalent source produces.
    inner <- bquote(exp(seq(.(str2lang(lo_txt)), log(.(str2lang(hi_txt))),
                             length.out = 100)))
    cl <- as.call(list(quote(data.frame), inner))
    names(cl) <- c("", time_var)
    return(cl)
  }

  # --- explicit DO list: DO MONTHS=1,2,3,6,12,24 TO 180 BY 12; --------------
  if (grepl(" DO ", body)) {
    m <- regmatches(body,
           regexec("DO ([A-Z_][A-Z0-9_]*) *= *([^;]+);", body))[[1L]]
    if (length(m) < 3L) return(NULL)
    var <- m[[2L]]
    parts <- trimws(strsplit(m[[3L]], ",", fixed = TRUE)[[1L]])
    elems <- list()
    for (part in parts) {
      rng <- regmatches(part,
               regexec("^([0-9.]+) +TO +([0-9.]+)( +BY +([0-9.]+))?$",
                       part))[[1L]]
      if (length(rng) >= 3L) {
        lo_txt <- rng[[2L]]
        hi_txt <- rng[[3L]]
        by_txt <- if (length(rng) >= 5L && nzchar(rng[[5L]])) rng[[5L]] else "1"
        lo <- suppressWarnings(as.numeric(lo_txt))
        hi <- suppressWarnings(as.numeric(hi_txt))
        by <- suppressWarnings(as.numeric(by_txt))
        # A range bound that failed to parse as numeric (rare, given the
        # digit-only capture group above, but not impossible -- e.g. a
        # malformed literal like "1.2.3"). Refuse the whole grid rather than
        # coerce silently to NA.
        if (is.na(lo) || is.na(hi) || is.na(by)) return(NULL)
        elems[[length(elems) + 1L]] <- bquote(
          seq(.(str2lang(lo_txt)), .(str2lang(hi_txt)), by = .(str2lang(by_txt)))
        )
      } else if (grepl("^[0-9.]+$", part)) {
        val <- suppressWarnings(as.numeric(part))
        if (is.na(val)) return(NULL)
        elems[[length(elems) + 1L]] <- str2lang(part)
      } else {
        # An element we cannot read (a SAS expression such as 1*DTY). Refuse
        # the whole grid rather than emit a partial one.
        return(NULL)
      }
    }
    inner <- as.call(c(quote(c), elems))
    cl <- as.call(list(quote(data.frame), inner))
    names(cl) <- c("", var)
    return(cl)
  }

  # SET-derived or anything else: not translatable.
  NULL
}

#' Parse a PROC HAZPRED block into predict() call(s).
#'
#' HAZPRED emits `_SURVIV`/`_CLLSURV`/`_CLUSURV` and
#' `_HAZARD`/`_CLLHAZ`/`_CLUHAZ`, so it maps to `predict()` call(s) with
#' `se.fit`. `conf.type` is deliberately NOT set: whether SAS's confidence
#' limits are log-log or logit is unresolved (spec section 8, open question
#' 2), and guessing would produce silently wrong bounds -- so this emits the
#' package default instead.
#'
#' `txt` is the whole normalised source, because HAZPRED's real input is the
#' `DATA=` prediction grid built by a preceding DATA step, not anything in
#' the PROC block itself. A grid that cannot be translated (`.hzr_parse_grid()`
#' returns `NULL`) is always recorded in `untranslated` here: a `predict()`
#' with no `newdata` is a hollow result, not "no grid needed".
#' @noRd
.hzr_parse_hazpred <- function(block, txt) {
  st <- strsplit(block$text, ";", fixed = TRUE)[[1L]]
  untr <- .hzr_untranslated_frame()
  seen <- 0L
  mapped <- 0L
  note <- function(kw, reason) {
    untr <<- rbind(untr, .hzr_untranslated_frame(NA_integer_, kw, reason))
  }

  toks <- strsplit(trimws(st[[1L]]), " ", fixed = TRUE)[[1L]]
  toks <- toks[nzchar(toks)]
  data_name <- NULL
  inhaz <- NULL
  want_surv <- TRUE
  want_haz <- TRUE
  want_cl <- TRUE

  for (tok in toks) {
    eqp <- .idx(tok, "=")
    key <- if (eqp > 0L) substring(tok, 1L, eqp - 1L) else tok
    val <- if (eqp > 0L) substring(tok, eqp + 1L) else ""
    token <- .hzr_sas_token(key, "HAZPRED", "HZPP")
    if (identical(token, "PROC") || identical(token, "HAZPRED")) next
    seen <- seen + 1L
    if (is.na(token)) {
      note(key, "unknown PROC HAZPRED option")
      next
    }
    mapped <- mapped + 1L
    switch(token,
      DATA    = data_name <- val,
      INHAZ   = inhaz <- val,
      OUT     = NULL,
      NOSURV  = want_surv <- FALSE,
      NOHAZ   = want_haz <- FALSE,
      NOCL    = want_cl <- FALSE,
      CLIMITS = want_cl <- TRUE,
      NOLOG = NULL, NONOTES = NULL,
      {
        mapped <- mapped - 1L
        note(key, "no predict() equivalent")
      }
    )
  }

  for (i in seq_along(st)[-1L]) {
    w <- strsplit(trimws(st[[i]]), " ", fixed = TRUE)[[1L]]
    w <- w[nzchar(w)]
    if (!length(w)) next
    kw <- w[[1L]]
    token <- .hzr_sas_token(kw, "HAZPRED", "STMT")
    seen <- seen + 1L
    if (is.na(token)) {
      note(kw, "unknown HAZPRED statement")
      next
    }
    mapped <- mapped + 1L
    if (!(token %in% c("TIME", "ID"))) {
      mapped <- mapped - 1L
      note(kw, "SAS listing control; no R effect")
    }
  }

  grid <- .hzr_parse_grid(txt, data_name)
  if (is.null(grid) && !is.null(data_name)) {
    note(paste0("DATA=", data_name),
         "prediction grid DATA step is not one of the translatable forms")
  }

  mk <- function(type) {
    args <- list(quote(fit))
    if (!is.null(data_name)) args$newdata <- as.name(data_name)
    args$type <- type
    args$se.fit <- want_cl
    as.call(c(quote(predict), args))
  }

  list(
    call = if (want_surv) mk("survival") else mk("hazard"),
    call_haz = if (want_surv && want_haz) mk("hazard") else NULL,
    inhaz = inhaz, grid = grid, untranslated = untr,
    tokens_seen = seen, tokens_mapped = mapped
  )
}
