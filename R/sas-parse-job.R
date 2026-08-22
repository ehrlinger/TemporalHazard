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
#' Builds an unevaluated `status` expression -- and, where needed, a
#' `time_lower` expression -- from the `EVENT`/`ICENSOR`/`LCENSOR`/`RCENSOR`
#' operands of a `PROC HAZARD` statements list, ready to drop into a
#' `hazard()` call.
#'
#' **`LCENSOR` is left-*truncation*, not left-censoring.** Its operand is the
#' counting-process *entry time* -- all four corpus uses are literally
#' `LCENSOR STARTTME`, paired with `TIME INT_TE` (see
#' `inst/dev/FIXTURE-GAP-LIST.md`, answered Q1). It never changes `status`:
#' HAZARD has no left-censoring in this package's sense, so `status = -1` is
#' never produced by this translator. A future reader must not "restore" an
#' `LCENSOR` -> `-1` branch.
#'
#' **`ICENSOR`'s first operand is an event *count* (OBS column 4, `C3`), not a
#' 0/1 flag.** Its grammar (`src/hazard/hazard_y.y`) is
#' `ICENSOR c3var '=' ctimevar;`, and the likelihood in `setlik.c` uses a
#' Nelson-type approximation `C3 * ln([CF(T) - CF(CT)] / (T - CT))`, so a row
#' is interval-censored where `C3 > 0`. The second operand (`CTIME`) is the
#' interval's *lower* bound -- the interval runs `CTIME` -> `TIME` (see
#' `FIXTURE-GAP-LIST.md` Q1). `EVENT` is optional: the reference `HAZARD`
#' program (`src/hazard/varterm.c`) terminates only when *both* `EVENT` and
#' `ICENSOR` are missing, so a job that specifies `ICENSOR` alone is
#' legitimate. In that case the base status is right-censored (`0`)
#' everywhere, overridden to `2` where `C3 > 0`. A job specifying neither is
#' rejected with an error, mirroring HAZARD's own termination.
#'
#' Because `C3` is a count, it also carries a *weight*. `setlik.c` multiplies
#' the whole log-likelihood contribution of an interval row by it
#' (`c3w = c3 * weight`, then both `-(c1w + c2 + c3w) * (cumhaz - cumhst)` and
#' `+= c3w * lct`), which is exactly what [hazard()]'s `weights` argument does
#' per row. So `weights_expr` is emitted as
#' `ifelse(<status> == 2, C3, 1)` whenever `ICENSOR` is present. Two details
#' there are not negotiable. Passing the bare `C3` instead would give every
#' event and right-censored row a weight of *zero* -- `C3` is 0 on those rows
#' -- deleting them from the likelihood without a word; `readc2.c` sets
#' `C2 = 1` on exactly those rows when no `RCENSOR` statement is given, so
#' their weight is 1. And the gate is on `status`, not on `C3 > 0`, so that a
#' row where `EVENT` and `C3` both fire keeps the event's weight, matching the
#' status mapping directly below.
#'
#' A `WEIGHT` statement, if also present, multiplies that expression, since
#' `c3w` is the product of the two.
#'
#' An event outranks interval censoring: a row that is an event is not
#' interval-censored, so `ICENSOR` overrides only where `EVENT == 0`. Treating
#' a malformed row where `EVENT` and `C3 > 0` both fire as an event is the
#' conservative reading; this is an assumption of this translation, not
#' something confirmed against a SAS reference run.
#'
#' `hazard()`'s `time_lower` carries a different meaning per row, selected by
#' `status` (see `R/hazard_api.R`): for status 2 (interval) it is the
#' interval's lower bound (`CTIME`), but for status 0/1 it is the
#' counting-process *entry time* (default `0`), not a censoring bound. `TIME`
#' is always the interval's upper bound, and that is exactly what
#' `hazard()`'s `time_upper` already defaults to when left `NULL` -- so
#' `time_upper` is never emitted by this translator; passing it would be
#' redundant, and omitting it cannot drift out of sync with `TIME`.
#' `time_lower` is therefore built as one of, depending on which statements
#' are present:
#' * `ICENSOR` only: `ifelse(status == 2, CTIME, 0)` (`0` is `hazard()`'s
#'   documented default entry time for status 0/1 rows).
#' * `LCENSOR` only: the `LCENSOR` variable directly -- it applies to every
#'   row, not just interval ones, so no `ifelse()` gating is needed.
#' * Neither: `time_lower` is omitted (`NULL`).
#'
#' **`LCENSOR` and `ICENSOR` together are refused, not translated.** There is
#' no expression that serves: `time_lower` is one column carrying two
#' meanings, so gating it on status hands the interval rows their lower bound
#' and thereby drops their entry time, fitting a left-truncated
#' interval-censored subject as at risk from time `0`. That converges and
#' reports plausibly -- the exact failure mode this package exists to refuse.
#' The reference implementation carries three distinct times (`TIME`, `CTIME`,
#' `STIME`) and subtracts `H(STIME)` for every row, interval rows included, so
#' translating it faithfully needs an entry-time argument [hazard()] does not
#' have. It is also rare: 0 of the 93 `PROC HAZARD` jobs across the production
#' studies use the combination. So the pair is recorded in `untranslated` and
#' `refused` is set, and the caller emits a `stop()` in place of the fit.
#'
#' The `ifelse()` forms are gated on a `status_name` placeholder (`.hzr_status`)
#' the caller assigns the `status_expr` to first -- never unconditionally,
#' which would trip `hazard()`'s finite-value check off the interval subset
#' and, where populated, silently redefine event/right-censored rows' risk-set
#' entry times.
#'
#' @param statements Named list of SAS statement operands, e.g.
#'   `list(EVENT = "DEAD", TIME = "T", ICENSOR = c("C3", "CTIME"))`.
#'   `ICENSOR` carries two operands (event-count variable, interval
#'   lower-bound time variable); `LCENSOR` and `RCENSOR` carry one. `EVENT`
#'   may be absent if `ICENSOR` is present.
#' @return `list(status_expr = <call>, status_name = <name|NULL>,
#'   time_lower = <call|NULL>, weights_expr = <call|NULL>,
#'   untranslated = <data.frame>, refused = <logical>)`. `status_expr` is
#'   a `bquote()`-built call, evaluable against an environment/list holding
#'   the named SAS variables. `time_lower` and `weights_expr`, when
#'   non-`NULL`, are `bquote()`-built calls; `status_name` is non-`NULL` only
#'   when `ICENSOR` is present, since only then does anything need to gate on
#'   it. `weights_expr` is likewise non-`NULL` only under `ICENSOR`, so a job
#'   without one emits no `weights` argument at all and fits with the unit
#'   weights it always did. `refused` is `TRUE` only for the `LCENSOR` +
#'   `ICENSOR` combination described above, and every other element is `NULL`
#'   when it is: there is nothing to emit. Both returns carry the same element
#'   names, so a caller reading one never meets an unexpected missing field.
#' @noRd
.hzr_censor_spec <- function(statements) {
  has_event <- !is.null(statements$EVENT)
  has_icensor <- !is.null(statements$ICENSOR)
  has_lcensor <- !is.null(statements$LCENSOR)

  if (!has_event && !has_icensor) {
    stop("The EVENT or ICENSOR variable must be specified.", call. = FALSE)
  }

  # hazard()'s time_lower has two meanings selected by status: entry time for
  # status 0/1, interval lower bound for status 2. One row cannot carry both,
  # so a left-truncated interval-censored subject would be fitted as at risk
  # from time 0 -- a converging, plausible, wrong answer. SAS carries three
  # distinct times (TIME, CTIME, STIME) and subtracts H(STIME) for every row.
  # Supporting this needs a new hazard() argument; refuse until it exists.
  if (has_lcensor && has_icensor) {
    return(list(
      status_expr = NULL,
      status_name = NULL,
      time_lower = NULL,
      weights_expr = NULL,
      untranslated = .hzr_untranslated_frame(
        NA_integer_, "LCENSOR + ICENSOR",
        paste("left truncation combined with interval censoring needs a",
              "separate entry-time argument hazard() does not have (#155);",
              "translate this job by hand")
      ),
      refused = TRUE
    ))
  }

  ev <- if (has_event) as.name(statements$EVENT)
  expr <- if (has_event) ev else 0

  # RCENSOR needs no branch: right-censoring is code 0, which is exactly what
  # EVENT already carries when the event did not occur. This is deliberate,
  # not an omission -- there is nothing for an RCENSOR flag to change.

  # LCENSOR is left-truncation, not left-censoring (see roxygen): it never
  # touches status, only time_lower below.

  if (has_icensor) {
    c3 <- as.name(statements$ICENSOR[[1L]])
    # C3 is an event count, not a 0/1 flag (see roxygen): a row is
    # interval-censored where C3 > 0. Interval censoring only applies where
    # the event did not occur; if it did, EVENT == 1 and the row stays coded
    # as an event, not an interval.
    expr <- if (has_event) {
      bquote(ifelse(.(ev) == 0 & .(c3) > 0, 2, .(expr)))
    } else {
      bquote(ifelse(.(c3) > 0, 2, .(expr)))
    }
  }

  status_name <- NULL
  time_lower <- NULL
  weights_expr <- NULL
  if (has_icensor) {
    status_name <- as.name(".hzr_status")
    ctime <- as.name(statements$ICENSOR[[2L]])
    # C3 is a count, and setlik.c multiplies the whole log-likelihood
    # contribution by it -- c3w = c3 * weight, then both
    # -(c1w + c2 + c3w) * (cumhaz - cumhst) and += c3w * lct. hazard()'s
    # `weights` multiplies each row's contribution the same way, so the count
    # belongs there; without it a row carrying 3 events is fitted as one
    # observation (#154). Status-gated, and emphatically not the bare count:
    # C3 is 0 on every non-interval row, and a weight of 0 deletes a row from
    # the likelihood outright. readc2.c sets C2 = 1 on exactly those rows
    # when no RCENSOR statement is given, so their weight is 1.
    weights_expr <- bquote(ifelse(.(status_name) == 2, .(c3), 1))
    # Status-gated, not unconditional (see roxygen): interval rows get the
    # interval's lower bound (CTIME); every other row falls back to
    # hazard()'s own entry-time default (0). LCENSOR cannot be present here
    # -- the pair is refused above -- so there is no other fallback to pick.
    time_lower <- bquote(ifelse(.(status_name) == 2, .(ctime), 0))
  } else if (has_lcensor) {
    # No ICENSOR, so no interval rows to gate around: LCENSOR's entry time
    # applies to every row, unconditionally.
    time_lower <- as.name(statements$LCENSOR)
  }

  list(
    status_expr = expr,
    status_name = status_name,
    time_lower = time_lower,
    weights_expr = weights_expr,
    untranslated = .hzr_untranslated_frame(),
    refused = FALSE
  )
}

#' Parse a PROC HAZARD block into a hazard() call, or a stop() call when the
#' block requests something this translator refuses (see .hzr_censor_spec()'s
#' LCENSOR + ICENSOR refusal and the SELECTION stepwise refusal below).
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
      ICENSOR    = {
        # ICENSOR's grammar (src/hazard/hazard_y.y) is
        # `ICENSOR c3var '=' ctimevar;` -- an event-COUNT variable (OBS
        # column 4, C3), not a 0/1 flag, and a second time variable (the
        # interval's lower bound), not two comma-separated bound variables.
        # Split on "=" and tolerate a stray trailing comma; anything else is
        # not this shape and must not be guessed at.
        parts <- trimws(sub(",$", "", strsplit(ops_text, "=", fixed = TRUE)[[1L]]))
        parts <- parts[nzchar(parts)]
        if (length(parts) == 2L) {
          statements$ICENSOR <- parts
        } else {
          mapped <- mapped - 1L
          note("ICENSOR", "expected 'ICENSOR count = timevar' operand shape")
        }
      },
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

  # The UNTRANSLATED callout alone is not enough here: a reader who renders
  # past it would still get a fit, and a fit over a mis-specified model is
  # this package's signature defect -- a result that looks like one and is
  # not. Emit a stop() in place of the hazard() call so the document fails
  # where the fit would have been, and say what has to be done by hand.
  # This return precedes SELECTION parsing, so a job carrying BOTH refusals
  # records only this one; it still stops loudly, and no corpus job does both.
  if (isTRUE(cens$refused)) {
    return(list(
      call = quote(stop(
        "This job combines LCENSOR (left truncation) with ICENSOR (interval ",
        "censoring). hazard()'s `time_lower` carries the entry time for ",
        "status 0/1 rows and the interval's lower bound for status 2 rows, ",
        "so one column cannot express both and any translation would fit ",
        "the interval-censored rows as at risk from time 0. Translate this ",
        "job by hand (#155)."
      )),
      status_call = NULL, outhaz = outhaz, untranslated = untr,
      tokens_seen = seen, tokens_mapped = mapped
    ))
  }

  # ICENSOR jobs get the status expression hoisted into its own chunk, named
  # .hzr_status so it cannot collide with a SAS variable, so that time_lower
  # can gate on it -- see .hzr_censor_spec() roxygen. LCENSOR-only jobs need
  # no gating (LCENSOR's entry time applies to every row unconditionally), so
  # status stays unhoisted there. time_upper is never emitted: the ICENSOR
  # interval's upper bound is TIME, and hazard()'s time_upper already
  # defaults to TIME when left NULL (see .hzr_censor_spec() roxygen).
  status_call <- NULL
  args <- list()
  if (!is.null(data_name)) args$data <- as.name(data_name)
  args$time <- as.name(statements$TIME)
  if (!is.null(cens$status_name)) {
    # The status chunk is evaluated outside hazard(), so hazard()'s data
    # masking does not reach it: wrap it in with() so its bare column names
    # resolve. A plain local binding (not a new column) keeps the caller's
    # data frame unmutated, and hazard()'s mask falls through to the caller's
    # frame to find it.
    status_expr <- if (is.null(data_name)) {
      cens$status_expr
    } else {
      call("with", as.name(data_name), cens$status_expr)
    }
    status_call <- call("<-", cens$status_name, status_expr)
    args$status <- cens$status_name
  } else {
    args$status <- cens$status_expr
  }
  if (!is.null(cens$time_lower)) args$time_lower <- cens$time_lower
  # Without fit = TRUE the emitted call returns an unfitted object: converged
  # is NA, objective is NA, and theta holds the SAS starting values, while
  # print.hazard() shows a populated summary that says none of that (#151).
  args$fit <- TRUE
  # A PARMS statement that actually specified shape parameters makes this a
  # multiphase job. hazard()'s `dist` defaults to "weibull", and its
  # `else if (!is.null(phases))` branch silently discards the entire phase
  # specification (with only a warning) when dist stays at that default --
  # so dist = "multiphase" must be emitted whenever phases were built. A job
  # with no PARMS statement at all (or one that specified no shape
  # parameters) has phase_calls = list(), i.e. parms$phases is the empty
  # `list()` call rather than NULL; omit both phases and dist on that path
  # so it stays a plain non-multiphase fit and does not trip that same
  # "'phases' is ignored" warning for a phases arg that was never meant to
  # carry anything.
  if (isTRUE(parms$has_phases)) {
    args$dist <- "multiphase"
    args$phases <- parms$phases
  }
  args$theta <- parms$theta
  # Two independent sources of a row weight, and they multiply, exactly as
  # setlik.c's c3w = c3 * weight does: the SAS WEIGHT statement's variable,
  # and ICENSOR's event count on interval rows (see .hzr_censor_spec()).
  wt_name <- if (!is.null(statements$WEIGHT)) as.name(statements$WEIGHT)
  if (!is.null(cens$weights_expr)) {
    args$weights <- if (is.null(wt_name)) {
      cens$weights_expr
    } else {
      bquote(.(cens$weights_expr) * .(wt_name))
    }
  } else if (!is.null(wt_name)) {
    args$weights <- wt_name
  }

  # Canonical control order, so the emitted call does not depend on the order
  # the options happened to appear in the SAS text.
  ctl <- ctl[intersect(c("maxit", "condition", "conserve", "method"),
                       names(ctl))]
  if (length(ctl)) args$control <- as.call(c(quote(list), ctl))

  head <- quote(hazard)
  if (!is.null(sel_ops)) {
    sel <- .hzr_selection_spec(sel_ops)
    untr <- rbind(untr, sel$untranslated)
    # A SELECTION statement that requests stepwise cannot be translated into
    # hzr_stepwise(): .hzr_refit_with_scope() (R/stepwise-refit.R) requires a
    # formula-interface base fit, and this translator always emits the
    # vector interface. Every candidate refit therefore errors, the forward
    # step downgrades those errors to warnings, and the run silently returns
    # zero steps -- indistinguishable from "nothing met slentry". That is
    # exactly the shipped-defect shape AGENTS.md warns about, so refuse
    # loudly instead, mirroring .hzr_censor_spec()'s LCENSOR + ICENSOR
    # refusal above: record an untranslated row and emit a stop() chunk in
    # place of the fit, rather than a hazard()/hzr_stepwise() call.
    if (isTRUE(sel$stepwise)) {
      untr <- rbind(untr, .hzr_untranslated_frame(
        NA_integer_, "SELECTION",
        paste("stepwise selection is not translated: hzr_stepwise()'s refit",
              "path needs a formula-interface base fit, and this translator",
              "emits the vector interface, so every candidate refit would",
              "error and the screen would silently report zero steps. Run",
              "this job's selection by hand (#160).")
      ))
      return(list(
        call = quote(stop(
          "This job's SELECTION statement requests a stepwise screen, which ",
          "is not translated: the stepwise-selection refit path needs a ",
          "formula-interface base fit, and this translator emits the ",
          "vector interface, so every candidate refit would error and the ",
          "screen would silently report zero steps -- indistinguishable ",
          "from \"nothing met slentry\". Run this job's selection by hand ",
          "(#160)."
        )),
        status_call = NULL, outhaz = outhaz, untranslated = untr,
        tokens_seen = seen, tokens_mapped = mapped
      ))
    }
  }

  list(call = as.call(c(head, args)), status_call = status_call,
       outhaz = outhaz, untranslated = untr, tokens_seen = seen,
       tokens_mapped = mapped)
}

#' Translate a SELECTION statement to hzr_stepwise() arguments.
#'
#' `hazard_y.y`'s `stepwisestmt` production is `STEPWISE stepwiseopts { setopt(33); }`,
#' and `stepwiseopts` can be empty -- the *statement* is what turns stepwise
#' on, not any particular direction keyword. So a bare `SELECTION;` (or one
#' carrying only `SLENTRY`/`SLSTAY`) legitimately enables stepwise; the
#' direction keywords only refine it. `ONEWAY` (aliases `NOSTEPWISE`/`NOSW`,
#' option 34) is the one option that turns stepwise back off.
#'
#' HAZARD's lexer also collapses FORWARD, FW, SW, SELECT and STEPWISE into
#' one token, which `hazard_y.y` maps to option 21; BACKWARD is option 22.
#' Option 21 is therefore two-way, so FORWARD maps to `direction = "both"`,
#' not `"forward"`. Mapping it to `"forward"` would silently change the
#' search on every stepwise job.
#'
#' Operands are resolved in `STEP` lexer context, matching HAZARD's own
#' `BEGIN STEP` start condition once inside a `SELECTION` statement. `SELECT`
#' is the one spelling that only resolves in `STMT` context (it is also an
#' alias for the statement keyword itself), so a `STEP`-context miss falls
#' back to `STMT` before being recorded as untranslated.
#'
#' When `ONEWAY`/`NOSTEPWISE`/`NOSW` disables stepwise, any `SLENTRY`/
#' `SLSTAY` also given are meaningless (there is no entry/stay search to
#' apply them to) and are moved into `untranslated` rather than silently
#' dropped -- discarding a parsed value with nothing recorded is exactly
#' the defect this package guards against.
#' @return `list(stepwise = <logical>, direction = <chr|NULL>,
#'   slentry = <dbl|NULL>, slstay = <dbl|NULL>, untranslated = <data.frame>)`.
#' @noRd
.hzr_selection_spec <- function(operands) {
  out <- list(stepwise = TRUE, direction = "both", slentry = NULL,
              slstay = NULL, untranslated = .hzr_untranslated_frame())
  for (op in operands) {
    eqp <- .idx(op, "=")
    key <- if (eqp > 0L) substring(op, 1L, eqp - 1L) else op
    val_txt <- if (eqp > 0L) substring(op, eqp + 1L) else NA_character_
    token <- .hzr_sas_token(key, "HAZARD", "STEP")
    if (is.na(token)) token <- .hzr_sas_token(key, "HAZARD", "STMT")
    if (is.na(token)) {
      out$untranslated <- rbind(out$untranslated, .hzr_untranslated_frame(
        NA_integer_, key, "unknown SELECTION option"
      ))
      next
    }
    switch(token,
      STEPWISE = out$direction <- "both",
      BACKWARD = out$direction <- "backward",
      ONEWAY   = {
        out$direction <- NULL
        out$stepwise <- FALSE
      },
      SLENTRY  = {
        val <- suppressWarnings(as.numeric(val_txt))
        if (is.na(val)) {
          out$untranslated <- rbind(out$untranslated, .hzr_untranslated_frame(
            NA_integer_, key, "non-numeric value for SLENTRY"
          ))
        } else {
          out$slentry <- val
        }
      },
      SLSTAY   = {
        val <- suppressWarnings(as.numeric(val_txt))
        if (is.na(val)) {
          out$untranslated <- rbind(out$untranslated, .hzr_untranslated_frame(
            NA_integer_, key, "non-numeric value for SLSTAY"
          ))
        } else {
          out$slstay <- val
        }
      },
      {
        out$untranslated <- rbind(out$untranslated, .hzr_untranslated_frame(
          NA_integer_, key, "no hzr_stepwise() equivalent"
        ))
      }
    )
  }
  if (!out$stepwise) {
    if (!is.null(out$slentry)) {
      out$untranslated <- rbind(out$untranslated, .hzr_untranslated_frame(
        NA_integer_, "SLENTRY",
        "stepwise disabled by ONEWAY/NOSTEPWISE/NOSW; SLENTRY has no effect"
      ))
      out$slentry <- NULL
    }
    if (!is.null(out$slstay)) {
      out$untranslated <- rbind(out$untranslated, .hzr_untranslated_frame(
        NA_integer_, "SLSTAY",
        "stepwise disabled by ONEWAY/NOSTEPWISE/NOSW; SLSTAY has no effect"
      ))
      out$slstay <- NULL
    }
  }
  out
}

#' Evaluate a SAS DATA-step numeric expression, using only known constants.
#'
#' Used to resolve explicit `DO` list elements such as `1*DTY`: `DTY` becomes
#' a plain number only when it was already folded from an earlier
#' `DTY=12/365.2425;` assignment in the *same* `DATA` step
#' (`.hzr_sas_data_constants()`). Anything else -- a function call, a
#' data-step variable, an unknown name -- must refuse rather than guess, so
#' `text` is checked against a strict whitelist regex (digits, the four
#' arithmetic operators, parentheses, whitespace, and known constant names)
#' *before* `parse()`/`eval()` ever see it; text that fails the whitelist is
#' never evaluated at all. `consts` doubles as the `eval()` environment, so
#' the only names that can resolve are exactly the ones the whitelist
#' allowed.
#' @param text A candidate numeric expression, e.g. `"1*DTY"` or `"24"`.
#' @param consts Named list of already-folded constants (name -> `double`).
#' @return A finite `double`, or `NA_real_` if `text` cannot be safely
#'   evaluated.
#' @noRd
.hzr_eval_sas_const <- function(text, consts) {
  text <- trimws(text)
  if (!nzchar(text)) return(NA_real_)
  nms <- names(consts)
  alt <- if (length(nms)) paste(nms[order(-nchar(nms))], collapse = "|") else "(?!)"
  whitelist <- sprintf("^(?:%s|[0-9.]+|[+*/()-]|[[:space:]]+)+$", alt)
  if (!grepl(whitelist, text, perl = TRUE)) return(NA_real_)

  # Safe by construction, not by luck: the whitelist above has already
  # rejected everything except digits/operators/parens/whitespace and
  # literal known-constant names, and `consts` (the only environment eval()
  # is given) holds nothing but doubles -- no function, including any base R
  # function, is reachable through it. There is no code path from here to
  # arbitrary execution.
  val <- tryCatch({
    expr <- parse(text = text, keep.source = FALSE)
    if (length(expr) != 1L) NA_real_ else eval(expr[[1L]], envir = consts)
  }, error = function(e) NA_real_)
  if (!is.numeric(val) || length(val) != 1L || !is.finite(val)) return(NA_real_)
  as.double(val)
}

#' Fold `NAME = <numeric expression>;` assignments into a constants map.
#'
#' Scoped to one `DATA` step's own preamble (the text before its `DO`
#' statement): collects assignments in order, so a later constant may
#' reference an earlier one (`INC=(5+LN_MAX)/99.9` after `LN_MAX=...`).
#' Only assignments `.hzr_eval_sas_const()` can actually evaluate end up in
#' the map -- an assignment whose right-hand side is not pure arithmetic
#' over already-known constants (a function call, an unresolved name) is
#' silently skipped here, not stored; it simply never becomes foldable.
#' @noRd
.hzr_sas_data_constants <- function(pre) {
  consts <- list()
  for (s in strsplit(pre, ";", fixed = TRUE)[[1L]]) {
    s <- trimws(s)
    if (!nzchar(s)) next
    m <- regmatches(s, regexec("^([A-Z_][A-Z0-9_]*) *= *(.+)$", s))[[1L]]
    if (length(m) < 3L) next
    val <- .hzr_eval_sas_const(trimws(m[[3L]]), consts)
    if (!is.na(val)) consts[[m[[2L]]]] <- val
  }
  consts
}

#' Read the step denominator out of a log-grid DATA step's own `INC=`.
#'
#' SAS writes the log-grid step as `INC = (<numerator>)/<denominator>`, and
#' the corpus uses three different denominators (`/49.9`, `/99.9`, `/999.9`).
#' The denominator is what sets both the step and the point count, so it is
#' read from the job rather than assumed.
#'
#' The numerator has to be the loop's own span, `log(hi) - lo`. Every corpus
#' job writes that span in one of two spellings -- `(5 + LN_MAX)` where the
#' loop starts at -5, or `(MAX - MIN)` -- and the two coincide only because
#' `lo = -5`. A numerator that is *not* the span means the emitted
#' `(log(hi) - lo)/denominator` step would not be the job's step, so this
#' returns `NA_real_` and the caller refuses the grid.
#'
#' @param pre The DATA step's text before its `DO` statement.
#' @param inc_var Name of the variable the `DO ... BY` clause steps by.
#' @param bound_var Name of the `DO ... TO` bound (e.g. `LN_MAX`).
#' @param ln_hi The bound's value, `log(hi)`.
#' @param span `log(hi) - lo`, the range the loop covers.
#' @return The denominator as a positive `double`, or `NA_real_` when the
#'   `INC=` assignment is absent or is not a form this can parse.
#' @noRd
.hzr_sas_grid_denom <- function(pre, inc_var, bound_var, ln_hi, span) {
  consts <- .hzr_sas_data_constants(pre)
  # LN_MAX = LOG(MAX) is a function call, so .hzr_sas_data_constants() never
  # folds it. The DO's TO bound is the log bound by construction, so bind it
  # here -- after the fold, so it wins over any earlier assignment.
  consts[[bound_var]] <- ln_hi

  rhs <- NA_character_
  for (s in strsplit(pre, ";", fixed = TRUE)[[1L]]) {
    s <- trimws(s)
    m <- regmatches(s, regexec("^([A-Z_][A-Z0-9_]*) *= *(.+)$", s))[[1L]]
    # Last assignment wins, as SAS's own sequential DATA-step semantics do.
    if (length(m) >= 3L && identical(m[[2L]], inc_var)) rhs <- trimws(m[[3L]])
  }
  if (is.na(rhs)) return(NA_real_)

  parts <- regmatches(rhs, regexec("^\\((.+)\\) */ *([0-9.]+)$", rhs))[[1L]]
  if (length(parts) < 3L) return(NA_real_)
  num <- .hzr_eval_sas_const(parts[[2L]], consts)
  den <- suppressWarnings(as.numeric(parts[[3L]]))
  if (is.na(num) || is.na(den) || den <= 0) return(NA_real_)
  if (!isTRUE(all.equal(num, span, tolerance = 1e-8))) return(NA_real_)
  den
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
    do_m <- regmatches(body, regexec(
      paste0("DO [A-Z_][A-Z0-9_]* *= *(-?[0-9.]+) TO ",
             "([A-Z_][A-Z0-9_]*) BY ([A-Z_][A-Z0-9_]*)"),
      body))[[1L]]
    # No `BY` at all means a SAS step of 1, which is not this stereotyped
    # form; refuse rather than reuse the step of the form it is not.
    if (length(do_m) < 4L) return(NULL)
    lo_txt <- do_m[[2L]]
    bound_var <- do_m[[3L]]
    inc_var <- do_m[[4L]]
    hi_txt <- local({
      # Anchor on a non-word character before MAX so LN_MAX=5.2 is not read
      # as MAX=5.2 -- that produced a grid ending at t = 5.2 instead of 181,
      # with no error (#153).
      m <- regmatches(body, regexec("(^|[^A-Z0-9_])MAX *= *([0-9.]+)", body))[[1L]]
      if (length(m) > 2L) m[[3L]] else NA_character_
    })
    if (is.na(hi_txt)) return(NULL)
    lo <- suppressWarnings(as.numeric(lo_txt))
    hi <- suppressWarnings(as.numeric(hi_txt))
    if (is.na(lo) || is.na(hi) || hi <= 0) return(NULL)
    span <- log(hi) - lo
    if (!is.finite(span) || span <= 0) return(NULL)
    do_at <- regexpr("DO [A-Z_][A-Z0-9_]* *= *[^;]+;", body)
    pre <- if (do_at > 0L) substring(body, 1L, as.integer(do_at) - 1L) else ""
    # The step is the job's own INC=, never an assumed one. Three
    # denominators appear across the corpus (/49.9, /99.9, /999.9), and
    # hardcoding /99.9 gave the /999.9 jobs a step ten times too large --
    # wrong times over a full progress bar, reported as fully translated.
    denom <- .hzr_sas_grid_denom(pre, inc_var, bound_var, log(hi), span)
    if (is.na(denom)) return(NULL)
    # SAS's DO LN_TIME = lo TO LN_MAX BY INC runs floor((LN_MAX - lo)/INC) + 1
    # times, and INC is span/denom, so the count follows the denominator:
    # /99.9 lands 100 points, /999.9 lands 1000. The last is
    # exp(lo + (n - 1) * INC), NOT LN_MAX -- seq(length.out = n) would imply
    # a /(n - 1) step, so only the first point would agree (#153).
    n <- floor(denom) + 1
    # Splice the matched text through str2lang(), not the numeric value: a
    # negative bound such as -5 parses (like source code) to a unary-minus
    # call, not a bare negative double, and only str2lang() reproduces that
    # so the emitted call matches what quote()ing the equivalent source
    # produces.
    lo_lang <- str2lang(lo_txt)
    inner <- bquote(
      exp(.(lo_lang) +
            seq(0, .(n - 1)) *
              ((log(.(str2lang(hi_txt))) - .(lo_lang)) / .(denom)))
    )
    cl <- as.call(list(quote(data.frame), inner))
    # predict.hazard() requires a column literally named `time`
    # (R/hazard_api.R:1006). Naming the grid column after the SAS DO variable
    # produced a newdata that predict() rejects outright, so the "grids
    # resolve" coverage figure counted grids that could not be used (#151).
    names(cl) <- c("", "time")
    return(cl)
  }

  # --- explicit DO list: DO MONTHS=1,2,3,6,12,24 TO 180 BY 12; --------------
  # or, with constants folded first: DO MONTHS=1*DTY,2*DTY,24 TO 180 BY 12;
  if (grepl(" DO ", body)) {
    do_at <- regexpr("DO [A-Z_][A-Z0-9_]* *= *[^;]+;", body)
    if (do_at == -1L) return(NULL)
    consts <- .hzr_sas_data_constants(substring(body, 1L, as.integer(do_at) - 1L))

    m <- regmatches(body,
           regexec("DO ([A-Z_][A-Z0-9_]*) *= *([^;]+);", body))[[1L]]
    if (length(m) < 3L) return(NULL)
    var <- m[[2L]]
    parts <- trimws(strsplit(m[[3L]], ",", fixed = TRUE)[[1L]])
    elems <- list()
    is_literal <- function(txt) grepl("^-?[0-9.]+$", txt)
    for (part in parts) {
      rng <- regmatches(part,
               regexec("^([A-Z0-9_.+*/()-]+) +TO +([A-Z0-9_.+*/()-]+)( +BY +([A-Z0-9_.+*/()-]+))?$",
                       part))[[1L]]
      if (length(rng) >= 3L) {
        lo_txt <- rng[[2L]]
        hi_txt <- rng[[3L]]
        by_txt <- if (length(rng) >= 5L && nzchar(rng[[5L]])) rng[[5L]] else "1"
        lo <- .hzr_eval_sas_const(lo_txt, consts)
        hi <- .hzr_eval_sas_const(hi_txt, consts)
        by <- .hzr_eval_sas_const(by_txt, consts)
        # A range bound this cannot be evaluated at all -- unknown name, a
        # function call, or a malformed literal. Refuse the whole grid
        # rather than coerce silently to NA or emit a partial grid.
        if (is.na(lo) || is.na(hi) || is.na(by)) return(NULL)
        # Splice bare numeric literals through str2lang(), not the double:
        # a negative bound such as -5 parses (like source code) to a
        # unary-minus call, not a bare negative double, and only str2lang()
        # reproduces that. A folded constant expression (e.g. 1*DTY) has no
        # such source form to preserve -- it becomes the evaluated number.
        elems[[length(elems) + 1L]] <- bquote(
          seq(.(if (is_literal(lo_txt)) str2lang(lo_txt) else lo),
              .(if (is_literal(hi_txt)) str2lang(hi_txt) else hi),
              by = .(if (is_literal(by_txt)) str2lang(by_txt) else by))
        )
      } else {
        val <- .hzr_eval_sas_const(part, consts)
        # An element that cannot be evaluated at all -- an unknown name, a
        # function call, or a data-step variable. Refuse the whole grid
        # rather than emit a partial one.
        if (is.na(val)) return(NULL)
        elems[[length(elems) + 1L]] <- if (is_literal(part)) str2lang(part) else val
      }
    }
    inner <- as.call(c(quote(c), elems))
    cl <- as.call(list(quote(data.frame), inner))
    # predict.hazard() requires a column literally named `time`
    # (R/hazard_api.R:1006). Naming the grid column after the SAS DO variable
    # produced a newdata that predict() rejects outright, so the "grids
    # resolve" coverage figure counted grids that could not be used (#151).
    names(cl) <- c("", "time")
    return(cl)
  }

  # SET-derived or anything else: not translatable.
  NULL
}

#' Parse a PROC HAZPRED block into predict() call(s).
#'
#' HAZPRED emits `_SURVIV`/`_CLLSURV`/`_CLUSURV` and
#' `_HAZARD`/`_CLLHAZ`/`_CLUHAZ`, so it maps to `predict()` call(s) with
#' `se.fit`.
#'
#' `conf.type = "logit"` is set on the survival call, because SAS HAZPRED's
#' survival confidence limits are on the logit scale (`hzp_calc_srv_CL.c`),
#' while `predict.hazard()` defaults to `"log-log"` (the survfit standard).
#' Emitting the default would produce bounds that silently disagree with the
#' job being reproduced. The hazard call sets nothing: hazard limits are on
#' the log scale in both engines, so the default already agrees.
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
    # SAS HAZPRED puts survival CLs on the logit scale; predict.hazard()
    # defaults to "log-log". Hazard CLs are log scale in both engines, so
    # only the survival call needs steering.
    if (identical(type, "survival") && isTRUE(want_cl)) {
      args$conf.type <- "logit"
    }
    as.call(c(quote(predict), args))
  }

  # NOSURV and NOHAZ together suppress both prediction families. The ternary
  # below tests only want_surv, so without this guard `call` would silently
  # fall back to a hazard predict() nobody asked for -- a populated result
  # over a request for nothing, this package's signature defect. Emit NULL
  # for both and record it, rather than guess which one the caller meant.
  if (!want_surv && !want_haz) {
    note("NOSURV NOHAZ", "both survival and hazard predictions suppressed")
  }

  list(
    call = if (want_surv) mk("survival") else if (want_haz) mk("hazard"),
    call_haz = if (want_surv && want_haz) mk("hazard") else NULL,
    inhaz = inhaz, grid = grid, untranslated = untr,
    tokens_seen = seen, tokens_mapped = mapped
  )
}
