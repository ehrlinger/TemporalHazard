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
