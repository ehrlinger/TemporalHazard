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
#' where the event did not occur), then `LCENSOR`, defaulting to the bare
#' `EVENT` variable, which is already `0`/`1` for right-censored/event and so
#' also covers `RCENSOR` with no extra branch (see below).
#'
#' @param statements Named list of SAS statement operands, e.g.
#'   `list(EVENT = "DEAD", TIME = "T", ICENSOR = c("LO", "HI"))`. `ICENSOR`
#'   carries two operands (lower, upper bound variables); `LCENSOR` and
#'   `RCENSOR` carry one.
#' @return `list(status_expr = <call>, time_lower = <name|NULL>,
#'   time_upper = <name|NULL>, untranslated = <data.frame>)`. `status_expr`
#'   is a `bquote()`-built call, evaluable against an environment/list
#'   holding the named SAS variables.
#' @noRd
.hzr_censor_spec <- function(statements) {
  ev <- as.name(statements$EVENT)
  expr <- ev

  # RCENSOR needs no branch: right-censoring is code 0, which is exactly what
  # EVENT already carries when the event did not occur. This is deliberate,
  # not an omission -- there is nothing for an RCENSOR flag to change.

  if (!is.null(statements$LCENSOR)) {
    flag <- as.name(statements$LCENSOR)
    expr <- bquote(ifelse(.(flag) == 1, -1, .(expr)))
  }

  if (!is.null(statements$ICENSOR)) {
    lo <- as.name(statements$ICENSOR[[1L]])
    # Interval censoring only applies where the event did not occur; if it
    # did, EVENT == 1 and the row stays coded as an event, not an interval.
    expr <- bquote(ifelse(.(ev) == 0 & !is.na(.(lo)), 2, .(expr)))
  }

  list(
    status_expr = expr,
    time_lower = if (!is.null(statements$ICENSOR)) as.name(statements$ICENSOR[[1L]]),
    time_upper = if (!is.null(statements$ICENSOR)) as.name(statements$ICENSOR[[2L]]),
    untranslated = .hzr_untranslated_frame()
  )
}
