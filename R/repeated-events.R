# repeated-events.R -- R reimplementation of the SAS %repeat macro
#
# PURPOSE
# -------
# Repeated-events HAZARD jobs build their fit input with the SAS macro
# %repeat (~/Documents/macro.library/repeat.sas, 161 lines, 2003-10-17).
# That input was never saved, so those jobs cannot be reproduced without
# rebuilding it.  This file rebuilds it in R.
#
# The macro is long-in, long-out: its input is already one row per
# candidate event per subject.  The transform is gap-filling and
# segmentation -- keep the event rows, guarantee at least one row per
# subject, append a terminal censored row at end of follow-up when the
# last event is earlier, then lag the event time within subject to turn
# a set of event *times* into a set of left-truncated *intervals*.
#
# SAS/C BRIDGE
# ------------
# Seven internal stage functions mirror the macro's seven DATA steps
# one-for-one, so the row ladder in the reference log is observable
# between stages.  See inst/dev/REPEATED-EVENTS-DESIGN.md.
#
# Three SAS subtleties are load bearing and are commented at their sites:
# "non-event" and "event" are not complements; stage 7 tests eventype=0
# exactly; and SAS sorts missing numerics first where order() sorts them
# last.

# Output columns the macro creates.  An input column of any of these
# names would be silently overwritten, so validation refuses it.
.hzr_re_output_cols <- c(
  "rcensor", "first", "last", "event", "event_no", "iv_start", "iv_seg", "renewal"
)

# SAS: (&eventype=0 or &eventype=.)
#
# NOT the complement of the event test, which is &eventype=1.  An
# eventype of 2 is neither a non-event nor an event.  Keep the two
# predicates separate.
.hzr_re_nonevent <- function(x) {
  is.na(x) | x == 0
}

# SAS: proc sort data=&out; by &id &iv_event;
#
# Missing numerics sort FIRST in SAS and LAST under order(), so missing
# times are mapped to -Inf for the sort key.  method = "radix" is forced
# because character-vector sort order for other methods depends on the
# locale's collating sequence, which would change which row is a subject's
# first row between machines.
.hzr_re_order <- function(data, id, time) {
  key <- data[[time]]
  key[is.na(key)] <- -Inf
  order(data[[id]], key, method = "radix")
}

# SAS: first.&id / last.&id, for data already sorted by id.
.hzr_re_flags <- function(ids) {
  n <- length(ids)
  if (n == 0L) {
    return(list(first = logical(0), last = logical(0)))
  }
  if (n == 1L) {
    return(list(first = TRUE, last = TRUE))
  }
  boundary <- ids[-1L] != ids[-n]
  list(first = c(TRUE, boundary), last = c(boundary, TRUE))
}

.hzr_re_validate <- function(data, id, time, followup, indicator) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  cols <- c(id = id, time = time, followup = followup, indicator = indicator)
  for (arg in names(cols)) {
    value <- cols[[arg]]
    if (!is.character(value) || length(value) != 1L || is.na(value)) {
      stop(sprintf("`%s` must be a single column name.", arg), call. = FALSE)
    }
    if (!value %in% names(data)) {
      stop(sprintf("Column \"%s\" (argument `%s`) not found in `data`.", value, arg), call. = FALSE)
    }
  }
  clash <- intersect(names(data), .hzr_re_output_cols)
  if (length(clash) > 0L) {
    stop(
      sprintf(
        "`data` already has column(s) %s, which this function creates and would overwrite.",
        paste(sprintf("\"%s\"", clash), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (anyNA(data[[id]])) {
    stop(sprintf("Column \"%s\" (argument `id`) has missing values.", id), call. = FALSE)
  }
  if (is.factor(data[[id]])) {
    stop(
      sprintf(
        paste(
          "Column \"%s\" (argument `id`) is a factor. Convert it to character first: a factor sorts",
          "by level order, not value, and can give a different subject ordering than SAS's `proc sort`."
        ),
        id
      ),
      call. = FALSE
    )
  }
  indicator_col <- data[[indicator]]
  if (!(is.numeric(indicator_col) || is.logical(indicator_col))) {
    stop(
      sprintf(
        "Column \"%s\" (argument `indicator`) must be numeric, integer or logical, not %s.",
        indicator, class(indicator_col)[1]
      ),
      call. = FALSE
    )
  }
  if (!any(!is.na(indicator_col) & indicator_col == 1)) {
    warning(
      sprintf(
        paste(
          "Column \"%s\" (argument `indicator`) has no value equal to 1: no events were found.",
          "This can be a legitimately all-censored cohort, but is often a data problem -- check that",
          "the indicator is coded as documented (1 = event, 0 or NA = no event)."
        ),
        indicator
      ),
      call. = FALSE
    )
  }
  if (anyNA(data[[followup]])) {
    stop(sprintf("Column \"%s\" (argument `followup`) has missing values.", followup), call. = FALSE)
  }
  invisible(NULL)
}

# Stage 1 -- SAS: data &out; set &in; &rcensor=0;
.hzr_re_stage1 <- function(data) {
  data$rcensor <- 0
  data
}

# Stage 2 -- SAS: if first.&id=0 and (&eventype=0 or &eventype=.) then delete;
#
# Keep the group's first row unconditionally, plus every row that is not a
# non-event.  Note the double negative is deliberate: "not a non-event" is
# wider than "is an event", and SAS's predicate is the non-event one.
.hzr_re_stage2 <- function(data, id, time, indicator) {
  data <- data[.hzr_re_order(data, id, time), , drop = FALSE]
  flags <- .hzr_re_flags(data[[id]])
  keep <- flags$first | !.hzr_re_nonevent(data[[indicator]])
  data[keep, , drop = FALSE]
}

# Stage 3 -- SAS:
#   if first.&id=1 and last.&id=1 and nonevent then do; &iv_event=&iv_end; &rcensor=1; end;
#   if first.&id>last.&id and nonevent then delete;
#
# Both statements read the SAME first./last. values, so the flags are
# computed once, before the padding mutation.  The two conditions are
# disjoint (first & last against first & !last), so order does not matter
# between them -- but recomputing the flags after the mutation would be
# wrong in principle and is avoided here on purpose.
.hzr_re_stage3 <- function(data, id, time, followup, indicator) {
  data <- data[.hzr_re_order(data, id, time), , drop = FALSE]
  flags <- .hzr_re_flags(data[[id]])
  nonevent <- .hzr_re_nonevent(data[[indicator]])

  solo <- flags$first & flags$last & nonevent
  data[[time]][solo] <- data[[followup]][solo]
  data$rcensor[solo] <- 1

  drop <- flags$first & !flags$last & nonevent
  data[!drop, , drop = FALSE]
}

# Stage 4 -- SAS:
#   first=first.&id; last=last.&id;
#   output;
#   if last.&id and &rcensor=0 then do;
#     &rcensor=1; &iv_event=&iv_end; &eventype=0; output;
#   end;
#
# `first` and `last` become ordinary data columns here.  Stage 6's renewal
# assignment reads THESE values, not the automatic variables of its own
# step, and stage 5 subsets the data in between.  They must be carried
# forward, never recomputed.  See REPEATED-EVENTS-DESIGN.md.
.hzr_re_stage4 <- function(data, id, time, followup, indicator) {
  data <- data[.hzr_re_order(data, id, time), , drop = FALSE]
  flags <- .hzr_re_flags(data[[id]])
  data$first <- as.numeric(flags$first)
  data$last <- as.numeric(flags$last)

  extra_rows <- which(flags$last & data$rcensor == 0)
  if (length(extra_rows) == 0L) {
    return(data)
  }
  extra <- data[extra_rows, , drop = FALSE]
  extra$rcensor <- 1
  extra[[time]] <- extra[[followup]]
  extra[[indicator]] <- 0

  # The appended row must land immediately after its source row, as the
  # second SAS `output` does.  Interleave rather than rbind-and-resort:
  # at a tie (last event exactly at &iv_end) a re-sort could place the
  # copy first, which would change which row stage 6 lags from.
  combined <- rbind(data, extra)
  place <- c(seq_len(nrow(data)), extra_rows + 0.5)
  combined[order(place, method = "radix"), , drop = FALSE]
}

# Stage 5 -- SAS:
#   if &rcensor=1 or &eventype=1;
#   if &eventype=1 then &event=1; else &event=0;
#
# A subsetting IF, then the event flag.  &eventype=1 is an exact test:
# a missing eventype fails it, and so does any other code.
.hzr_re_stage5 <- function(data, indicator) {
  is_event <- !is.na(data[[indicator]]) & data[[indicator]] == 1
  data <- data[data$rcensor == 1 | is_event, , drop = FALSE]
  data$event <- as.numeric(!is.na(data[[indicator]]) & data[[indicator]] == 1)
  data
}

# Stage 6 -- SAS:
#   retain lag_iv 0 number 0;
#   &iv_start=0; &event_no=0;
#   if first.&id=1 then do; number=&event; &event_no=&event; end;
#   else do; &iv_start=lag_iv; &event_no=number + &event; number=&event_no; end;
#   lag_iv=&iv_event;
#   if &iv_event=&iv_end then &rcensor=1;
#   &iv_seg=&iv_event-&iv_start;
#   &renewal=&event_no;
#   if (first and &event_no=0) or (last and &rcensor=1 and &event=0)
#      then &renewal=&event_no+1;
#
# The loop is row-sequential because SAS's RETAIN is.  lag_iv and number
# are DATA-step loop state and are deliberately NOT returned -- see
# REPEATED-EVENTS-DESIGN.md, "Column-count discrepancy".
.hzr_re_stage6 <- function(data, id, time, followup) {
  data <- data[.hzr_re_order(data, id, time), , drop = FALSE]
  flags <- .hzr_re_flags(data[[id]])
  n <- nrow(data)

  event_time <- data[[time]]
  event <- data$event
  iv_start <- numeric(n)
  event_no <- numeric(n)
  lag_iv <- 0
  number <- 0

  for (i in seq_len(n)) {
    if (flags$first[i]) {
      number <- event[i]
      event_no[i] <- event[i]
    } else {
      iv_start[i] <- lag_iv
      event_no[i] <- number + event[i]
      number <- event_no[i]
    }
    lag_iv <- event_time[i]
  }

  data$iv_start <- iv_start
  data$event_no <- event_no

  # SAS numeric missing is an ordered value that compares equal to itself, so
  # `if . = . then` is TRUE: a row with both &iv_event and &iv_end missing is
  # at end of follow-up. Only the one-side-missing case is false (`. = 9` is
  # false), so that case, and only that case, is excluded below.
  both_missing <- is.na(event_time) & is.na(data[[followup]])
  at_end <- both_missing |
    (!is.na(event_time) & !is.na(data[[followup]]) & event_time == data[[followup]])
  data$rcensor[at_end] <- 1

  data$iv_seg <- event_time - iv_start

  # `first` and `last` here are the stage-4 columns, deliberately stale.
  # `rcensor` is the value just updated on the line above, not stage 5's.
  data$renewal <- event_no
  bump <- (data$first == 1 & event_no == 0) |
    (data$last == 1 & data$rcensor == 1 & data$event == 0)
  data$renewal[bump] <- event_no[bump] + 1

  data
}

# Stage 7 -- SAS:
#   if &iv_seg=0 and &eventype=0 and (first.&id NE 1) then delete;
#
# Note &eventype=0 EXACTLY.  Unlike stages 2 and 3 this is not the
# "0 or missing" non-event test, so a zero-duration row with a missing
# indicator survives here.
#
# iv_seg can legitimately be NA (it is time - iv_start, and a missing event
# time propagates). SAS's `if &iv_seg=0 ...` is FALSE when iv_seg is missing
# (a missing never equals 0), so the row is kept -- guard the comparison so a
# missing iv_seg yields FALSE, not NA, otherwise `drop` is NA for that row and
# `data[!drop, ]` neither drops nor keeps it: it injects an all-NA phantom row.
.hzr_re_stage7 <- function(data, id, time, indicator) {
  data <- data[.hzr_re_order(data, id, time), , drop = FALSE]
  flags <- .hzr_re_flags(data[[id]])
  drop <- !is.na(data$iv_seg) & data$iv_seg == 0 &
    !is.na(data[[indicator]]) & data[[indicator]] == 0 &
    !flags$first
  data[!drop, , drop = FALSE]
}

#' Build repeated-event segments, reproducing the SAS `%repeat` macro
#'
#' Converts a long data set of candidate event times into one row per
#' inter-event segment, ready to fit as a repeated-events model. This is a
#' native R implementation of the SAS macro `%repeat` used by the
#' repeated-events HAZARD jobs, whose fit input was never saved.
#'
#' @details
#' The input holds one row per candidate event per subject, with an indicator
#' naming the event of interest: a value of 1 marks an event, and 0 or `NA`
#' marks its absence. For each subject the function keeps the event rows,
#' guarantees at least one row, appends a censored row at the end of follow-up
#' when the last event happened earlier, and then lags the event time within
#' subject so that each row describes the interval since the previous event.
#'
#' Two details are inherited from the macro and are not R conventions. An
#' indicator code that is neither 0, 1 nor missing is treated as neither an
#' event nor an absence, and rows are ordered with missing times first, as SAS
#' sorts them.
#'
#' @param data A data frame with one row per candidate event per subject.
#' @param id Name of the column identifying the subject.
#' @param time Name of the column giving the interval from time zero to each
#'   event.
#' @param followup Name of the column giving the interval from time zero to the
#'   end of follow-up.
#' @param indicator Name of the event indicator column: 1 marks an event, 0 or
#'   `NA` marks its absence.
#'
#' @return A data frame with one row per inter-event segment. Every input
#'   column is retained, with `time` and `indicator` altered where the macro
#'   alters them, and the following columns added:
#'   \describe{
#'     \item{`event`}{1 when the row is an event of interest, otherwise 0.}
#'     \item{`event_no`}{Running count of event repeats within the subject.}
#'     \item{`rcensor`}{1 when the row is right censored, otherwise 0.}
#'     \item{`iv_start`}{Interval from time zero to the start of the segment.}
#'     \item{`iv_seg`}{Duration of the segment, `time` minus `iv_start`.}
#'     \item{`renewal`}{Segment number under the modulated renewal
#'       formulation.}
#'     \item{`first`, `last`}{1 when the row was the subject's first or last
#'       before censored rows were appended, otherwise 0.}
#'   }
#'
#' @examples
#' # One row per candidate event per subject: s1 has two events, s2 one,
#' # and s3 none (its only row is a non-event, so it contributes no rows
#' # of its own before the function pads it with a censored one below).
#' events <- data.frame(
#'   id = c("s1", "s1", "s2", "s3"),
#'   t = c(1, 3, 2, NA),
#'   fu = c(10, 10, 10, 10),
#'   ev = c(1, 1, 1, 0)
#' )
#' # The result appends a censored row at end of follow-up for any subject
#' # whose last event happened before `fu`, so every subject gets at least
#' # one row even when, like s3, they never had an event.
#' hzr_repeated_events(events, id = "id", time = "t", followup = "fu", indicator = "ev")
#'
#' @export
hzr_repeated_events <- function(data, id, time, followup, indicator) {
  .hzr_re_validate(data, id, time, followup, indicator)

  data <- .hzr_re_stage1(data)
  data <- .hzr_re_stage2(data, id, time, indicator)
  data <- .hzr_re_stage3(data, id, time, followup, indicator)
  data <- .hzr_re_stage4(data, id, time, followup, indicator)
  data <- .hzr_re_stage5(data, indicator)
  data <- .hzr_re_stage6(data, id, time, followup)
  data <- .hzr_re_stage7(data, id, time, indicator)

  row.names(data) <- NULL
  data
}
