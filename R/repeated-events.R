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

  at_end <- !is.na(event_time) & !is.na(data[[followup]]) & event_time == data[[followup]]
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
