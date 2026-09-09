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
