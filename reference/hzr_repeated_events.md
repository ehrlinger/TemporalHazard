# Build repeated-event segments, reproducing the SAS `%repeat` macro

Converts a long data set of candidate event times into one row per
inter-event segment, ready to fit as a repeated-events model. This is a
native R implementation of the SAS macro `%repeat` used by the
repeated-events HAZARD jobs, whose fit input was seldom saved.

## Usage

``` r
hzr_repeated_events(data, id, time, followup, indicator)
```

## Arguments

- data:

  A data frame with one row per candidate event per subject.

- id:

  Name of the column identifying the subject.

- time:

  Name of the column giving the interval from time zero to each event.

- followup:

  Name of the column giving the interval from time zero to the end of
  follow-up.

- indicator:

  Name of the event indicator column: 1 marks an event, 0 or `NA` marks
  its absence.

## Value

A data frame with one row per inter-event segment. Every input column is
retained, with `time` and `indicator` altered where the macro alters
them, and the following columns added:

- `event`:

  1 when the row is an event of interest, otherwise 0.

- `event_no`:

  Running count of event repeats within the subject.

- `rcensor`:

  1 when the row is right censored at the end of follow-up, otherwise 0.
  This is set both for a genuine censoring row and for an event row
  whose event time equals the end of follow-up. `rcensor` and `event`
  are therefore not mutually exclusive; see the Note below.

- `iv_start`:

  Interval from time zero to the start of the segment.

- `iv_seg`:

  Duration of the segment, `time` minus `iv_start`.

- `renewal`:

  Segment number under the modulated renewal formulation.

- `first`, `last`:

  1 when the row was the subject's first or last row at stage 4,
  otherwise 0. The row stage 4 appends to pad a subject's censored tail
  inherits its source row's flags, so a subject can have two rows with
  `last == 1` together, or two with `first == 1` together.

## Details

The input holds one row per candidate event per subject, with an
indicator naming the event of interest: a value of 1 marks an event, and
0 or `NA` marks its absence. For each subject the function keeps the
event rows, guarantees at least one row, appends a censored row at the
end of follow-up when the last event happened earlier, and then lags the
event time within subject so that each row describes the interval since
the previous event.

Two details are inherited from the macro and are not R conventions. An
indicator code that is neither 0, 1 nor missing is treated as neither an
event nor an absence, and rows are ordered with missing times first, as
SAS sorts them.

The result can contain segments of zero length, where `iv_start` is not
less than `time`: for example, an event that shares its time with the
subject's previous row. They are kept, as the macro keeps them. A model
that treats `iv_start` as the time a subject enters the risk set needs
each segment to end after it starts, so adjust or remove these rows
before such a fit. The SAS jobs this function reproduces moved `time`
forward by about one hour on each.

`data` may not already contain a column named `rcensor`, `first`,
`last`, `event`, `event_no`, `iv_start`, `iv_seg` or `renewal`: this
function creates all eight and would silently overwrite an existing one
of the same name.

## Note

`rcensor` is not a plain censoring indicator: a row can have
`event == 1` and `rcensor == 1` together, when a subject's last event
happens to fall exactly at the end of follow-up. This is faithful to the
SAS macro, whose own header warns that this combination is not
compatible with a `HAZARD`-style fit. Do not pass `rcensor` directly as
a censoring indicator without first accounting for the overlap with
`event`.

## Examples

``` r
# One row per candidate event per subject: s1 has two events, s2 one,
# and s3 none (its only row is a non-event, so it contributes no rows
# of its own before the function pads it with a censored one below).
events <- data.frame(
  id = c("s1", "s1", "s2", "s3"),
  t = c(1, 3, 2, NA),
  fu = c(10, 10, 10, 10),
  ev = c(1, 1, 1, 0)
)
# The result appends a censored row at end of follow-up for any subject
# whose last event happened before `fu`, so every subject gets at least
# one row even when, like s3, they never had an event.
hzr_repeated_events(events, id = "id", time = "t", followup = "fu", indicator = "ev")
#>   id  t fu ev rcensor first last event iv_start event_no iv_seg renewal
#> 1 s1  1 10  1       0     1    0     1        0        1      1       1
#> 2 s1  3 10  1       0     0    1     1        1        2      2       2
#> 3 s1 10 10  0       1     0    1     0        3        2      7       3
#> 4 s2  2 10  1       0     1    1     1        0        1      2       1
#> 5 s2 10 10  0       1     1    1     0        2        1      8       2
#> 6 s3 10 10  0       1     1    1     0        0        0     10       1
```
