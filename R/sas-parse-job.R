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
#' Builds an unevaluated `status` expression (and, where needed, a
#' `time_lower` expression) from the `EVENT`/`ICENSOR`/`LCENSOR`/`RCENSOR`
#' operands of a `PROC HAZARD` statements list, ready to drop into a
#' `hazard()` call.
#'
#' **`LCENSOR` is left-*truncation*, not left-censoring.** Its operand is the
#' counting-process *entry time*: all four corpus uses are literally
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
#' interval's *lower* bound; the interval runs `CTIME` -> `TIME` (see
#' `FIXTURE-GAP-LIST.md` Q1). `EVENT` is optional: the reference `HAZARD`
#' program (`src/hazard/varterm.c`) terminates only when *both* `EVENT` and
#' `ICENSOR` are missing, so a job that specifies `ICENSOR` alone is
#' legitimate. In that case the base status is right-censored (`0`)
#' everywhere, overridden to `2` where `C3 > 0`. A job specifying neither is
#' rejected with an error, mirroring HAZARD's own termination.
#'
#' **`EVENT`, `ICENSOR` and `WEIGHT` are all counts**, and `setlik.c` combines
#' them for one OBS record as `c1c2c3 = c1w + c2 + c3w` (with `c1w = C1 * WT`
#' and `c3w = C3 * WT`), then
#' `llike = -(c1c2c3) * (CH(T) - CH(ST)) + c1w * log h(T) + c3w * lct`. That
#' decomposes into at most three independent contributions, each of which
#' [hazard()] expresses as one row's status and weight:
#' * `C1 > 0`: an event row of weight `C1 * WT` (status `1`);
#' * `C2 > 0`: a right-censored row of weight `C2` (status `0`);
#' * `C3 > 0`: an interval row of weight `C3 * WT` (status `2`).
#'
#' `C2` is the one **not** multiplied by `WT`, and that is not an oversight in
#' the reference: `readc2.c` sets `C2 = ONE` on exactly the rows where neither
#' `C1` nor `C3` fires, so SAS weights a right-censored row `1` whatever the
#' `WEIGHT` variable says. Emitting a bare `WT` there instead over- or
#' under-weights every censored row, and a `WEIGHT` that happens to be `0`
#' on censored rows deletes them from the fit outright, silently (#158).
#'
#' So `weights_expr` is emitted whenever `EVENT` or `ICENSOR` is present
#' (which is always, since a job with neither is rejected) as
#' `ifelse(EVENT > 0, EVENT * WT, ifelse(C3 > 0, C3 * WT, 1))`, dropping
#' whichever branches the job does not have and dropping `* WT` when there is
#' no `WEIGHT` statement. For the 0/1 `EVENT` and absent `WEIGHT` that most
#' jobs carry this evaluates to exactly the unit weights those jobs already
#' fitted with; it is the repeat-event and weighted rows that change.
#'
#' `status` therefore derives from `EVENT > 0`, never from `EVENT` itself
#' (#157). This package codes status `-1` left, `0` right, `1` event, `2`
#' interval, so a row recording two events (`EVENT = 2`) mapped straight onto
#' `status` becomes *interval-censored* (a different likelihood branch, not
#' an under-count), and `EVENT = 3` lands outside the coding altogether,
#' which also disables Conservation of Events (`coe_supported_data`).
#' `readc1.c` accepts any `C1 >= 0` and keeps a fractional count (`notintg`
#' only raises a report flag), so `> 0` is the right test and a fractional
#' count carries through as a fractional weight. A *negative* count SAS
#' deletes (`c1del`, `del = 1`) and a *missing* one likewise (`mc1del`,
#' `mdel = 1`); `readobs.c` then skips `setobs()` and subtracts the row from
#' `Nobs`, so it contributes nothing at all; it is emphatically not a
#' right-censored row of weight 1. This translator has no way to drop a row,
#' and dropping one silently would change `n` behind the reader's back, so the
#' hoisted status chunk opens with a guard per named count that stops and
#' names the variable, the SAS rule and the remedy (filter the rows out).
#'
#' **A row where `EVENT` and `ICENSOR` both fire is refused at fit time.**
#' `setlik.c` *sums* the two contributions, so such a row is simultaneously an
#' event of weight `C1 * WT` and an interval observation of weight `C3 * WT`;
#' [hazard()] carries one status and one weight per row and cannot express
#' both. Whether any row does this is a property of the *data*, not of the job
#' text, so it cannot be refused at translation time the way `LCENSOR` +
#' `ICENSOR` is; instead the hoisted status chunk opens with a guard that
#' stops before the fit and names the by-hand remedy (split the row in two).
#' Picking the event branch and discarding `c3w` (what this translator did
#' before) converges and reports plausibly, which is the shape this package
#' exists to refuse. The same guard covers `EVENT` + `RCENSOR` and `ICENSOR`
#' + `RCENSOR`; see the `RCENSOR` paragraph below.
#'
#' **`RCENSOR` names `C2` itself** (`hazard_y.y`:
#' `rcensorstmt : RCENSOR NAME { setvar(13,$2); }`; `rcnsprc.c` sets `c2name`
#' from statement field 13), and `C2` is likewise a *count*:
#' "COUNT OF CENSORED INDIVIDUALS AT TIME=T" in `setlik.c`'s header, and
#' `OBS(3)`, "NUMBER OF CENSORED OBSERVATIONS AT T". The two branches of
#' `readc2.c` are what decides the translation: when `c2name` is non-blank
#' `C2` is read straight from the data, and the `C2 = ONE` derivation is
#' skipped *entirely*; only a blank `c2name` derives it. So the censored
#' branch of `weights_expr` is the literal `1` for a job with no `RCENSOR`
#' and the `C2` variable itself for a job with one (#162), and it is never
#' multiplied by `WT` either way. Four censored individuals fitted as one
#' observation is the same under-weighting as #154/#157, arriving by the
#' third of the three counts.
#'
#' No corpus job exercises it: `hz.te123.OMC.sas` and `hz.tm123.OMC.sas` are
#' the only two of the 110 with an `RCENSOR` statement, and both set a 0/1
#' `CENSORED` that is mutually exclusive with `EVENT`, for which `C2` and the
#' derived `1` agree on every row. That is why it went unseen, not evidence
#' that it does not bite.
#'
#' Three edge cases follow from reading `C2` verbatim, and none is silent.
#' `readobs.c`'s all-zero deletion is narrower than it first looks: both of
#' its rules are guarded by a blank-name test, `ic10 && ic20` only when
#' `c3name` is blank and `ic30 && ic20` only when `c1name` is blank. So a row
#' with every count zero is deleted for `EVENT` + `RCENSOR` and for `ICENSOR`
#' + `RCENSOR`, and *kept* (contributing `c1c2c3 = 0`) when all three are
#' named. There is no such rule at all for `EVENT` + `ICENSOR`, which is
#' exactly the pairing with no `c2name`. This translator has no way to drop a
#' row, so an all-zero row arrives with weight `0`: no contribution to the
#' likelihood, which is what deletion means for the fit and what the
#' all-three case does anyway, though the row still counts toward `n`.
#' A *negative* `C2` SAS deletes (`c2del`) and a *missing* one it deletes too
#' (`mc2del`), the same rule `readc1.c` and `readc3.c` apply to their own
#' counts (guarded only by the name being non-blank), so both are caught by
#' the missing/negative guard above rather than by [hazard()]'s incidental
#' "non-negative and finite" check.
#'
#' **A row where two named counts both fire is refused at fit time.** With
#' `RCENSOR` present this is no longer only the `EVENT` + `ICENSOR` pair:
#' `readobs.c` deletes a row only when *both* of a pair are zero, so
#' `C1 > 0 & C2 > 0` and `C3 > 0 & C2 > 0` reach `setlik.c` too and are
#' summed there. Each such row is two observations at once (an event of
#' weight `C1 * WT` *and* a right-censored observation of weight `C2`, say),
#' and [hazard()] carries one status and one weight per row. A guard is
#' emitted for every pair of counts the job named, and only for named ones: a
#' *derived* `C2` fires on exactly the rows where no other count does, so it
#' cannot collide. Any job carrying a guard has its status hoisted into its
#' own chunk, so the guard runs, and is seen to run, before the fit.
#'
#' `hazard()`'s `time_lower` carries a different meaning per row, selected by
#' `status` (see `R/hazard_api.R`): for status 2 (interval) it is the
#' interval's lower bound (`CTIME`), but for status 0/1 it is the
#' counting-process *entry time* (default `0`), not a censoring bound. `TIME`
#' is always the interval's upper bound, and that is exactly what
#' `hazard()`'s `time_upper` already defaults to when left `NULL`, so
#' `time_upper` is never emitted by this translator; passing it would be
#' redundant, and omitting it cannot drift out of sync with `TIME`.
#' `time_lower` is therefore built as one of, depending on which statements
#' are present:
#' * `ICENSOR` only: `ifelse(status == 2, CTIME, 0)` (`0` is `hazard()`'s
#'   documented default entry time for status 0/1 rows).
#' * `LCENSOR` only: the `LCENSOR` variable directly; it applies to every
#'   row, not just interval ones, so no `ifelse()` gating is needed.
#' * Neither: `time_lower` is omitted (`NULL`).
#'
#' **`LCENSOR` and `ICENSOR` together are refused, not translated.** There is
#' no expression that serves: `time_lower` is one column carrying two
#' meanings, so gating it on status hands the interval rows their lower bound
#' and thereby drops their entry time, fitting a left-truncated
#' interval-censored subject as at risk from time `0`. That fit converges
#' and looks plausible, which is the failure mode this package exists to
#' refuse.
#' The reference implementation carries three distinct times (`TIME`, `CTIME`,
#' `STIME`) and subtracts `H(STIME)` for every row, interval rows included, so
#' translating it faithfully needs an entry-time argument [hazard()] does not
#' have. It is also rare: 0 of the 93 `PROC HAZARD` jobs across the production
#' studies use the combination. So the pair is recorded in `untranslated` and
#' `refused` is set, and the caller emits a `stop()` in place of the fit.
#'
#' The `time_lower` `ifelse()` is gated on a `status_name` placeholder
#' (`.hzr_status`) the caller assigns the `status_expr` to first. It is never
#' unconditional, which would trip `hazard()`'s finite-value check off the
#' interval subset and, where populated, silently redefine event/right-censored
#' rows' risk-set entry times. `weights_expr` gates on the count variables
#' themselves, not on `status`, so it stays valid on the un-hoisted paths.
#'
#' @param statements Named list of SAS statement operands, e.g.
#'   `list(EVENT = "DEAD", TIME = "T", ICENSOR = c("C3", "CTIME"))`.
#'   `ICENSOR` carries two operands (event-count variable, interval
#'   lower-bound time variable); `LCENSOR` and `RCENSOR` carry one. `EVENT`
#'   may be absent if `ICENSOR` is present.
#' @return `list(status_expr = <call>, status_name = <name|NULL>,
#'   time_lower = <call|NULL>, weights_expr = <call|NULL>,
#'   keep_expr = <call|NULL>, degenerate_expr = <call|NULL>,
#'   untranslated = <data.frame>, refused = <logical>)`. `keep_expr` and
#'   `degenerate_expr` are set only for an `ICENSOR` job: the rows
#'   `readct.c` keeps, and those it turns into exact events (#543).
#'   `status_expr` is
#'   a `bquote()`-built call, evaluable against an environment/list holding
#'   the named SAS variables. `time_lower` and `weights_expr`, when
#'   non-`NULL`, are `bquote()`-built calls; `status_name` is non-`NULL` on
#'   every non-refused path, because every such job names at least one count
#'   and every named count carries a missing/negative guard that has to run in
#'   its own chunk ahead of the fit.
#'   `weights_expr` is non-`NULL` on every non-refused path, because
#'   `EVENT` and `ICENSOR` are both counts and one of the two is always
#'   present. `refused` is `TRUE` only for the `LCENSOR` +
#'   `ICENSOR` combination described above, and every other element is `NULL`
#'   when it is: there is nothing to emit. Both returns carry the same element
#'   names, so a caller reading one never meets an unexpected missing field.
#' @noRd
.hzr_censor_spec <- function(statements) {
  has_event <- !is.null(statements$EVENT)
  has_icensor <- !is.null(statements$ICENSOR)
  has_lcensor <- !is.null(statements$LCENSOR)
  has_rcensor <- !is.null(statements$RCENSOR)

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
      keep_expr = NULL,
      degenerate_expr = NULL,
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
  c3 <- if (has_icensor) as.name(statements$ICENSOR[[1L]])
  c2 <- if (has_rcensor) as.name(statements$RCENSOR)
  wt <- if (!is.null(statements$WEIGHT)) as.name(statements$WEIGHT)

  # LCENSOR is left-truncation, not left-censoring (see roxygen): it never
  # touches status, only time_lower below.

  # EVENT and C3 are counts, not flags (see roxygen). status comes from
  # `> 0`, never from the count itself: EVENT = 2 read as a status is
  # interval-censored, a different likelihood branch (#157). RCENSOR adds no
  # branch of its own -- C2 > 0 is right-censoring, code 0, which is already
  # this expression's fallback. It changes the row's WEIGHT, not its status.
  # readct.c resolves a degenerate interval before the fit, on a row with
  # C3 > 0 (#543): CTIME == TIME makes it an exact event, C1 = C1 + C3 and
  # C3 = CT = 0 (readct.c:18-23), and a missing, negative or greater-than-TIME
  # CTIME deletes the row (readct.c:9-17). The status below mirrors the
  # first; `keep_expr` mirrors the second, and the fit reads only the rows it
  # keeps. hazard() itself refuses both shapes under objective = "sas".
  interval_code <- 2
  keep_expr <- NULL
  if (has_icensor) {
    ctime <- as.name(statements$ICENSOR[[2L]])
    tt <- as.name(statements$TIME)
    # NA-safe on TIME: a missing TIME leaves the row as it was (interval,
    # kept), for hazard() to reject with its own message, as before #543.
    interval_code <- bquote(ifelse(!is.na(.(ctime)) & !is.na(.(tt)) &
                                     .(ctime) == .(tt), 1, 2))
    keep_expr <- bquote(!(.(c3) > 0 & (is.na(.(ctime)) | .(ctime) < 0 |
                                         (!is.na(.(tt)) & .(ctime) > .(tt)))))
  }
  expr <- if (has_event && has_icensor) {
    bquote(ifelse(.(ev) > 0, 1, ifelse(.(c3) > 0, .(interval_code), 0)))
  } else if (has_event) {
    bquote(ifelse(.(ev) > 0, 1, 0))
  } else {
    bquote(ifelse(.(c3) > 0, .(interval_code), 0))
  }

  # The weight is the row's own count times the WEIGHT variable on event and
  # interval rows (c1w = C1 * WT, c3w = C3 * WT). C2 is the term that enters
  # the sum UNWEIGHTED (#158), and it is a count in its own right whenever
  # RCENSOR named it: readc2.c reads c2name straight from the data and only
  # DERIVES C2 = 1 -- on precisely the rows where neither other count fires
  # -- when that name is blank (#162). So the censored branch is the literal
  # 1 for a job with no RCENSOR and the C2 variable for a job with one.
  ev_w <- if (is.null(wt)) ev else bquote(.(ev) * .(wt))
  c3_w <- if (is.null(wt)) c3 else bquote(.(c3) * .(wt))
  c2_w <- if (has_rcensor) c2 else 1
  weights_expr <- if (has_event && has_icensor) {
    bquote(ifelse(.(ev) > 0, .(ev_w), ifelse(.(c3) > 0, .(c3_w), .(c2_w))))
  } else if (has_event) {
    bquote(ifelse(.(ev) > 0, .(ev_w), .(c2_w)))
  } else {
    bquote(ifelse(.(c3) > 0, .(c3_w), .(c2_w)))
  }

  # Every count the job NAMED, with the status code and the weight setlik.c
  # gives its contribution. A derived C2 is deliberately absent from this
  # table: readc2.c sets it on exactly the rows where no other count fires,
  # so it cannot collide with one. A named C2 can, and does (#162).
  counts <- list()
  if (has_event) {
    counts$EVENT <- list(v = ev, w = ev_w, code = 1L, readc = 1L,
                         what = "an event")
  }
  if (has_rcensor) {
    counts$RCENSOR <- list(v = c2, w = c2, code = 0L, readc = 2L,
                           what = "a right-censored observation")
  }
  if (has_icensor) {
    counts$ICENSOR <- list(v = c3, w = c3_w, code = 2L, readc = 3L,
                           what = "an interval observation")
  }

  # setlik.c SUMS the contributions (c1c2c3 = c1w + c2 + c3w), so a row with
  # two of these counts non-zero is two observations at once, and readobs.c
  # deletes a row only when BOTH of a pair are zero -- so such a row does
  # reach the fit. hazard() carries one status and one weight per row and
  # cannot say that. Whether any row does it is a property of the DATA, not
  # of the job text, so the refusal happens at fit time rather than at
  # translation time, the way the LCENSOR + ICENSOR one does.
  nms <- names(counts)
  for (i in seq_along(nms)) {
    for (j in seq_len(i - 1L)) {
      a <- counts[[nms[[j]]]]
      b <- counts[[nms[[i]]]]
      expr <- bquote({
        if (any(.(a$v) > 0 & .(b$v) > 0, na.rm = TRUE)) {
          stop(.(sprintf(paste(
            "This job has rows where the %s count (%s) and the %s count",
            "(%s) are both non-zero. SAS sums the two contributions",
            "(setlik.c: c1c2c3 = c1w + c2 + c3w), so such a row is at once",
            "%s of weight %s and %s of weight %s. hazard() carries one",
            "status and one weight per row and cannot express both. Split",
            "each such row into two -- a status %d row and a status %d row,",
            "same times, those two weights -- and fit by hand (%s)."),
            nms[[j]], deparse(a$v), nms[[i]], deparse(b$v),
            a$what, deparse(a$w), b$what, deparse(b$w), a$code, b$code,
            if ("RCENSOR" %in% c(nms[[j]], nms[[i]])) "#162" else "#157")),
            call. = FALSE)
        }
        .(expr)
      })
    }
  }

  # readc1.c, readc2.c and readc3.c apply one rule to every count the job
  # NAMED, and it is the same rule in all three: a missing value sets mdel, a
  # negative one sets del, and readobs.c then skips setobs() for that row and
  # subtracts it from Nobs. A deleted row is NOT a right-censored observation
  # of weight 1 -- it contributes nothing at all, and SAS reports the two
  # tallies (mcNdel, cNdel) in its listing. `> 0` alone cannot say that: it
  # folds a negative count into the censored fallback, and it propagates NA
  # out of ifelse() into both status and weights, where hazard() rejects it
  # with a message about weights that names neither the variable nor the
  # cause. Deleting the rows here instead would change n silently, so stop
  # and say which variable, what SAS does, and what to do about it.
  # Reversed so the guards read EVENT, RCENSOR, ICENSOR top to bottom in the
  # rendered document: each wrap goes outside the previous one.
  for (nm in rev(names(counts))) {
    a <- counts[[nm]]
    expr <- bquote({
      if (any(is.na(.(a$v)) | .(a$v) < 0)) {
        stop(.(sprintf(paste(
          "This job has rows where the %s count (%s) is missing or negative.",
          "SAS deletes such rows outright: readc%d.c sets mdel on a missing",
          "count and del on a negative one, and readobs.c then skips",
          "setobs() and subtracts the row from Nobs, so it contributes",
          "nothing to the likelihood -- it is not %s of weight 1. hazard()",
          "has no equivalent, and translating the row as censored would",
          "change the fit without saying so. Drop those rows before fitting",
          "-- subset to !is.na(%s) & %s >= 0 -- and compare your row count",
          "against the deletion tallies SAS prints for this job."),
          nm, deparse(a$v), a$readc, a$what,
          deparse(a$v), deparse(a$v))),
          call. = FALSE)
      }
      .(expr)
    })
  }

  # Every job carrying a guard gets its status hoisted into a chunk of its
  # own, not only the ICENSOR jobs whose time_lower has to gate on it: the
  # guard must run and be seen to run BEFORE the fit, and one buried inside
  # hazard(status = ...) is neither readable in the rendered document nor
  # separable from the fit that follows it. Since every non-refused job names
  # at least one count, and every named count carries the missing/negative
  # guard above, that is now every job. It also settles the evaluation order:
  # a missing count reaches the guard before it can reach hazard()'s
  # `weights` check, so the explanatory message wins over the incidental one.
  status_name <- as.name(".hzr_status")
  time_lower <- NULL
  if (has_icensor) {
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
    keep_expr = keep_expr,
    degenerate_expr = if (has_icensor) {
      bquote(.(c3) > 0 & !is.na(.(ctime)) & !is.na(.(tt)) & .(ctime) == .(tt))
    },
    untranslated = .hzr_untranslated_frame(),
    refused = FALSE
  )
}

# The U1 class for a job PROC HAZARD runs on a model this translation does
# not emit (#358). One lead, so every warning of the class reads the same.
.hzr_sas_not_mirrored_lead <- paste0(
  "This translation cannot emit PROC HAZARD's model for this job, so the ",
  "fit below stands in for a model PROC HAZARD does not fit: ")

#' Does a later `(` clear PROC HAZARD's syntax-error flag for this job?
#'
#' `hazard_l.l:56` is `\(  { BEGIN HZRP; yysynerr = 0; yylnctr = 1; }`. It
#' has no start condition, so it fires on every `(` the lexer reads, and it
#' clears the latch that `initprz.c:75-77` tests: a syntax error raised
#' before the job's last `(` no longer stops the job (#461). After that `(`
#' the lexer stays in its PROC-line state up to the `;`, where anything but
#' whitespace, `)` (`hazard_l.l:32`) and `= NUMBER` (`:53`, `:55`) is
#' unexpected text (`:176-179`) or a parse error, and sets the latch again.
#' A `(` in a `* ... ;` comment never reaches that rule (`:51`), and this
#' translation has stripped comments before it gets here.
#'
#' Judged by statement against the HAZARD binary on the package's `avc`
#' (`tests/testthat/fixtures/paren-reset-oracle.csv`):
#' * `"cleared"`: every statement carrying a syntax refusal comes before the
#'   last statement with a `(`, whose text after its last `(` is clean.
#'   PROC HAZARD does not refuse the job for them.
#' * `"same"`: as `"cleared"`, except that a refusal sits in that statement
#'   itself. PROC HAZARD's parser recovery then decides, and the binary
#'   fits `EARLY AGE*SEX, LOG();` but stops `EARLY 1AGE, SEX, LOG();` and
#'   `EARLY AGE=ABC, LOG();` before fitting.
#' * `"none"`: no `(` clears the refusals, so the refusal stands.
#'
#' The PROC line (statement 1) is never taken as the clearing statement: the
#' parser meets a PROC-line error at a token after it, often the `;`, so a
#' `(` in the same line clears too early (`MAXITER=5. ()` is refused). A
#' statement this translation cannot judge (an unknown keyword, a macro
#' reference) at or after the last `(` makes the verdict `"none"`: both
#' leave a refusal in place rather than clear one PROC HAZARD keeps.
#' @param st The job's statements, split at `;`, the PROC line first.
#' @param err Indices of statements carrying a syntax refusal.
#' @param blind Indices of statements this translation cannot judge.
#' @return `"cleared"`, `"same"` or `"none"`.
#' @noRd
.hzr_sas_paren_reset <- function(st, err, blind) {
  if (!length(err)) return("none")
  has <- grepl("(", st, fixed = TRUE)
  has[[1L]] <- FALSE
  if (!any(has)) return("none")
  last <- max(which(has))
  s <- st[[last]]
  tail <- substring(s, max(gregexpr("(", s, fixed = TRUE)[[1L]]) + 1L)
  tail <- trimws(gsub(")", " ", tail, fixed = TRUE))
  clean <- !nzchar(tail) ||
    (startsWith(tail, "=") &&
       .hzr_sas_lexer_number(trimws(substring(tail, 2L))))
  if (!clean || any(err > last) || any(blind >= last)) return("none")
  if (any(err == last)) "same" else "cleared"
}

# The verdict for a job whose syntax refusals a later `(` cleared (#461).
# It says the job is not refused FOR THEM, and no more: a job cleared here
# can still be refused later, at fit time, by SETG1 or SETG3 (measured:
# `... FIXTAU TAU=0 ... FIXG1; LATE LOG();` exits SEMANTIC), and that
# refusal carries its own warning.
.hzr_sas_paren_cleared <- function(what) {
  paste0(
    "PROC HAZARD does not stop this job for the syntax error in ",
    paste(what, collapse = "; "), ", because a `(` in a later statement ",
    "clears its syntax-error flag (hazard_l.l:56) before initprz.c:75-77 ",
    "tests it. Past its parse, the job is what its parser's error recovery ",
    "leaves, which this translation does not reproduce: ")
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
  toks <- .hzr_sas_join_spaced(toks[nzchar(toks)])
  # A bare `=` survives the joiner only when the grammar has nothing to pair
  # it with: `DATA = MAXITER = 50` is DATA='MAXITER' (a valid `DATA '=' NAME`,
  # because <HZRP>DATA switches the lexer to DSNM where MAXITER lexes as a
  # NAME, hazard_l.l:59,80) followed by a STRAY `=` and an orphan value. SAS
  # reaches `hazardopt : error` (hazard_y.y:76) there and does not run the
  # job. Recorded once, as the syntax error it is, rather than as one blank
  # "unknown option" row per leftover token (#433 review 2).
  stray <- which(toks == "=")
  if (length(stray)) {
    drop <- unique(c(stray, stray[stray < length(toks)] + 1L))
    leftover <- paste(toks[drop], collapse = " ")
    toks <- toks[-drop]
    proc_syntax_error <- paste0(
      "a stray `=` on the PROC HAZARD line (", leftover, "): the option ",
      "before it already took its value, so PROC HAZARD reaches ",
      "`hazardopt : error` (hazard_y.y:76) and rejects this job with a ",
      "syntax error")
    proc_syntax_what <- paste0("the stray `=` in `", leftover, "`")
  } else {
    proc_syntax_error <- NULL
    proc_syntax_what <- NULL
  }
  ctl <- list()
  data_name <- NULL
  # DATA= written, even with no usable name: the %HAZARD macro then finds it,
  # and the refusal below names PROC HAZARD's syntax error instead (#497).
  data_given <- FALSE
  data_raw <- NULL
  outhaz <- NULL
  # A PROC-line value the lexer does not read as a NUMBER (hazard_l.l:33-38,
  # the HZRP state at :53) is a syntax error: PROC HAZARD does not run the
  # job (U1, #403). as.numeric() reads 1E5 and 5., which the lexer does not.
  proc_rejected <- character(0)
  # The same refusals by construct alone, for a verdict that names them
  # without their reasons (#461).
  proc_what <- character(0)
  # A TIME or EVENT with no operand (#431), by token; fatal below only when
  # nothing else in the job supplies the variable.
  stmt_empty <- list()
  # Returns TRUE when it rejected the option, so the caller can stop rather
  # than add a second, sometimes contradictory row. `CONDITION=5.` used to say
  # both that PROC HAZARD's lexer rejects the number AND what its optimizer
  # does with the value, although a rejected job never runs (#433 review).
  # `DATA '=' dsfield` and `OUTHAZ '=' dsfield` (hazard_y.y:61-62), where
  # dsfield : NAME | LIBMEM (:80-81), have no form without a name, so an
  # empty value is the grammar refusing the option -- the same shape as
  # `MAXITER '=' NUMBER` with no number, and it belongs in the same presence
  # check. Before this, OUTHAZ= was dropped with no row and the job fitted,
  # and DATA= left data_name as "" and surfaced as an internal
  # "attempt to use zero-length variable name" (#433 review 2).
  check_name <- function(key, val) {
    if (.hzr_sas_is_macro(val)) return(FALSE)
    if (nzchar(val)) return(FALSE)
    proc_rejected <<- c(proc_rejected, paste0(
      key, ": no value, and PROC HAZARD has no form of this option without a ",
      "dataset name (hazard_y.y:61-62, :80-81), so it rejects this job with ",
      "a syntax error"))
    proc_what <<- c(proc_what, key)
    note(key, paste0("no value; PROC HAZARD has no form of this option ",
                     "without a dataset name (hazard_y.y:61-62, :80-81)"))
    TRUE
  }
  check_number <- function(key, val) {
    # A macro carries no verdict: SAS expands it before PROC HAZARD reads
    # the statement, so whether a NUMBER arrives is not knowable here.
    if (.hzr_sas_is_macro(val)) return(FALSE)
    if (!nzchar(val)) {
      # `MAXITER '=' NUMBER` and `CONDITION '=' NUMBER` (hazard_y.y:63-64)
      # have no form without a NUMBER, so `MAXITER=`, `MAXITER =` and a bare
      # `MAXITER` all fall to `hazardopt : error` (:76). This is the grammar
      # refusing the option, not the lexer refusing a value, hence the
      # different citation.
      proc_rejected <<- c(proc_rejected, paste0(
        key, ": no value, and PROC HAZARD has no form of this option without ",
        "one (hazard_y.y:63-64), so it rejects this job with a syntax error"))
      proc_what <<- c(proc_what, key)
      # A refused construct is listed as well as warned about: the warning
      # is read once at render, the row is what a reader greps afterwards.
      note(key, paste0("no value; PROC HAZARD has no form of this option ",
                       "without one (hazard_y.y:63-64)"))
      return(TRUE)
    }
    if (!.hzr_sas_lexer_number(val)) {
      proc_rejected <<- c(proc_rejected, paste0(
        key, "=", val, ": not a number PROC HAZARD's lexer reads ",
        "(hazard_l.l:33-38), so PROC HAZARD rejects this job with a syntax ",
        "error"))
      proc_what <<- c(proc_what, paste0(key, "=", val))
      note(paste0(key, "=", val),
           paste0("not a number PROC HAZARD's lexer reads ",
                  "(hazard_l.l:33-38)"))
      return(TRUE)
    }
    FALSE
  }

  for (tok in toks) {
    eqp <- .idx(tok, "=")
    key <- if (eqp > 0L) substring(tok, 1L, eqp - 1L) else tok
    val <- if (eqp > 0L) substring(tok, eqp + 1L) else ""
    token <- .hzr_sas_token(key, "HAZARD", "HZRP")
    if (identical(token, "PROC") || identical(token, "HAZARD")) next
    seen <- seen + 1L
    if (is.na(token)) {
      # An unknown word here IS a lexer catch-all in PROC HAZARD
      # (hazard_l.l:177-179), but this parser's block text is not guaranteed
      # to hold only PROC HAZARD statements: a %repeat call brings a DATA
      # step's own keywords through here. Claiming a syntax error on them
      # refused jobs that run, so it is recorded, not refused (see the
      # leftovers issue).
      note(key, "unknown PROC HAZARD option")
      next
    }
    # Eleven PROC HAZARD options are bare tokens with no `'=' value` form
    # (hazard_y.y:65-75). A value after one reaches `hazardopt : error`
    # (:76), yyerror latches yysynerr (yyerror.c:19) and initprz.c:75-77
    # stops the job with SYNTAX. The key resolved to a real option, so this
    # cannot be a %repeat DATA step's keyword (#431). A macro value does not
    # change the verdict: the `=` is written, whatever `&X` expands to.
    if (eqp > 0L && token %in% c("CONSERVE", "NOCONSERVE", "QUASINEWTON",
                                 "STEEPEST", "PRINTIT", "NOPRINT", "NOCOR",
                                 "NOCOV", "NOLOG", "NONOTES", "NUMERIC")) {
      proc_rejected <- c(proc_rejected, paste0(
        tok, ": ", key, " takes no value in PROC HAZARD (hazard_y.y:65-76), ",
        "so it rejects this job with a syntax error"))
      proc_what <- c(proc_what, tok)
      note(tok, "a value on an option that takes none (hazard_y.y:65-76)")
      next
    }
    mapped <- mapped + 1L
    switch(token,
      # A WORK. libref names the same dataset as the bare name. Dropping it
      # here makes the emitted data =, the status chunk and the guard use the
      # bare name, which is also how a %repeat OUT= is recorded.
      DATA        = {
        # Strip the WORK libref BEFORE the presence check: `WORK.` is neither
        # NAME nor LIBMEM (hazard_l.l:39-40), so SAS rejects it, and checking
        # the unstripped "WORK." let it through to leave an empty name and an
        # internal R error (#433 review 3).
        stripped <- sub("^WORK[.]", "", val)
        data_given <- TRUE
        data_raw <- val
        if (check_name(key, stripped)) {
          mapped <- mapped - 1L
        } else {
          data_name <- stripped
        }
      },
      OUTHAZ      = {
        if (check_name(key, val)) {
          mapped <- mapped - 1L
        } else {
          outhaz <- val
        }
      },
      MAXITER     = {
        if (check_number(key, val)) {
          # Rejected: one row, already recorded by check_number().
          mapped <- mapped - 1L
        } else {
          val_num <- suppressWarnings(as.numeric(val))
          if (is.na(val_num)) {
            mapped <- mapped - 1L
            note("MAXITER", "non-numeric value for MAXITER")
          } else if (val_num < 0) {
            # hazpprc.c:23-24: a negative value never reaches the iteration
            # limit, which keeps its default (stmtprc.c:73), so PROC HAZARD
            # fits as though MAXITER were absent. Emitting it handed
            # hazard() a negative maxit, and it returned the starting values
            # with converged = TRUE (#496).
            mapped <- mapped - 1L
            # A repeated option overwrites the one before it, so this also
            # clears an earlier MAXITER (measured: `MAXITER=0 MAXITER=-1`
            # fits, `MAXITER=-1 MAXITER=0` evaluates; Copilot on #539).
            ctl$maxit <- NULL
            note(paste0("MAXITER=", val), paste0(
              "negative, so PROC HAZARD keeps its default iteration limit ",
              "(hazpprc.c:23-24); not emitted"))
          } else {
            ctl$maxit <- val_num
          }
        }
      },
      # Recorded, never emitted: hazard() reads no `condition` (#384).
      # CONDITION=n (3 to 14, hazpprc.c:48-56) stops PROC HAZARD's optimizer
      # as ill-conditioned once log10 of its Hessian approximation's
      # condition estimate exceeds n (setopt.c:452-456). hazard()'s optimizer
      # has no such stop; it warns about the final Hessian after the fit.
      CONDITION   = {
        mapped <- mapped - 1L
        val_num <- if (check_number(key, val)) NULL else
          suppressWarnings(as.numeric(val))
        if (is.null(val_num)) {
          # Rejected: one row, already recorded by check_number().
        } else if (is.na(val_num)) {
          note("CONDITION", "non-numeric value for CONDITION")
        } else if (val_num < 3 || val_num > 14) {
          # hazpprc.c:48-56 stores only 3..14; otherwise the limit stays at
          # stmtprc.c:74's 0 and setopt.c:454 applies the built-in test.
          note("CONDITION", paste0(
            "CONDITION=", val, " is outside the 3 to 14 PROC HAZARD accepts ",
            "(hazpprc.c:48-56), so PROC HAZARD ignores it and applies its ",
            "built-in conditioning limits (setopt.c:458-466). hazard() has ",
            "no conditioning stop either way"))
        } else {
          note("CONDITION", paste0(
            "CONDITION=", val, " stops PROC HAZARD's optimizer as ",
            "ill-conditioned once, after its first iteration, log10 of the ",
            "Hessian approximation's condition estimate exceeds it ",
            "(setopt.c:452-456). hazard() has no such stop: it fits, then ",
            "warns if the final Hessian is ill-conditioned"))
        }
      },
      CONSERVE    = ctl$conserve <- TRUE,
      NOCONSERVE  = ctl$conserve <- FALSE,
      # Recorded, never emitted: hazard() reads no `method` (#384). QUASI
      # chooses PROC HAZARD's optimizer; hazard() has no choice to make.
      QUASINEWTON = {
        mapped <- mapped - 1L
        note("QUASINEWTON", paste(
          "QUASI chooses PROC HAZARD's quasi-Newton optimizer. hazard() has",
          "no optimizer choice to make: it fits by BFGS, a quasi-Newton",
          "method (a multiphase fit with fixed shapes may run a Nelder-Mead",
          "warm-up first), continued with stats::nlm() when SAS's gradient",
          "test fails. The search path can differ, and on a multimodal",
          "likelihood so can the optimum"))
      },
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
  # Which statements carry a syntax refusal, and which this translation
  # cannot judge, so .hzr_sas_paren_reset() can place them against the
  # job's last `(` (#461). The PROC line is statement 1.
  err_stmt <- if (length(proc_rejected) || !is.null(proc_syntax_error)) 1L else
    integer(0)
  blind_stmt <- integer(0)
  parms_stmt <- integer(0)
  phase_what <- character(0)
  semantic_seen <- FALSE
  statements <- list()
  icensor_unread <- NULL
  parms_ops <- character(0)
  sel_ops <- NULL
  sel_bad <- character(0)
  saw_restrict <- FALSE
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
    # A macro expands before PROC HAZARD's lexer runs, into anything,
    # including a `(` or a syntax error.
    if (.hzr_sas_is_macro(stmt_text)) blind_stmt <- c(blind_stmt, i)
    if (is.na(token)) {
      # Recorded, not refused, for the same reason as an unknown PROC option
      # above: the block text can carry another step's keywords. PROC
      # HAZARD's lexer rejects it all the same, so it can set the syntax
      # flag again after a `(` (#461).
      blind_stmt <- c(blind_stmt, i)
      note(kw, "unknown HAZARD statement")
      next
    }
    # Only PARMS and the phase statements are checked here for every syntax
    # error PROC HAZARD raises. Any other statement can hide one this
    # translation does not see (`RESTRICT A*B`, `SELECTION SLE=ABC`,
    # `WEIGHT 2W`), which sets the flag again after a `(` and was measured
    # to be refused, so it cannot follow a `(` that is to clear the job.
    # SELECTION's operands are now checked by .hzr_selection_syntax() (N3,
    # #504 review), but that check is not shown to catch every error PROC
    # HAZARD raises there, so SELECTION stays blind here.
    if (!token %in% c("PARAMETERS", "EARLY", "CONSTANT", "LATE")) {
      blind_stmt <- c(blind_stmt, i)
    }
    if (token %in% c("EARLY", "CONSTANT", "LATE")) {
      pc <- .hzr_parse_phase_covars(ops_text)
      if (length(pc$semantic)) semantic_seen <- TRUE
      if (length(pc$rejected) > length(pc$semantic) ||
          length(pc$not_a_name)) {
        err_stmt <- c(err_stmt, i)
        phase_what <- c(phase_what, paste(
          token, c(setdiff(pc$rejected_what, pc$semantic),
                   pc$not_a_name_what)))
      }
    }
    # TIME, EVENT, RCENSOR, LCENSOR and WEIGHT each take exactly one NAME
    # (hazard_y.y:106, :109, :112, :124, :127). Any other count falls to
    # `otherstmt : error` (:102) and initprz.c:75-77 stops the job with
    # SYNTAX, so taking ops[[1L]] and dropping the rest fitted a model the
    # job did not describe, silently (#431). A macro operand can expand to
    # any number of names, so the count carries no verdict when one is there.
    if (token %in% c("TIME", "EVENT", "RCENSOR", "LCENSOR", "WEIGHT") &&
        length(ops) != 1L && !any(.hzr_sas_is_macro(ops))) {
      why <- paste0(
        if (length(ops)) "more than one operand" else "no operand",
        "; PROC HAZARD's ", kw, " takes exactly one variable name ",
        "(hazard_y.y:106-127)")
      proc_rejected <- c(proc_rejected, paste0(
        stmt_text, ": ", why, ", so it rejects this job with a syntax error"))
      proc_what <- c(proc_what, stmt_text)
      err_stmt <- c(err_stmt, i)
      note(stmt_text, why)
      if (!length(ops)) {
        if (token %in% c("TIME", "EVENT")) {
          stmt_empty[[token]] <- paste0(stmt_text, ": ", why)
        }
        next
      }
      # The extra operands are dropped from the fit below; the warning and
      # the row above say so. Refused, so not counted as mapped.
      mapped <- mapped - 1L
    }
    # `parmsopts` needs at least one operand (hazard_y.y:133-134), so a bare
    # `PARMS;` falls to `otherstmt : error` (:102). The binary refuses it
    # with SYNTAX; it was dropped without a word (#461 review).
    if (token == "PARAMETERS" && !length(ops)) {
      why <- paste0("no operand; PROC HAZARD's ", kw, " takes at least one ",
                    "(hazard_y.y:133-134)")
      proc_rejected <- c(proc_rejected, paste0(
        stmt_text, ": ", why, ", so it rejects this job with a syntax error"))
      proc_what <- c(proc_what, stmt_text)
      err_stmt <- c(err_stmt, i)
      note(stmt_text, why)
      mapped <- mapped - 1L
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
        # That is the whole grammar (hazard_y.y:115-122). The ICNS lexer
        # state returns NAME and `=` only (hazard_l.l:55, :84, :174-175);
        # any other text or character is reported and dropped
        # (hazard_l.l:176-179), so a comma anywhere, a missing `=` or an
        # extra name is a syntax error and initprz.c:75-77 stops the job.
        # The translation used to strip a trailing comma and fit, silently
        # (#495). It now reads the operand as that lexer does: runs of
        # `[.-_A-Z0-9]`, `=`, and single other characters, where a run is a
        # NAME only if the whole run is one (the longer rule wins). `)` is
        # whitespace to the lexer (hazard_l.l:32). A `(` switches it out of
        # ICNS into the PROC-line state and clears the flag (hazard_l.l:56),
        # where anything but whitespace, `)` or another `(` sets it again
        # (measured: `C3=TL()` fits, `C3=TL(X)`, `C3=TL()=` and
        # `C3=TL() = 1` exit SYNTAX), so the text after it is checked for
        # that and the text before it is read as ICENSOR. A macro reference
        # is joined to the run it touches (`C&I`, `TL&I.`, `&&C3`), as SAS
        # resolves it into one name, and counts as a name. A macro CALL
        # (`%TRIM(TL)`) owns its `(`: it stands in for one macro token, so
        # the text around it is still read (`C3=TL, %TRIM(X)` keeps its
        # comma whatever the call expands to; Copilot on #546), but no
        # `count = timevar` is taken from a statement that holds one.
        call_re <- "%[A-Z_][A-Z0-9_]*[[:space:]]*[(][^()]*[)]"
        macro_call <- grepl(call_re, ops_text)
        icns_all <- gsub(")", " ", gsub(call_re, " &MACROCALL ", ops_text),
                         fixed = TRUE)
        icns <- sub("[(].*$", "", icns_all)
        tail <- if (!grepl("(", icns_all, fixed = TRUE)) "" else
          sub("^[^(]*[(]", "", icns_all)
        # Only the macro references themselves carry no verdict: text beside
        # them (`TL(X,&M)` keeps `X,`) sets the flag whatever they expand to
        # (Copilot on #546).
        tail_bad <- nzchar(gsub("[()[:space:]]", "",
                                gsub("(&+|%)[A-Z_][A-Z0-9_]*[.]?", "", tail)))
        toks <- regmatches(icns, gregexpr(
          "[-._A-Z0-9&%]+|=|[^[:space:]]", icns))[[1L]]
        is_macro_tok <- .hzr_sas_is_macro(toks)
        is_name_tok <- .hzr_sas_is_name(toks) | is_macro_tok
        err_tok <- !is_name_tok & toks != "="
        kept <- toks[!err_tok]
        kept_name <- is_name_tok[!err_tok]
        well_formed <- !macro_call && !tail_bad && !any(err_tok) &&
          length(kept) == 3L && kept_name[[1L]] &&
          identical(kept[[2L]], "=") && kept_name[[3L]]
        # A macro can expand to any number of names, so a count or order
        # mismatch carries no verdict when one is present. A stray character
        # outside the macro does: PROC HAZARD meets it whatever the macro
        # expands to (r-reviewer on #546).
        refused <- !well_formed &&
          (tail_bad || any(err_tok) || !any(is_macro_tok))
        if (refused) {
          why <- paste0("PROC HAZARD's ICENSOR is `ICENSOR count = timevar`, ",
                        "two names and nothing else (hazard_y.y:115-122; ",
                        "hazard_l.l:174-179)")
          proc_rejected <- c(proc_rejected, paste0(
            stmt_text, ": ", why, ", so it rejects this job with a syntax ",
            "error"))
          proc_what <- c(proc_what, stmt_text)
          err_stmt <- c(err_stmt, i)
          note(stmt_text, why)
          # Refused, so not counted as mapped, as for the #431 statements.
          mapped <- mapped - 1L
        }
        # What PROC HAZARD's parser reads once its lexer has dropped the
        # errors: the first NAME '=' NAME, whose actions fire before any
        # later token (setvar(14), setvar(15)). Measured after a clearing
        # `(`: `C3=TL,AGE` and `C3=TL AGE` use C3 and TL, and `C,3=TL`
        # uses C. Anything else leaves no ICENSOR to emit.
        if (!macro_call && length(kept) >= 3L && kept_name[[1L]] &&
              identical(kept[[2L]], "=") && kept_name[[3L]]) {
          statements$ICENSOR <- kept[c(1L, 3L)]
          if (!well_formed && !refused) {
            note(stmt_text, paste0(
              "a macro in this ICENSOR statement hides its shape; the fit ",
              "takes ", kept[[1L]], " = ", kept[[3L]], " and cannot tell ",
              "what PROC HAZARD reads once the macro expands"))
          }
        } else if (!refused) {
          # A macro hides the shape and no `count = timevar` can be read, so
          # the fit below omits interval censoring. That is a different
          # model, and it warns rather than leaving only a row (r-reviewer
          # pass 2 on #546).
          mapped <- mapped - 1L
          icensor_unread <- stmt_text
          note(stmt_text, paste0(
            "a macro in this ICENSOR statement hides its shape and no ",
            "`count = timevar` can be read, so the fit omits interval ",
            "censoring"))
        }
      },
      LCENSOR    = statements$LCENSOR <- ops[[1L]],
      RCENSOR    = statements$RCENSOR <- ops[[1L]],
      WEIGHT     = statements$WEIGHT <- ops[[1L]],
      # Accumulate rather than overwrite: PROC HAZARD keeps its PARMS fields
      # in one table that parmprc() reads once, after every statement has been
      # processed (stmtprc.c), so a second PARMS statement adds to the first
      # rather than replacing it. Overwriting also made the no-phase refusal
      # below fire on `PARMS MUE=0.2 THALF=1; PARMS FIXNU;` -- a job the
      # reference runs -- because only the trailing statement survived.
      PARAMETERS = {
        parms_ops <- c(parms_ops, ops)
        parms_stmt <- c(parms_stmt, i)
      },
      # SAS accepts `SLE = 0.2`. Splitting that on whitespace left three
      # tokens, recorded as untranslated, and the screen ran at the DEFAULT
      # threshold instead, so close the spaces around `=` first.
      # PROC HAZARD accumulates SELECTION statements across the job, and a
      # repeated option is last-wins (#505, measured: `SLE=0.05; SELECTION
      # SLS=0.1;` screens at 0.05, and BACKWARD in either statement makes
      # it backward). Assigning kept only the last statement.
      STEPWISE   = {
        new_ops <- strsplit(gsub("\\s*=\\s*", "=", ops_text), "\\s+")[[1L]]
        new_ops <- new_ops[nzchar(new_ops)]
        # The operands PROC HAZARD rejects with a syntax error (N3 and the
        # #504 review): the job warns under U1, as for MAXITER, and the
        # screen still runs on what .hzr_selection_syntax() keeps.
        chk <- .hzr_selection_syntax(new_ops)
        sel_ops <- c(sel_ops, chk$keep)
        if (length(chk$bad)) {
          proc_rejected <- c(proc_rejected, paste0(
            names(chk$bad), ": ", chk$bad, ", so PROC HAZARD rejects this ",
            "job with a syntax error"))
          proc_what <- c(proc_what, names(chk$bad))
          sel_bad <- c(sel_bad, chk$bad)
          err_stmt <- c(err_stmt, i)
        }
      },
      # RESTRICT constrains which variables the screen may select
      # (hazrd4.c's rsttbl). It is recorded here and refused below when the
      # job also has a SELECTION: a screen that ignored it would select by a
      # different rule than the job asked for.
      RESTRICT   = {
        saw_restrict <- TRUE
        # Without a SELECTION there is no screen for it to constrain, and
        # nothing here implements it, so it stays a recorded gap rather than
        # a mapped token. With one, the refusal below names it.
        mapped <- mapped - 1L
        note(kw, "no R equivalent")
      },
      # A second statement for a phase adds to its list (hazard_y.y appends
      # every phasevar); assigning replaced it and dropped the first (#342).
      EARLY      = covars$early <- c(covars$early, ops_text),
      CONSTANT   = covars$constant <- c(covars$constant, ops_text),
      LATE       = covars$late <- c(covars$late, ops_text),
      {
        mapped <- mapped - 1L
        note(kw, "no R equivalent")
      }
    )
  }

  # --- assemble -----------------------------------------------------------
  # The SELECTION spec is read BEFORE the PARMS block, because it decides
  # what a bare phase variable means: a candidate outside the model under
  # SELECTION, an ordinary covariate without it (setstat.c).
  sel <- if (is.null(sel_ops)) NULL else .hzr_selection_spec(sel_ops)
  parms <- .hzr_parse_parms(
    parms_ops, covars = covars,
    selection = if (is.null(sel)) FALSE else
      if (identical(sel$direction, "backward")) "backward" else "screen")
  untr <- rbind(untr, parms$untranslated)
  # PROC HAZARD itself refuses this job at parse: a syntax error in a phase
  # statement (hazard_l.l:176 -> initprz.c:75-77) or ORDER= with /E, /I or /S
  # (przconc.c:45-53 -> hazard.c:249-251). It produces no estimates, so a fit
  # here would answer a job the reference never runs (#340). Checked first:
  # SAS stops at parse, before anything the other refusals read -- including
  # the censoring spec, which throws on a job with no EVENT (#396 review).
  # Split by provenance, not by message text. A PHASE statement PROC HAZARD
  # refuses at parse has always stopped the document (#340) and still does.
  # A PARMS operand or PROC-line value it refuses used to emit a fit with an
  # untranslated row; since 2026-09-22 it emits the fit, the row AND a loud
  # warning, so a rendered document completes and carries the reason rather
  # than halting on it.
  # A TIME or EVENT with no operand (#431) takes this route too when no other
  # statement supplies the variable (a second TIME, or ICENSOR for EVENT):
  # there is then nothing to fit, so no fit to emit with a warning above it.
  # Otherwise it warns below like the other #431 refusals.
  stmt_fatal <- c(
    if (is.null(statements$TIME)) stmt_empty$TIME,
    if (is.null(statements$EVENT) && is.null(statements$ICENSOR))
      stmt_empty$EVENT)
  # A later `(` clears PROC HAZARD's syntax-error flag (hazard_l.l:56), so a
  # syntax refusal before the job's last `(` does not stop the job (#461).
  # A SEMANTIC refusal is not a syntax error and nothing clears it.
  # The PARMS operands are parsed together, so a PARMS refusal is placed at
  # the LAST PARMS statement. That can only keep a refusal PROC HAZARD would
  # have cleared, never clear one it keeps.
  if (length(parms$rejected_parms)) {
    err_stmt <- c(err_stmt, max(parms_stmt))
  }
  reset <- if (semantic_seen) "none" else
    .hzr_sas_paren_reset(st, err_stmt, blind_stmt)
  # The SELECTION rows say the job is refused only when no `(` clears it;
  # otherwise the verdict is the #461 warning's below.
  for (k in seq_along(sel_bad)) {
    note(names(sel_bad)[[k]], paste0(
      sel_bad[[k]],
      if (identical(reset, "none")) {
        ", so PROC HAZARD rejects this job with a syntax error"
      }))
  }
  what <- c(proc_syntax_what, proc_what, parms$rejected_parms_what, phase_what)
  if (length(parms$rejected_phase) || length(stmt_fatal)) {
    # With no TIME or EVENT there is nothing to fit whatever the flag says,
    # and a refusal in the `(`'s own statement was measured to stop before
    # fitting, so only "cleared" changes this verdict.
    msg <- if (identical(reset, "cleared") && !length(stmt_fatal)) {
      paste0("This translation cannot emit PROC HAZARD's model for this job: ",
             .hzr_sas_paren_cleared(what), "what it keeps may be text ",
             "this translation cannot place in a model (for example, after ",
             "`EARLY AGE=ABC;` it keeps AGE). Correct the statement(s) named ",
             "here and translate the job again.")
    } else {
      paste0(
        "PROC HAZARD does not run this job: ",
        paste(c(parms$rejected_phase, stmt_fatal), collapse = "; "),
        ". Correct the ",
        "statement(s) named here and translate the job again.")
    }
    return(list(
      call = as.call(list(quote(stop), msg, call. = FALSE)),
      status_call = NULL, outhaz = outhaz, untranslated = untr,
      tokens_seen = seen, tokens_mapped = mapped
    ))
  }
  refusal_warnings <- character(0)
  if (!is.null(icensor_unread)) {
    refusal_warnings <- c(refusal_warnings, paste0(
      .hzr_sas_not_mirrored_lead, "`", icensor_unread, "` carries a macro ",
      "this translation cannot resolve into `count = timevar`, so the fit ",
      "below has no interval-censored rows. Resolve the macro and translate ",
      "the job again."))
  }
  # A phase variable that is not a NAME is refused at parse too, but it
  # translated on main, so it warns here instead of stopping above (#440).
  rejected <- c(proc_rejected, parms$rejected_parms, parms$rejected_name)
  if (!is.null(proc_syntax_error)) {
    rejected <- c(proc_syntax_error, rejected)
    note("PROC HAZARD", proc_syntax_error)
  }
  # The verdict is the job's, not the construct's: each construct named is a
  # syntax error, and whether PROC HAZARD runs the job depends on where the
  # job's last `(` falls (#461). The rows keep their reasons either way. A
  # job PROC HAZARD runs on a model this translation cannot emit is the U1
  # class the FIXMNU1 warning below belongs to, so it takes that wording.
  if (length(rejected) && identical(reset, "cleared")) {
    refusal_warnings <- c(refusal_warnings, paste0(
      .hzr_sas_not_mirrored_lead, .hzr_sas_paren_cleared(what),
      "the recovery can keep text this fit leaves out (after ",
      "`EARLY AGE*SEX;` it keeps AGE) and drop text this fit keeps (after ",
      "`NU=ABC`, the rest of that PARMS statement). Correct the ",
      "statement(s) named here and translate the job again."))
  } else if (length(rejected) && identical(reset, "same")) {
    refusal_warnings <- c(refusal_warnings, paste0(
      .hzr_sas_not_mirrored_lead, "the syntax error in ",
      paste(what, collapse = "; "), " is followed by a `(` in the same ",
      "statement, and PROC HAZARD's lexer clears its syntax-error flag at ",
      "every `(` (hazard_l.l:56). Whether PROC HAZARD then fits the job ",
      "depends on its parser's error recovery, which this translation does ",
      "not reproduce, so it cannot tell whether PROC HAZARD refuses this job ",
      "or which model it fits: `EARLY AGE*SEX, LOG();` fits AGE alone, and ",
      "`EARLY 1AGE, SEX, LOG();` stops before fitting. Correct the ",
      "statement(s) named here and translate the job again."))
  } else if (length(rejected)) {
    refusal_warnings <- c(refusal_warnings, paste0(
      "PROC HAZARD does not run this job: ",
      paste(rejected, collapse = "; "), ". The fit below is this ",
      "translation's, not one PROC HAZARD would produce. Correct the ",
      "statement(s) named here and translate the job again."))
  }

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

  # The second refusal, and for the same reason as the LCENSOR + ICENSOR one
  # above: PROC HAZARD does not run this job at all. With no active phase
  # modterm.c:18-22 raises ERROR 1001 ("No phase selected") and hazard.c:299-302
  # exits via hzfxit("SEMANTIC") BEFORE results(), so the reference produces no
  # estimates -- there is no fit here for a translation to be faithful to. An
  # $untranslated row alone leaves the fit chunk in place, and a reader who
  # renders past the callout gets a converged single-distribution fit with a
  # populated summary standing in for a job that produced nothing.
  #
  # Placed AFTER the censoring refusal so a job carrying both keeps that one's
  # more specific message, matching the precedence already established there.
  if (isTRUE(parms$refused)) {
    return(list(
      call = quote(stop(
        "This PROC HAZARD job selects no phase: no PARMS statement named a ",
        "positive MUE, MUC or MUL. PROC HAZARD refuses such a job outright -- ",
        "modterm.c raises ERROR 1001, \"No phase selected\", and the ",
        "procedure exits before computing any results -- so there is no fit ",
        "to translate. Add the PARMS statement naming the phase(s) this ",
        "model needs, or fit it by hand as a single-distribution model if ",
        "that is what you intend.",
        call. = FALSE
      )),
      status_call = NULL, outhaz = outhaz, untranslated = untr,
      tokens_seen = seen, tokens_mapped = mapped
    ))
  }

  # A job SETG3 refuses. PROC HAZARD sets the error in shape() and exits
  # before results(), so nothing is fitted and there is no fit to translate.
  # The row alone left the hazard() chunk in place, and a reader who rendered
  # past the callout got a converged fit standing in for a job that produced
  # nothing -- the shape this package calls its signature defect. The code and
  # the operands travel in the message, so the reader knows what to change
  # (#359).
  if (!is.null(parms$refusal_reason) && !is.na(parms$refusal_reason)) {
    refusal_warnings <- c(refusal_warnings, paste0(
      "This PROC HAZARD job is refused before any fit is computed: ",
      sub("^PROC HAZARD refuses this job: ", "", parms$refusal_reason),
      ". ", if (grepl("SETG1 raises", parms$refusal_reason, fixed = TRUE))
        "SETG1" else "SETG3",
      " sets the error in shape() and the procedure exits before ",
      "results(), so PROC HAZARD produces nothing for this job and the fit ",
      "below stands in for no SAS result at all. Correct the PARMS ",
      "operand(s) named here, or fit the model by hand."))
  }
  # A job PROC HAZARD accepts and starts to fit, on a case its fit cannot
  # evaluate on the data measured (#424). Not a refusal, and not a different
  # model. "may_not_fit" is the data-dependent half: SAS stopped on one
  # dataset and fitted on another, so the tail says it MAY have no result
  # (#468 review).
  if (!is.null(parms$no_result_reason) && !is.na(parms$no_result_reason)) {
    refusal_warnings <- c(refusal_warnings, paste0(
      parms$no_result_reason, ". ",
      if (identical(parms$no_result_kind, "may_not_fit")) {
        paste0("On your data PROC HAZARD may have printed no estimates, so ",
               "check its listing before comparing the fit below with it.")
      } else {
        paste0("The fit below may stand in for no SAS result at all. ",
               "Correct the PARMS operand(s) named here, or fit the model ",
               "by hand.")
      }))
  }

  # A job PROC HAZARD runs, but on a model this translation does not emit:
  # FIXMNU1 constrains |M*NU| = 1 on the early phase (setg1.c:363-575,
  # hzd_early_t2p.c:65-77), and the emitted phase would estimate M and NU
  # without it. A fit here would be a different model standing in for the
  # job's, so the job warns and still fits (U1, #358). Mirroring the constraint is
  # new modelling, left out of 1.3.0.
  if (length(parms$not_mirrored)) {
    refusal_warnings <- c(refusal_warnings, paste0(
      .hzr_sas_not_mirrored_lead,
      paste(parms$not_mirrored, collapse = "; "), "."))
  }

  # A PARMS statement that builds no phase and is NOT refused -- operands this
  # parser could not read (a template's `MUE=?`). A MUE or MUL with no shape
  # operand no longer lands here: it builds on PROC HAZARD's own shape
  # defaults (#345). This is not a claim about PROC HAZARD, which is why it
  # is kept apart from `refused` above. Theta blocks are built only alongside
  # phases, so the call below would carry theta = c() under hazard()'s
  # default Weibull: a model PROC HAZARD never fits, with no starting values.
  # hazard() refuses that call now; before, it returned an unfitted object and
  # the chunk rendered as if it held a result. Stop here instead, where the
  # message can name the real cause.
  if (!isTRUE(parms$has_phases)) {
    return(list(
      call = quote(stop(
        "This job's PARMS statement builds no phase this translator could ",
        "use, so there are no starting values to fit from. The operands it ",
        "could not read or use are listed in $untranslated -- for example ",
        "`?` placeholders left for the reader to fill in, an operand ",
        "written with spaces around `=` (`MUE = 0.2`, which this translator ",
        "splits apart). This is a ",
        "limit of the translation, not a PROC HAZARD refusal. Correct those ",
        "operands and translate again, or fit the model by hand.",
        call. = FALSE
      )),
      status_call = NULL, outhaz = outhaz, untranslated = untr,
      tokens_seen = seen, tokens_mapped = mapped
    ))
  }

  # Why this job's SELECTION cannot be run faithfully, or NULL when it can.
  # A function because two paths need it: the screen below, and the no-DATA=
  # refusal, which returns first and must still record this reason so a reader
  # who adds DATA= as told does not meet a refusal the document never named.
  selection_refusal <- function(untr) {
    if (is.null(sel)) return(NULL)
    # Constructs with no faithful translation. A per-variable MOVE=
    # or ORDER= has no hzr_stepwise() equivalent at all (its max_move is
    # per run), and ORDER= drives entry order.
    per_var_opts <- grep("/(MOVE|ORDER)", untr$construct, value = TRUE)
    # force_in has no phase, so a variable held by /I in one phase and
    # movable in another would be pinned in BOTH. That is a wrong model,
    # not a path difference, which is where this draws the refuse line.
    # Both sides restricted to BUILT phases: a /I naming a phase this job
    # does not select has no phase to conflict with, and the discarded
    # phase's covariates are already recorded.
    movable_all <- unique(unlist(parms$selection$movable %||% list()))
    in_model_all <- unique(unlist(parms$selection$in_model %||% list()))
    cross_pinned <- intersect(
      intersect(parms$selection$force_in %||% character(0), in_model_all),
      movable_all)
    # A name PROC HAZARD accepts but R does not (`_X1`, `TRUE`) is NOT
    # refused. #411 refused it, citing a `/I` pin that never matched and a
    # score criterion that could not test the candidate. Measured through
    # this translator (#459), #437 had fixed the pin before the refusal
    # reached main, and #449 (PR #455) fixed the score path, so the job is
    # screened like any other.
    refusals <- c(sel$refuse,
                  if (saw_restrict) "RESTRICT",
                  if (length(per_var_opts)) per_var_opts,
                  if (length(cross_pinned)) {
                    paste0(cross_pinned, " (/I in one phase, movable in another)")
                  })
    if (!length(refusals)) return(NULL)
    # Name ONLY what fired. One boilerplate string listing every refusable
    # construct made the reason unreadable and, worse, made a test
    # asserting "FAST" pass for a MAXVARS job.
    why <- c(
      FAST = "FAST is a different search (H->f), not a stepwise run",
      MAXVARS = paste("MAXVARS caps the selected set and hzr_stepwise()",
                      "has no equivalent"),
      RESTRICT = paste("RESTRICT constrains which variables may be",
                       "selected (hazrd4.c's rsttbl)"))
    explain <- function(item) {
      if (!is.na(why[item])) return(unname(why[item]))
      if (startsWith(item, "MAXSTEPS=")) {
        return(paste0(item, ": PROC HAZARD refuses a negative MAXSTEPS ",
                      "(stpwprc.c:76-79), so there is no run to translate"))
      }
      if (grepl("/(MOVE|ORDER)", item)) {
        return(paste0(item, ": a per-variable MOVE= or ORDER= has no ",
                      "hzr_stepwise() equivalent (its max_move is per run)"))
      }
      paste0(item, ": force_in is keyed by variable name across phases, ",
             "so it would be pinned in the phase SAS leaves movable")
    }
    paste0(
      "SELECTION carries ", paste(refusals, collapse = ", "),
      ", which this translator cannot run faithfully. ",
      paste(vapply(refusals, explain, character(1)), collapse = "; ")
    )
  }

  # A phase covariate is evaluated only in `data`, and with no DATA= there is
  # none: hazard() refuses the call (#299) and tells the reader to pass
  # `data =`, which the SAS job never had. Read the emitted phase calls rather
  # than the SAS text, so this tracks what the chunk would carry (#311).
  # The emitted phases are NOT the whole question, though: SELECTION withholds
  # its candidates from them, and /E variables are never in them, yet the
  # screen's scope and the listwise guard read both from `data` too. Reading
  # only the phase calls made this refusal blind to exactly what SELECTION
  # withholds (#160), so every phase-statement variable outside the base
  # model (listwise_only) counts as well.
  phase_has_formula <- vapply(
    as.list(parms$phases)[-1L],
    function(ph) is.call(ph) && !is.null(ph[["formula"]]),
    logical(1)
  )
  # A SELECTION job always needs `data`, even with no phase variable at all:
  # hzr_stepwise() refits each candidate from it and stops without it.
  # And a job with no variable to read still names no DATA=, which the
  # %HAZARD macro refuses before PROC HAZARD runs (hazard.sas:11-15,
  # :145-152), so there is no fit to translate whatever the phases carry
  # (#497). Such a job used to fit from whatever the session held.
  if (is.null(data_name)) {
    # Say what is actually true of this job. "A phase has covariates" is
    # false when every named variable sits outside the emitted phases: a
    # SELECTION candidate, an /E variable, or a covariate of a phase that is
    # not built. Those are read by the screen or by the missing-row guard.
    why_data <- if (any(phase_has_formula)) {
      paste0("but a phase has covariates, and hazard() evaluates a phase's ",
             "covariates only in `data`.")
    } else if (!is.null(sel)) {
      paste0("but it carries a SELECTION statement, and hzr_stepwise() needs ",
             "`data` to refit each candidate",
             if (length(parms$listwise_only)) paste0(
               " (this job's are ", paste(parms$listwise_only, collapse = ", "),
               ")"),
             ".")
    } else if (length(parms$listwise_only)) {
      paste0("but its phase statements name ",
             paste(parms$listwise_only, collapse = ", "),
             ", which are outside the fitted model. With no dataset the ",
             "translation cannot say where to read them, and PROC HAZARD ",
             "deletes rows where any is missing.")
    }
    lead <- if (is.null(why_data)) {
      "names no DATA= dataset."
    } else {
      paste("names no DATA= dataset,", why_data)
    }
    # A block whose first statement is not the PROC statement (a %repeat
    # call it encloses, say) was read for options from the wrong statement,
    # so this cannot say SAS finds no DATA=.
    # A macro on the PROC statement is expanded before %HAZARD reads it, so
    # it may supply the DATA= this cannot see. And %HAZARD's test for DATA=
    # is a substring test (hazard.sas:11), so `OUTHAZ=MYDATA` passes it and
    # the macro reads whatever follows the next `=` (:13-17) as its dataset.
    # Claim the macro's refusal only where neither applies (#497 review).
    proc_macros <- toks[vapply(toks, .hzr_sas_is_macro, logical(1))]
    macro <- if (!startsWith(trimws(st[[1L]]), "PROC HAZARD")) {
      paste(
        "This translation read the PROC HAZARD options from the block's",
        "first statement, which is not the PROC HAZARD statement, so it may",
        "have missed a DATA= that SAS reads.")
    } else if (length(proc_macros)) {
      paste0(
        "The PROC HAZARD statement carries ",
        paste(proc_macros, collapse = ", "), ", which SAS expands before ",
        "the %HAZARD macro reads the statement, so it may supply a DATA= ",
        "this translation cannot see.")
    } else if (identical(data_raw, "WORK.")) {
      paste(
        "SAS does not run it either: PROC HAZARD reads WORK as the dataset",
        "name and the `.` after it as unexpected text (hazard_l.l:79-80,",
        ":176), a syntax error.")
    } else if (data_given) {
      paste(
        "SAS does not run it either: PROC HAZARD has no form of DATA=",
        "without a dataset name (hazard_y.y:61-62, :80-81), so it rejects",
        "the job with a syntax error.")
    } else if (grepl("DATA", st[[1L]], fixed = TRUE)) {
      paste(
        "SAS does not fit it from a dataset the job names either: the",
        "%HAZARD macro's test for DATA= is a substring test (hazard.sas:11),",
        "which the DATA elsewhere on the PROC statement satisfies, so the",
        "macro takes whatever follows the next `=` as its dataset",
        "(hazard.sas:13-17).")
    } else {
      paste(
        "SAS does not run it either: the %HAZARD macro finds its dataset only",
        "through DATA= on the PROC statement, and without it stops with",
        "\"HAZARD not attempted\" (hazard.sas:11-15, :145-152).")
    }
    untr <- rbind(untr, .hzr_untranslated_frame(
      NA_integer_, "DATA=",
      paste("the job", lead, macro, "Add DATA= to the job and translate",
            "again, or fit it by hand (#311, #497).")
    ))
    if (!is.null(sel)) {
      untr <- rbind(untr, sel$untranslated)
      reason <- selection_refusal(untr)
      if (!is.null(reason)) {
        untr <- rbind(untr, .hzr_untranslated_frame(NA_integer_, "SELECTION",
                                                    reason))
      }
    }
    return(list(
      call = as.call(list(quote(stop), paste(
        "This PROC HAZARD job", lead, macro,
        "Name the dataset with DATA= and translate the job again, or fit the",
        "model by hand."
      ), call. = FALSE)),
      status_call = NULL, outhaz = outhaz, untranslated = untr,
      tokens_seen = seen, tokens_mapped = mapped
    ))
  }

  # The status expression is hoisted into its own chunk, named .hzr_status so
  # it cannot collide with a SAS variable, because every job carries at least
  # one guard that must run, and be seen to run, ahead of the fit: the
  # missing/negative guard on each named count, plus the both-fire guard when
  # the job named more than one. ICENSOR jobs additionally need the name so
  # time_lower can gate on it. See .hzr_censor_spec() roxygen. The refused
  # path returned above, so status_name is non-NULL from here on.
  # time_upper is never emitted: the ICENSOR interval's upper bound is TIME,
  # and hazard()'s time_upper already defaults to TIME when left NULL (see
  # .hzr_censor_spec() roxygen).
  args <- list()
  if (!is.null(data_name)) args$data <- as.name(data_name)
  args$time <- as.name(statements$TIME)
  # The status chunk is evaluated outside hazard(), so it has to resolve
  # its own bare SAS column names. When the job reads a DATA= dataset the
  # hoisted status is written back INTO that data frame rather than bound
  # locally: hazard()'s vector path gives data columns precedence over the
  # calling frame, so a local binding is shadowed by any column of the same
  # name, and both `status =` and the `time_lower` gate would then read the
  # user's column instead of the computed censoring classification -- a
  # different censoring structure, from a fit that raises nothing but an
  # ambiguity warning. Derived as a column, the mask resolves to the value
  # this chunk just wrote, so it cannot be shadowed at all. `.hzr_` is this
  # package's reserved prefix, so overwriting a column of that name is the
  # intended consequence, not collateral damage. transform() masks exactly
  # as with() did, so the expression's column names still resolve against
  # the dataset. A job with no DATA= was refused above (#497), so there is
  # always a dataset here.
  derive <- as.call(list(quote(transform), as.name(data_name),
                         cens$status_expr))
  names(derive) <- c("", "", as.character(cens$status_name))
  # Degenerate ICENSOR bounds (readct.c, see .hzr_censor_spec()): the fit
  # reads only the rows PROC HAZARD keeps (`args$data` below), and this
  # column says how many rows each rule touched. The count is raised inside
  # the transform(), so the chunk keeps its `<data> <- transform(<data>,
  # ...)` shape, and local() binds nothing in the reader's session. The
  # caller's data frame keeps every row.
  if (!is.null(cens$keep_expr)) {
    derive$.hzr_keep <- bquote(local({
      .keep <- .(cens$keep_expr)
      .n_event <- sum(.(cens$degenerate_expr))
      .n_drop <- sum(!.keep)
      if (.n_event > 0 || .n_drop > 0) {
        warning("Degenerate ICENSOR intervals, resolved as PROC HAZARD ",
                "does (readct.c): ", .n_event,
                if (.n_event == 1) " row" else " rows",
                " with CTIME equal to TIME fitted as exact events ",
                "(readct.c:18-23), and ", .n_drop,
                if (.n_drop == 1) " row" else " rows",
                " with CTIME missing, negative or after TIME dropped from ",
                "the fit (readct.c:9-17).", call. = FALSE)
      }
      .keep
    }))
    derive$.hzr_icensor_event <- cens$degenerate_expr
  }
  status_call <- call("<-", as.name(data_name), derive)
  # Variables in some fitted phase formula: hazard() drops their missing rows
  # itself. Read by the listwise guard and the column check below.
  modelled <- unique(unlist(lapply(as.list(parms$phases)[-1L], function(ph) {
    if (is.call(ph) && !is.null(ph[["formula"]])) all.vars(ph[["formula"]])
  })))
  # Phase variables outside every fitted formula still delete their missing
  # rows in PROC HAZARD (see .hzr_parse_parms()), and hazard() cannot see
  # them. Stop in the status chunk, ahead of the fit, rather than fit more
  # rows than SAS did.
  if (length(parms$listwise_only)) {
    lw <- lapply(parms$listwise_only, as.name)
    any_na <- Reduce(function(a, b) call("|", a, b),
                     lapply(lw, function(v) call("is.na", v)))
    # Only rows hazard() would KEEP matter: where a modelled phase variable
    # is also missing, hazard() drops the row itself, which is what SAS does
    # too, so there is nothing to stop for (#340 item 8).
    if (length(modelled)) {
      kept <- Reduce(function(a, b) call("&", a, b),
                     lapply(lapply(modelled, as.name),
                            function(v) call("!", call("is.na", v))))
      any_na <- call("&", call("(", any_na), call("(", kept))
    }
    msg <- paste(
      "This job has rows where a phase variable outside the fitted base",
      "model --", paste(parms$listwise_only, collapse = ", "), "-- is",
      "missing. Such a variable is excluded (/E), belongs to a phase this",
      "job does not select, or is a SELECTION candidate the screen offers",
      "rather than the base model carrying it. PROC HAZARD still deletes",
      "those rows (getrisk.c lists every phase variable, and readobs.c drops",
      "a row where any is missing), and hazard() cannot, because the base",
      "model it fits does not contain the variable. Drop the rows before",
      "fitting, subsetting to complete values of those variables."
    )
    guarded <- bquote({
      if (any(.(any_na))) stop(.(msg), call. = FALSE)
      .(cens$status_expr)
    })
    status_call[[3L]][[3L]] <- guarded
  }
  # A phase variable the dataset does not contain failed deep inside the
  # chunk as "object 'ZZ' not found", naming neither the statement nor the
  # dataset (#340 item 9). Check the names first and say so.
  phase_vars <- union(modelled, parms$listwise_only)
  if (!is.null(data_name) && length(phase_vars)) {
    dsym <- as.name(data_name)
    # No assignment: every symbol in emitted code is read as a data column by
    # the test oracle's synthetic data, and a new name there changed it.
    present <- bquote({
      if (!all(.(phase_vars) %in% names(.(dsym)))) {
        stop("This job's phase statements name ",
             paste(setdiff(.(phase_vars), names(.(dsym))), collapse = ", "),
             ", which are not columns of ", .(data_name), ". Add them to ",
             .(data_name), " or remove them from the phase statements.",
             call. = FALSE)
      }
      # PROC HAZARD refuses a non-numeric phase variable (vfynvar.c:22-26,
      # "VARIABLE NOT NUMERIC", sets semerr; hazard.c:249-251 exits), where
      # hazard() would dummy-code it and fit. FUN is passed as a string so
      # no new symbol reaches the emitted code.
      if (!all(vapply(.(dsym)[.(phase_vars)], "is.numeric", NA))) {
        stop("PROC HAZARD refuses a phase variable that is not numeric ",
             "(vfynvar.c:22-26). These phase variables are not numeric: ",
             paste(names(which(!vapply(.(dsym)[.(phase_vars)], "is.numeric",
                                       NA))), collapse = ", "),
             ". Convert them to numeric codes, as the SAS dataset holds them.",
             call. = FALSE)
      }
    })
    # Ahead of transform(), not inside it: inside, a column named like the
    # dataset masks it, names() of that column is NULL, and the check refused
    # a job whose variables were all present (#396 review).
    status_call <- as.call(c(as.name("{"), as.list(present)[-1L],
                             list(status_call)))
  }
  # Degenerate ICENSOR bounds (readct.c, see .hzr_censor_spec()): the fit
  # reads only the rows PROC HAZARD keeps, and the status chunk says how
  # many rows each rule touched. The caller's data frame keeps every row.
  if (!is.null(cens$keep_expr)) {
    dsym <- as.name(data_name)
    args$data <- bquote(.(dsym)[.(dsym)$.hzr_keep, , drop = FALSE])
  }
  args$status <- cens$status_name
  if (!is.null(cens$time_lower)) args$time_lower <- cens$time_lower
  # An ICENSOR job is fitted on the interval term PROC HAZARD accumulates,
  # C3 * log([CF(T) - CF(CT)] / (T - CT)) (setlik.c; see the roxygen above),
  # not on hazard()'s default interval probability. Measured on the binary
  # over a constructed grid (18 optimised fits), the default moved MUC by up
  # to 21% and the objective by up to 291 units; "sas" reproduced both to
  # printed precision (#543, maintainer's ruling 2026-09-29). The value
  # reported is then SAS's objective, not a log-likelihood, and the
  # translated document says so above the fit.
  if (!is.null(statements$ICENSOR)) args$objective <- "sas"
  # Without fit = TRUE the emitted call returns an unfitted object: converged
  # is NA, objective is NA, and theta holds the SAS starting values, while
  # print.hazard() shows a populated summary that says none of that (#151).
  args$fit <- TRUE
  # A PARMS statement that builds a phase makes this a multiphase job.
  # hazard()'s `dist` defaults to "weibull", and its
  # `else if (!is.null(phases))` branch silently discards the entire phase
  # specification (with only a warning) when dist stays at that default --
  # so dist = "multiphase" must be emitted whenever phases were built. A job
  # with no PARMS statement at all (or one that builds no phase) has
  # phase_calls = list(), i.e. parms$phases is the empty
  # `list()` call rather than NULL; omit both phases and dist on that path
  # so it stays a plain non-multiphase fit and does not trip that same
  # "'phases' is ignored" warning for a phases arg that was never meant to
  # carry anything.
  if (isTRUE(parms$has_phases)) {
    args$dist <- "multiphase"
    args$phases <- parms$phases
  }
  args$theta <- parms$theta
  # One expression, built in .hzr_censor_spec() from every count the job
  # named -- EVENT's C1, ICENSOR's C3, RCENSOR's C2 -- together with the
  # WEIGHT variable, which is not a count but multiplies two of the three.
  # It cannot be composed by multiplying a WEIGHT variable onto a count-only
  # expression, because setlik.c's c2 term -- the right-censored row --
  # enters the sum unweighted (#158, #162).
  if (!is.null(cens$weights_expr)) args$weights <- cens$weights_expr

  # Canonical control order, so the emitted call does not depend on the order
  # the options happened to appear in the SAS text.
  ctl <- ctl[intersect(c("maxit", "conserve"),
                       names(ctl))]
  if (length(ctl)) args$control <- as.call(c(quote(list), ctl))

  head <- quote(hazard)
  stepwise_call <- NULL
  screen_check_call <- NULL
  if (!is.null(sel)) {
    untr <- rbind(untr, sel$untranslated)
    reason <- selection_refusal(untr)
    if (!is.null(reason)) {
      untr <- rbind(untr, .hzr_untranslated_frame(NA_integer_, "SELECTION",
                                                  reason))
      return(list(
        call = as.call(c(quote(stop), list(paste0(reason, ".")),
                        list(call. = FALSE))),
        status_call = NULL, outhaz = outhaz, untranslated = untr,
        tokens_seen = seen, tokens_mapped = mapped
      ))
    }

    sw_args <- list(as.name("fit_base"))
    # A backward screen takes no scope: hzr_stepwise() ignores it there
    # (#343), and under BACKWARD a bare variable starts IN the model
    # (setstat.c), so the base already carries the full candidate set.
    if (!identical(sel$direction, "backward")) {
      # ALWAYS emit scope, with an explicit NULL for a phase that offers no
      # candidate. hzr_stepwise(scope = NULL) enumerates every data-frame
      # column not already in the model, so omitting the argument handed the
      # screen the whole SAS dataset -- the EVENT count and the time
      # variable included.
      scope <- parms$selection$scope %||% list()
      # Built from symbols, not pasted text: a SAS name may begin with an
      # underscore and an R symbol may not (#411). See
      # .hzr_sas_covar_formula().
      sw_args$scope <- as.call(c(quote(list), lapply(scope, function(v) {
        if (length(v)) .hzr_sas_covar_formula(v) else NULL
      })))
    }
    sw_args$data <- args$data
    # hzr_stepwise() forwards `...` to every refit. Without the job's control
    # there, the selected model is refitted at hazard()'s defaults: a
    # NOCONSERVE job's base ran without Conservation of Events and its
    # selected model with it.
    if (!is.null(args$control)) sw_args$control <- args$control
    sw_args$direction <- sel$direction
    # SAS enters on the score statistic and removes on Wald, which is what
    # criterion = "score" does here.
    sw_args$criterion <- "score"
    sw_args$slentry <- sel$slentry
    sw_args$slstay <- sel$slstay
    # MOVE= is NOT mapped onto max_move: they count different things.
    # PROC HAZARD counts DELETIONS only (hazrd4.c:361-377, an addition counts
    # only under NOSTEPWISE) and keys the counter per (variable, PHASE) theta
    # slot (setstat.c:15). hzr_stepwise() counts entries AND exits, keyed by
    # variable NAME across phases, and a frozen variable is then pinned both
    # in and out. Emitting SAS's MOVE = 1 therefore froze a variable that
    # simply entered two phases, which PROC HAZARD leaves movable (#160
    # review). Recorded rather than translated into a value that is not the
    # same quantity.
    # MAXSTEPS: PROC HAZARD's default is INT_MAX (stpwprc.c:73-74), not
    # hzr_stepwise()'s 50, so the default is written out like the others.
    sw_args$max_steps <- if (is.null(sel$max_steps)) .Machine$integer.max else
      sel$max_steps
    # Only a variable a BUILT phase carries: sel_force_in is collected
    # across all three phases, and hzr_stepwise() ignores a name no phase
    # has, silently.
    built_vars <- unique(c(unlist(parms$selection$movable %||% list()),
                           unlist(parms$selection$in_model %||% list())))
    force_in <- intersect(parms$selection$force_in %||% character(0),
                          built_vars)
    if (length(force_in)) {
      sw_args$force_in <- as.call(c(quote(c), as.list(force_in)))
    }
    stepwise_call <- as.call(c(quote(hzr_stepwise),
                               Filter(Negate(is.null), sw_args)))
    # A screen can stop because no candidate could be SCORED, which reads
    # exactly like "nothing met slentry" (#159). hzr_stepwise() now says so
    # itself when it stops; the check below covers a screen that completed.
    # `fit` is substituted for this block's own slot name by the caller, in
    # the message string as well as the code: a second SELECTION block used
    # to tell the reader to look at `fit`, which in that document is a
    # different object.
    #
    # AUTHORING RULE for anything emitted from here: it runs in the READER'S
    # session, so it may use only base R at the version DESCRIPTION declares
    # plus this package's exports -- never this package's internals. %||% is
    # the worked example (#160 review): this package defines its own, so the
    # 62 internal uses are fine, but base R gained it in 4.4 while Depends
    # says 4.1, and an emitted chunk using it errored for a reader on 4.1 to
    # 4.3. A grep for the construct mostly finds those false positives; the
    # question is always which side of the namespace boundary the code runs
    # on. Nothing checks this automatically.
    # Keyed on the REASONS, not n_uncomputable_scores: that is an attempt
    # count which, since #399, also counts removals and entries whose WALD
    # test had no variance (`wald_no_variance`). Those are not uncomputable
    # scores, and hzr_stepwise() now names each such variable itself, as it
    # does when a screen STOPPED on uncomputable candidates. This check says
    # only what the screen has not already said (#400).
    # Every local is dot-prefixed: the chunk runs in the reader's session, so
    # a bare name would overwrite the reader's object of that name (#400).
    screen_check_call <- bquote({
      .reasons <- fit$criteria$uncomputable_reasons
      if (is.null(.reasons)) .reasons <- integer(0)
      .reasons <- .reasons[names(.reasons) != "wald_no_variance"]
      .n_unscored <- sum(.reasons)
      if (.n_unscored > 0L && !isTRUE(fit$criteria$stopped_uncomputable)) {
        warning(.n_unscored, " candidate score(s) could not be computed in ",
                "this screen (", paste0(names(.reasons), " = ", .reasons,
                                  collapse = ", "),
                "); see ", .(quote(fit_label)),
                "$criteria$uncomputable_reasons. A candidate the screen ",
                "could not score was not tested at that step.", call. = FALSE)
      }
      invisible(.n_unscored)
    })
  }

  # MAXITER=0 (#496). PROC HAZARD sets its iteration limit to 0
  # (hazpprc.c:20-22), as it does for any value below 1, which the
  # assignment to an int truncates (hazpprc.c:27, common.h:27; measured on
  # MAXITER=0.5 and .9). A job with more than one free parameter then skips
  # the optimizer: NOOPTIM() prints the log-likelihood at the starting values
  # (hazrd2.c:71-74, :133-145). Emitting `control = list(maxit = 0)` handed
  # hazard() a job to optimise, and it did, reporting converged = TRUE at a
  # likelihood PROC HAZARD never printed. The same evaluation is emitted
  # instead, as hzr_evaluate() on the job's model at its starting values.
  #
  # Unless NOCONSERVE is given, Conservation of Events has run before that
  # (setcoe(), shape.c:52; stmtprc.c:64 and :123-127 set the mode, and the
  # listing names it, cmpmeth.c:32-38): every MU is scaled by one factor so
  # that the predicted events equal the observed. Along that common scaling
  # the log-likelihood is E * s - exp(s) * S plus a constant, so the factor is
  # also where the log-likelihood peaks, and a one-dimensional maximisation
  # over a shift of every log(MU) finds it. Measured on the binary with and
  # without WEIGHT, LCENSOR and ICENSOR (tests/testthat/fixtures/
  # maxiter-zero-oracle.csv). An ICENSOR job's spec carries
  # objective = "sas" (#543), so its evaluation reproduces PROC HAZARD's
  # parameters and printed value too.
  code_body <- as.call(c(head, args))
  if (isTRUE(ctl$maxit < 1)) {
    maxit_label <- paste0("MAXITER=", format(ctl$maxit))
    note(maxit_label, paste0(
      "PROC HAZARD evaluates the log-likelihood at the starting values ",
      "without optimising (hazpprc.c:20-27, hazrd2.c:71-74); emitted as ",
      "hzr_evaluate(), which is not a fit"))
    refusal_warnings <- c(refusal_warnings, paste0(
      maxit_label, ": PROC HAZARD does not fit this job. It evaluates the ",
      "log-likelihood at the starting values (hazpprc.c:20-27, ",
      "hazrd2.c:71-74, :133-145)",
      if (!isFALSE(ctl$conserve)) {
        paste0(", after Conservation of Events has scaled every MU by one ",
               "factor (setcoe(), shape.c:52)")
      },
      ", and the chunk below does the same with hzr_evaluate(). Its result ",
      "is not a fit: it carries no standard errors, which PROC HAZARD may ",
      "print at those values, and predict() cannot use it",
      if (!is.null(stepwise_call)) {
        paste0(". The SELECTION screen is not run: PROC HAZARD still steps ",
               "through it, evaluating each step without optimising ",
               "(hazrd2.c:65-90), and the chunk evaluates the starting ",
               "model only")
      },
      "."))
    if (!is.null(stepwise_call)) {
      note("SELECTION", paste0(
        "not run under MAXITER=0: PROC HAZARD steps through the screen ",
        "without optimising (hazrd2.c:65-90); the emitted chunk evaluates ",
        "the starting model only"))
      stepwise_call <- NULL
      screen_check_call <- NULL
    }
    args$fit <- FALSE
    args$control <- NULL
    spec_call <- as.call(c(head, args))
    # PROC HAZARD refuses a job with more free parameters than events before
    # it evaluates anything (hazrd2.c:68-69, HAZ2TRM 7104; the count is C1 +
    # C3 over the rows, setobs.c:18). An evaluation always returns a number,
    # so without this the chunk reported a log-likelihood for a job PROC
    # HAZARD stops (r-reviewer pass 2 on #496). No events is one case of it.
    # Counted over the rows readobs() keeps, as tally is (Copilot on #539):
    # it drops a row whose TIME is missing or not positive (readt.c:9-15),
    # whose count is missing or negative (readc1.c:11-16; for
    # C3, readc3.c:10-15), or whose phase variable is missing
    # (readobs.c:128-134).
    counts <- Filter(Negate(is.null), list(
      if (!is.null(statements$EVENT)) as.name(statements$EVENT),
      if (!is.null(statements$ICENSOR)) as.name(statements$ICENSOR[[1L]])))
    keep <- c(list(bquote(!is.na(.(args$time)) & .(args$time) > 0)),
              lapply(counts, function(v) bquote(!is.na(.(v)) & .(v) >= 0)),
              lapply(phase_vars, function(v) bquote(!is.na(.(as.name(v))))))
    keep_expr <- Reduce(function(x, y) call("&", x, y), keep)
    events_expr <- Reduce(function(x, y) call("+", x, y), lapply(counts,
      function(v) bquote(sum(.(v)[.(keep_expr)]))))
    events_call <- call("with", args$data, events_expr)
    n_free <- parms$n_free
    guard <- bquote(if (.(events_call) < .(n_free)) {
      stop("PROC HAZARD stops this job before evaluating it: it has ",
           .(n_free), " free parameters and only ", .(events_call),
           " events (hazrd2.c:68-69, termination 7104). There is no ",
           "MAXITER=0 evaluation to report.", call. = FALSE)
    })
    if (isFALSE(ctl$conserve)) {
      code_body <- bquote(local({
        .(guard)
        hzr_evaluate(.(spec_call), theta = .(args$theta))
      }))
    } else {
      # For an ICENSOR job the spec carries objective = "sas" (#543), whose
      # argmax along the shift reproduces the MUE PROC HAZARD prints
      # (setcoe_obs_loop.c:114 counts C1 + C3 against the cumulative
      # hazard); the default interval likelihood peaks elsewhere
      # (r-reviewer on #496).
      # The search re-centres its bracket until the peak is inside it: a
      # start MU far from the CoE value (1e-15, or a late phase on a long
      # time scale) needs a shift past any fixed bracket, and optimize()
      # returns an edge without saying so. Ten re-centrings reach a factor
      # of about exp(300); a peak still at an edge stops the chunk rather
      # than report the likelihood there.
      code_body <- bquote(local({
        .(guard)
        .spec <- .(spec_call)
        .theta <- .(args$theta)
        .log_mu <- .(parms$log_mu_mask)
        .ll <- function(s) {
          suppressWarnings(hzr_evaluate(.spec, theta = .theta + s * .log_mu))$logLik
        }
        .shift <- 0
        for (.i in 1:10) {
          .at <- stats::optimize(.ll, .shift + c(-30, 30), maximum = TRUE,
                                 tol = 1e-10)$maximum
          .edge <- abs(.at - .shift) > 29.9
          .shift <- .at
          if (!.edge) break
        }
        if (!is.finite(.ll(.shift))) {
          stop("hzr_evaluate() cannot evaluate this model's likelihood at ",
               "the scaling of MU its search reached, so there is no ",
               "MAXITER=0 evaluation to report.", call. = FALSE)
        }
        if (.edge) {
          stop("No Conservation of Events scaling of MU was found within a ",
               "factor of exp(300) of the starting values, so there is no ",
               "MAXITER=0 evaluation to report.", call. = FALSE)
        }
        hzr_evaluate(.spec, theta = .theta + .shift * .log_mu)
      }))
    }
  }

  # John's 2026-09-22 decision, as amended at 19:51: a refusal warns and
  # emits the fit. That holds for the refusals that reach this return, not
  # for every refusal. Six paths above return a stop() in place of the fit:
  # a phase statement PROC HAZARD refuses at parse (#340), or a TIME or EVENT
  # with no operand and nothing else to supply it (#431), where a later `(`
  # changes only the stop's wording (#461); LCENSOR with ICENSOR (#155); no
  # phase selected (modterm.c ERROR 1001); a PARMS statement that builds no
  # phase this translation can use; no DATA= (#311); and a SELECTION
  # construct hzr_stepwise() cannot run (FAST, MAXVARS, RESTRICT, a negative
  # MAXSTEPS, a per-variable MOVE= or ORDER=, a cross-phase /I). A job with
  # no EVENT or ICENSOR statement at all never gets that far:
  # .hzr_censor_spec() raises during translation. Of the refusals that do
  # reach here, three (SETG3910, SETG3920, SETG3930) still halt, because SAS refuses them for a shape value that is
  # out of range (setg3.c:269-284) and hzr_phase() will not build a phase
  # from that same value. The warning is emitted in its own chunk ABOVE the
  # fit so that the real cause -- the SETG3 code and the operand -- is
  # RECORDED above it. Be clear about what that does and does not buy: under
  # Quarto the reader does NOT see it. knitr collects warnings INTO the
  # document, the chunk error then stops the render before any document is
  # written, and the console shows only hzr_phase()'s own
  # "gamma must be a positive scalar". The warning reaches a reader who runs
  # the chunks interactively, and the $untranslated row reaches anyone who
  # greps the job afterwards. An earlier version of this comment claimed the
  # cause was raised BEFORE the halt, which contradicted NEWS and was wrong
  # (#433 review). An earlier revision of the branch kept a stop() for those
  # three; it was replaced by this.

  list(call = code_body, status_call = status_call,
       stepwise_call = stepwise_call, screen_check_call = screen_check_call,
       outhaz = outhaz, untranslated = untr, tokens_seen = seen,
       tokens_mapped = mapped,
       sas_objective = identical(args$objective, "sas"),
       # Each is a reason PROC HAZARD would refuse this job, or would fit a
       # different model from the one emitted. They are carried out rather
       # than raised here: the point is that the RENDERED document warns, so
       # translate-sas.R emits them as a chunk immediately above the fit.
       refusal_warnings = refusal_warnings)
}

#' The SELECTION operands PROC HAZARD rejects with a syntax error.
#'
#' `stepwiseopt` (hazard_y.y:169-181) is `SLENTRY`, `SLSTAY`, `MOVE`,
#' `MAXSTEPS` or `MAXVARS` followed by `'=' NUMBER`, or a bare keyword. In the
#' STEP state a value is a NUMBER (hazard_l.l:33-38, :53) or unexpected text,
#' and a word the state has no rule for is unexpected text too
#' (hazard_l.l:176-179). Each of these sets the syntax-error flag, and the
#' binary refuses the job (measured on avc: `SLE=1E-3`, `BOGUS=1`, `BOGUS`,
#' `NOPRINTS=1`, `SLE 0.2` and `MOVE=ABC` exit SYNTAX; `NOPRINTS`, `SLE=0.2`
#' fit). A macro operand is SAS's to expand and carries no verdict.
#'
#' What is kept is what the screen still runs on: a value `as.numeric()`
#' reads (`1E-3`, as N3's warning says), and the keyword alone for a value
#' written on a bare keyword, so `BACKWARD=1` still screens backward. An
#' unknown option, a numeric option with no value and an unreadable value
#' are dropped, so .hzr_selection_spec() does not add a second row for them.
#' @param ops One statement's operands, spaces around `=` already closed.
#' @return `list(keep = <chr>, bad = <named chr>)`: `bad` maps each rejected
#'   construct to the reason, without the verdict.
#' @noRd
.hzr_selection_syntax <- function(ops) {
  numeric_opts <- c("SLENTRY", "SLSTAY", "MOVE", "MAXSTEPS", "MAXVARS")
  keep <- character(0)
  bad <- character(0)
  i <- 1L
  while (i <= length(ops)) {
    op <- ops[[i]]
    i <- i + 1L
    if (.hzr_sas_is_macro(op)) {
      keep <- c(keep, op)
      next
    }
    eqp <- .idx(op, "=")
    key <- if (eqp > 0L) substring(op, 1L, eqp - 1L) else op
    val <- if (eqp > 0L) substring(op, eqp + 1L) else ""
    # STEP context only. `SELECT` and the other statement keywords resolve
    # only in STMT context (hazard_l.l:112-114), so inside SELECTION they
    # are unexpected text: measured, `SELECTION SELECT;`, `SELECTION TIME;`
    # and `SELECTION SELECTION SLE=0.05;` exit SYNTAX.
    token <- .hzr_sas_token(key, "HAZARD", "STEP")
    if (is.na(token)) {
      bad[[op]] <- paste0("unknown SELECTION option, which PROC HAZARD's ",
                          "lexer reads as unexpected text (hazard_l.l:176-179)")
      next
    }
    if (token %in% numeric_opts) {
      if (!nzchar(val)) {
        # `SLE 0.2`: the number written without `=` is the same error.
        what <- op
        if (eqp == 0L && i <= length(ops) &&
              .hzr_sas_lexer_number(ops[[i]])) {
          what <- paste(op, ops[[i]])
          i <- i + 1L
        }
        bad[[what]] <- paste0("no value; PROC HAZARD has no form of this ",
                              "option without `= NUMBER` (hazard_y.y:169-173)")
        next
      }
      if (!.hzr_sas_lexer_number(val)) {
        bad[[op]] <- "not a number PROC HAZARD's lexer reads (hazard_l.l:33-38)"
        if (is.na(suppressWarnings(as.numeric(val)))) next
      }
      keep <- c(keep, op)
      next
    }
    if (eqp > 0L) {
      bad[[op]] <- "a value on an option that takes none (hazard_y.y:174-181)"
      keep <- c(keep, key)
      next
    }
    keep <- c(keep, op)
  }
  list(keep = keep, bad = bad)
}

#' Translate a SELECTION statement to hzr_stepwise() arguments.
#'
#' `hazard_y.y`'s `stepwisestmt` production is `STEPWISE stepwiseopts { setopt(33); }`,
#' and `stepwiseopts` can be empty: the *statement* is what turns stepwise
#' on, not any particular direction keyword. So a bare `SELECTION;` (or one
#' carrying only `SLENTRY`/`SLSTAY`) legitimately enables stepwise; the
#' direction keywords only refine it. `NOSTEPWISE`/`NOSW` (token `ONEWAY`,
#' option 34) does NOT turn it off: `stpwprc.c` leaves `sw = 1` and sets only
#' `nosw`, which caps each variable at one move, so it is a forward-only
#' screen (`direction = "forward"`). Every SELECTION statement therefore runs
#' a screen.
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
#' `SLENTRY` and `SLSTAY` are kept under `NOSTEPWISE` too: the forward
#' screen still applies its entry threshold.
#' @return `list(direction = <chr>, slentry = <dbl|NULL>, slstay = <dbl|NULL>,
#'   untranslated = <data.frame>, ...)`. A SELECTION statement always turns
#'   the screen on (hazard_y.y's `stepwisestmt`), so there is no on/off field.
#' @noRd
.hzr_selection_spec <- function(operands) {
  out <- list(direction = "both", slentry = NULL,
              slstay = NULL, max_steps = NULL, max_move = NULL,
              refuse = character(0), robust = character(0),
              untranslated = .hzr_untranslated_frame())
  saw_stepwise <- FALSE
  move_written <- FALSE
  saw_backward <- FALSE
  saw_oneway <- FALSE
  num_opt <- function(val_txt, key, what) {
    val <- suppressWarnings(as.numeric(val_txt))
    if (is.na(val)) {
      out$untranslated <<- rbind(out$untranslated, .hzr_untranslated_frame(
        NA_integer_, key, paste0("non-numeric value for ", what)
      ))
      return(NULL)
    }
    val
  }
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
      # Direction is resolved AFTER the loop, because stpwprc.c:16-22 is
      # order-independent: STEPWISE sets sw = 1 and BACKWARD then sets
      # sw = 0 and bw = 1 whatever order they appear in. A last-wins switch
      # ran a two-way screen for "SELECTION BACKWARD STEPWISE;".
      STEPWISE = saw_stepwise <- TRUE,
      BACKWARD = saw_backward <- TRUE,
      # NOSTEPWISE still screens: the SELECTION statement sets sw = 1
      # (hazard_y.y setopt(33), stpwprc.c) and NOSTEPWISE only sets nosw,
      # capping each variable at one move -- forward only (#342 review).
      ONEWAY   = saw_oneway <- TRUE,
      SLENTRY  = out$slentry <- num_opt(val_txt, key, "SLENTRY") %||% out$slentry,
      SLSTAY   = out$slstay <- num_opt(val_txt, key, "SLSTAY") %||% out$slstay,
      MAXSTEPS = {
        v <- num_opt(val_txt, key, "MAXSTEPS")
        # stpwprc.c:76-79 exits the job on a negative MAXSTEPS, so there is
        # no run to translate.
        if (!is.null(v) && v < 0) {
          out$refuse <- c(out$refuse, paste0("MAXSTEPS=", format(v)))
        } else if (!is.null(v)) {
          # stpwprc.c:82 casts to int, so a fractional MAXSTEPS truncates.
          out$max_steps <- trunc(v)
        }
      },
      MOVE     = {
        out$max_move <- num_opt(val_txt, key, "MOVE") %||% out$max_move
        move_written <- TRUE
      },
      # Printing only (H->nps / H->npq), so the fit and the screen are the
      # same with or without them: recorded, not refused.
      NOPRINTS = out$untranslated <- rbind(out$untranslated,
        .hzr_untranslated_frame(NA_integer_, key,
                                "suppresses a PROC HAZARD printout only")),
      NOPRINTQ = out$untranslated <- rbind(out$untranslated,
        .hzr_untranslated_frame(NA_integer_, key,
                                "suppresses a PROC HAZARD printout only")),
      # Refused: FAST is a different search (H->f) and MAXVARS caps the
      # selected set, neither of which hzr_stepwise() can express.
      FAST       = out$refuse <- c(out$refuse, "FAST"),
      MAXVARS    = out$refuse <- c(out$refuse, "MAXVARS"),
      # Recorded, not refused. ROBUST and SEMIROBUST choose PROC HAZARD's
      # OPTIMIZER for the stepwise step, nothing else: stpwprc.c:60-69 sets
      # swnewt/swhess, hazrd2.c:31-34 swaps them into H->newton/H->truhes
      # around the step and :87-88 restores them, and cmpmeth.c shows those
      # flags pick Newton vs quasi-Newton and a Hessian vs steepest-descent
      # start. They do not touch the variance, the converged estimates or
      # the Wald tests that drive selection. Both earlier rulings on ROBUST
      # (refuse, then translate "because it changes the variance") rested
      # on the false premise that they did (#160 review, pass 3).
      ROBUST     = out$robust <- c(out$robust, "ROBUST"),
      SEMIROBUST = out$robust <- c(out$robust, "SEMIROBUST"),
      {
        out$untranslated <- rbind(out$untranslated, .hzr_untranslated_frame(
          NA_integer_, key, "no hzr_stepwise() equivalent"
        ))
      }
    )
  }
  for (kw in out$robust) {
    out$untranslated <- rbind(out$untranslated, .hzr_untranslated_frame(
      NA_integer_, kw, paste0(
        "chooses PROC HAZARD's optimizer for the stepwise step (quasi-Newton, ",
        "started ", if (identical(kw, "ROBUST")) "by steepest descent" else
          "from the Hessian", "; stpwprc.c:60-69, hazrd2.c:31-34). ",
        "hzr_stepwise() uses its own optimizer. This changes the path to ",
        "convergence, and on a multimodal likelihood it can change the ",
        "optimum reached, and with it the selection")))
  }

  # BACKWARD wins over any other direction keyword, whatever the order
  # (stpwprc.c:16-22). NOSTEPWISE only caps moves, so it is forward.
  out$direction <- if (saw_backward) "backward" else
    if (saw_oneway) "forward" else "both"

  # PROC HAZARD's own defaults and its clamping, so the emitted call never
  # inherits hzr_stepwise()'s different defaults and never carries a value
  # PROC HAZARD would have thrown away (stpwprc.c:27-56, przconc.c:25).
  # An SLE outside (0, 1) is not an aggressive screen, it is a value SAS
  # replaces: slstay = 0 removes nothing and slentry = 5 enters everything.
  clamp <- function(val, key, lo_open, hi_open, default, why) {
    if (is.null(val)) return(default)
    if (val <= lo_open || val >= hi_open) {
      out$untranslated <<- rbind(out$untranslated, .hzr_untranslated_frame(
        NA_integer_, paste0(key, "=", format(val)), why))
      return(default)
    }
    val
  }
  out$slentry <- clamp(out$slentry, "SLENTRY", 0, 1, 0.3,
                       paste("PROC HAZARD ignores an SLENTRY outside (0, 1)",
                             "and uses 0.3 (stpwprc.c:27-35)"))
  sls_default <- if (identical(out$direction, "backward")) 0.05 else 0.2
  out$slstay <- clamp(out$slstay, "SLSTAY", 0, 1, sls_default,
                      paste("PROC HAZARD ignores an SLSTAY outside (0, 1)",
                            "and uses its default (stpwprc.c:36-44)"))
  if (is.null(out$max_move)) {
    out$max_move <- 1
  } else if (saw_oneway && out$max_move > 1) {
    out$untranslated <- rbind(out$untranslated, .hzr_untranslated_frame(
      NA_integer_, paste0("MOVE=", format(out$max_move)),
      paste("NOSTEPWISE caps every variable at one move, so PROC HAZARD",
            "resets MOVE to 1 (stpwprc.c:54-58, przconc.c:25)")))
    out$max_move <- 1
  }
  # Recorded for EVERY screen, not only a job that wrote MOVE=: the two
  # counters differ whatever the value, and PROC HAZARD's own default of 1
  # does not map onto hzr_stepwise()'s 4 either.
  out$untranslated <- rbind(out$untranslated, .hzr_untranslated_frame(
    NA_integer_,
    if (move_written) paste0("MOVE=", format(out$max_move)) else
      "MOVE (PROC HAZARD default 1)",
    paste(
      "PROC HAZARD's MOVE limit counts a variable's deletions, separately",
      "for each phase (hazrd4.c:361-362, setstat.c:15): at the default of 1",
      "a variable removed from a phase can never re-enter it",
      "(swvari.c:159-160). hzr_stepwise()'s max_move counts entries and",
      "exits together across every phase, so the two cannot be mapped onto",
      "each other. The screen runs with hzr_stepwise()'s own oscillation",
      "guard, so it can re-enter variables PROC HAZARD would keep out")))
  out
}

# ---------------------------------------------------------------------------
# PROC HAZPRED prediction grids: the job's own DATA steps, translated
# ---------------------------------------------------------------------------
#
# HAZPRED predicts at every row of its DATA= dataset, reading the time from
# the variable its TIME statement names (hazpred/timeprc.c:16-25) and every
# covariate of the model from the column of the same name
# (hazpred/hazpred.c:153-173). So the grid is whatever the job's DATA steps
# left in that dataset by the time the procedure ran: the covariates a
# `SET DESIGN` brings in, the rows a later `DATA PREDICT; SET PREDICT
# DIGITAL;` appends, and the time a `YEARS = MONTHS/12` derives. Reading
# only the first DO loop of the first `DATA <name>;` step dropped all three,
# and the prediction ran at covariates of zero and the wrong times with
# nothing recorded (#494).
#
# The steps are translated statement by statement into base R that builds
# the same data frame. What cannot be translated reliably is split two ways.
# A statement whose effect on the grid is known but whose value is not (a
# conditional assignment, a function this does not carry) is recorded, and
# the variable it sets is NA in the grid, so predictions that use it are NA
# rather than evaluated at a value SAS did not use. A statement that decides
# which rows the grid has (an input file, a MERGE, a loop whose bounds are
# data) refuses the whole grid, because no grid can be emitted that is not
# a guess.

#' Split normalised SAS source into statements, with their offsets.
#'
#' A `;` inside a quoted string does not end a statement.
#' @return A data frame with `start` (offset of the statement's first
#'   non-blank character in `txt`), `end` and `text` (trimmed).
#' @noRd
.hzr_sas_statements <- function(txt) {
  m <- gregexpr("(?:'[^']*'|\"[^\"]*\"|[^;'\"])+", txt, perl = TRUE)[[1L]]
  if (m[1L] == -1L) {
    return(data.frame(start = integer(0), end = integer(0), text = character(0),
                      stringsAsFactors = FALSE))
  }
  raw <- regmatches(txt, list(m))[[1L]]
  lead <- nchar(raw) - nchar(sub("^ +", "", raw))
  data.frame(start = as.integer(m) + lead,
             end = as.integer(m) + attr(m, "match.length") - 1L,
             text = trimws(raw), stringsAsFactors = FALSE)
}

#' Read one `KEY=value` dataset name out of a statement.
#' @noRd
.hzr_sas_opt_name <- function(stmt, key) {
  m <- regmatches(stmt, regexec(
    paste0("(^|[^A-Z0-9_])", key, " *= *([A-Z_][A-Z0-9_.]*)"), stmt))[[1L]]
  if (length(m) >= 3L) m[[3L]] else NULL
}

#' A dataset name as SAS resolves it: a one-level name lives in WORK, so
#' `WORK.PREDICT` is `PREDICT`. Any other libref names another dataset and
#' is kept.
#' @noRd
.hzr_sas_ds <- function(x) sub("^WORK[.]", "", x)

#' Every dataset a non-DATA statement writes.
#'
#' `OUT=` and its relatives (`OUTEST=`, `OUTSTAT=`, ...), `BASE=` (PROC
#' APPEND), and in PROC SQL `CREATE TABLE`, `CREATE VIEW`, `INSERT INTO`,
#' `DELETE FROM`, `UPDATE` and `ALTER TABLE`. In a PROC DATASETS
#' step every name counts, because `CHANGE`, `DELETE` and `MODIFY` all take
#' dataset names.
#' @noRd
.hzr_sas_written_names <- function(stmt, datasets = FALSE) {
  if (datasets) {
    nms <- regmatches(stmt, gregexpr("[A-Z_][A-Z0-9_.]*", stmt))[[1L]]
    return(setdiff(nms, c("PROC", "DATASETS", "CHANGE", "DELETE", "MODIFY",
                          "EXCHANGE", "AGE", "SAVE", "APPEND", "BASE", "DATA",
                          "LIBRARY", "LIB", "NOLIST", "KILL", "MEMTYPE")))
  }
  m <- regmatches(stmt, gregexpr(
    paste0("(^|[^A-Z0-9_])(OUT[A-Z]*|BASE) *= *[A-Z_][A-Z0-9_.]*|",
           "(CREATE +(TABLE|VIEW)|INSERT +INTO|DELETE +FROM|ALTER +TABLE|^UPDATE) +[A-Z_][A-Z0-9_.]*"),
    stmt))[[1L]]
  unique(sub("^.*[= ]", "", m))
}

#' The macro a statement calls, or `NULL`.
#'
#' `%INCLUDE` (and `%INC`) counts as a call: it runs code this cannot see.
#' The macro language's own statements and functions do not.
#' @noRd
.hzr_sas_macro_call <- function(stmt) {
  m <- regmatches(stmt, regexec("^%([A-Z_][A-Z0-9_]*)", stmt))[[1L]]
  if (!length(m)) return(NULL)
  nm <- m[[2L]]
  if (nm == "INC") nm <- "INCLUDE"
  lang <- c("LET", "PUT", "IF", "THEN", "ELSE", "DO", "END", "TO", "BY",
            "WHILE", "UNTIL", "GLOBAL", "LOCAL", "MACRO", "MEND", "GOTO",
            "RETURN", "ABORT", "SYMDEL", "SYSCALL", "SYSEXEC", "SYSLPUT",
            "SYSRPUT", "SYSMACDELETE", "WINDOW", "DISPLAY", "INPUT", "COPY",
            "STR", "NRSTR", "QUOTE", "NRQUOTE", "BQUOTE", "NRBQUOTE", "SUPERQ",
            "UNQUOTE", "EVAL", "SYSEVALF", "SYSFUNC", "QSYSFUNC")
  if (nm %in% lang) NULL else nm
}

#' The `%MACRO ... %MEND` definitions among a job's statements.
#' @return `list(scope, bodies)`: `scope[i]` numbers the outermost
#'   definition statement `i` belongs to, `%MACRO` and `%MEND` included, and
#'   is 0 outside every definition; `bodies` holds each macro's statements,
#'   named by the macro.
#' @noRd
.hzr_sas_macro_defs <- function(stmts) {
  scope <- integer(length(stmts))
  bodies <- list()
  depth <- 0L
  n_def <- 0L
  open <- character(0)
  for (i in seq_along(stmts)) {
    t <- stmts[[i]]
    if (grepl("^%MACRO ", t)) {
      if (depth == 0L) n_def <- n_def + 1L
      depth <- depth + 1L
      open <- c(open, sub("^%MACRO +([A-Z_][A-Z0-9_]*).*$", "\\1", t))
      bodies[[open[[depth]]]] <- character(0)
      scope[[i]] <- n_def
      next
    }
    if (depth == 0L) next
    scope[[i]] <- n_def
    if (grepl("^%MEND( |$)", t)) {
      open <- open[-depth]
      depth <- depth - 1L
      next
    }
    for (nm in open) bodies[[nm]] <- c(bodies[[nm]], t)
  }
  list(scope = scope, bodies = bodies)
}

#' The definition the statement holding offset `pos` stands in, or 0.
#' @noRd
.hzr_sas_scope_at <- function(st, defs, pos) {
  i <- which(st$start <= pos)
  if (length(i)) defs$scope[[max(i)]] else 0L
}

#' The datasets a macro body names as written, and whether it may write
#' others: through a macro variable (`DATA &DS;`), or by calling a macro or
#' an `%INCLUDE` of its own.
#' @noRd
.hzr_sas_body_writes <- function(body) {
  nms <- character(0)
  unknown <- FALSE
  for (t in body) {
    if (grepl("^DATA( |$)", t) && !grepl("^DATA *=", t)) {
      rest <- strsplit(gsub("[(][^)]*[)]", " ", trimws(sub("^DATA", "", t))), " +")[[1L]]
      nms <- c(nms, rest)
    } else {
      nms <- c(nms, .hzr_sas_written_names(t))
      if (!is.null(.hzr_sas_macro_call(t))) unknown <- TRUE
    }
    if (grepl("&", t) && grepl("^DATA |(OUT[A-Z]*|BASE) *= *&", t)) unknown <- TRUE
  }
  nms <- nms[grepl("^[A-Z_][A-Z0-9_.]*$", nms) & nms != "_NULL_"]
  list(names = unique(nms), unknown = unknown)
}

#' Every point in a job where a dataset is (re)defined, in file order.
#'
#' Five kinds. `data`: a `DATA <name>;` step, carrying its statements.
#' `hazpred`: a `PROC HAZPRED ... OUT=`, whose output has one row per row of
#' its `DATA=` dataset (`hazpred/obsloop.c:17-26`, `:67`). `sort`: a
#' `PROC SORT ... ; BY ...;`, which reorders its input. `opaque`: anything
#' else that writes a dataset (another procedure's or a macro's `OUT=`, a
#' multi-dataset or optioned `DATA` statement), which this cannot read and so
#' refuses if a grid depends on it.
#'
#' A `%MACRO ... %MEND` body is a definition: SAS runs its statements where
#' the macro is called, not where they stand. Each event carries `scope`,
#' the definition it stands in (0 outside every one), and
#' `.hzr_parse_grid()` reads a definition's events only for a PROC HAZPRED
#' in the same definition. A call of a macro defined in the file is an
#' `opaque` event for each dataset its body names as written. The fifth
#' kind, `call`, has no dataset: an `%INCLUDE`, a call of a macro the file does not define, or of
#' one whose body writes through a macro variable or runs code of its own.
#' It may rewrite any dataset and this cannot say which, so
#' `.hzr_parse_grid()` records it rather than refusing.
#' @noRd
.hzr_sas_dataset_events <- function(txt, blocks) {
  st <- .hzr_sas_statements(txt)
  defs <- .hzr_sas_macro_defs(st$text)
  scope_at <- function(pos) .hzr_sas_scope_at(st, defs, pos)
  ev <- list()
  add <- function(pos, name, kind, ...) {
    extra <- list(...)
    if (!is.null(extra$from)) extra$from <- .hzr_sas_ds(extra$from)
    ev[[length(ev) + 1L]] <<- c(list(pos = pos, name = .hzr_sas_ds(name), kind = kind,
                                     scope = scope_at(pos)), extra)
  }
  for (b in blocks) {
    head <- sub(";.*$", "", b$text)
    if (identical(b$proc, "HAZPRED")) {
      out <- .hzr_sas_opt_name(head, "OUT")
      from <- .hzr_sas_opt_name(head, "DATA")
      if (!is.null(out)) {
        if (is.null(from)) {
          add(b$start, out, "opaque", what = "PROC HAZPRED with no DATA=")
        } else {
          add(b$start, out, "hazpred", from = from)
        }
      }
    } else {
      out <- .hzr_sas_opt_name(b$text, "OUT")
      if (identical(b$proc, "REPEAT") && is.null(out)) out <- "EVENTS"
      if (!is.null(out)) add(b$start, out, "opaque", what = paste0("%", b$proc))
    }
  }

  in_block <- function(s, e) {
    any(vapply(blocks, function(b) e >= b$start && s <= b$end, logical(1L)))
  }
  cur <- NULL
  proc <- NULL
  last_name <- function(p) {
    before <- Filter(function(x) x$pos < p && !is.na(x$name), ev)
    if (!length(before)) return(NULL)
    before[[which.max(vapply(before, function(x) x$pos, numeric(1L)))]]$name
  }
  close_all <- function() {
    if (!is.null(cur) && !is.na(cur$name)) {
      add(cur$pos, cur$name, "data", stmts = cur$stmts)
    }
    cur <<- NULL
    if (!is.null(proc) && identical(proc$kind, "SORT")) {
      from <- if (is.null(proc$data)) last_name(proc$pos) else proc$data
      to <- if (is.null(proc$out)) from else proc$out
      if (!is.null(to)) {
        if (is.null(from) || proc$bad || !length(proc$by)) {
          add(proc$pos, to, "opaque", what = "PROC SORT")
        } else {
          add(proc$pos, to, "sort", from = from, by = proc$by)
        }
      }
    }
    proc <<- NULL
  }

  for (i in seq_len(nrow(st))) {
    t <- st$text[[i]]
    p <- st$start[[i]]
    if (!nzchar(t)) next
    if (in_block(p, st$end[[i]])) {
      close_all()
      next
    }
    if (grepl("^%(MACRO |MEND( |$))", t)) {
      close_all()
      next
    }
    mac <- .hzr_sas_macro_call(t)
    if (!is.null(mac)) {
      body <- defs$bodies[[mac]]
      if (is.null(body)) {
        add(p, NA_character_, "call", what = t)
      } else {
        w <- .hzr_sas_body_writes(body)
        for (nm in w$names) add(p, nm, "opaque", what = paste0("%", mac))
        if (w$unknown) add(p, NA_character_, "call", what = t)
      }
    }
    if (grepl("^DATA( |$)", t) && !grepl("^DATA *=", t)) {
      close_all()
      rest <- .hzr_sas_ds(trimws(sub("^DATA", "", t)))
      if (grepl("^[A-Z_][A-Z0-9_]*$", rest) && !identical(rest, "_NULL_")) {
        cur <- list(pos = p, name = rest, stmts = character(0))
      } else {
        # `DATA A B;`, `DATA A(KEEP=...)`, `DATA LIB.A;`: every name it writes
        # is recorded as unreadable, and the step's statements are swallowed.
        nms <- strsplit(gsub("[(][^)]*[)]", " ", rest), " +")[[1L]]
        for (nm in nms[grepl("^[A-Z_][A-Z0-9_.]*$", nms) & nms != "_NULL_"]) {
          add(p, nm, "opaque", what = paste0("DATA ", rest))
        }
        cur <- list(pos = p, name = NA_character_, stmts = character(0))
      }
      next
    }
    if (grepl("^PROC ", t)) {
      close_all()
      if (grepl("^PROC SORT( |$)", t)) {
        opts <- trimws(sub("^PROC SORT", "", t))
        opts <- gsub("(DATA|OUT) *= *[A-Z_][A-Z0-9_.]*", "", opts)
        proc <- list(kind = "SORT", pos = p, data = .hzr_sas_opt_name(t, "DATA"),
                     out = .hzr_sas_opt_name(t, "OUT"), by = character(0),
                     bad = nzchar(trimws(opts)))
      } else {
        proc <- list(kind = "OTHER", pos = p, what = sub("^(PROC [A-Z0-9_]+).*$", "\\1", t))
        for (out in .hzr_sas_written_names(t)) add(p, out, "opaque", what = proc$what)
      }
      next
    }
    if (grepl("^(RUN|QUIT)$", t)) {
      close_all()
      next
    }
    # A macro call or a procedure statement (OUTPUT OUT=, CREATE TABLE)
    # writes a dataset this cannot read. Inside a DATA step only a macro call
    # can: `OUT=X` there is an assignment.
    if (is.null(cur) || startsWith(t, "%")) {
      in_datasets <- !is.null(proc) && identical(proc$what, "PROC DATASETS")
      what <- if (is.null(proc) || startsWith(t, "%")) sub("[ (].*$", "", t) else proc$what
      for (out in .hzr_sas_written_names(t, in_datasets && !startsWith(t, "%"))) {
        add(p, out, "opaque", what = what)
      }
    }
    if (!is.null(cur)) {
      cur$stmts <- c(cur$stmts, t)
      next
    }
    if (!is.null(proc) && identical(proc$kind, "SORT") && grepl("^BY ", t)) {
      by <- strsplit(trimws(sub("^BY", "", t)), " +")[[1L]]
      if (!all(grepl("^[A-Z_][A-Z0-9_]*$", by)) || "DESCENDING" %in% by) proc$bad <- TRUE
      proc$by <- by
    }
  }
  close_all()
  ev[order(vapply(ev, function(x) x$pos, numeric(1L)))]
}

#' The last definition of dataset `name` that starts before offset `before`.
#' @return An index into `events`, or `NA_integer_`.
#' @noRd
.hzr_sas_resolve <- function(events, name, before) {
  hit <- which(vapply(events, function(e) identical(e$name, name) && e$pos < before,
                      logical(1L)))
  if (length(hit)) hit[[length(hit)]] else NA_integer_
}

#' Translate one SAS DATA-step arithmetic expression into an R call.
#'
#' Numbers, variable names, `+ - * / **`, parentheses, quoted strings, a lone
#' `.` (a missing value) and the functions `LOG` and `EXP`: the whole
#' vocabulary of the corpus's grid arithmetic. SAS and R agree on the
#' precedence of these, including `-2**2` (`-4` in both). Anything else (a
#' comparison, `AND`, another function) is declined rather than guessed at.
#' Built from a token whitelist, so the returned call can only reference the
#' variables listed in `refs`, `log`, `exp` and arithmetic.
#' @return `list(ok = TRUE, call, refs)` or `list(ok = FALSE, why)`.
#' @noRd
.hzr_sas_expr <- function(text) {
  fail <- function(why) list(ok = FALSE, why = why)
  s <- trimws(text)
  if (!nzchar(s)) return(fail("an empty expression"))
  fns <- c(LOG = "log", EXP = "exp")
  out <- character(0)
  refs <- character(0)
  take <- function(pattern) {
    m <- regmatches(s, regexpr(pattern, s, perl = TRUE))
    if (length(m)) m else NULL
  }
  while (nzchar(s)) {
    if (startsWith(s, " ")) {
      s <- sub("^ +", "", s)
      next
    }
    tok <- take("^(?:[0-9]+[.]?[0-9]*|[.][0-9]+)(?:E[+-]?[0-9]+)?")
    if (!is.null(tok)) {
      out <- c(out, tok)
      s <- substring(s, nchar(tok) + 1L)
      next
    }
    tok <- take("^(?:'[^']*'|\"[^\"]*\")")
    if (!is.null(tok)) {
      out <- c(out, encodeString(substr(tok, 2L, nchar(tok) - 1L), quote = "\""))
      s <- substring(s, nchar(tok) + 1L)
      next
    }
    tok <- take("^[A-Z_][A-Z0-9_]*")
    if (!is.null(tok)) {
      s <- substring(s, nchar(tok) + 1L)
      if (grepl("^ *[(]", s)) {
        if (!tok %in% names(fns)) {
          return(fail(paste0("calls ", tok, "(), which this translation does not carry")))
        }
        out <- c(out, fns[[tok]])
      } else {
        out <- c(out, paste0("`", tok, "`"))
        refs <- c(refs, tok)
      }
      next
    }
    if (startsWith(s, "**")) {
      out <- c(out, "^")
      s <- substring(s, 3L)
      next
    }
    ch <- substr(s, 1L, 1L)
    if (ch %in% c("+", "-", "*", "/", "(", ")", ",")) {
      out <- c(out, ch)
    } else if (ch == ".") {
      out <- c(out, "NA_real_")
    } else {
      return(fail(paste0("uses `", ch, "`, which this translation does not carry")))
    }
    s <- substring(s, 2L)
  }
  cl <- tryCatch(str2lang(paste(out, collapse = " ")), error = function(e) NULL)
  if (is.null(cl)) return(fail("does not read as arithmetic"))
  list(ok = TRUE, call = cl, refs = unique(refs))
}

#' The variables a DATA-step statement (or right-hand side) reads.
#'
#' Every name that is not quoted, not a function called, and not one of the
#' words SAS's own syntax uses around it.
#' @noRd
.hzr_sas_names_read <- function(text) {
  s <- gsub("'[^']*'|\"[^\"]*\"", " ", text)
  nms <- regmatches(s, gregexpr("(?<![A-Z0-9_.])[A-Z_][A-Z0-9_]*(?![A-Z0-9_]| *[(])",
                                s, perl = TRUE))[[1L]]
  setdiff(unique(nms), c("AND", "OR", "NOT", "EQ", "NE", "GT", "LT", "GE", "LE",
                         "IN", "IF", "THEN", "ELSE", "DO", "DELETE", "OUTPUT"))
}

#' Evaluate a translated expression over already-folded constants.
#' @return A length-one value, or `NULL` when it does not fold.
#' @noRd
.hzr_sas_fold <- function(call, const) {
  # Safe by construction: `call` comes only from .hzr_sas_expr(), whose token
  # whitelist admits numbers, strings, backquoted variable names, arithmetic,
  # log() and exp(). Every name was checked against `const` by the caller,
  # and base is the only other scope, so nothing else is reachable.
  val <- tryCatch(eval(call, envir = const, enclos = baseenv()),
                  error = function(e) NULL, warning = function(w) NULL)
  if (length(val) != 1L || !(is.numeric(val) || is.character(val))) return(NULL)
  if (is.numeric(val) && !is.finite(val)) return(NULL)
  val
}

#' Split `text` at commas that are outside parentheses and quotes.
#' @noRd
.hzr_sas_split_commas <- function(text) {
  chars <- strsplit(text, "", fixed = TRUE)[[1L]]
  depth <- 0L
  quote <- ""
  cut <- integer(0)
  for (i in seq_along(chars)) {
    ch <- chars[[i]]
    if (nzchar(quote)) {
      if (ch == quote) quote <- ""
    } else if (ch %in% c("'", "\"")) {
      quote <- ch
    } else if (ch == "(") {
      depth <- depth + 1L
    } else if (ch == ")") {
      depth <- depth - 1L
    } else if (ch == "," && depth == 0L) {
      cut <- c(cut, i)
    }
  }
  starts <- c(1L, cut + 1L)
  ends <- c(cut - 1L, length(chars))
  trimws(substring(text, starts, ends))
}

#' Translate the value list of `DO var = <spec>;` into an R vector call.
#'
#' SAS's list is comma separated; each item is a value or `a TO b [BY c]`,
#' and the loop takes every item in turn (so `DO T = a TO b BY c, b;` runs the
#' range and then takes `b` once more). Every bound must fold to a constant
#' from assignments earlier in the same step, as every corpus grid's does,
#' so the loop has the same values for every row it expands. A bound read
#' from data would give each row its own loop, which is refused.
#' @return `list(call, refs)` or `list(refuse = <why>)`.
#' @noRd
.hzr_sas_do_values <- function(spec, const) {
  refuse <- function(why) list(refuse = why)
  one <- function(txt) {
    e <- .hzr_sas_expr(txt)
    if (!e$ok) return(list(why = paste0("`", txt, "` ", e$why)))
    if (!all(e$refs %in% names(const))) {
      return(list(why = paste0("`", txt, "` is not a constant set earlier in the step")))
    }
    val <- .hzr_sas_fold(e$call, const)
    if (!is.numeric(val)) return(list(why = paste0("`", txt, "` does not evaluate to a number")))
    list(call = e$call, refs = e$refs, val = val)
  }
  elems <- list()
  refs <- character(0)
  for (part in .hzr_sas_split_commas(spec)) {
    rng <- regmatches(part, regexec("^(.+?) +TO +(.+?)(?: +BY +(.+))?$", part, perl = TRUE))[[1L]]
    if (length(rng)) {
      lo <- one(rng[[2L]])
      hi <- one(rng[[3L]])
      by <- if (nzchar(rng[[4L]])) one(rng[[4L]]) else list(call = 1, refs = character(0), val = 1)
      for (x in list(lo, hi, by)) if (!is.null(x$why)) return(refuse(x$why))
      # SAS runs a TO b BY c while the value has not passed b, counting down
      # for a negative c, as seq() does. A range SAS runs zero times (or
      # forever) would leave seq() to fail at render; refuse it here.
      steps <- (hi$val - lo$val) / by$val
      if (by$val == 0 || steps < 0 || steps > 1e6) {
        return(refuse(paste0("`", part, "` is a range SAS runs no times, or without end")))
      }
      elems[[length(elems) + 1L]] <- bquote(seq(.(lo$call), .(hi$call), by = .(by$call)))
      refs <- c(refs, lo$refs, hi$refs, by$refs)
    } else {
      x <- one(part)
      if (!is.null(x$why)) return(refuse(x$why))
      elems[[length(elems) + 1L]] <- x$call
      refs <- c(refs, x$refs)
    }
  }
  call <- if (length(elems) == 1L) elems[[1L]] else as.call(c(quote(c), elems))
  list(call = call, refs = unique(refs))
}

#' Classify one DATA-step statement.
#' @noRd
.hzr_sas_stmt_kind <- function(t) {
  if (grepl("^[A-Z_][A-Z0-9_]* *=", t)) return("assign")
  if (identical(t, "DO")) return("group")
  if (grepl("^DO +[A-Z_][A-Z0-9_]* *=", t)) return("do")
  if (grepl("^DO( |[(]|$)", t)) return("doother")
  if (grepl("^(IF|ELSE)( |$)", t)) {
    return(if (grepl("(^ELSE| THEN) DO$", t)) "ifdo" else "if")
  }
  if (identical(t, "OUTPUT")) return("output")
  if (identical(t, "DELETE")) return("delete")
  if (grepl("^DROP( |$)", t)) return("drop")
  if (grepl("^KEEP( |$)", t)) return("keep")
  if (grepl("^SET( |$)", t)) return("set")
  ignorable <- paste0("^(LIBNAME|FILENAME|TITLE[0-9]*|FOOTNOTE[0-9]*|OPTIONS?|",
                      "LENGTH|LABEL|FORMAT|INFORMAT|ATTRIB|PUT|FILE)( |$)")
  if (grepl(ignorable, t)) return("ignore")
  "other"
}

#' Parse a step's statements into a tree, matching each DO to its END.
#' @return `list(items, error)`; each item is `list(kind, text[, body])`.
#' @noRd
.hzr_sas_step_items <- function(stmts) {
  i <- 1L
  err <- NULL
  walk <- function(in_block) {
    items <- list()
    while (i <= length(stmts)) {
      t <- stmts[[i]]
      i <<- i + 1L
      if (identical(t, "END")) {
        if (in_block) return(items)
        err <<- "an END with no DO"
        return(items)
      }
      item <- list(kind = .hzr_sas_stmt_kind(t), text = t)
      if (item$kind %in% c("do", "group", "doother", "ifdo")) item$body <- walk(TRUE)
      items[[length(items) + 1L]] <- item
    }
    if (in_block && is.null(err)) err <<- "a DO with no END"
    items
  }
  items <- walk(FALSE)
  # An unconditional `DO; ... END;` group is only brackets: splice it in.
  flat <- function(items) {
    out <- list()
    for (it in items) {
      if (!is.null(it$body)) it$body <- flat(it$body)
      if (identical(it$kind, "group")) out <- c(out, it$body) else out[[length(out) + 1L]] <- it
    }
    out
  }
  list(items = flat(items), error = err)
}

#' The kinds, and the assignment targets, anywhere in an item tree.
#' @noRd
.hzr_sas_tree_kinds <- function(items) {
  unlist(lapply(items, function(it) c(it$kind, .hzr_sas_tree_kinds(it$body))))
}

#' What a conditional statement (or block) would change, if it ran.
#' @return `list(targets, deletes, bad)`: variables it assigns, whether it can
#'   delete rows, and the text of any statement inside it that decides rows
#'   in a way this cannot mark (an OUTPUT, a loop, SET).
#' @noRd
.hzr_sas_cond_effects <- function(item) {
  targets <- character(0)
  deletes <- FALSE
  bad <- NULL
  reads <- character(0)
  inner <- function(txt) {
    k <- .hzr_sas_stmt_kind(txt)
    if (identical(k, "if")) {
      visit(list(kind = "if", text = txt)) # ELSE IF ... THEN ...
    } else if (identical(k, "assign")) {
      targets <<- c(targets, sub(" *=.*$", "", txt))
    } else if (identical(k, "delete")) {
      deletes <<- TRUE
    } else if (!identical(k, "ignore")) {
      bad <<- c(bad, txt)
    }
  }
  visit <- function(it) {
    reads <<- c(reads, .hzr_sas_names_read(sub("^(IF|ELSE)( |$)", "", it$text)))
    if (identical(it$kind, "if")) {
      then <- regmatches(it$text, regexec("^(?:IF .*? THEN|ELSE) +(.+)$", it$text, perl = TRUE))[[1L]]
      if (length(then)) inner(then[[2L]]) else deletes <<- TRUE # a subsetting IF
    } else if (identical(it$kind, "ifdo")) {
      for (b in it$body) visit(b)
    } else {
      inner(it$text)
    }
  }
  visit(item)
  list(targets = unique(targets), deletes = deletes, bad = bad,
       reads = unique(setdiff(reads, targets)))
}

#' Translate one `DATA <name>;` step into R statements.
#'
#' @param ev The step's event from `.hzr_sas_dataset_events()`.
#' @param input A function of a dataset name returning `list(cols = ...)` for
#'   the definition this step's `SET` reads, or `list(refuse = <why>)`.
#' @return `list(code, cols, untr, warn)` or `list(refuse = <why>)`.
#' @noRd
.hzr_sas_translate_step <- function(ev, input) {
  refuse <- function(why) list(refuse = paste0("its DATA ", ev$name, " step has ", why))
  parsed <- .hzr_sas_step_items(ev$stmts)
  if (!is.null(parsed$error)) return(refuse(parsed$error))
  items <- Filter(function(it) !identical(it$kind, "ignore"), parsed$items)

  W <- as.name(ev$name)
  code <- list()
  untr <- .hzr_untranslated_frame()
  warn <- character(0)
  emit <- function(cl) code[[length(code) + 1L]] <<- cl
  set_col <- function(v, value) emit(call("<-", call("$", W, as.name(v)), value))
  u1 <- function(stmt, why) {
    untr <<- rbind(untr, .hzr_untranslated_frame(NA_integer_, stmt, why))
    warn <<- c(warn, paste0("`", stmt, "`: ", why))
  }

  # --- SET: the step's input rows ------------------------------------------
  set_names <- character(0)
  if (length(items) && identical(items[[1L]]$kind, "set")) {
    spec <- trimws(sub("^SET", "", items[[1L]]$text))
    set_names <- .hzr_sas_ds(strsplit(spec, " ", fixed = TRUE)[[1L]])
    items <- items[-1L]
  }

  kinds <- .hzr_sas_tree_kinds(items)
  stop_kinds <- c("set", "other", "doother")
  if (any(kinds %in% stop_kinds)) {
    hit <- NULL
    find <- function(its) {
      for (it in its) {
        if (is.null(hit) && it$kind %in% stop_kinds) hit <<- it$text
        find(it$body)
      }
    }
    find(items)
    return(refuse(paste0("`", hit, "`, which decides the grid's rows or values in a way ",
                         "this translation does not read")))
  }
  top <- vapply(items, function(it) it$kind, character(1L))
  loops <- which(top == "do")
  if (any(top == "delete")) return(refuse("an unconditional DELETE"))
  # OUTPUT is read at the top level, or as the last statement of one DO
  # loop. Anywhere else (under IF, or twice in a loop) it decides which rows
  # exist, and no grid emitted without it would be SAS's.
  nested_out <- any(vapply(items, function(it) {
    if (identical(it$kind, "do")) {
      b <- vapply(it$body, function(x) x$kind, character(1L))
      !length(b) || sum(.hzr_sas_tree_kinds(it$body) == "output") != 1L ||
        !identical(b[[length(b)]], "output") || any(b %in% c("do", "delete"))
    } else {
      "output" %in% .hzr_sas_tree_kinds(it$body)
    }
  }, logical(1L)))
  if (nested_out) return(refuse("an OUTPUT or DELETE inside a loop or a condition"))
  if (length(loops) > 1L || (length(loops) && any(top == "output"))) {
    return(refuse("more than one place that writes rows"))
  }
  n_out <- sum(top == "output")

  if (length(set_names)) {
    ins <- lapply(set_names, input)
    for (x in ins) if (!is.null(x$refuse)) return(x)
    each <- lapply(ins, function(x) x$cols)
    known <- unique(unlist(each))
    if (length(set_names) == 1L) {
      emit(call("<-", W, as.name(set_names)))
    } else {
      # SET A B stacks A's rows over B's; a variable only one of them has is
      # missing on the other's rows.
      parts <- lapply(seq_along(set_names), function(j) {
        src <- as.name(set_names[[j]])
        miss <- setdiff(known, each[[j]])
        if (length(miss)) {
          src <- as.call(c(list(quote(cbind), src),
                           stats::setNames(rep(list(NA_real_), length(miss)), miss)))
        }
        call("[", src, known)
      })
      emit(call("<-", W, as.call(c(quote(rbind), parts))))
    }
  } else {
    known <- character(0)
    emit(call("<-", W, quote(data.frame(row.names = 1L))))
  }
  const <- list()
  touch <- function(nms) {
    for (n in setdiff(nms, known)) set_col(n, NA_real_)
    known <<- union(known, nms)
  }

  assign_stmt <- function(text, carried = character(0)) {
    v <- sub(" *=.*$", "", text)
    e <- .hzr_sas_expr(sub("^[A-Z_][A-Z0-9_]* *= *", "", text))
    why <- if (!e$ok) {
      e$why
    } else if (length(intersect(e$refs, carried))) {
      paste0("reads ", paste(intersect(e$refs, carried), collapse = ", "),
             " as the previous loop pass left it")
    } else if (length(setdiff(e$refs, known))) {
      paste0("reads ", paste(setdiff(e$refs, known), collapse = ", "),
             ", which the grid does not carry at this point")
    }
    # Every variable the statement names is in SAS's program data vector,
    # missing until set, and so in the dataset it writes.
    touch(.hzr_sas_names_read(sub("^[A-Z_][A-Z0-9_]* *= *", "", text)))
    if (!is.null(why)) {
      set_col(v, NA_real_)
      u1(text, paste0("not translated: it ", why, ". ", v, " is NA in the grid"))
      const[[v]] <<- NULL
    } else {
      set_col(v, if (length(e$refs)) bquote(with(.(W), .(e$call))) else e$call)
      # SAS gives a missing value where R gives -Inf, Inf or NaN: LOG of 0 or
      # a negative, an EXP that overflows, a division by zero.
      if (is.call(e$call) && !grepl("\"", deparse1(e$call), fixed = TRUE)) {
        col <- call("$", W, as.name(v))
        emit(call("<-", call("[", col, call("!", call("is.finite", col))), NA))
      }
      val <- if (all(e$refs %in% names(const))) .hzr_sas_fold(e$call, const)
      if (is.null(val)) const[[v]] <<- NULL else const[[v]] <<- val
    }
    known <<- union(known, v)
  }
  cond_stmt <- function(it) {
    fx <- .hzr_sas_cond_effects(it)
    if (length(fx$bad)) return(paste0("`", fx$bad[[1L]], "` under a condition"))
    # A subsetting IF or an IF ... THEN DELETE decides which rows exist.
    if (fx$deletes) return(paste0("`", it$text, "`, which decides which rows exist"))
    touch(fx$reads)
    for (v in fx$targets) {
      set_col(v, NA_real_)
      const[[v]] <<- NULL
    }
    known <<- union(known, fx$targets)
    if (length(fx$targets)) {
      shown <- if (identical(it$kind, "ifdo")) paste0(it$text, "; ... END") else it$text
      u1(shown, paste0("a conditional statement is not translated: it sets ",
                       paste(fx$targets, collapse = ", "), ", NA in the grid"))
    }
    NULL
  }

  # Two OUTPUTs write the program data vector twice per input row, so every
  # variable exists at both, missing where not yet set.
  if (n_out > 1L) {
    all_targets <- unique(unlist(lapply(items, function(it) {
      if (identical(it$kind, "assign")) sub(" *=.*$", "", it$text) else .hzr_sas_cond_effects(it)$targets
    })))
    for (v in setdiff(all_targets, known)) set_col(v, NA_real_)
    known <- union(known, all_targets)
    emit(quote(.out <- list()))
  }

  drop <- character(0)
  keep <- NULL
  done <- FALSE
  k_out <- 0L
  for (it in items) {
    if (it$kind %in% c("drop", "keep")) {
      nms <- strsplit(trimws(sub("^(DROP|KEEP)", "", it$text)), " +")[[1L]]
      if (!all(grepl("^[A-Z_][A-Z0-9_]*$", nms))) {
        return(refuse(paste0("`", it$text, "`, which is not a plain list of variables")))
      }
      if (identical(it$kind, "drop")) drop <- c(drop, nms) else keep <- c(keep, nms)
      next
    }
    # After the step's only OUTPUT, or its output loop, nothing reaches the
    # grid: statements there change the program data vector and no row is
    # written from it.
    if (done) next
    if (identical(it$kind, "assign")) {
      assign_stmt(it$text)
    } else if (it$kind %in% c("if", "ifdo")) {
      bad <- cond_stmt(it)
      if (!is.null(bad)) return(refuse(bad))
    } else if (identical(it$kind, "output")) {
      if (n_out == 1L) {
        done <- TRUE
      } else {
        k_out <- k_out + 1L
        emit(bquote(.out[[.(k_out)]] <- .(W)))
      }
    } else if (identical(it$kind, "do")) {
      m <- regmatches(it$text, regexec("^DO +([A-Z_][A-Z0-9_]*) *= *(.+)$", it$text))[[1L]]
      var <- m[[2L]]
      vals <- .hzr_sas_do_values(m[[3L]], const)
      if (!is.null(vals$refuse)) return(refuse(paste0("`", it$text, "`: ", vals$refuse)))
      vals_call <- if (length(vals$refs)) bquote(with(.(W)[1L, , drop = FALSE], .(vals$call))) else vals$call
      # Each input row is repeated once per loop value, rows kept together,
      # as the DATA step writes them.
      emit(call("<-", quote(.v), vals_call))
      emit(bquote(.n <- nrow(.(W))))
      emit(bquote(.(W) <- .(W)[rep(seq_len(.n), each = length(.v)), , drop = FALSE]))
      set_col(var, quote(rep(.v, times = .n)))
      known <- union(known, var)
      const[[var]] <- NULL
      body <- it$body[-length(it$body)]
      assigned <- vapply(body, function(b) {
        if (identical(b$kind, "assign")) sub(" *=.*$", "", b$text) else NA_character_
      }, character(1L))
      # A variable an IF sets in the body may carry its value into the next
      # pass, whichever statement reads it: unresolvable for the whole loop.
      cond_set <- unique(unlist(lapply(body, function(b) {
        if (identical(b$kind, "assign")) NULL else .hzr_sas_cond_effects(b)$targets
      })))
      for (j in seq_along(body)) {
        b <- body[[j]]
        if (identical(b$kind, "assign")) {
          assign_stmt(b$text, carried = union(assigned[j:length(assigned)], cond_set))
        } else {
          bad <- cond_stmt(b)
          if (!is.null(bad)) return(refuse(bad))
        }
      }
      done <- TRUE
    }
  }
  if (n_out > 1L) {
    emit(bquote(.(W) <- do.call(rbind, .out)))
    emit(bquote(.(W) <- .(W)[order(rep(seq_len(nrow(.out[[1L]])), times = length(.out))), , drop = FALSE]))
  }
  if (length(drop)) {
    emit(bquote(.(W) <- .(W)[setdiff(names(.(W)), .(unique(drop)))]))
    known <- setdiff(known, drop)
  }
  if (!is.null(keep)) {
    emit(bquote(.(W) <- .(W)[intersect(names(.(W)), .(unique(keep)))]))
    known <- intersect(known, keep)
  }
  emit(bquote(rownames(.(W)) <- NULL))
  list(code = code, cols = known, untr = untr, warn = warn)
}

#' Translate the DATA steps that build a HAZPRED prediction grid.
#'
#' Resolves the grid as SAS does: the last definition of `name` before the
#' PROC HAZPRED block, each of its `SET` inputs the last definition before
#' that step, and so on back. The returned call builds every step it needs,
#' in file order, inside one `local()`, and ends by copying the column the
#' HAZPRED `TIME` statement names into the `time` column `predict()` reads.
#'
#' @param txt The whole normalised source.
#' @param name The HAZPRED `DATA=` dataset.
#' @param time_var The variable the HAZPRED `TIME` statement names.
#' @param before Offset of the PROC HAZPRED block; only definitions that
#'   start before it count.
#' @param blocks `.hzr_sas_blocks(txt)`.
#' @return `list(call, untranslated, reason)`. `call` is `NULL` when the grid
#'   cannot be built, and `reason` then says why. `NULL` means untranslated,
#'   never "no grid": the caller must record it.
#' @noRd
.hzr_parse_grid <- function(txt, name, time_var, before = nchar(txt) + 1L,
                            blocks = .hzr_sas_blocks(txt)) {
  empty <- .hzr_untranslated_frame()
  refuse <- function(why) list(call = NULL, untranslated = empty, reason = why)
  if (is.null(name) || !nzchar(name)) return(refuse("no DATA= dataset was named"))
  name <- .hzr_sas_ds(name)
  if (is.null(time_var) || is.na(time_var) || !nzchar(time_var)) {
    return(refuse(paste(
      "the block has no TIME statement, which PROC HAZPRED requires",
      "(hazpred/timeprc.c:10-14), so there is no time to predict at")))
  }
  events <- .hzr_sas_dataset_events(txt, blocks)
  # A step inside a %MACRO definition runs where the macro is called, so it
  # builds nothing where it stands, except for a PROC HAZPRED in the same
  # definition, which runs with it.
  st <- .hzr_sas_statements(txt)
  here <- .hzr_sas_scope_at(st, .hzr_sas_macro_defs(st$text), before)
  events <- Filter(function(e) e$scope %in% c(0L, here), events)
  memo <- list()
  built <- integer(0)
  # The offset at which each definition is last read: a call after the
  # definition and before that offset may have rewritten it.
  until <- rep(-Inf, length(events))
  build <- function(k) {
    key <- as.character(k)
    if (!is.null(memo[[key]])) return(memo[[key]])
    ev <- events[[k]]
    input <- function(nm) {
      j <- .hzr_sas_resolve(events, nm, ev$pos)
      if (is.na(j)) {
        return(list(refuse = paste0(ev$name, " reads ", nm, ", which no DATA step in this job ",
                                    "builds before it")))
      }
      until[[j]] <<- max(until[[j]], ev$pos)
      build(j)
    }
    res <- switch(ev$kind,
      data = .hzr_sas_translate_step(ev, input),
      hazpred = {
        r <- input(ev$from)
        if (!is.null(r$refuse)) r else
          list(code = list(call("<-", as.name(ev$name), as.name(ev$from))),
               cols = r$cols, untr = empty, warn = character(0))
      },
      sort = {
        r <- input(ev$from)
        if (!is.null(r$refuse)) {
          r
        } else if (!all(ev$by %in% r$cols)) {
          list(refuse = paste0("PROC SORT of ", ev$from, " is by a variable it does not carry"))
        } else {
          from <- as.name(ev$from)
          # SAS orders a numeric missing value below every number, so an
          # ascending sort puts it first. DESCENDING is refused as a sort.
          ord <- as.call(c(quote(order), lapply(ev$by, function(b) call("$", from, as.name(b))),
                           list(na.last = FALSE)))
          list(code = list(bquote(.(as.name(ev$name)) <- .(from)[.(ord), , drop = FALSE])),
               cols = r$cols, untr = empty, warn = character(0))
        }
      },
      list(refuse = paste0(ev$name, " is written by ", ev$what, ", which this translation ",
                           "does not read")))
    if (is.null(res$refuse)) built <<- c(built, k)
    memo[[key]] <<- res
    res
  }

  k <- .hzr_sas_resolve(events, name, before)
  if (is.na(k)) {
    return(refuse(paste0("no DATA step in this job builds ", name, " before this PROC HAZPRED")))
  }
  until[[k]] <- before
  res <- build(k)
  if (!is.null(res$refuse)) return(refuse(res$refuse))
  if (!time_var %in% res$cols) {
    return(refuse(paste0(
      "the TIME variable ", time_var, " is not a variable of ", name,
      ", so PROC HAZPRED stops (hazpred/timeprc.c:16-20)")))
  }
  ks <- sort(unique(built))
  parts <- lapply(ks, function(j) memo[[as.character(j)]])
  code <- unlist(lapply(parts, function(x) x$code), recursive = FALSE)
  untr <- do.call(rbind, c(list(empty), lapply(parts, function(x) x$untr)))
  # An %INCLUDE or a macro call this cannot read, between a definition the
  # grid uses and the step that reads it, may rewrite it. A dataset it is
  # known to write is an `opaque` event, and refuses the grid above, as a
  # macro's OUT= does. This one may write nothing, so the grid the job shows
  # is emitted and the call recorded, as the PARMS macros of #601 are.
  for (cl in Filter(function(e) identical(e$kind, "call"), events)) {
    hit <- vapply(ks, function(j) events[[j]]$pos < cl$pos && cl$pos < until[[j]], logical(1L))
    if (!any(hit)) next
    nms <- unique(vapply(ks[hit], function(j) events[[j]]$name, ""))
    nms <- paste(nms, collapse = " and ")
    untr <- rbind(untr, .hzr_untranslated_frame(NA_integer_, cl$what, paste0(
      "This call runs after the step that builds ", nms, " and before PROC HAZPRED reads ",
      "it, so it may rewrite ", nms, ", and hzr_translate_sas() cannot read what it ",
      "does. The grid emitted for DATA=", name, " is the one the job's DATA steps ",
      "show, without anything this call changes. Check that it leaves ", nms,
      " as it is, or build the grid by hand.")))
  }
  warn <- unlist(lapply(parts, function(x) x$warn))
  W <- as.name(name)
  body <- c(
    if (length(warn)) {
      list(call("warning", paste0(
        "This grid was built by SAS DATA step statements that hzr_translate_sas() does not ",
        "translate, and the variables they set are NA in this grid:\n",
        paste0("  ", warn, collapse = "\n")), call. = FALSE))
    },
    code,
    # PROC HAZPRED reads its time from the variable TIME names
    # (hazpred/timeprc.c:16-25); predict() reads a column named `time`.
    list(call("<-", call("$", W, as.name("time")), call("$", W, as.name(time_var))), W)
  )
  list(call = call("local", as.call(c(as.name("{"), body))), untranslated = untr, reason = NULL)
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
#' Hazard limits are on the log scale in both engines, so the hazard call
#' leaves `conf.type` at its default.
#'
#' `level` is set on both calls. PROC HAZPRED's default `CLIMITS` is 0
#' (`hazpred/stmtprc.c:14`), and any `CLIMITS` outside `(0, 1)` gives a
#' multiplier of one (`hazpred/hzpp.c:8-9`), so its default band is one
#' standard error, where `predict.hazard()` defaults to 95%. Until #493 the
#' level was never emitted, so every band was 1.96 times too wide, and
#' `CLIMITS=` was read and discarded.
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
  toks <- .hzr_sas_join_spaced(toks[nzchar(toks)])
  # `INHAZ= OUT=P`: PROC HAZPRED reads OUT as INHAZ's dataset name
  # (hazpred_l.l:36, :48) and meets a syntax error at the `=` after it. The
  # joiner pairs it the same way, which named OUT= as the missing option.
  # Split it back so the empty option is the one named (#498 review).
  ds_keys <- c("DATA", "INHAZ", "OUT")
  valueless <- character(0)
  i <- 1L
  while (i < length(toks) - 1L) {
    kv <- strsplit(toks[[i]], "=", fixed = TRUE)[[1L]]
    if (length(kv) == 2L && kv[[1L]] %in% ds_keys &&
          kv[[2L]] %in% c(ds_keys, "CL", "CLIMITS") &&
          identical(toks[[i + 1L]], "=") && !identical(toks[[i + 2L]], "=")) {
      valueless <- c(valueless, kv[[1L]])
      toks <- c(toks[seq_len(i - 1L)], paste0(kv[[1L]], "="),
                paste0(kv[[2L]], "=", toks[[i + 2L]]),
                toks[-seq_len(i + 2L)])
    }
    i <- i + 1L
  }
  # Same stray-`=` rule as the PROC HAZARD line. This caller shares the
  # joiner but had neither this nor a presence check, so a stray `=`
  # recorded a BLANK-keyword "unknown option" row and the prediction calls
  # were emitted for a job SAS rejects (#433 review 3).
  pred_syntax_error <- NULL
  stray <- which(toks == "=")
  if (length(stray)) {
    drop <- unique(c(stray, stray[stray < length(toks)] + 1L))
    leftover <- paste(toks[drop], collapse = " ")
    toks <- toks[-drop]
    pred_syntax_error <- paste0(
      "a stray `=` on the PROC HAZPRED line (", leftover, "): the option ",
      "before it already took its value, so PROC HAZPRED reaches a syntax ",
      "error (hazard_y.y:102) and rejects this job")
  }
  data_name <- NULL
  inhaz <- NULL
  out_given <- FALSE
  # `KEY '=' dsfield` (hazpred_y.y:50-52), and dsfield is a NAME or a
  # LIB.MEMBER (:62-64, hazpred_l.l:17-18). Any other value, `""` or `WORK.`
  # included, is a syntax error, and PROC HAZPRED stops (initprz.c:53-55).
  # A macro is exempt: SAS expands it first.
  given <- c(DATA = FALSE, INHAZ = FALSE, OUT = FALSE)
  bad_ds <- character(0)
  ds_re <- "^[A-Z_][A-Z0-9_]*([.][A-Z_][A-Z0-9_]*)?$"
  # A macro value passes only if it can expand to a name: each reference
  # (`%F(...)`, `&X.`, `&X`) stands in as a name, and the result must still
  # match. `.&X` and `1&X` begin with text PROC HAZPRED cannot read
  # (hazpred_l.l:56) whatever `&X` holds. A reference this cannot place
  # (nested parentheses) keeps the old blanket exemption.
  macro_can_name <- function(v) {
    s <- gsub("%[A-Z_][A-Z0-9_]*[(][^()]*[)]", "M", v)
    s <- gsub("&[A-Z_][A-Z0-9_]*[.]?", "M", s)
    grepl("[&%]", s) || grepl(ds_re, s)
  }
  raw_ops <- strsplit(trimws(st[[1L]]), " ", fixed = TRUE)[[1L]]
  check_ds <- function(key, val) {
    given[[key]] <<- TRUE
    if (grepl(ds_re, val) ||
          (.hzr_sas_is_macro(val) && macro_can_name(val))) {
      return(TRUE)
    }
    # Quote the value as written: the joiner splits `PRED(WHERE=(...))` at
    # its `=`, and `val` holds only the first piece.
    raw <- raw_ops[startsWith(raw_ops, paste0(key, "="))]
    if (length(raw) == 1L && nzchar(val)) val <- substring(raw, nchar(key) + 2L)
    why <- if (key %in% valueless) {
      paste0(key, "= has no dataset name: PROC HAZPRED reads the option ",
             "keyword after it as the name (hazpred_l.l:35-37, :48), and the ",
             "`=` that follows is a syntax error")
    } else if (!nzchar(val)) {
      paste0(key, "= has no dataset name, and PROC HAZPRED has no form of ",
             "it without one (hazpred_y.y:50-52, :62-64), a syntax error")
    } else {
      paste0(key, "=", val, " is not a NAME or a LIB.MEMBER ",
             "(hazpred_l.l:17-18, :47-48), and PROC HAZPRED rejects the text ",
             "it cannot read (:56-57), a syntax error")
    }
    why <- paste0(why, "; PROC HAZPRED then stops (initprz.c:53-55)")
    bad_ds <<- c(bad_ds, why)
    note(paste0(key, "="), why)
    FALSE
  }
  want_surv <- TRUE
  want_haz <- TRUE
  want_cl <- TRUE
  climit <- NULL
  time_var <- NULL

  if (!is.null(pred_syntax_error)) note("PROC HAZPRED", pred_syntax_error)
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
      DATA    = if (check_ds(token, val)) data_name <- val,
      INHAZ   = if (check_ds(token, val)) inhaz <- val,
      OUT     = if (check_ds(token, val)) out_given <- TRUE,
      NOSURV  = want_surv <- FALSE,
      NOHAZ   = want_haz <- FALSE,
      NOCL    = want_cl <- FALSE,
      # CLIMITS= sets only the level (hazpred/hazpprc.c:15-16); NOCL is a
      # separate flag that wins whatever the order (hzpp.c:5-6). Setting
      # want_cl here turned `NOCL CLIMITS=0.9` back into a banded job.
      # The value is the lexer's NUMBER (hazpred_l.l:13-16, unsigned); any
      # other value is a syntax error at hazpred_y.y:53.
      # A value that is not translated is not mapped, as MAXITER= and
      # CONDITION= count it, so the coverage figure does not claim it.
      CLIMITS = if (grepl("^([0-9]+|[0-9]*[.][0-9]+(E[+-]?[0-9]+)?)$", val)) {
        climit <- as.numeric(val)
      } else if (.hzr_sas_is_macro(val)) {
        # SAS expands a macro before PROC HAZPRED lexes the option, so it
        # may well be a valid number: no syntax-error verdict here.
        mapped <- mapped - 1L
        note(key, paste0(
          "CLIMITS=", val, " is a SAS macro reference, which SAS resolves ",
          "before PROC HAZPRED reads the option, so this translation cannot ",
          "tell what level it names; the bands are drawn at the one-SE ",
          "default"))
      } else {
        mapped <- mapped - 1L
        note(key, paste0(
          "CLIMITS= takes an unsigned number, not `", val, "`: PROC ",
          "HAZPRED reaches a syntax error (hazpred_y.y:53) and rejects this ",
          "job; the bands are drawn at the one-SE default"))
      },
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
    # `TIME NAME` (hazpred_y.y:80): the grid variable predictions are made
    # at. The last one written is the one setvar(11, ...) keeps.
    if (identical(token, "TIME")) time_var <- if (length(w) >= 2L) w[[2L]] else NA_character_
    if (!(token %in% c("TIME", "ID"))) {
      mapped <- mapped - 1L
      note(kw, "SAS listing control; no R effect")
    }
  }

  grid <- NULL
  grid_refused <- FALSE
  if (!is.null(data_name)) {
    g <- .hzr_parse_grid(txt, data_name, time_var,
                         before = if (is.null(block$start)) nchar(txt) + 1L else block$start)
    grid <- g$call
    untr <- rbind(untr, g$untranslated)
    grid_refused <- is.null(grid)
    if (grid_refused) {
      note(paste0("DATA=", data_name), paste0("prediction grid not translated: ", g$reason))
    }
  }

  # %HAZPRED requires DATA=, INHAZ= and OUT= on the PROC statement and stops
  # with "HAZPRED not attempted" when any is missing (hazpred.sas:13-33,
  # :153-163), so SAS predicts nothing. The block used to emit predict()
  # anyway: over the fitting rows with no DATA=, from `fit` with no INHAZ=
  # (#498). Each missing option is its own row; the stop() names them all.
  # sprintf, not paste0: paste0(character(0), "=") is "=", not empty.
  absent <- sprintf("%s=", names(given)[!given])
  or_list <- function(x) {
    if (length(x) > 1L) {
      paste(paste(x[-length(x)], collapse = ", "), "or", x[length(x)])
    } else {
      x
    }
  }
  # The refusal is claimed only where the macro makes it. A macro reference
  # on the statement is expanded before %HAZPRED reads &syspbuff, so it may
  # supply the option. And each of the macro's tests is a substring test
  # (%index for DATA, INHAZ and " OUT", hazpred.sas:13, :20, :27), so
  # `INHAZ=HAZDATA` passes the DATA= test; the macro then reads the value
  # after the next `=` in its place. The block stops either way (#498 review).
  pred_macros <- toks[vapply(toks, .hzr_sas_is_macro, logical(1))]
  probe <- c(DATA = "DATA", INHAZ = "INHAZ", OUT = " OUT")
  passes <- absent[vapply(sub("=$", "", absent), function(k) {
    grepl(probe[[k]], st[[1L]], fixed = TRUE)
  }, logical(1))]
  refused <- setdiff(absent, passes)
  macro_refusal <- if (length(pred_macros)) {
    paste0(
      "This PROC HAZPRED block names no ", or_list(absent), ", but its PROC ",
      "statement carries ", paste(pred_macros, collapse = ", "), ", which ",
      "SAS expands before the %HAZPRED macro reads the statement, so it may ",
      "supply the missing option(s). This translation cannot see what it ",
      "expands to. Write the option(s) out and translate the job again.")
  } else {
    paste0(
      if (length(refused)) paste0(
        "This PROC HAZPRED block names no ", or_list(refused), ", and the ",
        "%HAZPRED macro requires DATA=, INHAZ= and OUT= on the PROC ",
        "statement: without one it stops with \"HAZPRED not attempted\" ",
        "(hazpred.sas:13-33, :153-163), so SAS predicts nothing. "),
      if (length(passes)) paste0(
        "This PROC HAZPRED block names no ", or_list(passes), ", but the ",
        "%HAZPRED macro's test for it is a substring test (hazpred.sas:13, ",
        ":20, :27), which other text on the PROC statement satisfies, so ",
        "the macro takes the value after the next `=` in its place, not a ",
        "dataset the job names. "),
      "Add the missing option(s) and translate the job again.")
  }
  for (a in absent) note(a, macro_refusal)
  # A stray `=` is a syntax error PROC HAZPRED stops on too; recording it
  # alone left predict() in the document (#498, Copilot).
  ds_refusal <- c(
    if (!is.null(pred_syntax_error)) paste0(
      "PROC HAZPRED rejects this block: ", pred_syntax_error, "."),
    if (length(absent)) macro_refusal,
    if (length(bad_ds)) paste0(
      "PROC HAZPRED rejects this block: ", paste(bad_ds, collapse = "; "),
      ". Correct the option(s) named here and translate the job again."))

  mk <- function(type) {
    if (length(ds_refusal)) {
      return(as.call(list(quote(stop), paste(ds_refusal, collapse = " "),
                          call. = FALSE)))
    }
    # A refused grid is a refusal, not an absent argument. Emitting
    # predict(fit, newdata = <name>) when no chunk builds <name> leaves the
    # document to fail on an unbound name -- or, if an object of that name
    # happens to exist in the rendering session, to predict over unrelated
    # data and report it. That is the whole reason .hzr_parse_grid() refuses
    # a partially-read grid, so refuse in the emitted document too, the way
    # the LCENSOR + ICENSOR and SELECTION refusals do.
    if (grid_refused) {
      return(as.call(list(quote(stop), paste0(
        "The ", type, " predictions of this PROC HAZPRED block read the ",
        "grid ", data_name, ", which hzr_translate_sas() does not translate: ",
        g$reason, ". Build ", data_name,
        " by hand (a data frame with a `time` column) and call predict() ",
        "yourself; rendering over whatever else is named ", data_name,
        " would report predictions over the wrong times."
      ))))
    }
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
    # PROC HAZPRED's default CLIMITS is 0 (hazpred/stmtprc.c:14), and any
    # CLIMITS outside (0, 1) gives a multiplier of exactly one
    # (hazpred/hzpp.c:8-9): a one-SE band. predict.hazard() defaults to
    # 0.95, so leaving `level` out drew every band 1.96 times too wide (#493).
    if (isTRUE(want_cl)) {
      args$level <- if (!is.null(climit) && climit > 0 && climit < 1) {
        climit
      } else {
        quote(2 * stats::pnorm(1) - 1)
      }
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

# ---------------------------------------------------------------------------
# %repeat -> hzr_repeated_events()
# ---------------------------------------------------------------------------

# %repeat's keyword parameters and defaults, as the macro declares them
# (~/Documents/macro.library/repeat.sas), upper-cased as
# .hzr_sas_normalise() leaves every name.
.hzr_repeat_defaults <- c(
  IN = "BUILT", OUT = "EVENTS", EVENTYPE = "EVENTYPE", IV_EVENT = "IV_EVENT",
  IV_END = "IV_END", ID = "ID", EVENT = "EVENT", EVENT_NO = "EVENT_NO",
  RCENSOR = "RCENSOR", IV_START = "IV_START", IV_SEG = "IV_SEG", RENEWAL = "RENEWAL"
)

# hzr_repeated_events()'s fixed output names, each mapped to the macro
# parameter that renames it. first and last have no parameter: the macro
# always writes them under those names.
.hzr_repeat_outputs <- c(
  event = "EVENT", event_no = "EVENT_NO", rcensor = "RCENSOR",
  iv_start = "IV_START", iv_seg = "IV_SEG", renewal = "RENEWAL"
)

#' Translate one `%repeat(...)` call into an `hzr_repeated_events()` chunk.
#'
#' The emitted chunk renames the function's lower-case outputs to the names the
#' job uses, upper-cased like every name this translator emits; without that the
#' fit reads a column that does not exist. Before the call it drops any input
#' column named like one of those outputs. SAS overwrites such a column, and the
#' macro assigns each output on every row before reading it, so dropping is
#' exactly what SAS does. Left in place, the rename would produce two columns of
#' one name, and `$` would return the stale one.
#'
#' A call this cannot express is refused rather than guessed at: the chunk is a
#' `stop()`, and each problem is an `$untranslated` row. That covers an unknown
#' keyword, a positional argument, a value that is not a plain SAS name (a
#' `&macro` reference, say), an input column that is also an output
#' (`EVENTYPE=RCENSOR`: SAS zeroes that indicator before reading it, so the job's
#' own answer is not the model it describes), and two outputs with one name. A
#' keyword given twice, which SAS rejects, is refused too, and so is an output
#' named `LAG_IV` or `NUMBER`, the macro's own loop counters. The body is split
#' on every comma; a value holding parentheses, where a comma could be nested, is
#' refused as not a plain name either way.
#' @noRd
.hzr_parse_repeat <- function(block) {
  parts <- if (nzchar(block$text)) trimws(strsplit(block$text, ",", fixed = TRUE)[[1L]]) else character(0)
  args <- .hzr_repeat_defaults
  problems <- character(0)
  seen_keys <- character(0)

  for (p in parts) {
    kv <- regmatches(p, regexec("^([A-Z_][A-Z0-9_]*) ?= ?(.*)$", p))[[1L]]
    if (!length(kv)) {
      problems <- c(problems, sprintf("`%s` is not a KEY=VALUE argument", p))
      next
    }
    key <- kv[2L]
    val <- trimws(kv[3L])
    if (!key %in% names(args)) {
      problems <- c(problems, sprintf("`%s=` is not a %%repeat parameter", key))
      next
    }
    if (key %in% seen_keys) {
      problems <- c(problems, sprintf("`%s=` is given more than once", key))
      next
    }
    seen_keys <- c(seen_keys, key)
    name_re <- if (key %in% c("IN", "OUT")) {
      "^[A-Z_][A-Z0-9_]*([.][A-Z_][A-Z0-9_]*)?$"
    } else {
      "^[A-Z_][A-Z0-9_]*$"
    }
    if (!grepl(name_re, val)) {
      problems <- c(problems, sprintf("`%s=%s` is not a plain SAS name", key, val))
      next
    }
    args[[key]] <- val
  }

  # A WORK. libref names the same dataset as the bare name, and the rewrite
  # scan compares bare names. A fit written DATA=WORK.X keeps its libref, so it
  # is not matched to this OUT= and gets the ordinary "Assign" guard instead.
  args[c("IN", "OUT")] <- sub("^WORK[.]", "", args[c("IN", "OUT")])

  inputs <- unname(args[c("ID", "EVENTYPE", "IV_EVENT", "IV_END")])
  targets <- c(unname(args[.hzr_repeat_outputs]), "FIRST", "LAST")
  for (hit in intersect(inputs, targets)) {
    problems <- c(problems, sprintf(paste(
      "input column %s is also an output %%repeat writes; SAS overwrites it before reading it,",
      "so the job's own result is not the model it describes"
    ), hit))
  }
  dup <- unique(targets[duplicated(targets)])
  if (length(dup)) {
    problems <- c(problems, sprintf("two outputs are both named %s", paste(dup, collapse = ", ")))
  }
  for (hit in intersect(targets, c("LAG_IV", "NUMBER"))) {
    problems <- c(problems, sprintf(paste(
      "output %s has the name of a %%repeat loop counter; SAS's RETAIN reads it,",
      "so the job's own result is not the model it describes"
    ), hit))
  }

  if (length(problems)) {
    msg <- paste0("hzr_translate_sas() did not translate this job's %repeat call: ",
                  paste(problems, collapse = "; "), ".")
    return(list(
      call = bquote(stop(.(msg))),
      untranslated = .hzr_untranslated_frame(rep(NA_integer_, length(problems)),
                                             rep("%repeat", length(problems)), problems),
      tokens_seen = length(parts), tokens_mapped = 0L, in_name = NULL, out_name = NULL
    ))
  }

  in_name <- args[["IN"]]
  out_name <- args[["OUT"]]
  from <- c(names(.hzr_repeat_outputs), "first", "last")
  overwrite_msg <- paste0(" in ", in_name, ": %repeat writes columns of these names, so they were dropped ",
                          "before the call, as SAS overwrites them.")
  loop_msg <- paste0(" in ", in_name, ": SAS's %repeat reads a LAG_IV or NUMBER column in place of its ",
                     "own loop counters, so for this job the SAS listing's segment starts and event ",
                     "counts are not comparable with these.")
  call <- bquote(.(as.name(out_name)) <- local({
    d <- .(as.name(in_name))
    drop <- intersect(names(d), .(targets))
    if (length(drop)) {
      warning(paste(drop, collapse = ", "), .(overwrite_msg), call. = FALSE)
      d <- d[setdiff(names(d), drop)]
    }
    loop <- intersect(names(d), c("LAG_IV", "NUMBER"))
    if (length(loop)) warning(paste(loop, collapse = ", "), .(loop_msg), call. = FALSE)
    out <- hzr_repeated_events(d, id = .(args[["ID"]]), time = .(args[["IV_EVENT"]]),
                               followup = .(args[["IV_END"]]), indicator = .(args[["EVENTYPE"]]))
    names(out)[match(.(from), names(out))] <- .(targets)
    out
  }))

  list(call = call, untranslated = .hzr_untranslated_frame(),
       tokens_seen = length(parts), tokens_mapped = length(parts),
       in_name = in_name, out_name = out_name)
}

#' Steps between a `%repeat` call and a fit that may change the macro's `OUT=`.
#'
#' `segment` is the normalised source between the two, and the scan fails
#' closed. A `DATA` step is a hit when its output list names `out`; a step that
#' only reads it (`SET out`) is not. Any other step or statement that names
#' `out` is a hit as well -- `PROC SQL`, `PROC APPEND`, `PROC DATASETS`, a sort
#' with `NODUPKEY`, `OUT=` or `WHERE=`, a macro call -- because a false stop
#' costs the reader one deleted chunk and a missed rewrite fits data SAS did
#' not fit. The one exception is the plain `PROC SORT DATA=out`, which only
#' reorders rows and so cannot change the likelihood. A `WORK.` prefix names
#' the same dataset as the bare name. Each hit is returned quoted: the step's
#' statements up to the next `DATA`, `PROC`, `%HAZ`, `%REPEAT`, `RUN` or `QUIT`.
#' A step that changes `out` without naming it (a macro that writes it
#' internally) cannot be seen from here. Any statement starting with `%` (a
#' macro call) is a step of its own, so it cannot hide inside an exempt one.
#' A sort counts as plain only when every statement after its PROC line is a
#' `BY`: a `WHERE` statement subsets the data. A step that uses a macro
#' variable (`&name`) is a hit, because the variable could name `out`.
#' @noRd
.hzr_repeat_rewrites <- function(segment, out) {
  out <- sub("^WORK[.]", "", out)
  stmts <- trimws(strsplit(segment, ";", fixed = TRUE)[[1L]])
  stmts <- stmts[nzchar(stmts)]
  boundary <- "^(DATA |PROC |%|RUN$|QUIT$)"
  names_in <- function(s) sub("^WORK[.]", "", strsplit(s, "[^A-Z0-9_.]+")[[1L]])
  plain_sort <- paste0("PROC SORT DATA=", c(out, paste0("WORK.", out)))
  hits <- character(0)
  i <- 1L
  while (i <= length(stmts)) {
    last <- i
    if (grepl("^(DATA |PROC )", stmts[i])) {
      while (last < length(stmts) && !grepl(boundary, stmts[last + 1L])) last <- last + 1L
    }
    step <- stmts[i:last]
    s <- stmts[i]
    writes <- if (startsWith(s, "DATA ")) {
      targets <- strsplit(trimws(gsub("[(][^)]*[)]", " ", substring(s, 6L))), " +")[[1L]]
      out %in% sub("^WORK[.]", "", targets)
    } else if (s %in% plain_sort && all(grepl("^BY ", step[-1L]))) {
      FALSE
    } else {
      any(vapply(step, function(x) out %in% names_in(x), logical(1L)))
    }
    # A macro variable (&DSN) could name OUT; its value is not known here, so
    # a step that uses one is treated as naming it -- the same call
    # .hzr_parse_repeat() makes when it refuses IN=&DSN.
    writes <- writes || any(grepl("&[A-Z_]", step))
    if (writes) hits <- c(hits, paste0(paste(step, collapse = "; "), ";"))
    i <- last + 1L
  }
  hits
}

#' The stop() chunks and $untranslated rows for steps that may change `out`.
#'
#' One chunk per step `.hzr_repeat_rewrites()` finds in `segment`, each quoting
#' the step and telling the reader to replace it with R code, or delete it if
#' the step leaves `out` unchanged.
#' @noRd
.hzr_rewrite_stops <- function(segment, out) {
  steps <- .hzr_repeat_rewrites(segment, out)
  stops <- lapply(steps, function(step) {
    bquote(stop(.(paste0(
      "This job may change ", out, " after %repeat, in a SAS step that ",
      "hzr_translate_sas() does not translate: ", step, " Replace this chunk ",
      "with R code that makes the same change to ", out, ", or delete it if ",
      "the step leaves ", out, " unchanged."
    ))))
  })
  list(
    calls = stops,
    untranslated = .hzr_untranslated_frame(
      rep(NA_integer_, length(steps)), rep(paste0(out, " changed after %repeat"), length(steps)), steps
    )
  )
}
