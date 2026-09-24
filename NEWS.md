# TemporalHazard 1.2.11

## Breaking changes

* **A `hzr_translate_sas()` job PROC HAZARD refuses now warns loudly, and
  says so in `$untranslated`** (#359). When `SETG3` sets an error, the
  procedure exits in `shape()` before `results()`, so the job produces
  nothing. The translation recorded that as an untranslated row and emitted a
  `hazard()` chunk with nothing to mark it, and a reader who rendered past the
  callout got a converged fit standing in for a job with no result. The
  emitted document now carries a `warning()` immediately above the fit,
  naming the `SETG3` code, its cause and the `PARMS` operands that produced
  it, and the row is recorded as before. The fit is still emitted: a rendered
  document completes, and the reader is told what it stands in for.

  **Some of these still fail further down, and the warning says which.**
  Several `SETG3` refusals fire precisely because a shape value is out of
  range, and the same value is out of range for `hzr_phase()`, which will not
  build the phase. Rather than name a list of codes here, which drifted once
  already, the warning itself is derived by trying to construct the phase: it
  says `hzr_phase()` accepts the shape only when it does, and otherwise says
  the document stops at that check. A test executes every `SETG3` class's
  emitted chunks and requires the message and the outcome to agree, so the
  two cannot diverge again. In every case the warning is emitted in its own
  chunk **above** the fit, naming the `SETG3` code and the operand, so the
  cause is stated before `hzr_phase()` refuses; the render then stops there
  with `hzr_phase()`'s own message.

  Be aware of where that warning does and does not appear. When a chunk
  errors, Quarto writes no output document, and `knitr` collects warnings
  **into** the document rather than printing them, so the warning does not
  reach the render console either. What a reader has in that case is the
  emitted `.qmd` itself, where the `warning()` naming `SETG3910` sits
  immediately above the failing fit, and the `$untranslated` row on the
  translated job.

  The refusal is raised only where it is PROC HAZARD's. With `FIXGE2` or
  `FIXGAE2` and no `WEIBULL`, SAS reaches `SETG3` down a path the `setg3.c`
  trace does not model, so the trace's verdict is not used there. Only
  `SETG3`'s entry refusals are raised as refusals on that path (see the next
  entry).

* **More `hzr_translate_sas()` jobs that PROC HAZARD refuses, or fits
  differently, now warn loudly** (#358, #403, #421). Each was already recorded
  as an untranslated row, but nothing in the rendered document said so, and a
  reader met a converged fit with no sign that PROC HAZARD would not have
  produced it. Each now emits the fit, a `warning()` above it naming the
  cause, and the row:
  - a `PARMS` operand PROC HAZARD rejects with a syntax error: a value its
    lexer does not read as a number (`NU=1E-3`, `NU=2.`), a value keyword
    with no `= NUMBER`, a spaced operand that is invalid even joined, or a
    keyword outside its grammar (`FIXG1`);
  - a `MAXITER=` or `CONDITION=` value that its lexer does not read as a
    number, **or no value at all**: `MAXITER '=' NUMBER` and
    `CONDITION '=' NUMBER` (`hazard_y.y:63-64`) have no form without a
    number, so `MAXITER=`, `MAXITER =` and a bare `MAXITER` are each a
    syntax error and the job does not run;
  - a template's `?` placeholder in `PARMS`, which PROC HAZARD's lexer also
    rejects. It was filled from SAS's default and fitted; it now asks to be
    filled in;
  - a model this translation cannot emit:
    - `FIXMNU1` on an active early phase, which PROC HAZARD fits with
      `|M*NU| = 1`; this translation does not mirror that constraint;
    - `DELTA` other than 0 on an active early phase;
    - `FIXTAU` with no `TAU` written, which PROC HAZARD fixes at 0.75 of the
      longest follow-up;
    - `FIXGE2` or `FIXGAE2` without `WEIBULL`. That path is not modelled
      here, so the warning says the translation cannot tell whether PROC
      HAZARD refuses the job or which model it fits;
  - `SETG3`'s entry refusals, on every path;
  - `SETG1`'s refusals for an early phase (#424): `THALF` fixed at a value
    that is not positive (`SETG1910`); `M` and `NU` both fixed on a case no
    model takes (`SETG1940`, `SETG1950`, `SETG1960`, and `SETG1920` and
    `SETG1930` under `FIXMNU1`); and `DELTA` fixed outside `[-1, 1]`
    (`SETG1900`, `SETG1901`). They were fitted with no row, or warned for the
    wrong reason, that the model was not mirrored. Most of these values are
    out of range for `hazard()` as well, so the fit still fails after the
    warning;
  - an early phase PROC HAZARD may not fit (#424). With `NU=0` and `M` free,
    `SETG1` selects its limiting case, and these jobs fitted with no row.
    What the binary does next was measured on two datasets, and the warning
    says only what was seen. For `M=0 NU=0 FIXNU` it produced no result on
    both. For `M=1 NU=0` and `M=-1 NU=0`, with or without `FIXNU`, it stopped
    on a domain error (`DLG1980`) on one dataset and fitted on the other, so
    the warning says PROC HAZARD **may** print no estimates on your data. In
    neither case does it say PROC HAZARD fits another model;
  - a phase variable that is not a name to PROC HAZARD's lexer (#440):
    `AGE*SEX`, `LOG(AGE)`, `B SEX`, `1AGE`. A phase variable must be a NAME,
    `[_A-Z][_A-Z0-9]*` (`hazard_l.l:39`, `phasevar : NAME` at
    `hazard_y.y:213`), so PROC HAZARD rejects such a job at parse. The phase
    parser passed the text through as a column name, with no row and no
    warning, so `EARLY AGE=0.1, AGE*SEX=0.2;` emitted a fit on a column
    called `AGE*SEX` and, under `SELECTION`, a screen offering it as a
    candidate. The operand is now left out of the model, and out of `theta`
    with it, so the emitted formula has one starting value per term. This
    holds with and without `SELECTION`, and for every spacing of the
    operand. A name PROC HAZARD accepts, including `_X1` and a word that is
    a keyword elsewhere (`E`, `EARLY`), does not warn, and neither does a
    macro reference, which SAS expands before its lexer runs.

    Parentheses follow the lexer rather than a rule of thumb, and the
    verdicts are checked against the HAZARD binary (C-Version 4.4.4). `)` is
    whitespace to it (`hazard_l.l:32`), and `(` returns no token but
    switches it to its PROC-line state (`hazard_l.l:56`). So a last item
    `LOG()` or `LOG() = 0.2` is the variable `LOG`, which PROC HAZARD fits
    and the translation now keeps. `LOG(X)`, `AGE(1)`, `LOG() /I`, and every
    item after a `(` in the same statement are rejected. A `(` also clears
    PROC HAZARD's syntax-error flag, so a job whose phase statements carry
    one may run despite an earlier error: the binary runs
    `EARLY AGE*SEX, LOG();` and fits `AGE` alone. The translation does not
    reproduce which variables survive, and says so in the warning.
  - a value on a `PROC HAZARD` option that takes none (#431): `NOCOV=1`,
    `CONSERVE=YES`, `PRINTIT=1`, `NOPRINT=0`, in any spacing. Eleven options
    are bare tokens (`hazard_y.y:65-75`), so the `=` is a syntax error and
    the job does not run. The value was ignored and the job fitted with no
    row;
  - a `TIME`, `EVENT`, `RCENSOR`, `LCENSOR` or `WEIGHT` statement with other
    than one operand (#431). Each takes exactly one name
    (`hazard_y.y:106-127`). `EVENT DEAD EXTRA` fitted on `DEAD` and dropped
    `EXTRA` with nothing said; it now warns and fits on the first operand.
    With no operand at all, `WEIGHT`, `RCENSOR` and `LCENSOR` are left out
    of the fit, and `TIME` or `EVENT` stops the job, since there is no
    variable to fit, unless another statement supplies one (a second
    `TIME`, or `ICENSOR` for `EVENT`). An operand that is a macro reference is not counted,
    since it can expand to any number of names. The HAZARD binary is the
    oracle for both shapes, and it runs either job when a later phase
    statement carries a `(`, as above; the warning says so there.

  Every class above for `SETG1` was checked against the HAZARD binary
  (C-Version 4.4.4), with the data staged as PROC HAZARD reads it, on the
  package's `avc` data and on an independent seeded dataset. Where
  `SETG1` moves a starting value and runs, the translation now starts there
  too, with a row and no warning (#421): a free `THALF` that is not positive
  starts at 1, where it was emitted as written and the document stopped at
  its logarithm; `M=0 NU=0` with both free starts at `M = NU = 1`; and
  `NU=0` with only `M` fixed starts `NU` at 1.

  A `PARMS` or `PROC` value that carries a macro reference (`&X`, `%CALL`) is
  not refused, because SAS expands it before PROC HAZARD reads the statement.
  An operand this translation could not read, for that reason or any other,
  warns on a job whose phases it did build: the unread operand may be the one
  that sets a shape, and the emitted phase would then carry SAS's default
  where the job wrote something else.

  **A refused job no longer stops the render.** The warning is per job: a
  file holding several jobs emits one fit chunk each, and only the refused
  job's fit is preceded by a `warning()` chunk, so every other job is written
  out and runs unchanged. A document carrying a refusal renders to completion
  and shows the warning in its output, where an earlier draft of this work
  made it a `stop()` and Quarto then exited 1 and produced no output document
  at all, including for the jobs before the refused one.

  The exception is the refusals above whose shape `hzr_phase()` will not
  build. They are warned about and emitted like everything else, but
  `hzr_phase()` then refuses the out-of-range value, so a file containing such
  a job still yields no rendered output until it is corrected or removed. A reader who wants the other jobs'
  results in the meantime can delete that job from the file.

  **Which jobs stop and which warn, in one place.** A job stops only where it
  did before this release: a phase statement `PROC HAZARD` refuses at parse
  (#340) other than a phase variable that is not a name (#440, above), a
  `PARMS` statement that builds no phase this translator can use, a job
  with no `DATA=` whose phases name covariates (#311), a `SELECTION`
  job that selects no phase, and a `TIME` or `EVENT` statement with no
  operand (#431) that leaves nothing to fit. Everything else newly
  recognised in this release warns and still fits.

  The risk this accepts, deliberately: a rendered document that shows a
  warning and then carries on to a fit **can** be read as a clean result by
  someone who does not read the warning. That is why the warning is raised
  in its own chunk immediately above the fit rather than folded into it, and
  why every such job also carries a row in `$untranslated` -- the warning is
  read once at render, the row is what a reader can search for afterwards.

* **An operand written with spaces around `=` is read, not split apart**
  (#421). SAS's lexer skips whitespace (`hazard_l.l:32`), so `THALF = 0.3` and
  `MAXITER = 50` are the same jobs as the same operands written without the
  spaces. This translator
  split them on whitespace: the `PARMS` pieces were recorded and the phase was
  built from `PROC HAZARD`'s default instead of the written value, and the
  `PROC` line reported its pieces as unknown options. Operands are joined
  before parsing, on both. A joined operand `PROC HAZARD` still rejects
  (`THALF = ABC`, or `FIXNU = 1`, which takes no value) is a syntax error,
  just as it is when written without the spaces.

  Joining now works the way `PROC HAZARD`'s own lexer does. Whitespace only
  separates tokens there (`hazard_l.l:32`) and `=` is a token in its own
  right (`:55`), so every spelling of one statement is the **same** token
  stream to SAS. The operands are normalised to that token stream first and
  then paired as `KEY = VALUE` by the grammar, so all spellings of a
  statement give one answer by construction rather than by matching
  particular spellings. Two earlier attempts did match spellings, and each
  left another spelling reading a following option as a value: `PROC HAZARD
  DATA = MAXITER = 50` fitted with `data` set to `MAXITER=50` and the
  iteration limit silently dropped.

  A stray `=` left over after that pairing is now recorded and warned about
  as the syntax error it is. `DATA = MAXITER = 50` is read as SAS reads it
  --- `DATA` switches the lexer to its dataset-name state, where `MAXITER`
  is a name (`hazard_l.l:59, :80`), so the dataset is `MAXITER` and the
  trailing `= 50` is a stray `=` that sends `PROC HAZARD` to
  `hazardopt : error` (`hazard_y.y:76`).

* **`DATA=` and `OUTHAZ=` with no value are refused** (#433). `DATA '='
  dsfield` and `OUTHAZ '=' dsfield` (`hazard_y.y:61-62`), where a `dsfield` is
  a name or a libref-qualified name (`:80-81`), have no form without one, so
  the job does not run. `OUTHAZ=` was previously dropped with no row at all
  and the job fitted; `DATA=` surfaced as an internal R error naming neither
  the option nor what was lost. Both now warn and record the construct,
  alongside the existing check on `MAXITER=` and `CONDITION=`.

  One spelling is **not** covered, and fails before the joining can happen:
  `DATA = X` with spaces, on a `PROC HAZARD` line that is not wrapped in a
  `%HAZARD(...)` call. The scanner that cuts a file into blocks treats the
  word `DATA ` as the start of a new block, so the job is truncated after
  `PROC HAZARD` and the translation fails with "The EVENT or ICENSOR
  variable must be specified" for a job that does have an `EVENT` statement.
  This is unchanged from earlier releases; write `DATA=X` without the spaces,
  or wrap the job in `%HAZARD(...)`.

* **`hazard()` refuses a function-valued element of `data` (#420).** `data`
  masks the calling frame while `hazard()` evaluates `time`, `status`,
  `time_lower`, `time_upper` and `weights`, and while it evaluates the
  formula's `Surv()` response. R's function lookup walks past every binding
  that is not a function, so an element such as `rep = function(...) ...`
  was called in place of `base::rep()` by an expression like
  `weights = rep(1, n)`. The fit changed and nothing warned; this has
  shipped since 1.2.2 (#151). `stats::lm()` refuses the same shape.
  Both interfaces were affected. The formula path looked immune only
  because the column-reading step replicates each column to `nrow` and dies
  on a function while doing it -- at one row there is nothing to replicate,
  and a 1-row frame carrying a `round` made `Surv(round(tt), ss)` read the
  masked value as the response.
  **What now errors:** any `data` carrying a *named* element that is a
  function, whether or not an expression calls it. That includes using the
  mask to reach a helper, as in
  `hazard(time = f(t), status = s, data = list(t = ..., s = ..., f = myfun))`,
  and it includes an S4 generic or a reference-class generator, which are
  functions for this purpose, as is an element whose name is
  `NA_character_`, which R binds under the symbol `` `NA` `` and a call can
  reach. Remove the element and pass `data` without it: for a vector
  argument, or the formula's `weights`, define the helper in the calling
  environment; for a helper used inside the `Surv()` response, compute the
  value into a `data` column first, since the response is evaluated without
  the formula's environment.
  **Unaffected:** a numeric element or column of the same name, which was
  never consulted; a data-frame list-column of functions, which is a list;
  and an element with no name, which no expression can look up.

* **`hzr_bootstrap()` no longer counts replicates that estimated nothing as
  successes (#373).** The optimizer stands in 1e10 for a negative
  log-likelihood it could not evaluate, so a fit that never had a
  likelihood reports `objective = -1e10`, which is finite. A Weibull fit
  started at `theta = c(1e10, 1e10)` sits there, and each replicate
  reproduced it: 5 of 5 replicates were counted as successes, every `sd`
  was 0, and nothing warned. A replicate at the sentinel is now a failed
  replicate, with its own reason in `failure_reasons`, so `n_success` can
  be lower than before and the summary is built from fewer replicates; from `theta = 20`,
  where the optimizer moves before the clamp stops it, that is four
  replicates of five. And a free parameter that does not move across the
  replicates that estimated it, to within rounding, is named in a warning:
  from `theta = 50` the objective is finite but every replicate stays at
  its start, with an `sd` of about 1e-14 around 50. Parameters the fit
  holds fixed (by `hzr_phase(..., fixed =)`, by a constraint, or by
  Conservation of Events) are identical by design and are not named. A run
  in which only some replicates stay at their start while their objective is
  finite is not caught; that rests on the optimizer's convergence test
  (#351). What a sentinel objective should mean for a single fit is tracked
  separately (#351, #374).

* **`hzr_translate_sas()` no longer fits a job `PROC HAZARD` rejects: if
  you hold estimates from such a translation, they have no SAS run behind
  them (#340).** A phase statement with an option written in a form SAS's
  lexer or grammar rejects, such as `AGE/EI` (options glued together),
  `AGE/E/I`, `Y/`, `/S`, `AGE/MOVE` with no value, or a value after `=`
  that its lexer does not read as a number (`AGE=abc`, `AGE/MOVE=1E5`,
  `AGE/MOVE=Inf`), stops the job with a syntax error in `PROC HAZARD`, and
  `ORDER=` with `/E`, `/I` or `/S` stops it with "mutually exclusive". The
  translator used to record the text and fit anyway, with the variable in
  the model (or, for `AGE=abc`, left out of it). Each now emits a `stop()`
  naming the source. To check a job: search its phase statements for a run
  of option letters after one `/` (`/EI`, `/SI`), a second `/`, `ORDER=`
  beside `/E`, `/I` or `/S`, or a value after `=` that is not a plain
  decimal number (an exponent needs a decimal point: `1.0E5`, not `1E5`).
  None of these forms occurs in the reference corpus. Options separated by
  spaces (`AGE/E I`) are valid SAS and still translate, with `/E` taking
  precedence.

* **A named `dist` could skip the check on phase-scoped formula terms
  (#405).** `hazard()` accepts a named scalar such as
  `dist = c(model = "multiphase")`, but tested it with `identical()`, which a
  named value never matches. `hazard()` does this in its own checks, and so
  does code that later reads `fit$spec$dist`, such as the stepwise refit and
  the diagnostics. So `hazard()` did not refuse a
  term such as `constant(age)` in the global formula (#275). If a function of
  that name was visible, the term became an ordinary covariate and entered
  every phase, a different model with no warning. `hazard()` now drops the
  names from `dist` before reading it, so `fit$spec$dist` is stored without
  them.

* **Standard errors were too small for a late (`"g3"`) phase with free
  shapes: re-run any you have reported (#332).** At realistic optima the
  standard errors this package reported for such a fit were **12 to 14
  times too small**, so confidence intervals were far too narrow and Wald
  p-values far too significant. Anyone who has published or acted on a
  standard error, a confidence interval or a Wald test from a fit with a
  free `"g3"` shape should re-run it. Near the exponential limit (a fitted
  `alpha` of 0.0021) they were about four times too small, and where the
  data leave a shape undetermined the fit reported finite standard errors
  for a direction the likelihood does not determine at all.

  **The fit itself is unchanged unless its search passed through a small
  `alpha`.** For every fit without a `"g3"` phase, every fit whose `"g3"`
  shapes are fixed, and every free-shape fit whose search never took `alpha`
  to `1e-5` or below, the estimates and the log-likelihood are identical.
  What moves is the Hessian and everything read from it -- standard errors,
  Wald statistics, confidence intervals, the condition warnings, and the
  score test `hzr_stepwise()` uses to enter a variable. Where the search did
  pass through that region the optimizer's own gradient was wrong, and the
  estimates can move too; see below.

  The cause was numerical: the G3 second derivatives stepped every shape by
  a fixed amount, about 1.2e-4, and the score stepped a small `alpha` by
  1e-5. That is far too small a fraction of a large `gamma`, where rounding
  takes over, and far too large a fraction of a small `eta` or `alpha`,
  where truncation does. Each shape is now stepped in proportion to itself.
  Nothing warned, and the individual Hessian entries were never off by more
  than 0.77% -- it is the inversion of an ill-conditioned matrix that turned
  that into an order of magnitude in the standard errors, which is why the
  entry error alone is not the number to judge this by.

  - **If a fit has a free `alpha`, refit it under this version and compare
    -- its final `alpha` does not tell you whether this applies.** At
    `0 < alpha <= 1e-5` the gradient in `alpha`, which the optimizer uses,
    was about 50% off at `alpha = 1e-5` and approached 100% as `alpha` fell.
    A search only has to pass through that region for the difference to
    steer it. On the package's own `avc` data a free-shape fit started at
    `gamma = 0.1` went below `alpha = 1e-5` on its way, finished with `alpha`
    between 0.015 and 0.02 -- far outside the region -- and still stopped at
    `gamma` 100.9 where it used to stop at 94.4.

    Both versions stopped short of an optimum there, and both said so: the
    fit warns that its estimates fail the relative-gradient test SAS/C
    HAZARD requires (0.121 before this change, 0.482 after, against a limit
    of 6.06e-06). **If your fit carries that warning, its estimates were not
    a reliable optimum before this release either** -- refit with more
    starts or other starting values rather than reading a change in them as
    an improvement. On one weakly identified data set, fits started at such
    an `alpha` stopped a full log-likelihood unit below what was attainable
    and reported `converged = TRUE`. Such
    likelihoods are often multimodal, so a corrected fit is not guaranteed
    to end higher from every start. At such an `alpha` the Hessian is now
    evaluated, and is usually too ill-conditioned to invert: standard
    errors are unavailable, with a warning.
* **A formula fit given `weights = <name>` could silently use the wrong
  weights: re-run any where that name was also a variable in your session
  (#392).** On the formula interface, `hazard()` did not look `weights` up
  in `data`. A name that was only a column of `data` failed with "object not
  found", but a name that was both a column and a variable in the calling
  frame (`wc <- d$wc`, a leftover from an earlier step) silently used the
  variable, not the column, with no error and no warning. On a 40-row
  example with a unit-weight `wc` beside the intended column, the
  log-likelihood was -34.23592435 instead of -33.53365357; the size of the
  error depends on how far the two vectors differ. `weights` is now looked
  up in `data` first, then the calling frame, as the vector interface (`time
  =`, `status =`) already did. `stats::lm()` also looks in `data` first, but
  then in the formula's environment, not the caller's. When a name is both a
  column and a visible variable, the column is used and `hazard()` now
  warns, naming it, on both interfaces. If you see that warning, the fit may
  differ from one made by an earlier version. A name that is only a column,
  or only a variable, does not warn. Nor does the namespace or function name
  in a qualified call such as `base::abs(w)`, which is never looked up in
  `data`; that false warning had been raised on the vector interface since
  1.2.2 (#151), and only the call's own arguments are now checked. The
  warning's advice now names the data frame the call passed, as in `d$w`.
  Since 1.2.2 the vector interface had advised `data$<name>`, which finds
  `utils::data()` rather than the data frame, so following it was an error.
  The warning reads the names written in the expression: a name chosen at
  run time, as in `get(nm)`, is resolved the same way, column first, as in
  `stats::lm()`, but is not checked. That is not new; the vector interface
  has behaved so since 1.2.2. When the name is part of a larger expression,
  such as `(function(a) 1)(w)` or `{ w <- 1; w }`, the warning no longer
  says the column was used, because the expression may never read the name,
  may rebind it first, or may evaluate it elsewhere, as `with()` does; it
  says instead that which value it read, if any, is not checked.

* **A masked argument with an empty index, such as `time = m[, 1]`, no
  longer fails when `data` is given.** Since 1.2.2, `hazard(time = m[, 1],
  status = s, data = d)` stopped with "invalid first argument": the check
  for names that are both a column and a caller variable tried to look up
  the empty index as a name. It now skips it, on the vector interface and on
  the formula interface's `weights = m[, 2]`.

* **A multiphase phase formula without an intercept no longer drops its
  first term (#303).** `hzr_phase(formula = ~ 0 + age)` or `~ age - 1`
  fitted the phase without `age`, and `~ 0 + age + mal` without `age`, with
  no warning and no message. `~ 0 + age + grp`, for a factor `grp`, lost
  `age` and kept a column for every level of `grp`. `predict()` repeated
  the same design, so its results agreed with the wrong fit. A phase has no
  free intercept of its own (its scale parameter plays that role), so such
  a formula now builds exactly the design of the same formula with an
  intercept: `~ 0 + age` fits as `~ age`, and factors are coded as they
  would be with the intercept present, not with a column per level.
  `hzr_phase()` now warns that the removal is ignored, once, when the phase
  is created. Refit a model whose phase formula had no intercept: its
  estimates will change. A formula that fits as before starts with an
  unordered factor, character or logical column under the default
  treatment contrasts, which already built the design of the formula with
  an intercept. An ordered
  factor, or any factor under other `contrasts`, is now coded as it would
  be with the intercept (for an ordered factor, `o.L` and `o.Q` rather
  than `om` and `oh`), so its coefficients change meaning. An interaction
  first, such as `~ 0 + g:age`, lost a column and now keeps it.

  The score test in `hzr_stepwise()` rebuilt the other phases of the model
  the old way, so once the fit was corrected they would have disagreed:
  a candidate's score was computed against a design the model was not
  fitted with, silently when the column counts matched, and otherwise every
  candidate in the other phases could not be scored. It now builds them as
  the fit does, and declines to score (a score of `NA`) when a phase the
  step does not change rebuilds with different columns than the fit
  stored, as a model saved by an earlier version with such a formula can.
* **A model formula without an intercept now builds the design of the same
  formula with one, with a warning (#337).** Every distribution carries its
  own intercept, its baseline parameter (the one the covariates add to).
  Without one in the formula,
  `Surv(time, dead) ~ 0 + grp` coded a column for every level of a factor
  `grp`, collinear with that parameter. A Weibull fit lost its standard errors
  (`Hessian not invertible`). A multiphase fit inheriting the design
  stopped 9.5 log-likelihood units below the fit with the intercept, with
  warnings but wrong estimates. The formula now fits as
  `Surv(time, dead) ~ grp`, and factors are coded as they would be with the
  intercept. The design therefore has one column fewer, from the first
  factor term: without an intercept R codes only the first factor with
  every level, and later factors keep their usual coding, so
  `~ 0 + f1 + f2` had `f1a f1b f1c f2y` and now has `f1b f1c f2y`. In a
  single-distribution fit that is one coefficient fewer in `theta`. In a
  multiphase fit, every phase that inherits the global design loses that
  column, so `theta` is one entry shorter for each such phase; a phase
  with its own `formula` is unaffected. To reuse a `theta` supplied for the
  old design, drop the first factor's first-level coefficient from the
  global design, and from each phase that inherits it. Left unchanged,
  the fit stops: a single-distribution fit with a `non-conformable
  arguments` error, and a multiphase fit with an error naming both
  lengths and each phase's count (#408). A numeric term fits as before
  (`~ 0 + age` is `~ age`), apart from the new warning, which a
  `hzr_stepwise()` refit does not repeat. The phase-formula counterpart is
  #303.

* **`hzr_stepwise()` and `hzr_bootstrap()` now refuse a `scope` under
  `direction = "backward"` (#343).** A backward screen only drops terms the
  base model already has, and it never read `scope`. From `~ age + mal`,
  `scope = ~ age` still dropped `mal`, the variable left out of the scope,
  and a scope variable the base lacked was never tested. The result was the
  same as with no scope, and nothing said so. Such a call is now an error:
  pass the full model as the base and protect terms with `force_in`. In
  `hzr_stepwise()`, leave `scope` unset or empty (`~ 1`), since an empty
  scope offers nothing to enter and so agrees with a backward screen. In
  `hzr_bootstrap()`, an unset `scope` means no screen at all, so pass an
  empty scope such as `~ 1` with `direction = "backward"` to run a backward
  screen on each replicate. Under `direction = "both"`, `scope` names what
  may enter; as in SAS, the drop half still considers every term in the
  model except those in `force_in` and those frozen by `max_move` before
  the iteration began. `hzr_bootstrap()` refuses the combination before
  seeding.

* **`hzr_bootstrap()` now refuses a selection argument passed without
  `scope` (#343).** Without `scope` there is no screen, and `direction`,
  `criterion`, `slentry`, `slstay`, `max_steps`, `max_move`, `force_in` and
  `force_out` were ignored: `direction = "backward", force_in = "age"`
  returned a fixed-model bootstrap with every term at `pct = 100` and no
  message. A value other than the argument's default is now an error, so a
  wrapper that passes the defaults on still works. Pass `scope` to screen,
  or omit the argument to refit the exact model. A value equal to the
  default is accepted whether or not it was passed: without `scope` it asks
  for nothing, and nothing reads it.

* **A two-sided `scope` formula is now an error in `hzr_stepwise()` and
  `hzr_bootstrap()` (#343).** Only the right-hand side was read, so
  `scope = com_iv ~ age + mal` screened `age` and `mal` and never tested
  `com_iv`, with no message. The error names the left-hand side. The same
  applies to each element of a multiphase `scope` list, and a list that
  names a phase twice, whose second entry was never read, is refused too.

* **`hzr_stepwise()` and `hzr_bootstrap()` now check `...` against what a
  candidate refit may take from it (#386).** Both pass `...` on to every
  candidate refit, and `hazard()` stores a name it does not declare without
  reading it. So a misspelled argument had no effect: `slentyr = 1e-6`,
  meant as `slentry`, ran the screen at the default `slentry = 0.30` and
  selected three variables where the intended threshold selects one, with
  no message. Other `hazard()` arguments were no safer: a refit takes the
  response and `dist` from the base model, so `time_lower` was ignored on a
  formula fit and `dist` failed every candidate, and `weights` or
  `time_windows` changed the likelihood of the candidates but not of the
  base model they were compared with. Each is now an error naming the
  argument and, for a misspelling, the argument it was probably meant to
  be. `...` forwards `control`, including by an unambiguous abbreviation
  such as `contr`, and an `objective` equal to the base fit's. A differing
  `objective`, a `control` that is not a list, a repeated argument and an
  ambiguous abbreviation are refused at entry rather than failing every
  candidate. `hzr_bootstrap()` checks before seeding, so a refusal leaves
  the random number stream untouched and names `hzr_bootstrap()`. `?hzr_stepwise` had also described `...`
  as unused, because the print method's entry replaced it.

* **`hazard()` now refuses a multiphase phase formula with covariates when
  no `data` is supplied (#299).** Such fits previously ignored the phase
  formula. On the vector interface (`time =`, `status =`) without `data`, a
  phase's own formula was never evaluated: the phase took the global `x`, or
  no covariates at all, so `hzr_phase(formula = ~ mal)` fitted a model
  without `mal`, with no warning and no message. The call now stops with an
  error naming the phase and its formula, under `fit = FALSE` as well. Pass
  `data =` with the phase's variables as columns, or use
  `hazard(Surv(...) ~ ..., data = ...)`. A phase with no formula is
  unaffected, and so is an intercept-only `~ 1` unless the call also has a
  global `x` with at least one column. Beside such an `x`, a `~ 1` phase
  is refused too: without `data` it silently took `x`, and with `data` it
  has no columns, so the same call gave two different models. Pass
  `data =`, use the formula interface, or drop `x`. A constant term such
  as `~ log(2)`, which builds a column only in `data`, is refused as well.

* **`hzr_stepwise()` now refuses a fit saved by an earlier version whose
  phase formula was ignored this way (#299).** Given `data`, every refit
  built the ignored formula's columns into a model whose base never had
  them: a forward screen reported `ENTER age` over a final model that also
  carried the ignored `mal`, with no warning. The error names the phase and
  its formula. Refit the base model with `data =` and retry.

* **`hzr_stepwise()` and `hzr_bootstrap()` now also refuse a multiphase fit
  saved before 1.1.0 when they cannot tell whether its phase formula was
  used (#324).** Such an object stores neither its data frame nor a record
  of which phases used their formulas. On the vector interface, a call that
  names `data = dd` reads the same whether `dd` was a data frame or `NULL`,
  and a `NULL` there meant the phase formula was ignored. A fit with a phase
  whose stored columns are the ones it would have inherited was let through,
  and a screen then credited the ignored formula's columns to the candidate
  it was testing, with no warning. Such fits are now refused, naming the
  phase and its formula. This also refuses a fit that did use a phase
  formula whose columns are exactly the inherited ones: without the record
  the two cannot be told apart. The check turns on the stored data frame,
  not the version: a fit that kept `data$frame` (every fit from 1.1.0 on,
  unless that element was removed) is unaffected, and so is a call without
  `time =`. A call with `time =` is read as the vector interface even when
  it also passes a formula, so a wrapper's `formula = fml` beside it is
  judged here too. Refit the base model with the current version, passing
  `data =`, and retry.

* **`hazard()` now stops when two design columns share a name** (#298). A
  factor's dummy columns are named `<factor><level>`, so a factor `g` with
  level `b` and a numeric column `gb` both produced a column `gb`. The fit
  ran without a word: `coef()` carried two `gb` names, and `predict()` on a
  one-row `newdata` returned two values. The check covers the global design,
  an `x` matrix passed to the vector interface, the design after
  `time_windows` expansion, and, when `fit = TRUE`, each phase formula of a
  multiphase fit. The error names the colliding columns. Unnamed columns of
  `x` are allowed, but not under `time_windows`: the expansion names each
  window's column `<name>_w<k>`, so two unnamed columns both became `_w1`.
  A fit that used to run now stops: rename the numeric column, or rename the
  factor or change its levels (`relevel()`, `levels<-`).

* **An `offset()` term in a formula is now an error.** `hazard()` used to
  drop it without a word: `model.matrix()` leaves offsets out of the design
  and nothing read them back, so `Surv(time, status) ~ age + offset(z)`
  gave the same log-likelihood and the same coefficients as
  `Surv(time, status) ~ age`. **Fits written with an offset ignored it**
  (#297).
  The global formula, every `hzr_phase(formula = )`, and a
  `hzr_stepwise()` formula or scope now stop and name the offending term.
  Offsets are not supported: how one should enter each phase of the
  additive multiphase hazard is an open modelling question, not yet
  decided. `stats::offset(z)` is not refused, because R reads it as an
  ordinary covariate, not an offset, and fits a coefficient for it.

* **`predict(newdata = )` now takes only the columns of the model's `data`
  from `newdata`.** A term that uses row-level values kept outside `data`
  (a vector, matrix, list or environment in the formula's environment, as
  in `~ zz` or `~ ext$z`) is refused, even when `newdata` supplies the
  object, with an error naming the term (`term 'zz' of the model ...`)
  (#409). Such a term cannot be rebuilt for new rows. Move the variable
  into `data` as a column and refit. When the term's values simply do not
  line up with `newdata`'s rows, the error says so and gives both causes,
  since a length-changing function of a `data` column, such as
  `I(unique(age))`, reaches the same check. An error raised by your own code
  inside a term reaches you unchanged, with its own class and message, even
  when a different term is row-mismatched -- unless it comes from a
  `model.frame()` or `model.matrix()` call inside that term, which is read
  as the design build's own failure and replaced by the naming error.
  Before, a supplied `zz` or matrix `M` was used, but a missing or
  list-held one was silently read from the fitting rows (see Bug fixes).
  Formula constants, such as `cutoff` in
  `I(age > cutoff)` and spline knots, are unaffected. A fit saved by an
  earlier version without its data still takes such a variable from
  `newdata`, since it cannot tell it from a column.

* **`hazard(fit = TRUE)` without `theta` is now an error for the
  single-distribution models.** For `dist = "weibull"`, `"exponential"`,
  `"loglogistic"` and `"lognormal"`, the optimizer ran only when `theta` was
  supplied, so a call that left it out returned an unfitted object -- `NULL`
  coefficients, an `NA` objective -- with no error and no warning. `print()`
  gave no sign of it. The call now stops and asks for starting values, as it
  does for a zero-length `theta`, which used to fail inside `optim()` with a
  message that did not name the cause.
  `dist = "multiphase"` is unaffected: it assembles its own start from
  `phases`. `fit = FALSE` without `theta` still builds an unfitted model.

  Code that relied on the old behaviour was getting no fit. Two tests in this
  package were: they compared `NULL` coefficients with `NULL` coefficients,
  so the counting-process equivalence and epoch-split invariance they claimed
  for the Weibull were never checked. Both now fit, and both pass.

* **`hzr_translate_sas()` no longer drops phase covariates that follow a
  `/` option (#342).** SAS attaches `/ options` to one covariate at a
  time, so `EARLY AGE, MAL/I, OPMOS;` is three covariates. The translator
  cut the list at the first `/`, fitted `~AGE + MAL`, and recorded only
  that phase options were deferred, not that `OPMOS` was gone. Each
  covariate now keeps its own options: `/E` (`EXCLUDE`) leaves it out of
  the model, as `PROC HAZARD` does without a `SELECTION` statement; `/I`
  and `/S` leave it in; and a per-variable `MOVE=` or `ORDER=`, or any other
  option, is recorded in `$untranslated` under the
  variable's name.

  Two related fixes. A second `EARLY`, `CONSTANT` or `LATE` statement for
  the same phase now adds to that phase's covariates rather than replacing
  them. A covariate named twice in a phase is one parameter, as in
  `PROC HAZARD`, whose last mention sets its starting value and options;
  within one statement it used to put two starting values in `theta` for
  one column, shifting every later one.

  A `SELECTION NOSTEPWISE` (or `NOSW`) job is no longer read as no screen
  at all. It was translated to a plain fit with every candidate in the
  model, but `PROC HAZARD` still screens, forward only, with each candidate
  starting out of the model; it now translates as a forward-only screen
  (see the `SELECTION` entry under New features).

  A phase variable that is not in the fitted model (an `/E`
  variable, or a covariate of a phase the job does not select) still
  deletes its missing rows in `PROC HAZARD`, which `hazard()` cannot do
  for a variable it never sees, so the translated status chunk now stops
  when such a variable is missing and asks for those rows to be dropped.

* **`hzr_translate_sas()` now emits a `stop()` in place of the fit when a
  `PARMS` statement builds no phase it could use.** Operands the translator
  could not read (a template's `MUE=?`, or `MUE = 0.2` written with spaces
  around `=`, which `PROC HAZARD` accepts) are recorded in `$untranslated`,
  but the fit chunk
  used to be emitted anyway, as `hazard(fit = TRUE, theta = c())` under the
  default Weibull. That chunk rendered an unfitted object, and would now fail
  on the error above with a message about `theta` that does not name the real
  cause. The emitted `stop()` names it. This is a limit of the translation,
  not a `PROC HAZARD` refusal, so it is kept apart from the existing
  "selects no phase" stop.

* **`hzr_translate_sas()` no longer writes `CONDITION=` or `QUASI` into
  `hazard()`'s `control` (#384).** They were emitted as `condition` and
  `method`, which `hazard()` never reads, so the translation counted two
  options as mapped while they did nothing. No fit changes. Both are now
  recorded in `$untranslated` with the reason. `CONDITION=` stops
  `PROC HAZARD`'s optimizer when its Hessian approximation becomes too
  ill-conditioned, and `hazard()` has no such stop; it warns about the final
  Hessian instead; a `CONDITION=` outside 3 to 14, which `PROC HAZARD`
  itself ignores, is recorded as such. `QUASI` chooses `PROC HAZARD`'s
  optimizer, and `hazard()` has no choice to make: it fits by BFGS, a
  quasi-Newton method, after a Nelder-Mead warm-up in some multiphase fits.

* **`hzr_translate_sas()` now emits a `stop()` in place of the fit when a
  job has no `DATA=` and a phase has covariates (#311).** A phase's
  covariates are evaluated only in `data`, and such a job's fit chunk has
  none, so `hazard()` stopped with advice to pass `data =`, an argument the
  SAS job never had. Before the refusal above (#299) the same chunk fitted
  the phase without its covariates. The job is now recorded in
  `$untranslated` with the reason, and the emitted `stop()` says to add
  `DATA=` and translate again. A job with no `DATA=` and no phase
  covariates still translates to a fit, unless it has a `SELECTION`
  statement (see the next entry).

* **The no-`DATA=` refusal (#311) counts every variable a phase statement
  names, not only the covariates of the base model (#160).** A `SELECTION`
  job withholds its candidates from the base model, so a refusal reading only
  the base model could not see them. When the job's `SELECTION` could not be
  run anyway (`FAST`, say), that reason is recorded in `$untranslated` beside
  the `DATA=` one, so adding `DATA=` does not reveal a second refusal. A
  `SELECTION` job with no `DATA=` is refused even if its phase statements
  name no variable, because a screen refits every candidate from `data`. **This
  widens #311 on purpose.** A job with no `DATA=` is now refused, with or
  without `SELECTION`, when its only phase variables are any of these:
  - excluded with `/E`;
  - on a `LATE` statement when `PARMS` has no `MUL`;
  - on a `CONSTANT` statement when `PARMS` has no `MUC`;
  - on an `EARLY` statement when `PARMS` has no `MUE`.

  **Such jobs used to translate.** Their missing-value guard then read those
  variables from whatever environment rendered the document. The refusal
  names the variables rather than claiming a phase has covariates.

* **`predict(newdata = )` matches covariates by name, so `newdata` with
  other names now stops.** A fit made through the vector interface with a
  named `x`, say `x = cbind(age = , mal = )`, needs `newdata` columns
  called `age` and `mal`. Before, `predict()` matched any names, or a bare
  matrix, to the coefficients by position. That was right only when the
  order happened to agree, and nothing said when it did not (#267). Code
  that passed such `newdata` now gets an error naming the missing columns:
  rename the columns to match `x`. A fit made with an unnamed `x` still
  matches by position. A formula fit saved by an earlier version stored no
  formula design; `predict()` rebuilds it (see below), and a fit whose
  design it cannot rebuild is matched on its design-matrix columns: a
  factor must be given
  as `grpyoung` and a transform as `log(age)`. Refit it to give the
  formula's variables instead.

* **A formula fit saved by an earlier version is matched by name, as a new
  fit is.** Such a fit stored no formula design, so `predict(newdata = )`
  rebuilds it from the fit's stored formula and data frame, and uses it
  only if it reproduces the fitted design matrix exactly: the same columns
  with the same values. The fit then takes the formula's variables
  (`grp = "old"`), ignores unused columns, and lets a variable win over a
  contradicting design column, all as a new fit does (#301). One value
  changes with it: `newdata` giving a numeric variable and a column built
  from it, such as `age` and `I(age^2)` at the design-column means, is now
  rebuilt from `age`, as for a new fit, where its columns used to be taken
  as given. The design is not rebuilt for a 1.0.3-era fit, which kept no
  data frame; for a formula that no longer reproduces the fit; for one
  computed in the call (`as.formula(...)`), which would be re-run; and for
  a formula that uses anything but the data's columns and a short list of
  R's own design functions (arithmetic and comparisons, `I()`, `log()`,
  `exp()`, `sqrt()`, `abs()`, `pmin()`, `pmax()`, `c()`, `factor()`,
  `relevel()`, `scale()`, `poly()`, `splines::ns()` and `splines::bs()`).
  A constant
  such as `k` in `I(age > k)`, a function of the user's, or even `pi`
  could have changed since the fit without changing the fitted rows, so it
  is not trusted; a number written into the formula, as in
  `I(age > 50)`, is. For these, `newdata` with columns other than the
  design columns and `time` is refused, with
  `This fit was saved by an earlier version of TemporalHazard, without a
  stored formula design, ...; it also has '...'. Refit the model with the
  current version, or pass only the design columns.` Without a design
  nothing can tell such a column from a formula variable that contradicts
  a design column, which would otherwise be ignored silently. Before, such
  a column made the positional match fail, so this is as loud as it was
  (#272).

* **`predict(newdata = )` stops when `newdata` gives some of the formula's
  variables beside the fitted design columns, with others missing.** Say
  a fit of `~ age + grp + sex`, with `newdata` holding `grp = "old"`,
  `grpyoung = 1` and `sexM = 0` but no `sex`. That used to return a value,
  and silently the wrong one: the design columns were used, and the
  `grp` it was given was ignored. It now stops with the error
  `'newdata' gives the formula variable(s) 'grp' but lacks 'sex', while
  carrying the fitted design columns. Give all of the formula's variables,
  so the design can be rebuilt from them.` Pass all of the formula's
  variables (here `age`, `grp` and `sex`), or the design columns alone
  (`age`, `grpyoung`, `sexM`). A numeric covariate such as `age` is both a
  variable and a column, so it makes a mix only when another column is
  built from it, as in `~ age * grp` or `I(age^2)`: a changed `age` would
  leave those columns stale. The
  wrong-answer fix itself is under Bug fixes (#272).

* **`predict(newdata = )` stops for the time-based predictions of a model
  with a covariate named `time`.** In `newdata` the column `time` is the
  prediction time for `"survival"`, `"cumulative_hazard"`, every multiphase
  type and any fit with `time_windows`, so such a covariate could not be
  given its own value. When it was the only covariate it was dropped, and
  the prediction silently came back at the baseline: a Weibull fit of
  `~ time` gave 0.1414 where the covariate made it 0.1420. Beside other
  covariates, the one column served as both, so the covariate was always
  set to the prediction time. A formula constant named `time`, as in
  `I(age > time)`, was replaced the same way, because a formula looks its
  symbols up in `newdata` first. None of this gave an error. These calls
  now stop and ask for the variable to be renamed and the model refitted,
  whether it is in the global formula, a named `x` or a multiphase phase
  formula. For a single-distribution model, `hzr_gof()` and
  `hzr_deciles()` pass the fitted design columns, which are used as they
  are, so they stop only for a design column named `time` itself, which
  follow-up time used to overwrite silently. For a multiphase fit they
  use the fitted per-phase designs and re-evaluate no formula, so a
  `time` variable does not stop them there. Nor does a phase formula that
  the fit did not use: a vector-interface fit ignores one. `log(time)`, a constant such as
  `I(age > time)` and a list element such as `cfg$time` do not stop them.
  In the global formula, neither does a value that `scale()` stored at
  fit time, in `predict()` either. A phase formula is re-evaluated as
  written, so a `time` constant there stops even inside `scale()`.
  `"linear_predictor"` and single-distribution `"hazard"` have no
  prediction time, so they read a `time` column as the covariate, now
  also when it is the only one (it used to stop with "Predictors are
  required"). They refuse only a `time` constant that a `time` column in
  `newdata` would mask. `predict()` without `newdata` is unaffected
  (#270).

* **A multiphase formula that names a phase as a function is now an error**
  (#275). `hazard(Surv(int_dead, dead) ~ constant(age), dist = "multiphase",
  phases = ...)` read `constant(age)` as a phase-scoped term, then replaced
  the whole right-hand side with `~ 1`, and nothing sent the term to its
  phase. The fit converged without an `age` coefficient, with no error and
  no warning. `hazard()` now stops, names the phase or phases it found, and points to
  `hzr_phase(..., formula = ~ var)`, which is where a phase's covariates
  belong. Code that relied on the old behaviour was fitting a model without
  those covariates; drop the terms from the formula to keep that model, or
  move them into `hzr_phase()` to get the one the formula described.
  Formulas that call an ordinary function such as `log()`, and plain global
  covariates, are unaffected.

* **`hzr_bootstrap()` now refuses a fit whose formula uses a per-row
  variable that is not a column of its `data`** (#278). Replicates resample
  the rows of `data`, so such a variable was held fixed while the rows moved
  under it. The interval was wrong, and every replicate still reported
  success, with no warning: on `avc`, a copy of `age` kept outside `data`
  gave an interval that excluded its own estimate. The check covers the
  response, the covariates of the global and phase formulas, and a
  select-mode `scope`, and reads variables from the terms, so `log(age)`
  needs only the column `age`. A constant outside `data`, such as `pi`, a
  cutoff or a knots vector, is still allowed. The error names the
  variables; add them to `data` and refit. A scope variable that only the
  scope formula's own frame can see, which the refits could never test, is
  refused the same way. So is a vector-interface fit whose design matrix
  was passed directly as `x`: it was re-evaluated without resampling in
  every replicate, or in select mode dropped from the candidate refits, so
  all of them succeeded and the interval was wrong.

* **`hzr_stepwise()` now refuses a base fit written as `Surv(...) ~ .`**
  (#279). It read the base model's terms without the data, which cannot
  expand `.`, and treated the failure as a model with no terms: the screen
  reported zero steps, which looks the same as finding nothing to drop. It
  now stops before printing anything and asks for the base model's terms to
  be written out. A `scope` of `~ .`, which used to give a screen with no
  candidates, stops the same way, and so does a multiphase base fit whose
  global formula uses `.` while a phase has no formula of its own and so
  inherits it. A screen whose base model has its terms written out is
  unchanged, as is a multiphase screen in which every phase has its own
  formula.

* **`hzr_stepwise()` now refuses a vector-interface multiphase fit whose
  phases inherit a design matrix passed directly as `x`** (#284). A phase
  with no formula of its own uses that `x`, and a candidate refit had no
  terms to rebuild it from, so every refit dropped its columns from every
  such phase: on `avc`, entering `mal` into one phase took `age` out of
  both, and the log-likelihood fell from -196.44 to -202.17 while the step
  table reported a plain "enter". The only warning blamed a refit that "did
  not converge". The screen now stops before printing anything. Give each
  phase its covariates with `hzr_phase(formula = ~ ...)`, with the columns
  in `data`; a fit whose phases all have their own formulas is unaffected.
  It refuses the same way when a phase the screen can step inherits
  time-varying coefficients (`time_windows`), which rebuilding the phase
  from its terms would fit as one constant effect, or a term that expands
  to more than one column, such as a factor with more than two levels or
  `poly()`, which a step cannot add or drop as one coefficient. Both used
  to change the phase's design without a word. A forward screen steps only
  the phases its `scope` names, so an inheriting phase outside the scope
  is left alone; a backward or two-way screen can drop from any phase.

* **An entry time after the exit time is now an error, and
  `time_lower = time` now means "no entry" in every family** (#253). On a
  row with status 0 or 1, `time_lower` is the counting-process entry time
  when `0 < time_lower < time`. `hazard()` used to warn when
  `time_lower >= time` on such a row and fit anyway, and each family then
  did something different. Multiphase returned a "log-likelihood" of
  +47915.76 with `converged = TRUE` on the AVC data. The Weibull read the
  rows as entering at time 0, and the other three families ignored
  `time_lower` altogether.
  - `time_lower > time` on a status 0/1 row now stops, as SAS HAZARD
    rejects a start time after the exit time (error `SETCOE960`).
  - `time_lower == time` on a status 0/1 row is read as no entry time, with
    no warning, in all five families. This is the mixed-interval layout,
    where exact and right-censored rows carry `time_lower = time` and only
    interval-censored rows (status 2) carry a real lower bound. It was
    already the Weibull rule; multiphase used to degenerate on it.
  - `time_lower = 0` still means no entry time.
  - Rows with `time_lower == time > 0` beside rows with a genuine entry
    time now stop. In counting-process data they are zero-length epochs,
    which `hzr_repeated_events()` can emit; read as "no entry", each would
    be charged its full cumulative hazard from time 0.

* **`"time"` is now a reserved phase name, like `"total"`** (#224). The
  decomposed output of `predict(decompose = TRUE)` starts with a `time`
  column holding the prediction times, then adds one column per phase under
  the phase's name. A phase called `time` therefore overwrote the requested
  times with its own cumulative hazard, and the times were lost with no
  error. `hazard()` and `hzr_theta_names()` now stop on a phase named
  `time`, before any fitting, and ask for a different name. Rename the
  phase; nothing else about the model changes.

* **`hazard()` now stops on zero observations** (#231). Given a `time` of
  length 0, or a formula whose `data` has no rows, every distribution
  returned a `hazard` object anyway, and with `fit = TRUE` it reported
  `converged = TRUE`. The warnings were about the Hessian (`rcond = 0`, not
  invertible), which read as a conditioning problem rather than as no data.
  The call now errors before any fitting, under `fit = FALSE` as well. The
  same holds when no row contributes to the likelihood: every row has
  weight 0, or is right-censored at time 0, where the cumulative hazard is
  0. Such a fit came back converged at its starting values with an
  objective of 0. Check that the data frame, or the subset passed to
  `data`, has rows, and that some row with positive weight is an event or
  is followed past time 0.

* **`hazard()` now stops on a status code other than -1, 0, 1 or 2**
  (#231). Every likelihood branches on those four codes, so a row coded
  anything else fell through all of them and contributed nothing, with no
  warning. The likeliest way in was `survival::Surv(type = "interval")`'s
  own codes passed as a plain vector, where 3 means interval-censored: those
  rows were silently dropped, and data coded only that way returned its
  starting values with `converged = TRUE`. The error names the rows. Pass
  a `Surv` object as the response, or as `status`, and it is translated.
  A character or factor `status` is refused too: it passed as text, and the
  exponential, Weibull, lognormal and log-logistic fits then returned their
  starting values as a converged fit. A logical `status` is still accepted.
  A classed numeric such as `bit64::integer64`, which `data.table::fread()`
  and `arrow` produce, is now read as its values in `time`, `status`,
  `time_lower`, `time_upper` and `weights`, and as a column of `data`,
  where `Surv()` and the model formulas read it. Before, those fits read
  its stored bits and returned their starting values as converged.

* **A forward candidate whose refit adds no design column of its own is now
  a recorded refit failure, not an error out of `hzr_stepwise()` (#442).**
  Under `"wald"` and `"aic"`, `hzr_stepwise()` used to stop with an error
  ("added no design-matrix column", or "does not add a column" when the
  refit only changes the parameterisation). It now issues a warning, adds the
  candidate to `$criteria$refit_failures` with the refusal as its reason in
  `$criteria$refit_failure_reasons`, and the screen goes on, as `"score"`
  already did. Code that caught the error with `tryCatch(..., error = )`
  will no longer see it; read `$criteria$refit_failures` instead.

* **`hzr_stepwise()`'s `$steps$variable` records the model's term label on
  every row (#449).** An entry row used to carry the name as the `scope`
  wrote it and a drop row the `terms()` label, so a non-syntactic column
  `_X1` entered as `_X1` and left as `` `_X1` ``, and a literal column
  `age:mal` entered under the interaction's spelling `age:mal`. Every
  row now uses the label, whatever form the `scope` took. For a syntactic
  name the label is the name, so nothing changes; code that matched an
  entry row of a non-syntactic column by its bare name should match the
  backquoted label instead.

## New features

* **`hzr_translate_sas()` now says when `PROC HAZARD` rewrote a shape operand
  before fitting, instead of emitting the rewritten value silently.** Under
  `FIXGE2` or `FIXGAE2` with `WEIBULL`, `SETG3` moves the late shape onto the
  constraint before the fit (`setg3.c:449-467, :827`), so a job written
  `ALPHA=2 GAMMA=5 ETA=1 FIXGAE2 WEIBULL` is fitted by `PROC HAZARD` at
  `ALPHA=2.5`, not at the 2 on the statement. The translation already emitted
  `alpha = 2.5`, the model `PROC HAZARD` fits, but said nothing, so a reader
  comparing the emitted call against the job saw a value they had not written
  and no reason for it. Such a rewrite is now recorded, naming the operand and
  both values (`ALPHA=2 -> 2.5`).

  **No fit changes.** The emitted call is the same on both sides; what is new
  is the row and the "untranslated construct(s)" warning that goes with it.
  A job that translated cleanly and reported no rows may now report one.

* **A fit now records why its gradient test was not run, not merely that it
  was not** (#351). SAS/C HAZARD accepts an optimum only when the relative
  gradient is small enough, and every fit reports that test in
  `fit$fit$rel_gradient`. `NA` there has always meant "not evaluated", never
  a pass -- but it did not say why, and "not evaluated at the estimates" on
  its own reads like a failure the fit is declining to name. It usually is
  not one. Under Conservation of Events the test is computed by
  differencing the log-likelihood, so a point the difference needs can fall
  outside the region where the likelihood is finite while the estimates
  themselves are sound. Reading that as a failure would condemn a good fit.
  The reason is now recorded in `fit$fit$rel_gradient_reason` --
  `NA_character_` when the test did run -- and `print()` and `summary()`
  append it, so the routes to a missing result are told apart from each
  other and from a test that ran and failed.
  A test that ran still reports "met" or "not met" exactly as before.

* **`hzr_translate_sas()` now translates a `SELECTION` statement into an
  `hzr_stepwise()` call** (#160). Such a job used to emit a `stop()`: the
  refit path needed a formula-interface base fit, so every candidate refit
  would have failed and the screen would have reported zero steps, which
  reads exactly like "nothing met `slentry`". The refit is phase-aware now,
  so the job translates into two chunks, the shape-fixed base fit and the
  screen, carrying the job's own candidates, per-variable flags and
  thresholds: a bare phase variable is a candidate offered through `scope`
  and withheld from the base model, `/S` starts in the model, `/I` becomes
  `force_in`, and `/E` appears nowhere. `PROC HAZARD`'s defaults are always
  written out (`SLE` 0.3, `SLS` 0.2, or 0.05 under `BACKWARD`), so
  the call never inherits a different default from `hzr_stepwise()`. A
  `BACKWARD` job gets no `scope` and a base carrying every candidate, which
  is where `PROC HAZARD` starts one.

  **The screen may select a different model than `PROC HAZARD` did**, and
  the rendered document says so in a callout above the chunk: `PROC HAZARD`
  uses approximate variances during selection, which the entry statistic
  here reproduces (except for a candidate refitted because its information
  is indefinite) but the Wald removal tests do not, and `force_in` is keyed by variable name across phases where
  SAS's `/I` holds a variable in one phase. Read the result as this
  package's screen of the job's candidates, not as a reproduction of the SAS
  run.

  **`ROBUST` and `SEMIROBUST` translate: they choose an optimizer, not a
  variance.** In `PROC HAZARD` they select the algorithm for the stepwise
  step (quasi-Newton, started by steepest descent or from the Hessian), not
  the variance. The option is recorded in `$untranslated`, and the screen
  uses this package's own optimizer. A different optimizer takes a different
  path, and on a multimodal likelihood it can reach a different optimum. They appear on 90.5% of the `SELECTION`
  statements in the production corpus, so this is the common case.

  **What is refused, so you can tell in advance which of your jobs are
  covered.** A `SELECTION` this translator cannot run faithfully emits a
  `stop()` rather than a screen: `FAST` (a different search), `MAXVARS`
  (caps the selected set), `RESTRICT` (constrains which variables may be
  selected), a per-variable `MOVE=` or `ORDER=`, and a variable held by
  `/I` in one phase but movable in another (`force_in` is not phase-keyed,
  so it would be pinned in both). On the reference corpus **2 of 4
  `SELECTION` jobs translate**; the two refusals are a cross-phase `/I` and
  a `RESTRICT` statement.

  **The screen can re-enter a variable `PROC HAZARD` would keep out.**
  `PROC HAZARD`'s `MOVE` limit counts a variable's *deletions*, separately
  for each phase, and at its default of 1 a variable removed from a phase
  can never return to it. `hzr_stepwise()`'s `max_move` counts entries and
  exits together across every phase and lets a removed variable re-enter.
  The two are not the same quantity, so the emitted call carries no
  `max_move`, `MOVE=` is recorded in `$untranslated`, and the callout names
  the difference. On the one reference job that translates, an unbounded
  screen re-entered five variables `PROC HAZARD` would have kept out.
* **`hzr_evaluate()` evaluates a model at parameters you supply (#144).**
  A parity check needs the likelihood at another program's converged
  estimates, evaluated by this package's own likelihood, and there was no
  way to ask for it. `hzr_evaluate(object, theta)` returns the
  log-likelihood of the model's data at `theta`, for a fitted model or one
  built with `fit = FALSE`, and with `times` (multiphase only) the hazard
  and cumulative hazard there for a covariate-free subject. The result is
  not a fit and does not read as one: it carries no standard errors, no
  convergence status and no covariance, and `print()` says so on its first
  line. A phase built with `hzr_phase(constraint = )` has its derived shape
  re-derived here, as the fit re-derives it, so a contradictory value passed
  in `theta` is replaced rather than used as given. At a fitted model's own
  estimates it returns that fit's objective, except where the fit reports an
  objective it is not at: under Conservation of Events the conserved scale is
  re-solved after the objective is recorded (#362), and the two then differ.

* **`hzr_phase()` can derive one late-phase shape from the others (#325).**
  The new `constraint` argument covers SAS/C's two late-phase constraints:
  - `"alpha_gamma_eta"` holds `alpha = gamma * eta / 2` (`FIXGAE2`);
  - `"eta_gamma"` holds `eta = 2 / gamma` (`FIXGE2`).

  The derived shape is recomputed from the others at every step of the fit.
  It is not estimated, and it cannot be fixed. A value supplied for it,
  through `hzr_phase()` or in `hazard(theta = )`, is replaced with a warning
  when it differs. Its standard error is the delta-method one, so `predict()`
  confidence limits carry its uncertainty; `summary()` shows it but does not
  test it against zero.

  Before this, `hzr_translate_sas()` recorded both flags as untranslated but
  still emitted a runnable fit, which estimated the derived shape freely.
  That is a different model. On a production `FIXGAE2` job it came out 2
  log-likelihood units above the SAS fit, with shapes 2 to 4 times SAS's and
  a Hessian that was not positive definite. The translator now maps both
  flags onto `constraint` for a `WEIBULL` late phase, following the rules in
  `setg3.c`. That job now reproduces PROC HAZARD's log-likelihood and
  estimates to about 2e-5. Other combinations of the flags stay recorded as
  untranslated, with the reason:
  - either flag without `WEIBULL`;
  - both flags together;
  - a combination PROC HAZARD refuses (`SETG3990`, `SETG31000`).

* **Every fit now says what it did not do** (#242, following #197). A
  `hazard` object carries `degraded`, the steps the fit did not perform, and
  `degraded_causes`, the reason for each. `print()` and `summary()` always
  show them as a "Not done in this run" block, and the block reads "none"
  when nothing was lost: a line that appears only on bad news cannot be told
  from one that was never written. Five steps are recorded:
  - fitting itself: `fit = FALSE`, or a fit imported from SAS output;
  - standard errors, naming whether numDeriv was missing, `numDeriv::hessian()`
    failed, or the Hessian was non-finite or singular, and, when a covariance
    was computed, which estimated parameters were left without a usable
    variance;
  - the variance of the phase that Conservation of Events conserves;
  - the weak-direction check;
  - Conservation of Events itself.

  The block replaces two notes that could name the wrong cause: "not examined
  for a weakly identified direction" and "standard errors unavailable; the
  Hessian could not be inverted". `fit$fit$weak` keeps its meaning, and is
  `NA` exactly when `"weak_direction_check"` is listed. An object saved by an
  earlier version prints "not recorded" rather than "none".

* **`hzr_translate_sas()` translates a `%repeat` call** into
  `hzr_repeated_events()` (#241), renaming its outputs to the names the job
  gives them so the job's fit reads them. The macro's input is still the
  reader's to supply. Any step between the macro and the fit that names the
  macro's output, or uses a macro variable that might, stops the document
  with the step quoted, rather than fitting data the job changed; a plain
  `PROC SORT` is the one step let through. A step that changes the output
  without naming it, such as a macro that writes it internally, is not
  detected.

## Bug fixes

* **A `"hazard"` phase fitted outside your data is now reported (#444).** The
  `"hazard"` phase type is −log(1 − G(t)), which grows without bound as G
  approaches 1, and nothing held its `t_half` inside the observed times. A fit
  could walk `t_half` below the first observation, evaluate the whole data
  range where G is essentially 1, and return a log-likelihood of **+290082**
  with `converged = TRUE` and no warning — a supremum reported as an
  interior optimum. On the shipped `cabgkul` data the fitted `t_half` was
  0.000352 against a first observed time of 0.0329.

  Such a fit now **warns** and records what was found in `fit$fit$boundary`.
  Nothing is bounded and no estimate moves: this reports, it does not
  constrain.

  The warning states a plain fact — `t_half` is below the observed support
  — with no tuned threshold, and the magnitude is in the record so you can
  judge it. Alongside the ratio, each entry carries **1 − G(t_min)**, the
  phase's remaining mass at the first observed time, which is the mechanism
  itself: a fit merely hugging the edge of its data measures 0.053, while the
  `cabgkul` fit above measures 3.8e-12.

  `fit$fit$boundary` is `NULL` when the check ran and found nothing, a list of
  records when it found something, and `NA` when it did not run — with the
  reason in `fit$degraded_causes`, so "nothing found" stays distinguishable
  from "never looked". Catch the warning with `hzr_unbounded_phase`, or the
  whole boundary family with `hzr_boundary`.
* **The G3 phase's `log_tau` derivative is now taken in `log_tau` (#352).**
  `.hzr_g3_phase_derivatives()` described itself as taking "central
  differences for log_tau" and stepped `tau` linearly instead, with an
  absolute floor of `1e-10`. Once `tau` fell below that floor the step was
  larger than `tau * h`, so the step stopped shrinking with `tau` and became
  a large *relative* step; below `tau = 1e-10` it also exceeded `tau` itself
  and the difference turned one-sided. The derivative the optimizer and the
  Hessian both use was **99.4% wrong at `tau = 1e-12`**, 1.4% wrong at
  `1e-9` and 0.012% wrong at `1e-8` — the second and third of those from the
  relative-step effect alone, with the branch still central — measured
  against an analytic derivative of the closed form. The step is now
  proportional to `tau` at every scale, and the one-sided fallback is removed
  because it can no longer be reached.

  Where `tau` is so small that multiplying it by `exp(1e-5)` returns the same
  number — below about `5e-319`, at the bottom of double precision — or where
  `tau` is infinite, the two evaluation points coincide. That is now reported as `NaN` rather than the
  plausible `0` a coincident difference quotient produces.

  **Some late-phase (`g3`) fits will move.** Where `tau` is small the
  optimizer now follows a more accurate gradient and can land somewhere
  measurably different: across six trial two-phase fits, two moved by more
  than `1e-8` relative, one of them by **21% on a parameter and 42% on a
  standard error**, with the objective **0.0126 log-likelihood units better**
  — a better optimum, not merely a different one. Fits whose shapes stay
  above about `1e-5` move by around `1e-10` relative, which is the precision
  the step change itself carries; between `1e-8` and `1e-5` the old
  derivative was wrong by between `1e-4` and `1e-10`, so fits there can move
  by more than that. **No fit in this package's own test suite
  moves**: its results are identical before and after, to every assertion.

  The `gamma` and `eta` steps keep their existing floors deliberately. `G3`
  is very nearly linear in each of them near zero, so the floor stays small
  relative to the scale on which the function varies even when it is 100% of
  the parameter, and both measure accurate to `1.1e-6` or better at the
  shapes where the `tau` derivative failed.

* **A single-distribution `theta` must have one entry per parameter, and a
  Weibull scale and shape must be positive, fitted or not (#375, #383).**
  `hazard()` compared a supplied `theta` only with the design's column
  count, as a lower bound, so a wrong length was caught only sometimes, and
  when it was not, the result could be wrong. With `fit = TRUE`, some wrong
  lengths failed with an unrelated error (`non-conformable arguments`), and
  some fitted silently: a `theta` holding only the shape parameters fitted
  the model with its covariates dropped, and on a one-covariate model a
  `theta` one entry too long returned its starting values unfitted. With
  `fit = FALSE` the object was built, and `predict()` then either failed
  with an unrelated message or, for a model with no covariates given an
  extra entry, applied it to a `newdata` column the model never had and
  returned a wrong prediction with no warning. A Weibull scale or shape at
  or below zero failed with `non-finite value supplied by optim`. Both are
  now refused, naming the lengths or the parameter, for example
  `'theta' has 2 entries, but this weibull model takes 3: 2 shape
  parameters, then one coefficient per column of the design (1 column).`
  The count is the likelihood's, so `control$shape_param_count`, which the
  likelihood ignores, does not change it. Unlike a multiphase model (#408),
  an unfitted single-distribution model is refused too: its parameter count
  is known without fitting, and an object of the wrong length could not be
  predicted from correctly. A multiphase specification may carry fewer
  entries until a fit resolves its phases' designs.

  The same check now runs in `predict()`, ahead of the type dispatch rather
  than inside one branch of it. A stored `theta` longer than the design
  allows, in a hand-edited or legacy object, was refused by
  `type = "hazard"` and `"linear_predictor"`, where the design is multiplied
  as a matrix, but `"survival"` and `"cumulative_hazard"` recycled the
  surplus coefficients into an outer product and returned two values per row
  with no error. They now refuse, naming both counts, as `hzr_evaluate()`
  already did.

  Separately, `predict(newdata = )` now **warns** when it matches `newdata`'s
  columns to a model's coefficients **by position**. That happens only for an
  object that stored no design matrix, where position is the only mapping
  left, and it means reordering or renaming `newdata`'s columns silently
  changes the predictions. The warning names how many coefficients are being
  matched, shows the columns it used, and says to refit so the design is
  stored and the mapping is by name. The behaviour is unchanged: `hazard()`
  already refuses to build such an object, so one can only arrive from an
  older version or by hand, and it still predicts.

* **A data defect reaching the score criterion is no longer reported as a
  numerical failure (#407).** The score path absorbs a Hessian it cannot
  build or invert and says so, which is right, but it absorbed *every*
  error, so a defect in the data raised inside the likelihood came back as
  "the current model's information matrix could not be inverted" and sent
  the reader looking at conditioning. Since the likelihood's data guards
  carry the class `hzr_data_error` (#426), the six sites that wrap the
  likelihood, its gradient and its Hessian now let that class through and
  go on absorbing everything else. Two of those six can actually receive
  one, the multiphase Hessian and the multiphase gradient; the other four
  wrap code that never calls the guards, so their narrowing is defensive
  rather than a behaviour change, and both live sites are pinned by a
  test. A genuinely numerical failure is
  unchanged: it is still swallowed and still reported as one. The two
  `tryCatch` calls that wrap linear algebra rather than the likelihood,
  `solve()` on the information block and `model.frame()` on a phase
  formula, are untouched, since a failure there really is what they exist
  to absorb.

* **A Conservation of Events fit now reports the log-likelihood of the
  estimates it returns (#362).** Under CoE the conserved phase's scale is
  re-derived after the optimizer finishes, and the reported objective was the
  optimizer's own value, taken before that step. The two could describe
  different parameter vectors: on one fit of the shipped `avc` data the fit
  reported `-71.934410193` while the likelihood of its own returned `theta`
  was `-78.1497959154`, a gap of 6.22. Anything reading the objective
  inherited the discrepancy, including `print()`, `summary()`, the `logLik`
  and `delta_logLik` columns of `hzr_stepwise()$steps`, and the log-likelihood
  the score criterion works from. The objective is now recomputed at the
  returned estimates, so `objective` and `theta` describe the same point. The
  one exception is loud: if the likelihood cannot be evaluated there, the
  optimizer's own value is kept and a warning says so. **No estimate
  changes**: `theta` is untouched and only the number
  reported beside it moves, and only for fits where the two had diverged. On
  the datasets measured that was 2 fits in 15.

  **This corrects the report, not the fit.** The largest gaps arose where the
  estimates are themselves unsound: standing on a discontinuity in the
  likelihood, where a change of one floating-point step in a parameter moves
  the log-likelihood by several units (see Known limitations, #448). A fit in
  that state reports `converged = TRUE`, and that flag does not mean the fit
  is sound there. Read the relative-gradient test beside it, which such fits
  fail.

  What a sentinel objective should mean for a single fit is tracked
  separately (#351, #374).
* **A ridge is no longer named from a covariance that is not a covariance
  (#416).** `summary()`'s weak-direction report reads the flat direction from
  the correlation of the estimates. When the Hessian was taken where it is not
  negative definite, typically short of the optimum, standardising it produced
  a matrix with "correlations" outside -1 to 1, and a ridge was reported from
  it: on the fits measured, correlations of 1.01, 1.31 and 11.3. Those are
  impossible, and every one of them leaves a negative eigenvalue, so the
  report is now declined for such a matrix, with `weak` set to `NA` and the
  reason `"covariance is not positive definite"` rather than a named set of
  parameters. A genuine ridge is unaffected: a real correlation matrix is
  positive semi-definite, so a true flat direction sits at or above zero.

  The eigenvalue test uses a tolerance scaled to the numerical error of the eigenvalue computation,
  `n * eps * max|lambda|`, taken from the correlation matrix. An earlier draft
  used a fixed `-sqrt(eps)`, about `-1.49e-08`, which still admitted matrices
  that are indefinite far beyond rounding error, so a ridge was named for one whose
  off-diagonal read `1.00000001`. The scale is taken from the correlation
  matrix and not the covariance deliberately, since the correlation matrix is
  scale-free and the decision must not depend on whether a time was recorded
  in days or years.

* **`hzr_stepwise()` and `hzr_bootstrap()` warn once about a `control`
  element the fit does not read, not once per candidate refit (#410).**
  Since #376 made `hazard()` warn about an element a fit ignores rather
  than accept it silently, both functions have handed `control` to every
  candidate refit, so one warning became **six** in a three-step screen and
  **three** in a select-mode bootstrap, one per candidate refit of the
  screen it runs on the real data before resampling. (The replicate screens
  run muffled, so the bootstrap's count did not grow with `n_boot`.) No
  result changed: the harm is to
  **other** warnings, because past 50 R prints only "There were 50 or more
  warnings", so the repeats can bury the ill-conditioned-Hessian and
  gradient-test warnings that say a fit is not to be trusted. The forwarded
  `control` is now validated once, at the call, and the refits are given
  what survives, so they have nothing left to warn about. An element the fit
  does read is still forwarded and still takes effect. One case gains a
  warning rather than losing repeats: a screen that never refits a candidate,
  such as `criterion = "score"` with a threshold nothing clears, reported an
  ignored `control` name not at all, because the warning came from the
  refits.

* **`force_in`, `force_out` and a character `scope` now name a variable by
  looking it up, and warn when a name matches nothing (#437).** `terms()`
  backquotes a label whose variable is not syntactic, so the column `_X1`
  appears as `` `_X1` `` among a model's terms, while the three arguments
  are documented as variables and the SAS translator emits bare names. The
  two spellings never met: a pinned variable was **dropped with no warning
  naming it**, a `force_out` one was still offered, and a character `scope`
  re-offered a variable the model already had. Names that reach R this way
  are ordinary in translated work: a leading underscore, a dot, a reserved
  word. Each name is now resolved once, when the screen starts, by lookup
  rather than by reading the string: a name that is exactly a column of
  `data` is that column, so `"_X1"` pins `_X1` and `"TRUE"` pins a column
  named `TRUE`; otherwise a name that is exactly a term label of the model
  or `scope` is that term, so `` "`_X1`" `` and `"age:mal"` work as well;
  and a name that is neither is ignored **with a warning naming it**, where
  before it was ignored in silence. The column is looked up first, so when
  `data` has a column literally named `age:mal`, `"age:mal"` resolves to
  that column, and the interaction can be named only in a formula `scope`.
  Distinct columns stay distinct: `age`, `age ` and `age # x` are three
  columns, and a column literally named `` `x` `` is not `x`. A string that
  only resembles a name is not read as one: `"age "` when there is no such
  column, or `` "`age`" ``, which `terms()` never writes, is warned about
  and ignored. The names ignored are also recorded on the result, in
  `$scope$unresolved` (and `$unresolved` of a select-mode `hzr_bootstrap()`),
  and `print()` shows them, so `suppressWarnings()` or a saved object does
  not lose them; a character `scope` emptied this way says so where the
  screen stops, rather than "no further action".

  This is about MATCHING: which variables are pinned, excluded or in the
  scope. How a matched candidate then enters the model is the #449 entry
  below.

* **A stepwise candidate is now scored and entered as the column it names
  (#449, #438, #441).** The refit wrote the candidate's name, as spelled,
  into the formula text, so a column whose name reads as a different term
  entered as that term, with no warning. A column `age ` or `age # x`
  beside `age` refit as `age`, through the default `scope = NULL`, a
  character `scope`, a multiphase default scope and the screens
  `hzr_bootstrap()` runs; a literal column `age:mal` refit as the
  interaction, while under `"score"` its entry p-value was the column's;
  and a multiphase default scope was built the same way, so a strong
  column `x2 ` read as the noise column `x2` and was never scored. The
  refit, the multiphase default scope and the multiphase score, which
  builds the candidate's phase formula from text as well, now write the
  label `terms()` gives the resolved column, which reads back as that
  column, and the score reads the values of that column and compares it
  with the model's terms by that label. So the column scored is the
  column entered, and `"wald"`, `"aic"` and `"score"` reach the same
  model: with the interaction `age:mal` in the model, a literal column
  `age:mal` is a variable of its own under all three, where `"score"`
  had declined it as the interaction. The same change fixes two loud
  failures: a bare non-syntactic name such as `"_X1"` in a character
  `scope` now enters under `"wald"` and `"aic"`, where its refit failed
  to parse (#441), and under `"score"` a non-syntactic candidate written
  as its label, as a formula `scope` writes it, is read from its column
  rather than reported as not found in `data` and skipped (#438). An
  interaction is no longer scored from a literal column that shares its
  spelling; the score declines it, as it does any term that is not a
  column, and its warning now says that instead of "not found in `data`",
  with the new reason `not_single_column` in `$criteria$uncomputable_reasons`
  where it read `non_numeric`.
  A column no formula can name, such as one called `.`, is not offered as
  a candidate, and the screen says so once.

  `$steps$variable`, `$scope$frozen` and `$criteria$wald_untested_entries`
  name such a variable by its label, as the breaking change above on
  `$steps$variable` sets out; `$criteria$refit_failures` still names a
  failed candidate as the scope wrote it.

* **`hzr_translate_sas()` no longer fails on a SAS covariate whose name begins
  with an underscore** (#411). `PROC HAZARD`'s lexer reads a name as
  `[_A-Z][_A-Z0-9]*` (`hazard_l.l:39`), so `_X1` is a legal phase-statement
  covariate. The translation built each phase formula by pasting the names
  into `str2lang()`, and an R symbol may not begin with an underscore unless
  it is backquoted, so the job stopped with R's own parser error
  (`unexpected symbol`) rather than anything about the job. Formulas are now
  built from symbols, at the phase statements and at the `SELECTION` scope
  alike, and `deparse()` backquotes such a name so the emitted document
  re-parses to the same call.

  The failure was loud, so no fit stood in for one: such a job produced
  nothing. On a sample of the studies share, 13 distinct `PROC HAZARD` steps
  fail this way and none of them carries a macro, making this a second and
  independent cause of an unreadable step.

  Note that the column must really be named `_X1` in the data frame. R's
  `data.frame()` renames it to `X_X1` unless you pass `check.names = FALSE`,
  and a renamed column is **refused** by name rather than quietly dropped, so
  a fit cannot come back short a covariate without saying so.

  **A `SELECTION` job carrying such a name is refused rather than screened.**
  `hzr_stepwise()` spells a non-syntactic name two ways at once: backquoted
  in the `terms()` labels its candidates are keyed on, bare in `force_in`.
  The two never match, so a `/I` pin is ignored and a `BACKWARD` screen can
  drop a variable `PROC HAZARD` holds in, with no warning naming it; and the
  score criterion, the only one this translator emits, indexes the data by
  the backquoted label and skips the candidate as "not found". Both are wrong
  models delivered as populated results, so such a job now stops and names
  the cause. For a name like `_X1` this costs nothing: the job stopped before
  this release too, one step earlier, in the phase formula.

  Text `PROC HAZARD` does not accept as a name is **not** refused here,
  because it no longer reaches this check. `AGE*SEX` is not a name
  (`hazard_l.l:39`, `hazard_y.y:213`), so `PROC HAZARD` rejects that job at
  parse; the phase parser now leaves such an operand out of the model, with
  a warning and an untranslated row, rather than stopping a job that
  translated before (#440, under Breaking changes).

* **`hzr_translate_sas()` builds a phase whose `PARMS` writes only its scale**
  (#345). An active `MUE` or `MUL` with no shape operand used to be recorded
  as untranslated and build no phase. PROC HAZARD runs that phase on its own
  shape defaults (early `THALF` 1, `NU` 2, `M` 1; late `GAMMA` 1, `ALPHA` 1,
  `ETA` 2), which do not depend on the data, so the translation now builds it
  the same way, provided it read the whole `PARMS` statement. If any operand
  could not be read (for example one written with spaces around `=`, which
  `PROC HAZARD` accepts), a shape may have been written that the translator
  did not see, so the phase is recorded and not built. A `PARMS` value that
  `PROC HAZARD`'s lexer does not read as a number (`1E-3`, `2.`, `+0.2`) is
  now recorded as a syntax error, not read by R and fitted. The late `TAU` start (0.75 of the longest follow-up) is the one
  value that depends on the data, and it is recorded, as it already was for a
  late phase written without `TAU`. That record now says what it means:
  because the multiphase likelihood is multimodal, a different start can
  change the estimates, not only the path to them; and with `FIXTAU` on an
  unwritten `TAU`, PROC HAZARD holds `TAU` at that data-dependent value while
  the translation holds it at 1, a different model.

* **Two `hzr_translate_sas()` rows now state their consequence** (#345 review).
  - `FIXMNU1` on an active early phase is a real PROC HAZARD constraint
    (`|M*NU| = 1`) that the translation does not apply. It was recorded as "PARMS
    token has no phase target", which read as a parsing gap; the row now says
    the constraint is not applied and the emitted phase is a different model.
  - A `PARMS` keyword that is not in PROC HAZARD's grammar (for example
    `FIXG1` or `FIXG3`, which are internal flags, not options) is one PROC
    HAZARD rejects, so its job does not run. The row keeps its "unresolved
    PARMS keyword" prefix and now says that.

* **A translated `SELECTION` job's check chunk no longer repeats
  `hzr_stepwise()`'s own warnings, or calls a failed Wald test a score it
  could not compute (#400).** Since #399, `hzr_stepwise()` warns when a
  screen stops on candidates it could not test, listing every reason. When
  the screen completed, it names each variable it kept in the model without
  a Wald test of its removal. The check chunk after the screen
  still warned twice more. First it warned that the model had "no usable
  standard error". Then it warned that
  `N candidate score(s) were uncomputable`, counting the Wald failures as
  scores. It now reads
  `$criteria$uncomputable_reasons` without `wald_no_variance`, names each
  reason, and stays quiet when the screen stopped, because that warning
  already lists every reason.

* **A multiphase fit stops when `theta` does not have one entry per
  parameter, instead of fitting with extra entries or failing obscurely
  (#408).** `hazard()` compared a supplied `theta` only with the global
  design's column count, and only as a lower bound. A multiphase `theta`
  that was too long therefore fitted with the extra entries carried along,
  and reported `converged = TRUE` with a `theta` longer than the model; one
  that was too short failed inside the fit with `'names' attribute [9]
  must be the same length as the vector [7]`. The fit now stops, naming
  both lengths and each phase's share: `'theta' has 11 entries, but this
  model takes 9 (early 6, constant 3)`. A phase with its own `formula` is
  counted from that formula, every other phase from the global design.
  Unfitted (`fit = FALSE`), `theta` is still returned as supplied.

* **A translated job's missing-value guard no longer stops on rows
  `hazard()` drops anyway, and an absent phase variable is named (#340).**
  The guard for a variable outside the fitted model (`/E`, say) stopped
  whenever it was missing, even on rows where a modelled variable was
  missing too; `hazard()` drops those rows itself, as SAS does. It now stops
  only on rows `hazard()` would keep. A phase-statement variable that the
  dataset lacks used to fail as "object ... not found"; the status chunk
  now names every such variable and the dataset. A phase-statement variable
  that is not numeric (a character or factor column) now stops the job too,
  as `PROC HAZARD` does ("VARIABLE NOT NUMERIC"); the translation used to
  dummy-code it and fit a model SAS never ran.

* **`hazard()` now warns about every `control` element it does not read
  (#376).** `control` used to accept any name silently, so a mistyped one,
  such as `n_startz` for `n_starts`, left the default in force and said
  nothing. `hazard()` accepts `maxit` and `reltol` for every model,
  `shape_param_count` for a single-distribution one, and `n_starts`,
  `conserve`, `phase_share_tol` and `start_seed` for a multiphase one. The
  fit reads all of them except `shape_param_count`, which the
  single-distribution stepwise refit and score test read back from the fit.
  Any other element now draws one warning that names it and says why it
  has no effect, and the fit proceeds unchanged, as `stats::optim()` does
  for unknown `control` names. The element is dropped before the fit: R's
  `$` matches a partial name, so `n_starts_extra` used to be read as
  `n_starts`, and `fit$spec$control` now keeps none of the ignored
  elements. That covers a misspelling, an unnamed
  element, five names `?hazard` used to document as accepted although no
  fit read them (`abstol`, read only by a bounded optimizer that no fit
  uses, and `method`, `condition`, `nocov` and `nocor`), `fix` and `quasi`,
  a multiphase element such as `n_starts` given to a single-distribution
  fit, and `shape_param_count` given to a multiphase fit, where nothing
  reads it. No name is an error (a bad value for an element the fit reads,
  such as `maxit = "a"`, still stops a fit that reads it): `hzr_stepwise()` and
  `hzr_bootstrap()` pass `control` to every candidate refit, and an error
  there counts as a failed candidate, so a screen would report success
  having tested nothing. The SAS options `CONDITION=`, `NOCOV`, `NOCOR` and
  `QUASI`, which the SAS-to-R migration vignette used to show as `control`
  elements, have no `control` equivalent.

  **If you used `control$fix`, your fit was not constrained.** It was never
  documented, and the fitting code never read it: a fit given
  `control = list(fix = ...)` was the unconstrained fit, with its "fixed"
  parameters free. It now draws a warning saying so. To hold a parameter at
  its starting value, use `hzr_phase(fixed = )` on a phase of a multiphase
  model and re-run; a single-distribution model has no mechanism for fixing
  a parameter. Unlike the other entries that ask you to re-check results,
  this one has no known affected user: nothing in this package or its tests
  passed `control$fix` to `hazard()`, so the exposure is limited to anyone
  who found and used the undocumented name. `control$quasi` was never read
  either, and warns the same way.

* **The `objective = "sas"` interval checks name the real defect (#340).**
  The objective's own guard let an interval row with an `NA` bound through:
  `which()` drops an `NA` comparison, so the row became `-Inf`, a value the
  optimizer walks away from, while the entry check stops on the same row.
  Both now stop on it. The entry check also reports a `time_lower` or
  `time_upper` shorter than `status` as a length mismatch, naming `time`
  when the bound was left to default to it, where it reported an `NA`
  bound. Both are reachable only by calling the internals
  directly or editing a fit's stored data, since `hazard()` checks lengths
  and missing bounds first.

* **A score-criterion `hzr_stepwise()` screen on a multiphase base that
  dropped rows with missing covariates now stops and says so (#372).** Such a
  fit drops every row whose phase covariate is missing, `NA` or `NaN` (in the
  data, or made so by a transform such as `sqrt()` or `log()` of a negative
  value), but keeps the full response in `$data`. The score test's row check
  counted that full response, so it passed; every candidate then failed to
  line up with the fit's design, and the screen stopped with no steps, blaming
  each candidate as `not_expandable`. The check now counts the rows the fit
  was estimated on and stops with an error naming how many rows the base
  dropped and the remedy: refit the base model on only the rows it used and
  pass that data frame. `hzr_bootstrap()` with `scope` and
  `criterion = "score"` on such a base changes the same way: it used to run
  replicates that could score nothing, and now stops before the first. This
  was never a silent wrong answer (the screen always warned that nothing could
  be scored), but it reported the wrong cause, and a screen that used to
  finish with zero steps now stops with an error. A base fitted on complete
  data is unaffected.

* **`hzr_stepwise()` now says when it could not run a Wald test for an entry
  or a removal (#389).** A Wald test needs the model's variance for the
  coefficient. When that variance was missing, the p-value was `NA`, and the
  variable was treated exactly as if it had been tested: an untested removal
  stayed in the model as if it met `slstay`, and an untested entry stayed out
  as if it missed `slentry`. There was no warning, and
  `$criteria$n_uncomputable_scores` stayed 0. The common case is an
  interval- or left-censored multiphase fit on an installation without
  `numDeriv`, where no coefficient has a variance. A backward screen there
  stopped after 0 steps and kept variables it drops when `numDeriv` is
  present; a forward Wald screen stopped after 0 steps and entered nothing.
  Such tests are now counted in `n_uncomputable_scores` under the reason
  `wald_no_variance`, and the variables are listed in the new
  `$criteria$wald_untested_removals` and `$criteria$wald_untested_entries`.
  A screen whose last iteration could test none of its candidates for entry,
  or none for removal, sets `stopped_uncomputable` and warns which. Any other
  variable decided without a test is named once in a warning.
  `hzr_bootstrap()` runs its replicates quietly, so it now warns with the
  number of replicates that decided a variable untested. `stopped_uncomputable`
  is now decided by the last iteration alone, so a two-way screen that could
  not test anything at one step and recovered at the next is no longer
  reported as stopped. A forced-in variable
  is never a removal candidate and is not counted. The stop warnings of
  `hzr_stepwise()` and `hzr_bootstrap()` now name both criteria, and
  `hzr_bootstrap()` no longer says such replicates contribute no selections:
  a replicate may stop after several steps, and an untested removal counts
  as selected.

* **A backward `hzr_stepwise()` drop's refusal now names the reduced design,
  and its pre-check no longer repeats parse-time warnings (#343).** The
  reason for refusing a drop that removes no column read "the refit's
  design", but on the single-distribution path that design is built before
  any refit, so it described something that did not exist; it now reads "the
  reduced design", on both paths. The same pre-check parsed the current and
  the reduced formula before the refit parsed the reduced one again, which
  doubled any warning raised while building the design: 8 per step instead
  of 4. It now parses quietly. When `data` is the frame the base was fitted
  on, the current formula's warnings already surfaced from the base fit, and
  the reduced formula's warnings surface again from the refit if the drop
  goes ahead. When `data` differs, the current formula's warnings on it are
  no longer shown; they cannot change the pre-check's decision, which
  compares column counts.

* **A classed numeric *matrix* column was read as its raw storage, when
  fitting as well as predicting (#371).** `hazard()` and
  `predict(newdata = )` read a classed numeric column as its values (#231,
  #347), but a column carrying a `dim` was left alone entirely, because a
  genuine matrix column (`I(cbind(p, q))`, a `Surv`) must not be flattened.
  A `bit64::integer64` matrix column therefore reached the model as its
  stored doubles. On 120 rows of `avc` with `age` as such a column, the fit
  ran on 9e-300 in every row: the covariate coefficient stayed at its
  starting 0.004 and the log-likelihood came out -93.05 against -92.45 for
  the same numbers as a plain column, with no error and no warning. The same
  column in `newdata` predicted 0.1, 0.1414 and 0.2236 where the values give
  0.1271, 0.2027 and 0.2521. Such a column is now read as its values and
  keeps its shape. Which classes need reading is decided by behaviour rather
  than by a list: the class's own `as.numeric()` is compared with the stored
  doubles, and the column is replaced only when they differ, so a `Surv`
  column, whose stored doubles are its values, keeps its class.

* **A fit made with `survival::Surv()`'s own status codes was wrong, not
  empty, and said nothing (#231).** `Surv()` codes interval-censored rows
  `3`, and this package codes them `2`. Passing survival's integers as a
  plain `status` vector -- what `unclass(sv)[, "status"]` or `sv[, 2]` gives
  -- left those rows out of the log-likelihood while the analytic gradient
  still counted them, so the fit converged to the optimum of neither model.
  On 200 rows with 50 interval-censored, the scale parameter came out 14.6%
  away from the same data coded correctly, with no error and no warning. If
  you have fitted interval- or left-censored data by passing `Surv()`'s
  codes through, re-run it: either pass the `Surv` object itself, which is
  translated, or use this package's codes (`-1` left, `0` right, `1` event,
  `2` interval). Such a `status` is now refused, naming the offending rows.
* **`hzr_translate_sas()` now mirrors PROC HAZARD when `FIXGE2` or `FIXGAE2`
  meets `SETG3_ignore_tau()`** (#328, #329 review). That branch runs when
  both flags are set, or when either is set with `ALPHA` fixed at 1. PROC
  HAZARD then fixes all four late shapes: `TAU` = 1, `ALPHA` = 1, and `GAMMA`
  and `ETA` at 2 and 1 (or 1 and 2 when the job wrote `ETA = 2`). The
  translation used to record these jobs as untranslated and still emit a
  phase with `GAMMA` free. It now emits the fixed phase, records any value
  the job wrote that neither program uses, and records `SETG3940` when
  `ALPHA` is fixed at anything other than 1.

* **`hzr_phase("g3")` now refuses an infinite `tau`, `gamma` or `eta`, and a
  derived shape that is not a finite positive number** (#329 review). An
  infinite shape was accepted and failed only inside the optimizer. Under
  `constraint`, finite sources could still overflow to an infinite derived
  shape, or underflow to `alpha = 0`, which would silently select the
  exponential limiting case. `hzr_translate_sas()` now records a job whose
  late shape is not finite as written (`GAMMA=1e400` reads as `Inf`) or after
  a `FIXGE2`/`FIXGAE2` rewrite, since its emitted `hzr_phase()` call can no
  longer be built.

* **`hazard(fit = FALSE)` applies a phase constraint to a supplied `theta`
  when phases carry covariates** (#328). It used to warn that the derived
  slot could not be located, even for a `theta` already on the constraint.
  The slot is now located the way the fit locates it, so an on-constraint
  `theta` passes silently and an off-constraint one is replaced with a
  warning.
* **A multiphase fit no longer stops with "t_half must be a positive
  scalar" when the optimizer steps a phase's time scale out of range
  (#262).** `t_half` and `tau` are carried on the log scale, and a step can
  take them past what `exp()` represents, to `0` or `Inf`. Where
  `exp(log_t_half)` came back as 0, the log-likelihood, its gradient and
  the Conservation of Events solve all raised that error rather than
  treating the point as infeasible, so a single-start fit
  (`control = list(n_starts = 1)`) failed outright, and with several starts
  that start was lost. Elsewhere the three disagreed: the score raised
  `missing value where TRUE/FALSE needed` at `t_half` or `tau` of `Inf`,
  where the log-likelihood returned `-Inf` (`t_half`) or a finite value
  (`tau`). A phase whose time scale is not finite and positive is now
  infeasible on all three paths, and the optimizer backs away from it.

  That last case is a **deliberate behaviour change**: at `tau = Inf` the
  late phase switches off and the log-likelihood used to take that limiting
  value, so a fit whose `tau` ran past `exp(709.78)` could return a finite
  objective there while its score could not be evaluated. Such a point is
  now infeasible. Fits whose scales stay in range are unchanged.
* **The saturated-phase warning no longer claims the likelihood is unchanged
  when the fit evaluates the phase elsewhere (#228).** The warning reports a
  phase whose contribution is constant across the event times, and said its
  shape parameters were unidentified and "the likelihood is unchanged whether
  they are pinned or fitted". An interval-censored or left-truncated
  likelihood also evaluates the phase at the interval bounds and the
  counting-process entry times, where a phase flat across the event times can
  still be climbing: on one such fit the likelihood moves 103 units in a
  parameter the message called unchanged. Where the fit has such points, the
  warning now says how many there are and that the shapes may be identified
  there. **The warning still fires on such a fit**: which phases are
  reported, and `fit$fit$phase_share`, are unchanged, and what those further
  points are worth is not measured. The diagnostic says less than it did,
  not something different.
* **`predict()` on a model built with `fit = FALSE` now says so (#144).**
  Its numbers come from the starting values the model was given, not from
  estimates, and it said nothing: predicting from non-estimates in silence
  is the defect a `fit = FALSE` object invites. It now warns, under the
  condition class `hzr_unfitted_prediction`, which
  `options(TemporalHazard.warn_unfitted_prediction = FALSE)` switches off
  for code that means it. Predicting from such a model remains supported
  for every distribution but `"multiphase"`.

* **`predict()` on an unfitted multiphase model now says what is missing
  (#144).** It failed with an unrelated internal error,
  `missing value where TRUE/FALSE needed`, from `predict()` reading the
  per-phase covariate counts that an unfitted object does not carry. The
  error now says that a multiphase model built with `fit = FALSE` has no
  per-phase design matrices, because they are resolved when the model is
  fitted, and points at `fit = TRUE` or the new `hzr_evaluate()`. Models of
  the other distributions built with `fit = FALSE` still predict from their
  supplied parameters, with the warning above; only multiphase ever failed.

* **`predict(newdata = )` now warns when `newdata` is evaluated differently
  from the fitting data (#331, #334, #335).** Predicted values are
  unchanged; the warning names the cause. It fires for a term that computes
  a statistic over the rows, such as `I(age - mean(age))`,
  `I(scale(age)^2)` or a `factor()` nested inside another call, and for a
  column whose type differs from the fitting data's, such as a numeric
  column given as character or a `difftime` in other units. It also fires
  when a fit saved by 1.2.10 or earlier has its design rebuilt under a
  contrasts function other than `contr.treatment` or `contr.poly`, which
  that fit did not record. The new section "How `newdata` is evaluated" in
  `?predict.hazard` describes these cases and two that are not detected: a
  formula-environment constant changed since the fit, and collation in
  string comparisons. Fits saved before 1.1.0 kept no fitting data, and
  their column types are not checked.
* **One row at time 0 no longer empties an exponential, Weibull or
  log-normal fit (#341).** A row right-censored at time 0, as
  `Surv(0, NA, type = "interval2")` gives, contributes nothing to the
  likelihood. With any left- or interval-censored row in the data, these
  three families refused the whole fit instead: the objective was clamped,
  and `hazard()` reported `converged = TRUE` at the starting values, with
  no warning about the data. The log-normal refused such a row on any data,
  and also refused an interval opening at 0, `(0, u]`, which is left
  censoring at `u`. Each is now evaluated as the row it is, matching the
  log-logistic and multiphase fits. Fits without such rows are unchanged.

* **The vignettes no longer skip every chunk in silence when the rendering
  session cannot see the installed package (#276).** Each vignette gates its
  chunks on `requireNamespace("TemporalHazard")`, so a render session that
  could not load the package produced a complete-looking document with no
  computed output and no error. Rebuilding the vignettes now puts the
  library that `R CMD build` installs the package into back on
  `.libPaths()`, and under continuous integration a package that will not
  load stops the build and says why.

* **A multiphase `hzr_stepwise()` `scope` whose element is not a formula now
  says so (#328).** `scope = list(late = c("age", "mal"))` stopped with
  "formula must be a `formula` object", which named neither `scope` nor the
  phase. It now names `scope$late` and says a multiphase `scope` is a named
  list of formulas keyed by phase. A character vector remains a valid scope
  for a single-distribution fit.

* **`predict(newdata = )` on a multiphase fit saved before this version no
  longer gets `scale()`, `poly()` or `ns()` in a phase formula silently wrong
  (#307).** Such a fit stored no phase design, so the phase was rebuilt from
  `newdata` alone, and those terms took their centering, scaling or basis
  from `newdata`'s own rows instead of the fitting data. At all of the
  fitting rows that reproduces the fit; at any other `newdata` it does not.
  At three of the fitting rows, `scale(age)` was off by up to 96%,
  `poly(age, 2)` by a factor of 3e5, and `ns(age, df = 3)` by 66%. One row
  of a `scale(age)` phase came back as a zero-length prediction.

  Such a fit's phase is now rebuilt only when that can be checked, and is
  otherwise refused with advice to give `newdata` the fitted design columns
  or to refit. A fit saved by 1.1.0 or later kept its fitting data, and its
  phase design is rebuilt from that data exactly as the fit built it.
  Reproducing the fitted rows is not enough, since a `cutoff` moved between
  two fitted ages changes no fitted row. So the phase formula must use only
  the kept data's columns and R's own design functions (a user's function
  of the same name is not one), hold no term coded by contrasts (a factor,
  character or logical column, `cut()`), whose coding the fit did not
  record, and rebuild the fitted columns exactly. Each function must be
  written as a plain name or as `pkg::fn` with both parts written as names;
  a quoted spelling such as `base::"log"(age)` is not recognised, and such a
  phase is refused at `newdata` rather than rebuilt. A fit saved by 1.0.3 or
  earlier kept neither design nor data, so it cannot say which of its
  formula's names were data columns: a constant `k` in `I(age * k)` that is
  gone at predict time would be taken from a `newdata` column named `k`.
  Its phase is refused at `newdata`; it still predicts without `newdata`,
  and from its design columns. A fit made before duplicated design column
  names were refused (#296) is refused at `newdata` too, since no selection
  by name can tell its columns apart. Current fits still evaluate a formula
  against `newdata` and their environment, as `lm()` does; that is #331.

* **`hzr_bootstrap()` now bootstraps a vector-interface fit made without
  `data =`** (#259, #312). It counted the rows to resample in the fit's data
  frame, and such a fit has none, so every one was refused with a message
  that named its vectors `'NA', 'NA'` and sent you to the formula interface.
  The stored `time`, `status`, `time_lower`, `time_upper` and `weights` are
  now resampled together, and the replicates match those of the same model
  fitted with a formula and `data =`. Select mode (`scope =`) on such a fit
  still stops, now saying why: its candidate columns have no data frame to
  be resampled with. The other refusals now give their real reason too:
  vectors that do not have one value per row of `data =`, an object missing
  a stored vector, which names the missing argument, and a `data =` that is
  a list rather than a data frame. Only the vector interface accepts a list,
  and its bootstrap already stopped with the same `'NA'` message, so a list
  `data =` still does not bootstrap; only the message is new. One kind of fit
  is still refused, for its real reason: a multiphase fit saved before
  `hazard()` refused a phase formula without `data =`, whose formula was
  ignored. Its stored call can no longer be refit, and resampling it
  returned no replicates and no error, so `hzr_bootstrap()` now refuses it
  with the message `hzr_stepwise()` gives. That refusal takes only the
  ignored-formula check: a fit `hzr_stepwise()` declines to step for other
  reasons, such as a phase inheriting a factor with more than two levels,
  still bootstraps.

* **`hzr_bootstrap()` now says why replicates failed, and warns when every
  one did.** Each replicate catches its own error so one bad resample cannot
  end the run, and it used to drop the message: a run could fail every
  replicate and return an empty `replicates` table with only `n_failed` to
  show for it. The result gains `failure_reasons`, a named integer vector
  counting each failure by its error message, or by
  `"non-finite objective (did not converge)"`, most common first. It sums to
  `n_failed`, and is empty but present when nothing failed. When no
  replicate succeeds, `hzr_bootstrap()` warns, naming the most common
  reason. Partial failure does not warn; its reasons are in
  `failure_reasons`.

* **`hzr_bootstrap()` no longer stops the whole run when a replicate's refit
  returns something other than a fit (#333).** Reading the objective off a
  bare vector was an error outside the replicate's own error handling. Such
  a replicate now counts as failed, under the reason
  `"refit returned a <class>, not a fit object"` (`an <class>` when the
  class begins with a vowel, as in `"an integer"`). `hazard()` never returns
  one; a stored call rewritten to another function can. The refusal of a
  `data =` that is not a data frame now names `data` in its message.

* **A single-distribution `hzr_bootstrap()` screen that selects nothing now
  warns.** The "selected no covariate" warning compared the replicates'
  parameters with `names(coef())`, which are `NULL` for a single-distribution
  fit, so its shape parameters `param_1` and `param_2` counted as selected
  covariates. Such a screen returned a summary of only those parameters, each
  at `pct = 100`, and said nothing. Multiphase screens already warned.

* **`hzr_bootstrap()` names two more kinds of refit that are not a fit
  (#343).** A refit returning a list with no `fit` was tallied as a
  convergence failure; it is now ``"refit returned a <class> with no `fit`, not a
  fit object"`` (`an <class>` before a vowel), in both modes. A refit whose fit held a finite
  objective but no estimates counted as a success and then ended the run
  building its replicate row; it is now a failed replicate,
  `"refit returned no parameter estimates"`. `hazard()` returns neither. A
  fit with no `data` frame whose call names a formula that was `NULL` when
  it ran, as a wrapper forwarding its own `formula` argument can leave, now
  stops with a message saying there are no rows to count, instead of
  `length(n) == 1L is not TRUE`.

* **The G3 late-phase shape is now accurate where `(t/tau)^gamma`
  underflows.** With a large `gamma`, event times well below `tau` take
  `(t/tau)^gamma` past double-precision underflow (about `exp(-708)`), and
  a very large `alpha` can underflow the same quantity divided by `alpha`.
  `hzr_decompos_g3()` then clamped the value to the smallest double, which
  froze `G3` below that time and put the log of `g3` wrong by more than
  100. The likelihood of those events was wrong, and the analytic Hessian,
  which differences the shape across that cliff, read a `log_tau` diagonal
  of 1.4e6 against a true 5.1e3, with no warning. Once the quantity falls
  below 1e-10, both logs are now computed in their limiting form, linear in
  `log(t/tau)`, which is accurate to about 1e-10. Fits with no time in that
  region are unchanged, and fits with one change only in the last digits
  unless they reached the old clamp. The SAS/C `HAZARD` code has a cliff
  here too: for `alpha > 0` its `ln(e^x + 1)` returns 0 below underflow,
  and for `alpha = 0` it already breaks down once `(t/tau)^gamma` is below
  about 1e-16. Fits that reach this region can differ from `HAZARD`; this
  package takes the accurate value.

* **`predict(newdata = )` on a model with no covariates now ignores
  `newdata`'s unused columns**, as it already did for models with
  covariates. For a `Surv(time, status) ~ 1` fit, or a vector-interface fit
  without `x`, every column other than `time` was taken as a covariate.
  The Weibull, log-logistic and log-normal fits then stopped with an error.
  The exponential fit, whose only baseline parameter is the log rate
  `log_lambda` and which has no shape parameter, used that log rate as the
  coefficient: it returned `age` times the log rate as the linear predictor
  with no error, -280 for `age = 70`, where the answer is 0 (#300).

* **A backward `hzr_stepwise()` drop now has to remove a column** (#320).
  Under treatment contrasts, `model.matrix()` codes an interaction whose main
  effect is absent with a full set of dummies: `~ z:f` gives `z:fa, z:fb`,
  which spans what `z, z:fb` spans. So dropping `z` from `~ z + z:f` removed
  no column, and the "reduced" model was the model it started from: the same
  coefficient count, the same column space, the same likelihood. The step
  accepted that drop and reported a p-value for it, while the fit still
  carried the variable. An uncapped run then stopped at the next step, where
  the interaction had become two columns, so the wrong step was masked by an
  unrelated error; a run that ended right after the drop (`max_steps`)
  returned it as a result, with `$steps` and the final model disagreeing.
  The post-drop refit is now checked against the model it came from, and a
  drop that does not reduce the design is refused with a reason in
  `$criteria$refit_failure_reasons`, as a failed refit already was. The
  forward step has refused the mirror of this, a candidate that adds no
  column, since #306. A single-distribution fit gets the same refusal and
  the same reason (#323). Its refit warm-starts from a `theta` one element
  shorter than such a design needs, so it used to fail to conform first and
  report "non-conformable arguments", which named the symptom and not the
  cause. Its reduced design is now decided before the refit.

* **A stepwise refit failure now says why.** `hzr_stepwise()` catches each
  candidate's refit error so one bad candidate cannot end the screen, and it
  used to drop the message with it: the warning read `candidate refit failed
  for gb.` and nothing more, even when `hazard()` had stopped with a message
  naming the problem. The four refit warnings (candidate, Wald fallback,
  post-entry, post-drop) now carry the refit's error message, or say that it
  did not converge, and `$criteria$refit_failure_reasons` keeps each one,
  named by its `refit_failures` token. `refit_failures` itself is unchanged.
  Under `criterion = "score"`, which builds each candidate's design without
  refitting, a candidate whose column name collides with a factor's dummy
  column (numeric `gb` beside factor `g` with level `b`) is now declined with
  the reason `duplicate_column`, and a run that completes anyway warns about
  it. The check uses the name the refit's `model.matrix()` would give the
  column, so a logical `flag` (column `flagTRUE`) beside a factor dummy `flag`
  is still scored. The single-distribution score path used to
  score it against a design with two `gb` columns, so it could win the step
  and fail only at the post-entry refit, which stopped the screen with every
  other candidate untested; the multiphase score path declined it as
  `not_expandable`.

* **Backward `hzr_stepwise()` now tests a dropped variable on its own
  coefficient** (#315). The drop test looked a variable's coefficient up
  by its bare name. `model.matrix()` names a logical `flag`'s column
  `flagTRUE`, while a factor `fla` with level `g` owns a column named
  `flag`, so in `~ fla + flag` the name `flag` found the factor's dummy.
  The run stopped with "expands to multiple coefficients" on `fla` before
  it decided anything, and that error is all that kept a wrong p-value out
  of the table: `flag` was being tested on the dummy (z = -1.40, against
  11.39 on its own column). The variable is now found by its term in the
  design the fit stored, for single-distribution and multiphase fits. A
  fit with no stored design (the `time =` / `x =` interface, or a formula
  fit saved by 1.2.10 or earlier) still uses the name, which is not safe
  from this collision; there the factor's own "expands to multiple
  coefficients" error still stops the run. Refit such a model with this
  version before a backward screen. On a single-distribution fit,
  a logical or two-level factor with no colliding column used to stop with
  "not found in the design matrix"; it is now tested.

* **Backward `hzr_stepwise()` no longer stops when `theta` is named after
  the covariates** (#304). A single-distribution fit started from
  `theta = c(mu = 0.1, nu = 1, age = 0, mal = 0)` failed its first drop
  test with "Unknown coefficient name(s): 'beta1'". The step names a
  covariate's coefficient by its position (`beta1`, `beta2`, ...), while
  the Wald test took the names from `theta` whenever all of them were set.
  Unnamed and `beta`-named starting values ran, so the naming a user is
  most likely to write was the one that failed. The Wald test now names
  single-distribution coefficients by position, whatever `theta` is called.
  That also matters for correctness: had the names been matched, a
  covariate called `nu` would have been tested as the shape parameter.
  A fit whose `theta` holds only the shape parameters (`c(mu = 0.2, nu = 1)`
  for `~ gb + z`) fits without its covariates, and position would then test
  `mu` and `nu` as `gb` and `z`. An unnamed `c(0.2, 1)` already did so,
  reporting both covariates as highly significant. The Wald test now stops
  with an error when `theta` does not hold one value per shape parameter and
  per covariate column. It also refuses a `time_windows` fit by name. There
  each covariate has a coefficient per window, and `beta1` tested only the
  first window's, while `beta2` could be another window's coefficient of the
  same covariate, a p-value for the wrong coefficient with no warning.

* **`hzr_stepwise()` now tests an entering candidate on its own
  coefficient** (#305). The Wald criterion, and the Wald fallback the score
  criterion uses for a candidate it cannot score, looked the new coefficient
  up by the variable's bare name. `model.matrix()` names a logical `flag`'s
  column `flagTRUE`, so when a factor already in the model had a dummy
  column named `flag` (a factor `fla` with level `g`), the step tested that
  dummy instead: a p-value for the wrong coefficient, with no error and no
  warning. On a simulated example the screen reported p = 0.16 for a
  candidate whose own p-value was below 1e-29, and did not enter it. The
  candidate is now found as the design-matrix column its refit added, the
  rule the score criterion already used. Single-distribution and multiphase
  fits were both affected. A candidate that adds no column, such as `z`
  added to `~ z:f` (the columns `z:fa, z:fb` become `z, z:fb`), leaves the
  likelihood unchanged; it was reported with a small p-value and entered,
  and is now an error.

  Single-distribution fits also now accept the one-column terms multiphase
  fits already did. A logical, two-level factor or character candidate used
  to stop under `criterion = "wald"` with "not found in the design matrix",
  although the score criterion's refusal of such a column names
  `criterion = "wald"` as the way to test it. It is now tested.

* **`predict()` on a multiphase fit now returns an unnamed vector.** For
  `type = "cumulative_hazard"`, `"survival"` and `"hazard"`, every element
  was named after a parameter -- `"constant.log_mu"`, say, on each row --
  because the phase parameters are named elements of `theta` and R carried
  a name onto the prediction: from any phase without covariates, and, for a
  single row of `newdata`, from almost any phase parameter. The values were
  right; the names were meaningless, and they followed the result into
  anything built from it. With one row, the `decompose = TRUE` and
  `se.fit = TRUE` data frames also took such a name as their row name.
  Single-distribution `predict()` carried names the same way; see the next
  item (#309).

* **`predict()` on a single-distribution fit now returns an unnamed vector
  too.** For `type = "survival"` and `"cumulative_hazard"`, a lognormal fit
  named every value `"mu"`, and a Weibull, exponential or log-logistic fit
  did the same for a single row of `newdata` (#309). As above, `mu` is a
  named element of `theta`, and R carried the name onto the prediction,
  through `rep()` for the lognormal and through any length-1 operand when
  there is one row. With one row, the `se.fit = TRUE` data frame also took
  `"mu"` as its row name. The values were right. `type = "hazard"` and
  `"linear_predictor"` were already unnamed.

* **`predict(newdata = )` no longer lets a design column override the
  formula variable it contradicts.** `newdata` may give a factor as its
  level (`grp = "old"`) or as the fit's design column (`grpyoung = 1`).
  When it held both and they disagreed, the design column silently won: a
  Weibull fit of `~ age + grp` returned the "young" cumulative hazard,
  0.362, for a row given `grp = "old"`, where the answer is 0.180. When
  the formula's variables are all present, `newdata` is now rebuilt from
  them and a design-named column is an unused extra. Design columns alone
  are used when no formula variable is given: a fit saved by an earlier
  version, or `newdata` given as design columns only. Some variables
  beside the design columns, with others missing, is now an error
  rather than a guess, because either route would ignore part of it. That
  includes a changed numeric `age` beside columns built from it, such as
  `I(age^2)` or `age:grpyoung`, which would otherwise be used stale: a
  copy of the fitted design with `age` edited gave 0.362 for
  `~ age * grp` where the answer is 0.660. `hzr_deciles()` and
  `hzr_gof()`, which evaluate at fitted design rows or their means, declare
  that themselves, so `hzr_gof()` still reports the curve at `mean(age^2)`
  for an `I(age^2)` term, not at `mean(age)^2` (#272).

* **`predict(newdata = )` no longer matches a single-distribution model's
  covariates by column position.** For `dist = "weibull"`,
  `"exponential"`, `"loglogistic"` and `"lognormal"`, the covariates in
  `newdata` were multiplied into the coefficients in the order they
  appeared, whatever their names. Reordered columns gave a wrong answer
  with no error: a Weibull fit of `~ age + mal` given `newdata` with `mal`
  before `age` returned a cumulative hazard of 4.85e20 in place of 0.25,
  and a survival of 0 in place of 0.78. A `newdata` missing a covariate
  could also return a value. Every prediction type was affected, as was
  the time-varying expansion.
  Covariates are now matched by name, through the same design
  reconstruction as the multiphase fix below, so a factor can be given as
  a level label. `newdata` may instead carry the fit's design-matrix
  columns by name (`grpyoung`), which is how `hzr_deciles()` and
  `hzr_gof()` call it. A column the model does not use is ignored, and one
  it needs but `newdata` lacks is an error that names it, even when an
  object of that name exists in the workspace. An unused column also
  stays unused when it shares its name with a constant in the formula,
  such as `cutoff` in `I(age > cutoff)`. It used to replace the constant
  silently. A fit made with an unnamed
  `x` matrix still matches by position, since there is nothing else to
  match on, and a `newdata` with only a `time` column still evaluates the
  baseline (#267). This rejects some `newdata` that was accepted before;
  see Breaking changes.

* **`predict(newdata = )` now evaluates a multiphase fit that has both a
  global covariate and phase-formula covariates.** A phase without its own
  formula inherits the global design, but at `newdata` it was built from
  every non-time column, so the global phase received the phase formulas'
  variables as well as its own. No `newdata` could satisfy both kinds of
  phase: `hazard(Surv(t, d) ~ age, phases = list(early = hzr_phase(...,
  formula = ~ mal), constant = hzr_phase("constant")))` stopped with
  "non-conformable arguments" for every prediction type, and a factor global
  covariate stopped with a different error. Extra or reordered columns
  failed the same way. Such a phase is now rebuilt from the global formula's
  own terms, factor levels and contrasts, as `predict.lm()` does. A factor
  can be given as a single label, and data-dependent terms such as
  `scale(x)` and `poly(x, 2)` reuse the fit's centre, scale and basis
  instead of recomputing them from the new rows. The global formula now also
  finds a non-column variable (`cutoff` in `I(x > cutoff)`) in the
  environment the formula was written in, as `model.frame()` does; before,
  only a global variable was found. `hazard()` stores these in
  `object$data$x_design`. Fits made through the vector interface select
  their columns by name, or by position when `x` was unnamed. Point
  predictions, `se.fit = TRUE` and `decompose = TRUE` are each checked
  against `exp(x beta_j) H0_j(t)` per phase (#266).

* **`predict(newdata = )` on a multiphase fit with `time_windows` and a
  global covariate now returns one value per row.** The fit expands the
  design that a phase without its own formula inherits into one column
  per window (`age_w1`, `age_w2`). At `newdata` that design was rebuilt
  without that expansion, and meeting the per-window coefficients it
  returned a
  flattened matrix: four numbers for two rows, with no error. It is now
  expanded at the prediction times, and matches `predict()` at the
  fitted data.

* **A multiphase fit made through the vector interface, given a phase
  formula it did not use for fitting, now predicts from the design it was
  fitted on.** Only the formula interface builds a phase from its own
  formula. With `hazard(time =, status =, x =)` the phase inherits the
  global `x`, and `hzr_phase(formula = )` is ignored. `predict(newdata = )`
  rebuilt that unused formula anyway. With `formula = ~ log(age)` and
  `x = cbind(age)` it multiplied `log(age)` by a coefficient fitted on
  `age`: 0.046 0.206 0.525 0.525 where the fit gives 0.023 0.122 0.082
  0.153. Under `time_windows` it returned eight values for four rows.
  Neither gave an error. `predict()` now routes each phase the way the fit
  built it, from the fit's own record, as `hzr_gof()` already did; the two
  share one rule. A fit saved before that record existed is routed by its
  stored columns, so an old fit keeps the phase formula it was fitted with.

* **`predict(newdata = )` no longer lets `newdata` stand in for what a
  formula takes from outside `data`.** A formula can use a variable that is
  not a column of `data`: a constant, such as `cutoff` in `I(age > cutoff)`
  or spline knots, or an object with one value per fitting row. Both global
  and phase designs were rebuilt from all of `newdata`, which caused silent
  errors:
  - An extra column masked a constant. A `cutoff = 0` column turned
    `I(30 > 50)` into `I(30 > 0)`: 0.425 for 0.191 on a Weibull fit, and
    0.074 for 0.373 on a phase formula.
  - A row-level object kept outside `data` and absent from `newdata` was
    read from the fitting rows, in fitting order. That covered a vector,
    and a list, environment or data frame read with `$`, as in `~ ext$z`.
    With a one-row `newdata`, the prediction came back with one value per
    fitting row.

  `newdata` now supplies only the columns of `data`, so every other formula
  symbol comes from the formula's environment, and a term that uses
  row-level values from outside `data` is refused (see Breaking changes).
  A rebuilt design whose row count differs from `newdata`'s is refused too.

* **`predict(newdata = )` on a multiphase fit now codes a phase formula's
  factors as the fit did.** It rebuilt a phase's design with a bare
  `model.matrix()` at `newdata`, with no stored levels or contrasts, so a
  factor given as one label (`grp = "young"`, including a one-row `newdata`
  with a single character label) stopped with "contrasts can be
  applied only to factors with 2 or more levels", and a factor whose levels
  were in another order was coded against the wrong level with no error: a
  wrong cumulative hazard. The fit now stores each phase formula's terms,
  factor levels and contrasts (`fit$x_design`), and `predict()` rebuilds the
  phase's columns from them, matched by name. A level the fit never saw is an
  error. So is a covariate the phase uses that `newdata` lacks, where it was
  silently taken from a same-named object in the workspace; a `newdata` with
  only a `time` column still evaluates the baseline, every covariate at 0.
  `newdata` carrying only the phase's design columns by name (`grpyoung`) is
  taken as it is; when it also carries the formula's variables, the variables
  win, so `grp = "old"` beside `grpyoung = 1` is the old value, not the young
  one (#272). Some of the variables beside the design columns, with others
  missing, is an error, and so is a changed variable that another design
  column is built from (`age` beside a stale `age:grpyoung`) when the others
  are missing. A fit saved by an earlier version rebuilds from its variables
  as before whenever they are all given, and refuses some of them beside its
  design columns, as it errored before. It also refuses a phase covariate
  missing from `newdata` rather than taking a same-named object, when it
  kept its fitting data (saved by 1.1.0 or later). A fit saved by 1.0.3 or
  earlier did not keep it, cannot tell a missing covariate from a formula
  constant, and still takes a same-named object, as before.

* **`predict(type = "survival", se.fit = TRUE)` now reports the standard
  error of the survival probability.** The `se.fit` column held the standard
  error of the cumulative hazard, `se(H)`, bit-identical to the column that
  `type = "cumulative_hazard"` returns, under a survival label. It now holds
  `S * se(H)`, the delta-method standard error of `S = exp(-H)`, which is
  what `summary.survfit()` reports as `std.err`. For a Weibull fit of
  `Surv(int_dead, dead) ~ age + mal` to `na.omit(avc)`, at `time = 5`,
  `age = 60`, `mal = 1`, where `S = 0.700`, the old column read 0.0706
  against the correct 0.0495.
  Every path was affected: all four single distributions, multiphase fits,
  and `hzr_read_outhaz()` objects. The confidence limits were already right
  and have not changed, so the `PROC HAZPRED` parity of `lower` and `upper`
  still holds. `PROC HAZPRED` prints no standard error, so there was no SAS
  value for this column to reproduce.

* **The multiphase gradient and Hessian are now right when an early phase's
  `m` is near 0.** Both differentiate in `m` by finite differences, and
  their stencils straddled 0: the gradient's (half-width about 6e-6)
  whenever `|m|` was smaller than that, the analytic Hessian's (half-width
  1.2e-4) whenever `nu > 0` and `|m|` was below 1.2e-4. `hzr_decompos()`
  changes formula at `m = 0`, and the `m < 0` family meets the `m >= 0` one
  in a cusp rather than continuing it, so each difference mixed two
  branches and returned neither side's derivative. The gradient gave +8.4
  where the true value was -20.2, on a 13-parameter fit that converged to
  `m = 2.9e-6`; the Hessian put -4.8e5 on the `m` diagonal where the value
  is about 0.7, so the standard errors of such fits were wrong too. Every
  other parameter's gradient was unaffected.

  Both stencils now keep the sign of `m`, one-sided on the `m >= 0` side,
  where the family is smooth: there the gradient and the analytic Hessian
  are now right. Below 0 the gradient's step is at most 1% of `|m|`, because
  the cusp varies on that scale, floored at 1e-10 so rounding stays bounded.
  The Hessian's steps below 0 are one-sided but still 1.2e-4 wide: they stop
  the branches mixing, but for `m` within about 1e-4 below 0 and `nu < 2`
  they understate the curvature, so the standard error of `m` there is too
  large (about five times, on the fit the tests use). At `nu = 0` that path
  used to stop with an error; it now returns these values. Fits whose
  standard errors come from `numDeriv` -- those with left- or
  interval-censored rows -- and the score test's information still
  difference across 0. Likelihood values are unchanged.

  The likelihood itself is still not differentiable at `m = 0`, so a fit
  whose optimum sits there reports a nonzero gradient. SAS/C never meets the
  point: it estimates `log|M|` with the sign fixed by the starting value, so
  `M` cannot reach or cross 0. `hazard()` estimates `m` directly and can.

* **A fit that reports convergence is now checked against SAS/C HAZARD's
  own test for it, and continued with `stats::nlm()` when it fails the
  test.** `hazard()`'s BFGS optimizer stops on the relative change in the
  log-likelihood (`control$reltol`, default 1e-5), which lets a flat ridge
  end short of the maximum with `converged = TRUE`: a 13-parameter
  early-CDF plus late-G3 model stopped 0.013 below the SAS listing's
  log-likelihood, and synthetic fits of the same shape up to 5 units below.
  SAS/C accepts an optimum only when the relative gradient,
  `max |g_i| * max(|x_i|, 1) / max(|f|, 1)`, is at most `eps^(1/3)`, about
  6e-6. When BFGS reports convergence and that test fails, every
  distribution's fit is now continued with `stats::nlm()`, the
  Dennis-Schnabel algorithm SAS/C's optimizer was ported from, at SAS's
  tolerances, and the continued point is kept only if the log-likelihood
  improves. The default `reltol` is unchanged: tightening it instead cost
  30% to 60% more time on the test suite and broke eight or nine tests.

  Every fit records the test in `fit$fit$rel_gradient` (`NA` when the test
  was not applied, because the optimizer did not report convergence, or the
  gradient cannot be evaluated) and, when the continuation improved the fit,
  `nlm()`'s termination code in `fit$fit$polish_code`, and `print()` and
  `summary()` show it. Under Conservation of Events the analytic score omits
  how the conserved scale moves with the other parameters, so there the test
  is computed from finite differences of the log-likelihood, as SAS/C does.
  The continuation keeps the analytic score, so a CoE fit can honestly end
  with the test not met. Only SAS/C's two hard failures warn: code 4, the
  iteration limit, and code 5, where the likelihood kept rising along some
  direction and may have no maximum. Codes 2 and 3, where SAS/C prints a
  caution and retries, are recorded without a warning. The test is relative
  to the size of the log-likelihood, so a fit that meets it is within SAS's
  tolerance of the maximum rather than exactly at it.

  Estimates of fits that used to stop short now change. One test depended
  on a detail of where BFGS stopped: it showed `gamma` and `eta`
  non-identified at `alpha = 1` by a large standard error. On that exactly
  flat ridge the Hessian is singular in theory, so whether a finite standard
  error comes out at all is numerical noise, and at the polished point it
  does not. The test now accepts either a missing or a 100-fold larger
  standard error, and also checks that both fits reach the same
  log-likelihood and the same `gamma * eta`.

* **Exponential, log-logistic and log-normal fits ignored left truncation**
  (#253). On status 0/1 rows these three families used `time_lower` only as
  a censoring bound, which applies to status 2, so a left-truncated fit was
  silently fitted as if every subject had been at risk from time 0. They
  now subtract the cumulative hazard at entry, H(time) - H(time_lower), as
  the Weibull and multiphase likelihoods already did. The log-likelihood,
  its gradient and the closed-form Hessian all carry the entry term.

* **`hzr_gof()` reported a conservation ratio that was not one** (#254).
  For a model with covariates it computed expected events from a single
  curve at the covariate means, then printed the total as the
  "Conservation ratio (E/O)". On the covariate model in the clinical
  walkthrough vignette that printed 0.606, while the fit conserved events exactly (68.000 expected
  against 68 observed). Expected events are now summed per subject, each
  subject's cumulative hazard at exit minus that at entry. For Weibull,
  exponential and multiphase fits with conservation of events, E/O is then
  the conservation-of-events identity; for log-logistic and log-normal fits
  it checks calibration in total. For a weighted fit, both observed and
  expected events now carry the case weights, since that is what a weighted
  fit conserves (the sum of w·H equals the sum of w·d). The `par_surv` and
  `par_cumhaz` columns are still the covariate-mean curve, for plotting
  against Kaplan-Meier, and the risk-set counts and Kaplan-Meier columns stay
  unweighted. Unweighted intercept-only fits without entry times are
  unchanged.

* **`hzr_gof()` drew the mean-patient curve at covariates of 0** for a
  multiphase fit whose covariates enter only through the phase formulas.
  Such a fit has no global design matrix, so the `par_surv` and
  `par_cumhaz` columns were predicted from time alone, which set every
  phase covariate to 0. They now use each phase's design-matrix column
  means, so a factor enters as the proportion of patients in each level.

* **`hzr_gof()` stopped on a multiphase fit with both a global covariate and
  phase-formula covariates, and every multiphase fit carried a
  `par_cumhaz_time` column** (#263, #264). With `Surv(...) ~ age` and a phase
  formula `~ mal`, the mean-patient curve was built from the global
  covariates alone, so `predict()` found no `mal` column and stopped before
  the expected-event tally. The curve is now evaluated at the column means
  of each phase's own design matrix, for global, phase-formula and mixed
  covariates alike. With `time_windows`, a multiphase fit's output had twice
  as many rows as grid times. The mean patient now carries the covariate
  means in the window that contains each time, for every phase built on the
  global covariates, including one whose formula the fit could not evaluate
  without `data`. A multiphase fit that dropped rows with a missing phase
  covariate is now refused: its design matrix is shorter than the data, and
  the per-subject tally recycled it and gave a wrong total with only a
  length warning. Separately, the `par_cumhaz_<phase>` columns were chosen
  by dropping `total` from the decomposition, which let its `time` column
  through as a phase; they are now chosen by phase name.

* **`hzr_gof()` had five smaller errors**, found by Copilot's and
  r-reviewer's reviews of #285.
  - `seq()` builds grid times that differ from the data times in the last
    binary digits, and they were matched with a fixed tolerance of 100
    machine epsilons. Above 128 that is smaller than the gap between
    adjacent doubles, so at times in days or months events fell off the
    grid: 180 of 197 were counted on a `seq(150, 300, by = 0.1)` grid. The
    tolerance now scales with the time.
  - With a custom `time_grid`, `n_risk` between Kaplan-Meier times carried
    the previous count forward, so it kept subjects who had left and
    missed ones who had entered. It now counts the risk set at each grid
    time.
  - An unsorted grid made the cumulative columns non-cumulative. The grid
    is now sorted, with repeated times dropped.
  - With an event at time 0, `km_surv` at 0 was averaged with 1.
  - For a fit with both `time_windows` and entry times, expected events
    took H(entry) in the entry-time covariate window, while the likelihood
    uses the exit-time window. E/O came out 1.052 on a Weibull fit that
    conserves events.

* **`hzr_gof()` places each subject at the Kaplan-Meier time `survfit()`
  gave it** (#286). `survfit()` merges exit times that are closer together
  than its tolerance. The per-subject tallies added above for #254 matched
  each subject's raw time instead, so a merged subject fell off the default
  grid and out of both tallies. On a fit with entry times, 55 of 68 events
  were counted and E/O read 1.22. This came in with the #254 change and
  never shipped. `hzr_gof()` now also warns when a subject cannot be placed
  on the default grid, rather than leaving it out silently.

* **A Weibull fit with one masked variance reported the others on the wrong
  scale.** When the Hessian inverse has a non-positive variance, its row and
  column are set to `NA`. The delta-method transform from the internal
  `(alpha, psi)` scale to the reported `(mu, nu)` scale was then skipped for
  the whole matrix, so the surviving standard errors stayed on the internal
  scale beside `(mu, nu)` estimates, with no error. A masked `alpha` printed
  `SE(psi)` as the standard error of `nu`, off by a factor of `nu`; a masked
  `psi` printed `SE(alpha)` as the standard error of `mu`, which depends on
  `psi` and has no valid standard error there. The transform now runs on the
  finite block and carries the mask through the Jacobian: a reported
  parameter is `NA` exactly when it depends on a masked one. Exponential,
  log-logistic and log-normal fits are unaffected; they report on the scale
  they are optimised on.

* **A `survival::Surv()` object passed as `status` is now translated, as the
  formula interface always did** (#226). `Surv()` codes censoring with
  different integers from this package, and the vector interface took the
  object's second column unchanged. Under `type = "left"` a left-censored row
  was fitted as right-censored; under `"interval"` and `"counting"` the
  second column is not the status at all, so the fit read `time2` or `stop`
  as status codes. There was no error and no warning, and
  `objective = "sas"` could not see a left-censored row to refuse it. Both
  interfaces now read the `Surv` through one internal helper, driven by its
  `type`, so they store the same status and bounds and give the same fit.
  The bounds a `Surv` carries are taken from it; a `time`, `time_lower` or
  `time_upper` that disagrees with them is an error rather than being
  silently replaced. `hzr_bootstrap()` resamples those bounds too, although
  they never appear in the stored call.

* **`.` in a `hazard()` formula no longer puts the response in the design**
  (#273). `Surv(int_dead, dead) ~ .` expanded `.` to every column of `data`,
  including `int_dead` and `dead`, so the outcome was fitted as a predictor.
  With starting values sized for those extra columns the fit converged, with
  no error and a log-likelihood far above the correct model's. With starting
  values sized for the real covariates it stopped with "non-conformable
  arguments", which did not
  name the cause, and `predict(newdata = )` demanded the response columns.
  `.` now means every column the `Surv()` term does not use, as in
  `survival::coxph()`, so a `~ .` fit gives the same design and estimates as
  the formula written out in full. A `data` with no other column gives a
  model with no covariates, and there `.` beside other terms is an error.
  **Estimates from an earlier `~ .` fit change**, and so does the length of
  `theta` it needs. A `.` in `hzr_phase(formula = )` is fixed separately
  (#277). A right-hand-side variable that is not a column of `data` is now
  looked up where the formula was written, so a
  variable local to the calling function resolves instead of failing with
  "object not found". `hzr_bootstrap()`, which resamples only the rows
  of `data`, refuses such a fit (#278).

* **`hzr_phase(formula = ~ .)` no longer puts the response in the phase
  design** (#277). A phase formula's `.` was expanded by `model.frame()`
  against every column of `data`, including the columns of the `Surv()`
  term, so the outcome was fitted as a phase covariate and the fit
  converged with no error. `hazard()` now writes `.` out once, before
  fitting, the same way as for the global formula (#273): every column the
  `Surv()` term does not use. The fitted object stores the written-out
  formula, so `predict(newdata = )` no longer needs the response columns.
  On the vector interface (`time =`, `status =`) no `Surv()` term says which
  columns hold the response, so a phase formula with `.` is now an error
  there; write the phase's terms out. **Estimates from an earlier fit with
  `.` in a phase formula change.**

* **`hzr_argument_mapping()` listed DELTA as implemented.** Its
  `implementation_status` was `"implemented"` and its `r_parameter` read
  "(absorbed by decompos)", while the row's own notes say a non-zero DELTA
  is refused or flagged and never fitted (#181). The row is now
  `"planned"` with `r_parameter` "(not implemented)", so
  `hzr_argument_mapping(include_planned = FALSE)` no longer includes it.

* **A multiphase stepwise step now adds to the covariates a phase inherits,
  instead of replacing them** (#284). A phase with no formula of its own
  uses the global formula's covariates. Entering a variable into such a
  phase built a formula holding only the new variable, so `early.age`
  became `early.mal`: the refitted model had lost a covariate, its
  log-likelihood fell (-196.44 to -204.43 on `avc`, with the control the
  tests use), and the step reported `mal` entering at p = 5.9e-05, a
  p-value from that smaller model. The
  model the step describes, `early ~ age + mal`, gives p = 5.4e-04. The
  step now starts from the inherited terms, so it fits that model and
  reports its p-value. Drops start from the inherited terms too, and the
  inherited covariates are now drop candidates. The default score
  criterion used to stop on this case, reporting that the candidate "could
  not be added to the model"; it now scores it. Two related defects are
  fixed with it. Dropping a phase's last covariate set its formula to
  NULL, which made the phase inherit the global covariates again, so a
  later "drop" could bring `age` back; it now leaves `~ 1`. And the score
  test pinned a candidate in the phase's last slot, which is wrong after an
  interaction: `model.matrix()` puts main effects first, so with
  `age * mal` it scored `age:mal` in the candidate's place. It now finds
  the candidate's column by name. That one affected phases with their own
  formula too. **Selected models and p-values change** for any multiphase
  screen that stepped a phase without a formula while the global formula
  had covariates, or scored a candidate for a phase with an interaction.

## Known limitations

* **A multiphase fit can come to rest on a discontinuity in the likelihood,
  and still report `converged = TRUE` (#448).** When a `cdf` phase's shape
  `nu` is driven towards zero, the phase's `(t_half/t)^(1/nu)` term acquires
  an exponent of order `1e15`, so the phase approaches a step at `t_half`. If
  `t_half` then comes to rest within a floating-point step of one or more
  observed event times, the log-likelihood is discontinuous there: on a fit of
  the shipped `avc` data a one-step change in `log(t_half)` moves the
  log-likelihood by 6 to 30 units, in no consistent direction, with five tied
  event times accounting for the whole of it. Such a fit reports
  `converged = TRUE` while failing the relative-gradient test by six orders of
  magnitude, so **read `rel_gradient` and the phase's `nu` before trusting a
  multiphase fit**, and treat a `nu` at the boundary as a warning that the
  estimates are not identified. This release does not change the behaviour:
  whether the reference `PROC HAZARD` reaches the same state on the same job
  has not been established, and that answer decides whether the fix is a
  parity break or a shared degeneracy.

* **In a two-way `hzr_stepwise()` screen, `$scope$frozen` can name a
  variable the final model excludes (#378).** With `direction = "both"`,
  each iteration makes a forward step and then a backward step, and the
  protected sets are fixed when the iteration starts, so a variable that the
  forward step freezes can still be dropped by the backward step that
  follows. Forward-only and backward-only screens are not affected: neither
  makes both steps in one iteration, so their `$scope$frozen` and final model
  agree. It is then reported as frozen while the selected model
  does not contain it, and nothing warns. When they disagree, **trust the
  final model and `$steps`**, which records both the `"frozen"` row and the
  `"drop"` after it; read `$scope$frozen` as the variables that reached the
  `max_move` cap, not as variables held in the model. This release does not
  change the behaviour: correcting the timing alone was measured to keep
  variables above `slstay` at the default `max_move`, so it is deferred to
  be fixed together with how moves are counted (#379).

# TemporalHazard 1.2.10

## New features

* New `hzr_repeated_events()` rebuilds the input to a repeated-events hazard
  model. From a long data set with one row per candidate event per subject, it
  returns one row per inter-event segment, with the segment's start time,
  duration and running event count. It reproduces the SAS macro `%repeat`,
  which built this input for the repeated-events `HAZARD` jobs and whose output
  was seldom saved, so those jobs can now be run again in R. It refuses input
  that would otherwise give a plausible but wrong result -- a non-numeric time,
  follow-up or indicator column, a factor `id`, a missing `followup` value, or
  an empty data frame -- and warns, naming the subjects, when an event falls
  after the end of follow-up, when `followup` varies within a subject, or when
  a missing time leaves a segment undefined. As in the macro, `rcensor` and
  `event` can both be 1 on the same row; see `?hzr_repeated_events`.

## Bug fixes

* **`hzr_translate_sas()` now starts an unspecified shape parameter where
  `PROC HAZARD` starts it, not where `hzr_phase()` does.** A `PARMS` statement
  that named only some of a phase's shape operands had the rest filled from
  `hzr_phase()`'s defaults, and three of them disagree with the reference C
  (`src/hazard/stmtprc.c`): `NU` starts at 2 rather than 1, `M` at 1 rather
  than 0, and `ETA` at 2 rather than 1. Since the multiphase likelihood is
  multimodal, a different starting vector can reach a different optimum, so
  this was a fidelity divergence rather than a cosmetic one. The defaults are
  now `PROC HAZARD`'s, and the emitted `hzr_phase()` call names every shape
  argument explicitly so that what is printed and what reaches the optimizer
  cannot disagree. `TAU` is the exception, because `PROC HAZARD` derives it
  from the data: an unspecified `TAU` becomes `0.75 * Tmax`
  (`src/hazard/readobs.c`) and one written as non-positive becomes
  `2 * Tmax / 3` (`SETG3()`). Neither can be reproduced at parse time, so such
  a phase is emitted at `tau = 1` and recorded in `$untranslated`, naming
  whichever rule applies -- unless `SETG3_ignore_tau()` does, which pins `TAU`
  at 1 anyway. No job in the *public corpus* is partially specified, so no corpus
  translation changes; the package's own end-to-end fits test does carry a
  partial block (`MUE THALF NU MUC`, no `M`), and its early phase now starts at
  `m = 1`.

* **`hzr_translate_sas()` now pins `TAU` where `SETG3_ignore_tau()` pins it.**
  When `ALPHA` is fixed at 1, `setg3.c:378-379` sets `TAU` to 1 *and* fixes
  it; the emitted `hzr_phase()` call mirrored neither, leaving `TAU` free. At
  `alpha = 1` the `G3` form collapses to `(t/tau)^(gamma*eta)`, so `log_mu`
  and `log_tau` are exactly aliased: the translated fit converged onto a flat
  ridge and returned no standard errors for either, with nothing in
  `$untranslated` to say so. The emitted call now carries `tau = 1` and
  `"tau"` in `fixed`, which identifies `log_mu` again. Jobs whose `PARMS`
  named a different `TAU` additionally record a row, since that value is used
  by neither `PROC HAZARD` nor the translation. No corpus translation changes:
  every late-phase block in the public corpus already writes `TAU=1 FIXTAU`,
  so the emitted calls are byte-identical before and after (checked by running
  the parser over all of them). What changes there is that those jobs no
  longer need a warning.

  The same branch's `ETA` fix (`setg3.c:405`) is mirrored too: with `GAMMA` and
  `ETA` both free at `alpha = 1` the pair is exactly singular, and
  `hzr_phase()` previously left both free and fitted the ridge. On the
  package's own fixture that moves `gamma`'s standard error from 7.5 to 0.02.
  This replaces an `$untranslated` row with a faithful translation.

* **`hzr_translate_sas()` now reports the whole of what `SETG3()` would do to
  a late phase, and mirrors only the part that has to be mirrored.**
  `PROC HAZARD` does not optimize from the operands `PARMS` supplies: `SETG3()`
  rewrites them first, and refuses some jobs outright. The translator now walks
  that function (`src/model/setg3.c`) and records both, splitting them on a
  single principle:

  - A rewrite that resolves an **exact non-identifiability** is *mirrored*,
    because the degeneracy is algebra and is just as real in R. There are two,
    both in `SETG3_ignore_tau()`: the `TAU` pin and the `ETA` fix, described
    above.
  - Every other rewrite keeps `PROC HAZARD` inside a numerical branch it can
    evaluate -- the role `g3flag` plays, which `hzr_decompos_g3()` does not
    need because it carries the general four-parameter `G3` form. Copying those
    would import a SAS limitation into R, and reaching a late shape the
    reference implementation cannot is a purpose of this package. They are
    *recorded* in `$untranslated`, so a SAS parity run knows why the starting
    values differ.

  In practice this covers `SETG3_verify_ge_2()`'s push of `GAMMA * ETA` clear
  of 2 (which the `PROC HAZARD` defaults `gamma = 1`, `eta = 2` trip exactly,
  so it applied to every defaulted non-`WEIBULL` late phase), the value
  substitutions in all eight sign branches, and `SETG3_alpha_fixup()` /
  `SETG3_alpha_gener()` deriving `ALPHA` from `GAMMA * ETA`.

  Nine refusal codes can now be recorded rather than emitted as runnable fits,
  where previously only `SETG3980` was -- and that one only for `alpha = 0`,
  not for the negative `ALPHA` that raises it too. (All sixteen are mirrored
  from the C, but seven guard conditions the entry checks at `setg3.c:269-284`
  have already refused, so no input reaches them; an exhaustive search in the
  tests pins which nine are live.)

  One correction to a rule this package had recorded wrongly: an **unspecified**
  `TAU` does not reach `SETG3()` as the value `stmtprc.c` starts it at, 0.
  `src/hazard/readobs.c:153-154` replaces it with `0.75 * Tmax` on an active
  late phase, and `readobs()` runs before `SETG3()` (`hazard.c:276` against
  `:292`). So an absent `TAU` starts at `0.75 * Tmax`, not `2 * Tmax / 3` --
  that rule (`setg3.c:317`) governs only a `TAU` the job wrote as non-positive
  -- and a bare `FIXTAU` cannot raise `SETG3900`. `ALPHA = 0` **with**
  `FIXALPHA` is the limiting exponential in both implementations and still
  translates cleanly; left free, `SETG3` derives an `ALPHA` instead, which is
  now reported.

  No corpus job is affected: every late-phase block in the public corpus
  carries `WEIBULL`, which returns at `setg3.c:347` before the sign dispatch,
  and all of them write `TAU=1 FIXTAU ALPHA=1 FIXALPHA` with a positive `GAMMA`
  and `ETA`. Verified by running the parser over all of them before and after.
* **`hzr_translate_sas()` no longer builds a phase that `PROC HAZARD` would
  not.** A `PARMS` statement names its phases with `MUE`, `MUC` and `MUL`; the
  shape operands (`THALF`/`NU`/`M` early, `TAU`/`GAMMA`/`ALPHA`/`ETA` late)
  only shape a phase that already exists. The translator had this the other way
  round and keyed on the shape operands, so
  `PARMS MUE=0.2 THALF=1 NU=1 TAU=2 GAMMA=1.5` emitted a two-phase model
  against `PROC HAZARD`'s one -- carrying an invented `mu` starting value of
  0.1 for the phase that should not have been there, and offering nothing in
  `$untranslated` to say so. In the reference C only `setparmno()` sets
  `C->phase[n]`, and only when that `MU` is greater than zero
  (`src/hazard/setparmno.c`); the seven shape operands are registered by
  `setprmf()`, which never touches it. A phase is now built only when its own
  `MU` was specified and positive -- so a `MUC` of exactly zero no longer
  builds a constant phase either, where before the translator tested only
  whether the keyword was present. Shape operands belonging to a phase that
  never activated are recorded in `$untranslated` rather than dropped, since
  `PROC HAZARD` zeroes them (`src/hazard/stmtprc.c`) and skips their covariates
  (`src/hazard/setstat.c`). A job that activates no phase at all is now
  **refused** rather than translated -- whether its `PARMS` named no positive
  `MU`, or it carried no `PARMS` statement whatsoever. `PROC HAZARD` does not
  run such a job: `src/hazard/modterm.c` raises `ERROR 1001: No phase
  selected` and the procedure exits before computing any results, so there is
  no fit for a translation to be faithful to. The translator previously
  emitted a runnable single-distribution fit and reported full token coverage;
  it now emits a `stop()` in place of the `hazard()` call, as it already did
  for `LCENSOR` combined with `ICENSOR`. The two cases are one state rather
  than two: `src/hazard/stmtprc.c` zeroes all three phases at initialization
  and only `setparmno()` turns one back on, so a job with no `PARMS` has no
  active phase for the same reason a `PARMS` naming `MUE=0` does. `modterm()`
  is reached on every job, not only once a multiphase model has been selected
  -- its one call site in `outmods()` is unconditional in the procedure's main
  sequence.

  The refusal fires only when every `PARMS` operand was understood. A
  statement this parser could not read is recorded, operand by operand, but
  never refused: `PARMS MUE = 0.2 THALF = 1` (spaces around `=`) parses to
  nothing here while `PROC HAZARD`'s own lexer discards whitespace and runs
  the job with an active early phase, so refusing it would stop a job the
  reference accepts. A second `PARMS` statement also now adds to the first
  rather than replacing it, matching the single field table `parmprc()` reads
  once after all statements are processed.

  One job in the public `hazard` corpus changes, and only to drop a claim that
  was wrong: the second `%HAZARD` block of
  `dist/examples/hm.dthar.TGA.sas` is a documentation template carrying
  literal `?` placeholders, which the reference would reject as a syntax
  error rather than as `ERROR 1001`. (That file's first block is a valid
  `PARMS` activating two phases and is unaffected.) Its per-operand rows are
  unchanged. No corpus job is refused.

* **`hzr_translate_sas()` no longer discards the `ALPHA` and `ETA` a `PARMS`
  statement specified alongside `WEIBULL`.** The translator read the bare
  `WEIBULL` keyword as a request to constrain the late phase to
  `alpha = eta = 1` and pinned both. `SETG3_weibull()` in the reference C
  (`src/model/setg3.c`) is the *generalized* Weibull -- it admits all positive
  parameter values, validates `gamma > 0`, `eta > 0` and `alpha >= 0`, and
  assigns nothing. A production job carrying
  `alpha = 2.501719 ... eta = 0.1365255 weibull` therefore translated to a
  seven-parameter fit where `PROC HAZARD` estimated nine, with both starting
  values replaced by 1 and no warning. `WEIBULL` now leaves `ALPHA` and `ETA`
  as `PARMS` gave them, free unless an explicit `FIXALPHA`/`FIXETA` pins them.
  The `G3`-collapses-to-Weibull identity at `alpha = eta = 1` is real and
  unchanged; it was simply not what the keyword requests.

# TemporalHazard 1.2.9

## Bug fixes

* **The phase-identifiability warning no longer names a cause that did not
  occur.** `variation` is the relative range of a phase's contribution across
  the observed times, so it collapses for every phase when those times carry
  nothing that separates them. The saturated scan fired anyway: with
  `time = rep(2, 20)` it told a `constant` phase that its shape parameters were
  unidentified, when a `constant` phase has none, and blamed a short half-life.
  A single-row fit produced the same text.

  That case is now reported in its own words, under a condition that takes two
  things together: the observed times' own relative range is below the
  tolerance, *and* no phase's contribution varies across them. Neither half is
  sufficient. The range alone would condemn times bunched far from the origin,
  where a steep phase can be fully identified; the phase measures alone would
  condemn well-spread times whenever a single phase saturates. Because the
  range is relative, near-ties reach the condition as exact ties do --
  `c(100, 100.000001)` produced the old wording verbatim. Per-phase saturation
  and absence over well-spread times are reported exactly as before, and the
  condition is now caught even when every phase carries covariates, where the
  per-phase measures are withheld and previously nothing was said at all.

  Two related corrections. A phase that has not started is reported as absent
  even in a left-truncated fit, since that test rests on the share and not on
  how the times are spread. And when the counting-process entry times or
  interval bounds supply enough evaluation points, the degenerate-times verdict
  is withheld -- and "enough" is counted rather than assumed. With the event
  times tied, each added point supplies one more functional of the parameters,
  so one more than the number of added points must reach the number of **free**
  parameters. The package's default two-phase model has five, and pinning a
  phase's shapes lowers the bar because those parameters are no longer free.
  Points that add nothing do not count: a zero entry time, a bound equal to an
  event time -- both of which `Surv(type = "interval")` synthesises for every
  row -- or a bound the objective never reads, since `time_upper` is live only
  for left-censored and interval rows. Counting any of those silenced this
  warning on data that needed it while changing no likelihood value.

  The message also no longer names shape parameters for a model that has none.
  A single `constant` phase is now silent, since one functional determines its
  `mu` exactly; two `constant` phases still warn, since only their sum is
  determined, but the wording says the phases cannot be told apart rather than
  naming shapes they do not have.

  Where the times are NOT degenerate, the per-phase saturation message is
  reported regardless of the bounds, so adding a `start` column no longer
  removes a correct diagnostic (#211).

  Known limitation, unchanged by this work and stated rather than implicit: the
  per-phase measures are taken over the event times alone, and the count of
  evaluation points is deliberately conservative -- it does not credit the
  extra functional a zero entry time makes separately observable, so a fit can
  be warned about while being identified. That cuts both ways.
  A phase can be flat across the event times while its shape is identified
  through the bounds, so the saturated message can overstate what is lost; and
  the counting rule above is a *local* identification argument for tied event
  times, so a fit can clear it and still have a flat direction the guard does
  not see, in which case nothing is said at all. Tracked in #228.

* **The SAS nomogram parser no longer deletes a whole-year time column.** The
  `PROC PRINT` observation counter is identified as a leading column running
  exactly `1..n`, which rested on no measurement column ever being a gapless
  1-based integer run -- true of the corpus in hand, not of SAS listings. A
  counter-suppressed nomogram on a whole-year grid prints `YEARS` as `1.0000,
  2.0000, 3.0000`, and the parser deleted the time key. Because the remaining
  columns then matched the header count, the shape check that would have caught
  it passed: the listing parsed, and the comparison ran against a table with no
  time column.

  The header now decides: a leading column the header names as a counter is
  dropped whatever values it carries, which also fixes the reverse case -- page
  two of a paginated listing starts at 41, and requiring a `1..n` run left that
  counter in the data. The labels are matched case-insensitively, after the
  leading underscore the parser strips, so SAS's `_N_` is recognised. A `1..n`
  run under a name the parser does not know is now kept and warned about rather
  than silently dropped: the listing cannot distinguish a counter spelled a new
  way from a measurement that happens to run `1..n`, and an extra column is
  visible where a deleted one is not (#212).

* **A phase named `total` is now rejected rather than silently losing its
  contribution.** `total` is the key the multiphase cumulative-hazard
  accumulator is stored under, so `phases = list(total = hzr_phase("cdf"), ...)`
  overwrote that phase's own contribution vector with the summed hazard. The fit
  returned normally; the phase then reported a contribution share of exactly 1
  and the identifiability diagnostic read it as dominating the model. The name is
  reserved across the decomposition and prediction output as well, where it
  labels the total column. `hazard()` and `hzr_theta_names()` -- which applies
  the same validation -- now stop with an explanatory message, and the
  reservation is documented on the `phases` argument (#214).

* **`objective = "sas"` data defects are reported as data defects.** The two
  conditions the SAS objective imposes -- no left-censored rows, and a positive
  width on every interval-censored row -- are pure functions of the data, but
  were checked inside the objective. The optimizer's per-start handler caught
  them and reported "produced no usable fit from N starts", which reads as a
  convergence problem and invites raising `n_starts`: a remedy that cannot work,
  because the condition is identical at every start. `hazard()` now checks both
  before any optimization, and the interval-width message reports the offending
  row in the data rather than its position among the interval rows. The guards
  inside the objective and gradient are unchanged, since the gradient is
  reachable without `hazard()` (#213).

* **Behaviour change:** because those preconditions are properties of the data
  and not of the fit, they are now checked whenever `objective = "sas"` is
  specified, including `fit = FALSE`. `hazard(..., fit = FALSE, objective =
  "sas")` on left-censored or zero-width-interval data previously returned a
  `hazard` object carrying an objective it could never be fitted with; it now
  stops. An `NA` in `status` under `objective = "sas"` is also reported by
  argument and row, rather than as a bare `missing value where TRUE/FALSE
  needed`.

  Note this reaches only the codes `hazard()` is given. On the **vector**
  interface a `survival::Surv()` object is unclassed without translating its
  codes, so they mean something else on arrival, and what goes wrong depends on
  `type`: under `type = "left"` a left-censored row is coded `0` and is fitted
  as *right-censored*; under `type = "interval"` it is coded `2` and is read as
  an *interval* of zero width, which this guard then rejects; and a genuine
  interval row, coded `3`, matches no branch of the likelihood and contributes
  nothing. The formula interface translates and is guarded. That asymmetry is a
  pre-existing defect of the vector path, tracked in #226.

## Testing

* **The SAS parity tests for `hz.te123.OMC` fit 1 and `hz.tm123.OMC` still
  described the P1 #6 Conservation-of-Events gap that PR #65 closed, and their
  tolerances were set from that gap rather than from the code's behaviour.**
  Fit 1 asserted the log-likelihood within `0.2` where R and SAS now agree to
  1.3e-04, `hz.tm123.OMC` within `0.5` against 4.6e-04, and neither asserted `MUE`
  at all, each carrying a comment quoting a value the fix had already moved. Restoring
  the pre-PR-#65 defect -- dropping the entry-time term from
  `.hzr_conserve_events()` -- left every one of those assertions passing, so
  they could not have caught the regression they were nominally about. The
  log-likelihood tolerances are now `1e-5` (relative), both
  Conservation-of-Events intercepts `MUE` and `MUL` are asserted, and each fit
  first checks `conserve_applied`, since CoE disables itself silently on
  unsupported data. The same mutation now fails six assertions.

  `MUL` is compared as a ratio to the SAS value rather than directly, and that
  distinction is the point rather than a detail. `expect_equal()` divides by
  `mean(abs(expected))` only when that exceeds `tolerance`; `MUL` is 2.1e-04,
  below any tolerance worth setting for it, so a direct comparison silently
  becomes an *absolute* one -- `MUL = 0` passes at `5e-04`, and so does the
  regression above. A first version of this change asserted `MUL` directly and
  reproduced, in the fix, the defect it was removing. Against an expected value
  of `1` the comparison is relative, as the tolerance implies. `hz.te123.OMC`
  fit 2's pre-existing `MUL` assertion sat 7% above that branch point and moves
  to the same form. `hz.te123.OMC` fit 2's log-likelihood tolerance
  goes from `1e-2` to `1e-5` for the same reason; it carried no stale claim, but
  `1e-2` admitted a drift of 3.1 in log-likelihood on a quantity matching to
  1.5e-04. No package code changed.

# TemporalHazard 1.2.8

## New features

* **New `hzr_theta_names()`** returns the names of a multiphase `theta`
  vector, in the order `theta` requires, before any fit runs. Use it to check
  a hand-written starting vector against the specification it belongs to:

  ```r
  phases <- list(early = hzr_phase("cdf"), late = hzr_phase("g3"))
  stopifnot(length(theta0) == length(hzr_theta_names(phases)))
  setNames(theta0, hzr_theta_names(phases))
  ```

  `theta` is positional and its entries are not on a common scale -- the late
  phase logs `mu` and `tau` but carries `gamma`, `alpha` and `eta` naturally --
  and **wrapping the wrong element in `log()` produces a fit, not an error**.
  A comment describing the order is therefore not enough for a template that
  ships to many studies, which is what prompted this: the order is a property
  of the phase specification and changes the moment an author adds, removes or
  retypes a phase.

  The function is a thin wrapper over the naming the optimizer itself uses,
  not a second implementation of it. `hazard()`'s optimizer, the score test's
  re-expansion and this function now all go through one internal helper, so
  the order documented here cannot drift from the order a fit produces --- a
  test asserts the two are identical for three different phase
  specifications. Phase validation is shared too, so unnamed phases get the
  same `phase_1` / `phase_2` labels a fit will give them.

* **`hzr_stepwise()`'s `$steps` frame gains a `stat_type` column**, saying
  what the `stat` on that row is and so which reference distribution
  recomputes its p-value: `"score_q"` (chi-square on `df`), `"wald_z"`
  (standard normal) or `"wald_chisq"` (chi-square on `df`).

  `df` could not tell these apart. A scalar Wald is reported as a *z*, not as
  its square, so it and a score Q are both recorded at `df = 1` while calling
  for different distributions. The gap was widest on a candidate rescued by
  the Wald fallback added in 1.2.7: that row carries a Wald z under
  `criterion = "score"`, where every neighbouring row carries a Q. The
  selection was right and `p_value` was right, but a reader recomputing a
  p-value from `stat` the way the neighbouring rows permit got an answer wrong
  by dozens of orders of magnitude. Under `criterion = "score"` the *entry*
  rows reading `"wald_z"` are the ones the fallback rescued --- filter on
  `action == "enter" & stat_type == "wald_z"`. Drop rows are always
  Wald-tested under that criterion, following SAS, so they read `"wald_z"`
  whether or not the fallback fired.

* **`hzr_bootstrap()` reports Wald fallbacks in select mode**, through
  `$n_wald_fallback_replicates` and `$n_wald_fallbacks`, with a warning when
  either is non-zero. Every replicate runs under `suppressWarnings()`, so a
  run in which the fallback fired throughout previously reported nothing at
  all -- and the bootstrap is where a wholesale substitution matters most,
  since those entries drive the pooled selection frequencies.

* **The nomogram parser accepts `lines =` as well as `path =`.** A multi-fit
  listing needs each nomogram attributed to the fit whose block contains it,
  not to the file. Previously a caller who had already split the listing by
  fit had to write each block back out to a temporary file to parse it, which
  is a workaround rather than an interface. Passing the lines directly now
  works, and is what makes the multiple-nomogram warning actionable.

## Bug fixes

* **`objective` is now recorded on the fit, so refits keep the estimand.**
  `hazard(objective = "sas")` was accepted and acted on by the optimizer, but
  the choice was never stored on the object. Everything that rebuilds a
  `hazard()` call from `fit$spec` therefore reverted to
  `objective = "likelihood"`: every `hzr_stepwise()` candidate refit, the
  Wald-fallback refit and each accepted move, plus the score test's numeric
  Hessian and gradient.

  Selecting on a `"sas"` base fit consequently differenced `delta_logLik`,
  `aic` and `delta_aic` across *two different estimands*, and the score
  statistic was computed against a likelihood the fit had not been fitted to.
  On the esophagectomy reference the two objectives differ by about 22
  log-likelihood units -- larger than most single-variable effects -- and the
  run produced a full `$steps` table with no error and no warning.

  `fit$spec$objective` now records it, and one accessor feeds every consumer
  so they cannot drift apart. A fit that carries no `objective` predates the
  argument and is read as `"likelihood"`, which is what it was.

  Note the asymmetry this had created: `hzr_bootstrap()` *refit* mode
  re-evaluates the stored call, which did carry `objective = "sas"`, so it was
  already consistent; *select* mode goes through the scope refit and was not.

  Following from this, **`objective` cannot be overridden per refit.** Passing it
  through `hzr_stepwise()`'s `...` is accepted when it restates the base fit's
  own objective and refused when it conflicts, with a message saying why:
  a candidate refit under a different objective would make `delta_logLik`,
  `aic` and `delta_aic` differences between two estimands rather than between
  two models. It is refused rather than quietly ignored, because a full
  `$steps` table that silently disregarded an explicit argument is the failure
  this entry is about.

* **A candidate that neither criterion could test is now reported as such.**
  `criterion = "score"` declines a candidate whose observed information is
  indefinite at `beta = 0` -- which happens when the effect is *large* -- and
  then refits it and tests it by Wald instead. That rescue can converge,
  producing a perfectly good point estimate, while its Hessian is singular:
  no standard error, so the Wald test cannot be computed either.

  The fallback dropped that case silently. The row kept the *score's* reason,
  `information_indefinite`, which describes the first of two independent
  failures and says nothing about the second, and no counter moved. A
  strongly predictive variable vanished and the screen rendered as an honest
  "nothing met `slentry`" -- the same failure this package has twice fixed
  elsewhere, where a clean-looking screen and an honest null result cannot be
  told apart.

  Such rows now report `fallback_no_variance`, distinct from
  `information_indefinite` (the rescue errored or did not converge, and is
  listed in `$criteria$refit_failures`). `hzr_stepwise()` and
  `hzr_bootstrap()` both warn on either, saying the candidate was tested by
  **neither** criterion and which mechanism applied. The
  `information_indefinite` prose claimed "that refit also failed", which was
  wrong for the new case and sent readers to an empty `refit_failures`; both
  it and the bootstrap warning now describe the two mechanisms separately.

  No behaviour changed in what gets selected: a candidate that could not be
  tested is still not entered. What changed is that the run says so.

* **`fit$fit$weak` now distinguishes "no ridge" from "not checked".** The
  weak-identification diagnostic introduced in 1.2.7 returned `NULL` both when
  a fit had been examined and found well identified and when it could not be
  examined at all -- most importantly when no Hessian was available, which
  happens on an install without the suggested numDeriv package and on fits whose
  rows are left- or interval-censored, where the analytic Hessian declines by
  design. `NEWS` offered `fit$fit$weak` as the programmatic check, so on those
  installs it certified as well identified a fit nothing had looked at.

  The field now takes three values: a list when a ridge was found, `NULL` when
  the fit was examined and is well identified, and `NA` when the check could
  not run. Test it with `is.list(fit$fit$weak)` rather than `!is.null()`.
  `summary()` prints a note in the `NA` case saying the check did not run, and
  a fit imported with `hzr_read_outhaz()` -- which has no R Hessian to examine
  -- now reports `NA` rather than reading as certified clean.

* **A fit with more than one flat direction says so.** The detector reported a
  single direction, which invited reading every parameter it did not name as
  identified. It now reports `n_directions`, the number of distinct
  near-flat directions found, and the warning says when there is more than
  one. Distinctness is counted over the parameters spanning each direction,
  not over eigenvectors: a ridge between two parameters clears the gate twice,
  once on the flat direction and once on its stiff partner, so counting
  eigenvectors would report one ridge as two.

* **A rescued candidate that goes on to win is no longer fitted twice.** The
  Wald fallback refits each candidate it rescues, and the acceptance step then
  refit the winner again with identical arguments. The two fits were
  bit-identical, so this was cost rather than incorrectness -- but it doubled
  the price of every accepted fallback entry, against a criterion whose whole
  advantage is that it does not refit per candidate. The rescuing fit is now
  kept and reused, as the Wald path already did with its candidate fits.

* **`DELTA` is not implemented, and two comments said it was absorbed.** The
  headers of `R/decomposition.R` and `R/argument_mapping.R` both stated that
  the C `DELTA` parameter's time transformation
  `B(t) = (exp(delta * t) - 1) / delta` is "absorbed by `decompos()`". It is
  not absorbed; it is unimplemented, and `delta = 0` is assumed. `DELTA`
  enters the C reference in three separate places -- it builds `rho` from
  `B(t_half)` rather than `t_half`, it replaces the time argument with
  `B(t)`, and `delta * t` enters the log-density additively so the density
  carries a factor of `exp(delta * t)` -- and R computes the `delta = 0`
  branch of all three.

  The comment was the harmful part. It made the omission look deliberate and
  safe, so a reader looking for exactly this discrepancy was told to stop
  looking, while a `PROC HAZARD` job with `DELTA != 0` was reproduced against a
  different function with no error. Both comments now say what is true.

  The SAS-facing paths now distinguish the two cases rather than treating
  `DELTA` as one unmapped keyword. `hzr_read_outhaz()` already stopped on a
  non-zero `DELTA`; `hzr_translate_sas()` now records `PARMS DELTA = <nonzero>`
  as untranslated with a reason saying the emitted call fits a *different*
  model, and the `.lst` natural-estimates parser warns when a listing carries
  one. `DELTA = 0` and a bare `FIXDELTA` are treated as faithful translations,
  because that is the branch R implements -- previously both the safe and the
  unsafe case produced the identical generic note "PARMS keyword has no phase
  target", which distinguished nothing.

* **A multiphase fit now records whether Conservation of Events was actually
  applied.** CoE counts exact events, so it is disabled whenever any `status`
  falls outside \{0, 1\} -- which interval or left censoring guarantees -- and
  whenever the model has fewer than two phases. That is deliberate and
  correct. What was missing is that **nothing on the returned object said it
  had happened**: a caller who passed `control = list(conserve = TRUE)` got a
  fit carrying `conserve = TRUE` over a computation that did not run.

  It is not an edge case. `ICENSOR` appears on 42 to 74 blocks per production
  study, and `ICENSOR` guarantees `status` leaves \{0, 1\}, so the auto-disable
  fires constantly. Production also writes `NOCONSERVE` on 16 to 68 blocks per
  study, so R and SAS usually agree on the *outcome* -- but for different
  reasons, and a job carrying both `CONSERVE` and `ICENSOR` is exactly where
  the two could diverge unobserved.

  A fitted multiphase object now carries two new fields, both alongside the
  requested `conserve` under `fit$spec$control`:

  * `fit$spec$control$conserve_applied` -- logical, whether CoE was actually
    applied;
  * `fit$spec$control$conserve_disabled_reason` -- one of `"not_requested"`,
    `"unsupported_censoring"`, `"single_phase"`, `"no_events"` or
    `"setup_failed"`, and `NA` when CoE was applied.

  Read `fit$spec$control$conserve_applied`, not `fit$spec$control$conserve`:
  the latter says only what you asked for. `conserve` is a
  `dist = "multiphase"` control; the single-distribution fits do not use it.

  The reason is recorded rather than a bare logical because the causes want
  different responses -- and a bare `FALSE` reads as "you turned it off" to a
  user who did the opposite.

* **The SAS `.lst` nomogram parser no longer returns the first of several
  tables silently.** `.hzr_parse_sas_nomogram()` matched every nomogram header
  in a listing and read only the first, with no warning and nothing in the
  return value to say a second existed. That is the same shape as the three
  layout defects fixed in 1.2.1 -- silent, and discoverable only by pointing
  the parser at a second study.

  It now warns when a listing holds more than one, and attaches `n_found` to
  the returned frame whatever the count, so a caller can distinguish "one
  nomogram" from "the first of several" without reading the source. The
  behaviour is otherwise unchanged: the first table is still what comes back.

  Every file in the corpus this was found against happens to print exactly one
  nomogram, so the parser got the right answer there -- by luck of the corpus
  rather than by construction.

## Documentation

* **Corrected the documented paths for two fit fields.** The weak-identification
  result is at `fit$fit$weak` and the phase shares at `fit$fit$phase_share`;
  `NEWS.md` and, for the shares, `?hazard` had both named them one level too
  high. Code following the documented recipe got `FALSE` from
  `is.list(fit$weak)` for every fit, ridge or not -- the same wrong answer the
  1.2.8 three-value change was made to prevent, reached by a different route.

* **`stat_type` no longer over-claims which rows the Wald fallback rescued.**
  Under `criterion = "score"` drop rows are always Wald-tested, following SAS,
  so they read `"wald_z"` whether or not the fallback fired. The rescued rows
  are the *entry* rows: filter on
  `action == "enter" & stat_type == "wald_z"`.

* **`predict()` and `hzr_nelson()` now say that SAS draws narrower bands.**
  SAS `%KAPLAN`, `%NELSONT` and `PROC HAZPRED` all take their band width from
  `CLEVEL`, whose default is `0.68268948` -- documented in the macro source as
  "(1 sd)". That makes the multiplier `1` to seven decimals -- the literal is
  truncated -- so the band is one standard error, 68.3%, not 95%.

  Nothing here computed the wrong thing: parity is tested and passing, and
  `hzr_kaplan()` already documented the convention. The gap was that the other
  three entry points did not, and they are the ones a reader meets when
  checking an R fit against an existing SAS figure. At the R default of
  `level = 0.95` the reproduced band is about 1.96 times wider than the one
  being checked against, with no error on either side -- so the two look like
  they disagree numerically when they do not.

  The defaults are unchanged. `0.95` is the right R-side default, and adopting
  SAS's silently would make `predict()` disagree with every other R modelling
  function. The help pages now carry the level to pass instead:
  `level = 2 * stats::pnorm(1) - 1`.

* `summary()`'s documentation and the *Inference and diagnostics* vignette
  both listed the notes the method prints and had not been updated for the
  ridge note. Both now include it, and the vignette says plainly that a flat
  direction means the point estimates along it are unreliable, not only their
  standard errors.

* **`hazard()` now states that SAS's `STEEPEST` has no equivalent.**
  `PROC HAZARD` jobs write `STEEPEST QUASI` together -- steepest descent, then
  quasi-Newton -- and `STEEPEST` appears 14 to 109 times per study across the
  corpus. `QUASI`/`QUASINEWTON` maps to `method = "bfgs"`; there is no
  steepest-descent option and no two-stage strategy. Since the multiphase
  likelihood is multimodal, a different descent path can land on a different
  optimum, so a fit translated from such a job may not reproduce SAS's
  estimates. `hzr_translate_sas()` already recorded the keyword as
  untranslated rather than dropping it; the `control$method` documentation now
  says why.

## Internal

* Two tests in `test-score-wald-fallback.R` did not catch the mutations their
  comments named. Widening `.hzr_score_fallback_reasons` to include `constant`
  and `collinear` left both green, because a degenerate candidate still fails
  to enter -- it merely costs a refit on the way out -- and the "noise stays
  out" assertion is guarded by `slentry` rather than by how narrow the
  fallback is (the fixture's own Wald p-values are 0.0997 and 0.149, so a
  fallback that refit everything would still decline both). Both now assert
  `n_wald_fallbacks`, which is the quantity that moves.

* A roxygen block in `R/score-test.R` bound to the character vector declared
  after it rather than to the function it documents. `@noRd`, so no Rd was
  affected; source readability only.

# TemporalHazard 1.2.7

## New features

* **A fit sitting on a likelihood ridge now says so, and names the
  parameters.** An ill-conditioned Hessian already warned that standard
  errors were unreliable. That understates the problem when the
  ill-conditioning is a ridge: the likelihood is near-flat along some
  combination of parameters, and there the individual *point estimates* are
  not determined by the data either -- only the combination is. The fit still
  reports `converged`, and the coefficient table still prints a number for
  every parameter, so nothing on the object signalled it.

  `hazard()` now warns, once per fit, naming the parameters that span the flat
  direction along with their correlation and the Hessian's `rcond`, and
  `summary()` prints the same note. The finding is recorded on the object as
  `fit$fit$weak`, so it can be checked programmatically rather than scraped from a
  warning. See the 1.2.8 notes below for the three values that field takes.

  The direction is read off the *correlation* of the estimates rather than
  their raw covariance. Parameters here sit on very different scales -- an `m`
  of 27 against a `nu` of 0.027 -- and in raw units a direction that moves
  both equally in statistical terms loads almost entirely on the larger one,
  which would report a two-parameter ridge as a single unidentified parameter.
  A parameter that is merely imprecise, without trading off against another,
  is deliberately not reported: that is ordinary low precision, and the
  existing `rcond` warning and the parameter's own standard error already
  cover it.

  The check is generic -- it runs for every distribution and knows nothing
  about phase shapes -- and is gated on the `rcond` threshold the package
  already uses, so it never fires where the ill-conditioning warning stays
  silent.

## Bug fixes

* **The score criterion no longer declines a candidate for being too
  predictive.** `criterion = "score"` computes SAS HAZARD's Q exactly --
  `Q = grad^2 * I22`, the reciprocal Schur complement of the *observed*
  information at `beta = 0` (`src/vars/q1.c`). When a candidate's true effect
  is far from zero the log-likelihood is convex there, the Schur complement
  turns negative, and Q is undefined. The criterion therefore declined
  candidates in proportion to how predictive they were: on a fixture with one
  planted effect (`beta = 0.9`, LR = 178) and two pure-noise columns, the
  screen entered both noise columns and never tested the real one (#130).

  This is not a deviation from the reference -- it is inherited from it. SAS
  documents the same failure in `q1.c` ("IT IS POSSIBLE THAT THE PROGRAM WILL
  RETURN A NEGATIVE Q VALUE ... THE USER SHOULD USE THE MORE EXPENSIVE Q2 AS
  AN ALTERNATIVE"), and `dqstat.c` declines the candidate with `p = 1`. `Q2`
  is named once in the C tree and never implemented.

  So Q itself is unchanged and stays bit-faithful; only the *handling*
  diverges, and only where SAS says its own answer is unusable. A candidate
  the score cannot test is now refit and tested by Wald -- the substitute the
  unbuilt `Q2` was for. This is a deliberate, documented divergence from the
  reference implementation.

  The fallback is deliberately narrow: it applies to `information_indefinite`
  and `coefficient_diverging`, the two causes that mean "the approximation at
  zero broke down". Collinear, constant and non-numeric candidates are still
  declined without a refit, so the screen keeps the speed advantage that the
  score criterion exists for -- the cost is paid only on the few candidates
  that trip it. The returned object's `$criteria` gains `n_wald_fallbacks`,
  so the substitution is reported rather than silent, and a fallback refit
  that fails is recorded in `$criteria$refit_failures` and warned about
  rather than leaving a row indistinguishable from one never refit.


# TemporalHazard 1.2.6

## Bug fixes

* `.hzr_parse_sas_nomogram()` no longer discards a nomogram whose PROC PRINT
  counter column is labelled `OBS` rather than `Obs` (#184). The counter was
  dropped by exact name, so on the other casing it stayed among the header
  names, the row-width guard rejected every data row, and the parser returned
  `NULL` -- indistinguishable from a listing that printed no nomogram at all.
  A corpus sweep therefore reported its own parse failures as gaps in the SAS
  output. The counter is now identified structurally, as a leading column
  running exactly `1..n`, so any label parses.

* `.hzr_parse_sas_nomogram()` now warns rather than returning `NULL` in
  silence when a nomogram header matched but no data rows could be read.
  `NULL` again means "no such table", and only that.

# TemporalHazard 1.2.5

## New features

* A multiphase fit now warns when a phase has effectively left the model. Such
  a phase is silent in every other way: the fit converges, reports no trouble,
  and the affected parameters simply drift.

  Two modes are distinguished, because the consequences differ. A phase that is
  **absent** -- contributing essentially none of the cumulative hazard at any
  observed time -- has not started by the end of follow-up, and neither its
  `mu` nor its shape is identified. A phase that is **saturated** -- one whose
  contribution is constant across the observed times, typically a `cdf` phase
  whose half-life is far shorter than the first observation -- has already
  finished, and then acts as a constant offset: its `mu` stays well identified
  while the shape parameters (`t_half`, `nu`, `m`) go exactly flat. Pinning
  those at any value leaves the log-likelihood unchanged.

  The distinction is the point. It is tempting to describe a phase that
  supplies no late hazard as one whose `mu` has stopped being identified;
  `mu` is in fact the one parameter that survives, through the offset the
  phase already contributed. The share is measured against the cumulative
  hazard rather than the instantaneous hazard for the same reason.

  The shares are recorded on the fit as `fit$fit$phase_share`, so the warning can
  be checked rather than taken on trust, and the threshold is
  `control$phase_share_tol` (default 1e-8).

# TemporalHazard 1.2.4

## New features

* `hazard()` gains an `objective` argument. The default, `"likelihood"`, is
  unchanged: interval-censored rows contribute the interval probability
  `log(S(l) - S(u))`. The new `"sas"` reproduces what `PROC HAZARD` actually
  accumulates for such a row -- the ordinary event-density term with the
  *instantaneous* hazard replaced by the *interval-mean* hazard over (l, u]:

  ```
    d * log[ S(u) * (Lambda(u) - Lambda(l)) / (u - l) ]
  ```

  which makes the three row types one family: right-censored contributes
  `log S(u)`, an exact event `log S(u) + log h(u)`, and an interval-censored
  row `log S(u) + log h_bar(l, u]`. Exact-event and right-censored rows are
  untouched by the switch, and it applies only to `dist = "multiphase"`.

  **This is a different estimator, not a reparameterization, and must not be
  used for new analyses.** It is a density, not a probability, so it is
  inconsistent for wide intervals -- on a 12-year-interval reference fit the
  two forms differ by 22 log-likelihood units. It exists to reproduce legacy
  `PROC HAZARD` runs, and it is deliberately an explicit top-level argument
  rather than a `control` element, because it changes the estimand.

  Interval-censored rows with `u <= l`, and any left-censored row, are errors
  under `"sas"` rather than silently-substituted values: `PROC HAZARD` has no
  left-censoring statement, so no SAS run corresponds to such a result.

* New dataset `uslife2023`: the NCHS United States life table for 2023 on a
  synthetic 100,000 radix, 124 rows, every one interval-censored and exactly
  one year wide. Published aggregate counts only. It is the reference fixture
  for `objective = "sas"`, which reproduces its SAS log-likelihood of -410414
  at the printed estimates and at three off-optimum points of SAS's own
  iteration trace.

## Internal

* The interval-censored contribution was written twice -- once in the
  log-likelihood and again in the finite-difference closure inside the
  gradient. Those copies had to agree or the optimizer would step by the
  gradient of a different objective than it evaluated. Both now delegate to a
  single `.hzr_logl_interval()`. Behavior under the default is unchanged and
  bit-identical, log-likelihood and gradient alike.


# TemporalHazard 1.2.3

## Bug fixes

* A logical column no longer stops a stepwise screen dead. `hzr_stepwise(scope
  = NULL)` enumerates its own candidates and counts logical columns among them,
  on the grounds that a 0/1 field arriving logical rather than numeric is a
  property of the reader that produced the frame, not of the variable. Both
  criteria then refused what the package had offered: the score criterion
  errored with "is not numeric (logical)", and the Wald criterion failed
  looking up `phase.var` when `model.matrix()` had named the column
  `phase.varTRUE`. Either way the screen stopped before its first step, on a
  column nobody had chosen by hand.

  A logical candidate is now modelled as the 0/1 predictor it is, and gives
  the same screen as the identical numeric column: same variables entered, in
  the same order, at the same p-values, with the same coefficients.

  The coefficient-name half of this also fixes a two-level factor, which
  expands the same way (`varb`). A candidate that expands to more than one
  column is still refused, unchanged.

  That makes `criterion = "wald"` a real answer for a two-level factor or
  character column named in an explicit `scope`, which the score criterion
  still cannot expand. Its refusal used to say switching criterion would not
  help, and now points at it instead.


# TemporalHazard 1.2.2

## New features

* `hazard()`'s vector interface now evaluates `time`, `status`, `time_lower`,
  `time_upper` and `weights` in `data`'s scope. `hazard(data = df, time = tt)`
  previously failed with `object 'tt' not found`: `data` was consulted only by
  the formula path, and the vector path accepted it and ignored it. The rule is
  `subset()`'s -- a column of `data` wins, and `df$col`, a local vector or a
  literal falls through to the calling frame unchanged -- so with `data = NULL`
  nothing changes and the formula path is untouched.

  Because a column winning can silently redirect a wrapper that forwards its
  own argument by name, `hazard()` now **warns**, once per call, when a symbol
  is both a column of `data` and visible from the calling frame -- that frame
  or a lexical parent of it, up to and including the global environment --
  naming every such symbol and the argument it appeared in. `data` must now be
  a data frame or a list: `hazard(data = <matrix>, ...)` errors, where a matrix
  was previously accepted and silently ignored along with everything else in
  `data`.

* `hzr_translate_sas()` translates a SAS `PROC HAZARD` / `PROC HAZPRED` job
  into a Quarto document of equivalent R calls. It parses the SAS statements,
  builds the calls, and renders them into `.qmd` chunks -- the model state is
  stored as unevaluated calls, so rendering is `deparse()`, not string
  templating.

  **This function is experimental.** A job that translates now renders: the
  emitted `hazard()` chunk binds its fit to a name and asks for an actual
  fit, and the `predict()` chunks have something to predict from. Measured
  on the public `hazard` corpus of 110 `.sas` files, 57 translate into 22
  distinct documents; the 11 of those that synthetic data can drive end to
  end evaluate every chunk and bind a converged fit, and the other 11 --
  `PROC HAZPRED`-only jobs with no local fit to bind -- are exercised up to
  their fit chunks. Read that as a measurement, not as "the translator
  works": the rest of the corpus is refusals or jobs whose external `INHAZ=`
  could not be resolved. It remains a translation aid rather than a turnkey
  reproduction, and the API, the `hzr_sas_job` field layout and the emitted
  document format may all still change.

  Two SAS constructs are refused outright rather than mistranslated into
  something that computes a wrong answer. Each records an `UNTRANSLATED` row
  and emits a `stop()` in place of the fit, so the document fails where the
  fit would have been:

  - a `SELECTION` statement requesting a stepwise screen. `hzr_stepwise()`'s
    refit path needs a formula-interface base fit and this translator emits
    the vector interface, so every candidate refit would error and the screen
    would report zero steps -- indistinguishable from "nothing met
    `slentry`" (#152, #160; the underlying `hzr_stepwise()` silent no-op is
    #159).
  - `LCENSOR` combined with `ICENSOR`. `hazard()`'s single `time_lower`
    argument carries the entry time for status 0/1 rows and the interval's
    lower bound for status 2 rows, so one column cannot express both (#155).

  Two gaps that made the emitted calls compute a different answer from the
  SAS job are closed. The log prediction grid now takes its step from the
  job's own `INC=` expression rather than a hardcoded one (#153): three
  denominators appear across the public corpus (`/49.9`, `/99.9`, `/999.9`)
  and the denominator sets both the step and the number of points SAS's
  `DO lo TO hi BY INC` lands, so reading every job as `/99.9` gave the
  `/999.9` jobs 100 points on a step ten times too large -- every time
  wrong, over an empty `untranslated` and full coverage. The span is the
  loop's own `log(hi) - lo`, which coincides with `5 + log(hi)` only because
  every corpus job starts at -5, and an `INC=` in a form the parser cannot
  read is refused with an `UNTRANSLATED` row rather than stepped by a guess.
  Every log grid in the public corpus also writes its `DO` statement with an
  explicit trailing element (`DO lo TO hi BY INC, hi`), which SAS's `DO` list
  syntax evaluates as the loop *plus* one final point at exactly `hi` -- the
  translator emitted only the loop, so every translated grid stopped short
  of the time the job actually asked for (about 8% short for a `/99.9`
  step). The trailing element is now read from the job's own `DO` statement
  and emitted as the grid's last point; a trailing element this cannot
  resolve to the loop's own bound is refused with an `UNTRANSLATED` row.
  `EVENT`, `ICENSOR` and `WEIGHT` are **counts** in the reference
  implementation, not flags, and `setlik.c` combines one record's
  contribution as `c1c2c3 = c1w + c2 + c3w` with `c1w = C1 * WT` and
  `c3w = C3 * WT`. All three now reach the fit that way. An `ICENSOR` event
  count is no longer discarded (#154); an `EVENT` count carries into
  `weights` and `status` derives from `EVENT > 0`, where `EVENT = 2` used to
  map straight onto `status = 2` and be fitted as **interval-censored** --
  a different likelihood branch, not an under-count (#157); and a `WEIGHT`
  variable no longer weights right-censored rows, because `c2` is the one
  term entering that sum unweighted and `readc2.c` sets it to `1` on exactly
  those rows -- a `WEIGHT` that was `0` there previously deleted them from
  the fit silently (#158). A row where the `EVENT` and `ICENSOR` counts both
  fire is two contributions at once, which one `status` and one `weight`
  cannot express, so the emitted status chunk now stops before the fit
  rather than picking the event branch and discarding the interval one.

  `RCENSOR` is the third of those counts and was being ignored outright
  (#162). It names `C2` -- "COUNT OF CENSORED INDIVIDUALS AT TIME=T" --
  and when a job names it, `readc2.c` reads the column straight from the
  data and skips the `C2 = 1` derivation that a job without `RCENSOR` gets.
  Four censored individuals were therefore fitted as one observation. The
  censored branch of `weights` is now that variable, still unmultiplied by
  `WEIGHT`, and the both-fire guard now covers `EVENT` + `RCENSOR` and
  `ICENSOR` + `RCENSOR` as well: `readobs.c` deletes an all-zero row only
  when `RCENSOR` is named *and* exactly one of the other two is, so a row
  with two counts positive always survives to be summed.

  A `0/1` `RCENSOR` flag that is exactly `1 - EVENT` -- which is what both
  corpus jobs carry, and what the statement is usually used for -- fits
  exactly as before. A `0/1` flag that is **not** its complement does not,
  and both ways it can differ are deliberate: a row with the event and the
  censoring flag both set now stops the document instead of being fitted as
  an event alone, and a row with neither set now carries weight `0` instead
  of a fabricated weight of `1`, matching the row SAS would have deleted.
  A count column that is **negative or missing** on any row now stops the
  document with a message naming the variable, rather than being folded into
  the censored branch or propagating `NA` into `weights` until `hazard()`
  refused it as "non-negative and finite". `readc1.c`, `readc2.c` and
  `readc3.c` apply the same rule to every count the job names -- a missing
  value sets `mdel`, a negative one sets `del` -- and `readobs.c` then skips
  `setobs()` for that row and subtracts it from `Nobs`. Such a row
  contributes nothing at all, so translating it as a right-censored
  observation of weight `1` adds survival mass at a time SAS had removed.
  This is the one case where a missing count is **not** interchangeable with
  a zero one: a zero count is kept and contributes, a missing one is deleted.
  The translator cannot drop rows without changing `n` behind the reader's
  back, so it stops and says to filter them.

  Every translated job now emits a `status` chunk ahead of its fit, where
  these guards live -- previously only jobs with `ICENSOR` or more than one
  named count did. The emitted document format remains experimental.

  Loading a fit from an external `INHAZ=` dataset returns a classed
  `hzr_outhaz` object with a `predict()` method (#151). That method takes
  the same arguments in the same order as `predict.hazard()` --
  `newdata`, `type`, `decompose`, `se.fit`, `level`, `conf.type` -- so a
  positional call means the same thing for both methods of the generic, and
  `conf.type` (the `PROC HAZPRED` parity switch the translator emits) is a
  real argument rather than one that a misspelling could drop into `...`,
  returning the log-log limits the SAS job did not ask for. Its *value* is
  checked only on the survival standard-error path that reads it, exactly as
  `predict.hazard()` does, so an ignored value does not fail a point or
  hazard prediction; a mistyped argument *name* still errors. `type` defaults to
  `"hazard"`, as in `predict.hazard()`, and `decompose = TRUE` is an error:
  an `OUTHAZ=` dataset carries fitted parameters, not a per-phase
  decomposition. Point predictions work; `se.fit = TRUE` is **refused** whenever the SAS fit *estimated* a
  late shape parameter that `PROC HAZARD` put on a composite scale --
  `log(GAMMA*ETA - 2)` and friends, which is the generic unconstrained
  three-phase case rather than an exotic one -- and likewise under `FIXMNU1`
  or where one late parameter is derived from another. A translated `PROC
  HAZPRED` block requests confidence limits unless the SAS job says `NOCL`,
  so such a job stops at its `predict()` chunks with an explicit message
  naming the parameter and its scale, rather than reporting standard errors
  built on the wrong one.

  Treat `job$coverage` as a measure of *parsing* -- tokens recognised --
  not of whether the result runs. The
  parameter translation itself is verified separately: refitting the
  `hz.death.AVC.sas` job's parameters through `hazard()` directly reproduces
  the SAS log-likelihood to the six significant figures the reference
  listing prints (`-210.501`).

  The keyword grammar behind the parser -- 122 keyword rules, 67 of them
  mapped to an R target -- is **generated from the reference
  `HAZARD`/`HAZPRED` C implementation's own lex sources**
  (`data-raw/hazard-grammar.R`), not hand-written. Only the extracted table
  ships; no GPL-2 source enters the tarball. A hand-written table would
  capture only the spellings a study happened to use, and the grammar has
  real context-dependent collisions --
  `M` means a phase shape parameter inside `PARMS` and `MOVE` inside a
  `PHOP`/`STEP` statement -- that a context-free lookup gets silently wrong.

  Constructs the translator does not cover are recorded on the returned
  `hzr_sas_job` object and rendered as visible `UNTRANSLATED` callouts in
  the `.qmd`, never dropped. Two limits are worth stating plainly:

  - **Prediction grids built from `SET`-derived values, function calls, or
    unknown names are not translated.** The parser resolves a `PROC
    HAZPRED` grid's `DO` loop bounds when they are literal numbers or
    DATA-step constants it can fold (e.g. `DO MONTHS = 1*DTY, 2*DTY, ...;`
    with `DTY` assigned earlier in the same DATA step) -- but a bound
    read from `SET`, computed by a function call, or naming something the
    parser can't resolve is refused whole rather than partially read: a
    partially read grid is a partial `newdata`, which is a hollow result.
    Such grids emit an explicit `UNTRANSLATED` block instead, and the
    `predict()` chunks that would have read the grid become a `stop()`
    naming it: emitting `predict(fit, newdata = <name>)` that nothing
    builds either fails on an unbound name or, if the rendering session
    happens to hold an object of that name, reports predictions over
    unrelated times. On the
    public corpus, grid resolution is 19 of 55 (35%), up from 10 of 55
    (18%) before constant folding.
  - **An unresolved `INHAZ=` fails the render, on purpose.** A `PROC
    HAZPRED` job whose fitted-model dataset can't be located -- neither
    from another translated job's `OUTHAZ=` nor from the `librefs`
    argument -- gets an `inhaz-unresolved` chunk, ahead of the grid and
    `predict()` chunks, whose whole body is a `stop()` naming the
    unresolved libref. The document fails to render rather than reporting
    predictions over a model it never loaded.

## Bug fixes

* `$se` on a fitted object is now one standard error per parameter, whatever
  the variance matrix looks like. A multiphase fit legitimately carries NA
  variance rows for the parameters it holds fixed, and a single NA anywhere in
  the matrix collapsed the whole vector to a length-1 `NA`. A five-parameter
  fit came back with a length-1 `$se`, so naming the standard errors against
  the parameters failed with "'names' attribute [5] must be the same length as
  the vector [1]".

  A scalar `NA` now carries the meaning it already has for `vcov()`, that there
  is no variance matrix at all. Where a matrix is present, `$se` is computed
  element by element: a parameter held fixed carries an NA variance row and
  earns an `NA` standard error, while every parameter whose variance *was*
  computed keeps its own. Sizing the vector to the matrix is not enough on its
  own -- filling it with `NA` throughout would leave `$se` conformable and
  empty, and contradicting `summary()`, which reads the same matrix and reports
  those standard errors. `summary()` was never affected either way, so this was
  a quiet inconsistency on the fit object rather than a visible break.

* `hzr_decompos()` no longer returns a wrong value for large `|m|`. The three
  branches with a nonzero `m` all formed `2^m` and the terms built from it
  all three lost the answer well inside the range a fit can reach.

  For `m > 0` the failure is overflow. `2^m` is `Inf` from `m = 1024`, but
  `bt^(-1/nu)` goes first: at `t/t_half = 0.5` with `m * nu = 3` it overflows by
  `m = 750`, and sooner for `nu > 1`. Either way `btnu` becomes `Inf`, and
  `Inf^(-1/m)` is `0`. So `hzr_decompos(0.5, t_half = 1, nu = 3/1000, m = 1000)`
  reported `G = 0` where the answer is `0.3969`, with `g` and `h` `NaN`. Nothing
  warned. `G = 0` is a perfectly ordinary probability, and a fit whose optimizer
  wandered into large `m` used it. The `nu < 0` branch collapsed the same way,
  to `G = 1`.

  For `m < 0` the failure is cancellation instead. `1 - 2^m` rounds to exactly
  `1` once `2^m` falls below machine epsilon, so `(1 - 2^m)^(-nu) - 1` is `0`,
  `rho` is `Inf`, and `G` is again `0`. That collapse is at `m = -53`, which an
  optimizer reaches much more easily than `m = 750`, and the accuracy decays
  before it: at `m = -20` the old code was already wrong in the tenth digit.

  All three branches now work on the log scale. The `m` in the `(2^m - 1)/m`
  factor of `rho` cancels the explicit multiplier, so `m * bt^(-1/nu)` is
  exactly `(t_half/t)^(1/nu) * (2^m - 1)`, and `log(btnu)` follows from
  `hzr_log1pexp()` applied to the log of that product. The `m < 0` branch takes
  `log(1 - 2^m)` from `hzr_log1mexp()` rather than forming the difference. Both
  primitives were already in the package. Checked against a reference computed
  at 100 or more decimal digits, `G` and `g` are now accurate to machine
  precision from `m = -1000` up to `m = 5000`, they track the analytic
  large-`m` limit at `m = 1e6`, and they are unchanged where the old code was
  already right. Small `m` improves too, by eight orders of magnitude or more:
  `log(2^m - 1)` is taken as `x + hzr_log1mexp(x)` for `x = m * log(2)`, which
  holds at both ends, where the direct `log1p(-2^(-m))` decays from about
  `m = 1e-3` down and reaches `-Inf` once `2^(-m)` rounds to `1`.

  One boundary remains, and it is now visible rather than silent. Below about
  `m = -1074` the term `2^m` underflows outright and no rearrangement recovers
  it in double precision; `hzr_decompos()` returns `NA` there.

  The multiphase log-likelihood is evaluable again over the same range. On a
  two-phase fixture it returned `-Inf` from about `m = 450`, for this same
  reason: the smallest observed time is what makes `log(t_half/t)` largest, so
  the overflow arrives earlier than the `m = 750` above. That is what made the
  likelihood surface along the ridge `m * nu = const` hard to characterize.

* `hzr_stepwise()` can no longer return a zero-step result that is silently
  empty. Every accepted move goes through a refit, and a refit that failed was
  downgraded to a warning and then dropped: the returned object carried no
  record of it, so a screen that could not fit a single candidate looked
  exactly like one that tested them all and liked none. `$criteria` now carries
  `refit_failures`, `n_refit_failures` and `stopped_refit_failed`; a run that
  ends on an iteration with failed refits warns that its candidates were never
  tested; and the trace names the cause instead of claiming "no further
  action". A base fit built with the vector interface (`time =` / `status =`)
  stores no formula for the refit to mutate, so every candidate would fail --
  `hzr_stepwise()` now rejects it up front with one message naming the remedy,
  through the same predicate the refit itself uses.

* A multiphase fit is now reproducible. `hazard(dist = "multiphase")` offsets
  the starting values for every optimization start after the first, and those
  offsets were drawn from the ambient RNG stream. The identical call run twice
  returned a different answer: on a 150-row two-phase fit the estimates moved
  by about 0.3 on the log scale and the objective by about 0.08, which is
  enough to change what the fit says. Fitting also advanced the caller's
  stream, so a later `sample()` or `rnorm()` depended on whether a model had
  been fitted first.

  Where the assembled starting values did not converge on their own, the draw
  decided whether there was a fit at all: the fit succeeded only from a
  perturbed start, and about a quarter of draws stopped outright. The same call
  could raise that error on one run and not the next. The reason those starting
  values failed is fixed below, so the draw no longer decides that; it decides
  only which optimum is reached.

  The offsets now come from an internally seeded stream, and the ambient
  `.Random.seed` is restored afterwards. The same data and the same control
  give the same fit, with no `set.seed()` needed, and fitting leaves the
  caller's stream where it found it. The new `control$start_seed` (default 3)
  selects a different ensemble of starts. That is worth reaching for when a fit
  looks like it settled in a local optimum: fit at a few seeds and compare the
  `objective` values. See the note on `starts` below for how to read two
  objectives that differ, which is not always a pair of rival optima.

  `start_seed` takes any whole number within integer range, negatives included
  -- `set.seed(-1)` is perfectly valid and deterministic, so restricting to
  non-negative values would discard half the seed space for no reason. A
  fractional value is rejected rather than truncated: `set.seed()` truncates,
  so `3.9` and `3` would select the same ensemble, and a sweep over
  `3.1 / 3.5 / 3.9` would report three fits having tried one set of starts.
  Coercing quietly would keep that aliasing and merely move it. A value out of
  integer range is rejected too, because `set.seed()` would otherwise fail with
  "supplied seed is not a valid integer" and name neither the argument nor the
  fit it came from.

  `hzr_bootstrap()` draws its own resample before each refit, so replicates are
  still distinct. Its numbers do shift, because the refits no longer advance
  the stream between resamples, and a run with `seed=` is now reproducible end
  to end.

* A multiphase optimization start no longer dies on an infeasible shape. The
  multiphase cumulative hazard short-circuits to an infinite hazard when a
  phase's `m` and `nu` are both negative, so the optimizer sees a penalty and
  backs out of the region. Asked for a per-phase decomposition it
  short-circuited in the wrong shape -- a bare vector where the caller expects
  a named list -- and the Conservation-of-Events adjustment, which runs inside
  the objective on every evaluation, raised `$ operator is invalid for atomic
  vectors`. BFGS steps into that region routinely, so the error came back out
  of `optim()` and the multi-start loop threw the whole start away.

  A discarded start was reported as a failure to converge, so a crash read as a
  numerical problem. On the two-phase fixture in `test-multiphase-gradient.R`
  it cost the fit its assembled starting values outright: `n_starts = 1`
  stopped with an error, the fit survived only on a perturbed start, and 12 of
  50 `start_seed` values failed. All 50 converge now, `n_starts = 1` converges
  on its own, and across 50 seeds every one of the 250 starts is usable.

* A multiphase fit now says which of its starts survived. `fit$fit$starts`
  gives one row per optimization start: its `status`, its `objective`, its
  `convergence` code from `optim()`, whether it was the `best` one and so the
  fit you are looking at, and the `message` of any error it raised. Worth a
  look when a fit is in doubt: on the fixture above the assembled start reaches
  -159.15 and a perturbed start -158.30, and start 1 wins 5 of 50 seeds at the
  default `n_starts = 5` (17 of 50 at `n_starts = 3` -- the rate depends on how
  many starts there are to lose to, so read it against your own setting).

  Read two such numbers as objectives, not as rival optima. That fixture has no
  interior maximum in `m`. Profiled, its objective climbs past -158.30 toward a
  finite limit of -157.88 that is reached only as `m` grows without bound, so
  the better number is a point on a flat ridge where the optimizer met its
  tolerance, and its standard error on `m` is 42.8 against an estimate of 27.2.
  Starts that disagree like that are telling you the shape is barely
  identified, which is the reading `starts` is there to support. On data that
  does identify the shape the picture is the ordinary one: the `avc`
  early+constant profile has an interior maximum near `m = 1` and falls away on
  either side.

  `status` separates the four ways a start can end, and in particular a start
  that stopped at `maxit` reads as `"nonconverged"`, not `"ok"`. That
  distinction is not cosmetic: `optim()` attaches a perfectly finite objective
  to a run it abandoned at the iteration limit, and such a start can carry a
  better objective than one that genuinely converged and so become the
  reported fit. Which start wins is unchanged -- it is still the best
  objective -- but you can now see whether it converged. `fit$fit$converged`
  continues to report that for the fit as a whole.

  A start that errors now also warns rather than being absorbed, and when every
  start fails the error names what was raised instead of calling it a
  convergence failure. An error thrown inside the objective used to be
  indistinguishable from a start that merely optimized badly, which is how the
  defect above stayed hidden.


# TemporalHazard 1.2.1

## Breaking changes

This release contains a breaking change but ships as a minor version. The
`1.x` line is the run-up to a first production release; the major digit is
reserved for that milestone rather than spent on a single changed default.
The change below is also closer to a correction than a redesign — the previous
default deviated from the SAS/C reference this package exists to reproduce.
Read the entry regardless: it can change which variables a stepwise run
selects.

* `hzr_stepwise()` now defaults to `criterion = "score"`, reproducing SAS/C
  HAZARD's `SELECTION` statistic. Previously it defaulted to `"wald"`, which
  refit the model once per candidate and used the refit's Wald chi-square --
  a deviation from the reference implementation this package exists to
  reproduce. **Re-running an existing stepwise analysis can now select a
  different variable set**, because the score and Wald paths take different
  step sequences. Pass `criterion = "wald"` to restore the previous behavior
  exactly.

  The score criterion also removes the per-candidate refit, which dominated
  runtime: a 92-variable two-phase screen fell from roughly 25 minutes per
  bootstrap replicate to seconds.

  Following SAS, the variance used during *selection* is approximate --
  shaping-parameter covariances are ignored. Final-model standard errors are
  unchanged and still use the full Hessian.

  Score is an *entry* criterion. The drop path never refit per candidate, so
  removals are still tested on the current model's Wald p-value against
  `slstay`, as SAS does; drop rows in `$steps` are labelled `"wald"`
  accordingly.

## New features

* The SAS `.lst` parsers now ship with the installed package, under
  `sas-parity/` (`inst/sas-parity/` in the source tree -- `R CMD INSTALL`
  strips the `inst/` prefix). They previously lived in `tests/testthat/`,
  which `R CMD INSTALL` skips unless `--install-tests` is passed -- so a plain
  `install.packages()` or `remotes::install_github()` left them unreachable,
  and a downstream analysis wanting to check its own SAS output against them
  had to clone the repository. Reach them with:

  ```r
  source(system.file("sas-parity", "helper-sas-parity.R",
                     package = "TemporalHazard"))
  ```

  The parsers themselves are unchanged; only their location is. The package's
  own parity tests load them through a shim at
  `tests/testthat/helper-sas-parity.R`, so testthat's helper auto-sourcing
  still applies and no test file changed.

  These functions remain internal (`.hzr_`-prefixed) and unexported. They
  parse a specific vintage of SAS HAZARD listing output and carry no API
  stability guarantee.

* `hzr_bootstrap(verbose = TRUE)` now shows a text progress bar over the
  bootstrap replicates (via `utils::txtProgressBar()`) instead of an
  every-50-replicates message.

* `hzr_bootstrap()` gains a `scope` argument for embedded stepwise variable
  selection during each bootstrap replicate -- the R equivalent of SAS's
  `%HAZBOOT` procedure. **This is experimental**: the selection arguments and
  the shape of what they return may change in a future release, and
  `?hzr_bootstrap` says why under "Selection mode is experimental". The
  fixed-formula bootstrap (`scope = NULL`) is unaffected and unchanged.
  The short version: the design is still being read off production runs, and
  a screen large enough to matter runs for hours while this function writes
  nothing until its last replicate, so splitting a run across processes is
  currently the caller's job. Each replicate runs a fresh `hzr_stepwise()`
  selection (starting from a fixed-shape refit of the base model) instead
  of a plain refit, so `summary$pct` reports the variable's selection
  frequency across resamples and `summary$mean`/`sd`/`ci_*` describe the
  coefficient distribution conditional on selection. `scope = NULL`
  (the default) preserves the original fixed-formula bootstrap unchanged.

* `hzr_read_outhaz()` reads a `PROC HAZARD` `outhaz=` estimate dataset,
  returning the estimates, each parameter's free/fixed status, the
  variance-covariance matrix over the free parameters, and the model-structure
  flags. `outhaz` stores its numbers at full double precision where the
  printed `.lst` carries about seven significant figures, so for any quantity
  it holds it is the better parity reference -- print precision stops being
  the binding constraint and optimizer convergence takes over. The
  log-likelihood is not among them; that still comes from the `.lst`.

## Bug fixes

* **`hzr_stepwise()` never checked that an accepted step improved the fit.**
  A forward step enters a model that *contains* the one it started from, so at
  the optimum the log-likelihood cannot fall. It was written into `$steps` at
  every step and compared at none, so a step whose refit failed to converge
  entered anyway and every later step was then scored against a model that was
  not at its own optimum. In the production screen that surfaced this, three of
  ten steps lowered the log-likelihood and the run still reported convergence,
  ten entries and `p = 0.000` throughout; the final 19-coefficient model fitted
  57 units worse than the nested 16-coefficient model from three steps earlier,
  which cannot happen at a maximum.

  `$steps` now carries `delta_logLik`, a forward step that lowers the objective
  warns and is counted in `$criteria$n_nonmonotone_entries`, and
  `hzr_bootstrap(scope = )` reports `$n_nonmonotone_replicates` — a replicate
  whose path went backwards still contributes its selections to the pooled
  frequencies, and each replicate runs under `suppressWarnings()` so the
  step-level warning cannot reach the user. The comparison carries a small
  tolerance so optimizer noise does not fire it. Reported as issue #134.
* **A score statistic could be finite, enormous and meaningless.** A
  production screen accepted a candidate with `stat` = 92,211 on 1 df and
  `p = 0.000`, after which the refit made the model worse. Neither existing
  guard reached it: the adjusted variance stayed positive and well above the
  collinearity floor, so the statistic was reported as evidence.

  Near-collinearity alone does not do this — as a candidate approaches
  collinearity its score shrinks along with its variance and `Q` stays small.
  `Q` explodes when the model being scored against is *not at its optimum*,
  because the reduced-model score is then no longer zero: the numerator is
  inflated while the denominator stays small. That is the state a failed refit
  leaves behind. `hzr_stepwise()` now declines a candidate whose implied
  coefficient exceeds ±50, reporting `coefficient_diverging`, which is what
  the SAS/C reference has always done (`dqstat.c` rejects `|QBETA| > 50` as
  "the model is going to infinity"). Measured on the bundled `avc` data, a
  model displaced 0.25 from its optimum produced `Q` = 6.5e7 with no reason
  reported at all; a legitimate candidate reaches an implied coefficient of
  about 14, so the threshold has real headroom. Reported as issue #134.

* **The multiphase gradient and Hessian disagreed with the log-likelihood on
  left-truncated data.** For a row with `status` in `{0, 1}` the log-likelihood
  subtracts `H(time_lower)` unconditionally, but the analytic derivatives
  defined the entry time with an extra `time_lower < time` filter. A subject
  entering the risk set at its own event or censoring time was therefore
  differentiated as though it had no entry time, while its weight was still
  applied -- so the derivative was taken of a different function from the one
  being evaluated, and the optimizer left any sensible region immediately.
  Measured on `avc` at fixed parameters, the analytic gradient was out by 382
  where every row entered at its exit time, and by 126 where only *some* did
  -- which is ordinary left-truncated data, not a pathological input. Both
  derivatives now define the entry time exactly as the likelihood does, and
  new tests assert agreement with `numDeriv` across five entry-time layouts
  rather than the one the old filter happened to admit.

* **`hazard()` documented `time_lower` incorrectly, and now warns when it is
  self-defeating.** The argument was described only as the lower bound of a
  censoring interval, "defaulting to `time` if NULL". For `status` in
  `{0, 1}` it is in fact the counting-process **entry time**, and leaving it
  `NULL` means entry at `0`, *not* at `time`. Read literally, the old wording
  said that passing `time_lower = time` changes nothing; it in fact states
  that every subject left the risk set at the instant it entered, which
  removes every such row from the likelihood and leaves the objective
  unbounded above. The documentation now gives both roles, and supplying
  `time_lower >= time` on a `status` 0 or 1 row warns, naming the count and
  the `NULL` default. Reported as issue #136.

* **`hzr_stepwise()` now says *why* a candidate could not be scored, and warns
  when the reason is that the candidate looks strong.** Under
  `criterion = "score"` a candidate whose Q statistic cannot be computed drops
  out of the step, and the run previously reported only a count of them. Two of
  the causes mean opposite things. A collinear column should be dropped. But
  the observed information at `beta = 0` is not positive definite away from a
  maximum, and when a candidate's effect is *large* the log-likelihood curves
  upward there, the adjusted variance goes negative, and the candidate is
  declined -- so the criterion is least able to score exactly the variables a
  screen most wants to find. The old warning attributed both to "a degenerate
  or collinear candidate column", which tells a user to discard their best
  variable.

  `$criteria$uncomputable_reasons` (and `$uncomputable_reasons` on a
  `mode = "select"` bootstrap) now counts the causes by name, `$all_scores`
  carries a `reason` column per candidate, and both warnings name them. The
  reference implementation separates these too, and the R side now matches its
  split: a candidate whose *own* observed information is not positive is
  reported apart from one that is unusable only given what is already in the
  model (`information_nonpositive` against `collinear` and
  `information_indefinite`). The first is reachable on a multiphase fit with a
  large share of interval-censored rows. A run
  that *completed* while declining a candidate for this reason now warns too:
  it previously returned a selection -- sometimes an empty one -- in complete
  silence, which is the case where the omission is least visible. The
  underlying limitation of the score criterion is unchanged and is tracked
  separately; `criterion = "wald"` tests these candidates.

  One behavior change comes with it: the guard on the adjusted variance is now
  a magnitude test rather than a sign test. A variance within rounding distance
  of zero is reported as collinear whichever side of zero it lands on, and only
  a materially negative one is reported as indefinite. The previous floor was
  signed and relative to `I_bb`, so where `I_bb` was itself negative a slightly
  negative variance passed through and produced a negative Q.

* `hzr_bootstrap()` now resamples fits built with the **vector interface**
  (`time =` / `status =` rather than a formula plus `data`). Previously it
  resampled `data` only, but a vector-interface call stores `time = d$col` as
  an *expression*, so every replicate re-evaluated it against the original
  data and returned the original fit. The result was `n_success = n_boot`,
  `n_failed = 0`, no warning, and `n_boot` **identical** replicates -- a
  summary table that looked complete and contained nothing, with `sd` exactly
  0 on every parameter. The evaluated `time`, `status`, `time_lower` and
  `time_upper` vectors are already stored on the fitted object, so they are
  now resampled by the same index as the rows and rewired into each
  replicate's call, exactly as `data` and `weights` already were. Both
  interfaces now produce identical bootstrap replicates for the same model,
  data and seed. Found running a 500-replicate production bagging job that
  completed in 9.5 minutes and produced no usable output.

* **The formula interface mistranslated left- and interval-censored
  `Surv()` objects.** `survival::Surv()` and this package use different
  integer codings for censoring status, and the parser passed `Surv()`'s
  through unchanged. `Surv(time, event, type = "left")` codes a left-censored
  row as `0`, which this package reads as *right*-censored: a wrong answer
  with no error, warning, or other outward sign. Under
  `type = "interval"` / `"interval2"`, `Surv()` codes rows `0`/`1`/`2`/`3`
  for right / event / left / interval against this package's `0`/`1`/`-1`/`2`,
  so left-censored rows were read as interval-censored and interval rows
  carried a status the likelihood does not recognise at all.

  Two related faults in the same branch: `Surv()` stores the status in its
  `time2` column for every non-interval row, and the parser read that
  sentinel as an upper bound; and it set `time_lower` for every row, which
  the likelihood treats as a counting-process *entry* time when status is
  `0` or `1`, cancelling each exact-event and right-censored row out of the
  likelihood. Together these made an interval-censored formula fit return
  the optimizer's failure sentinel rather than a fit.

  Status codes are now translated, an upper bound is taken only from a
  genuine interval row, and `time_lower` left-truncates only interval rows.
  A regression test asserts that a `Surv(type = "interval")` fit reproduces
  the equivalent vector-interface fit to 1e-8 in log-likelihood.
  Found when a production study's three interval-censored records could only
  be expressed through the vector interface.

* A fit that cannot compute a Hessian now says so. The analytic Hessian
  declines for left- and interval-censored rows by design, the optimizer falls
  back to `numDeriv::hessian()`, and `numDeriv` is a `Suggests` -- so on a
  machine installed without Suggests, an interval-censored multiphase fit
  produced no standard errors, `rcond = NA`, `pd = NA` and a `vcov()` of bare
  `logical`, with nothing naming the cause. The user-visible symptom was
  `diag(vcov(fit))` reporting an invalid `'nrow'`, which is unrecognisable
  from the cause. Three paths now warn: `numDeriv` absent (naming it and the
  install command), `numDeriv::hessian()` failing (carrying its message), and
  no Hessian available at all. A `hessian_fn` hook that *errors* is also no
  longer swallowed into silence, so a broken analytic hook is distinguishable
  from one that deliberately declines. Behavior is unchanged -- the
  diagnostics are still `NA` -- but the reason is now stated. Found while
  fitting a production interval-censored study.

* **The score criterion could not test a single candidate on an interval- or
  left-censored multiphase fit.** The analytic multiphase Hessian declines by
  design for `status` in `{-1, 2}`, and the score path had no fallback on that
  branch -- the single-distribution branch has had one all along. The `NULL`
  propagated into the step's reusable nuisance block, every candidate scored
  `NA`, and `hzr_stepwise()` stopped having tested nothing, reporting it in the
  language of a degenerate candidate. Both halves became reachable in this
  release and only together: the `Surv()` translation fix above made left- and
  interval-censored rows expressible through the formula interface, and
  `criterion = "score"` became the default. No test exercised the two at once.

  The observed information is now computed numerically where the analytic form
  declines, as the single-distribution path already did. It agrees with the
  analytic Hessian to 1e-4 on the equivalent right-censored fit, which is what
  licenses using it in place of one. The cost is a numeric Hessian per
  candidate -- the per-candidate work the score criterion exists to avoid --
  but it is paid only where there would otherwise be no information matrix at
  all, and slower is the right trade against selecting nothing. `numDeriv` is a
  `Suggests` here as elsewhere: when it is absent this now stops and names both
  it and `criterion = "wald"`, rather than returning a screen that tested
  nothing.


* **`hzr_stepwise(scope = NULL)` still failed on a formula passed by
  variable.** The fix for that defect reached `.hzr_refit_with_scope()` but
  not three sibling sites, so the default-scope path still raised
  `invalid formula "f": not a call` -- the very string the entry below says
  no longer occurs. All four sites now resolve the stored formula through one
  internal helper, so a fifth cannot drift: `match.call()` records `formula`
  unevaluated, and `deparse(quote(f))` is `"f"`, which `as.formula()` rejects.

  Two consequences of that path becoming reachable, both fixed here.
  `scope = NULL` now skips columns it cannot model instead of erroring on
  them -- numeric and logical columns are kept, since whether a 0/1 field
  arrives logical or numeric depends on the reader that built the frame
  rather than on the variable:
  under an explicit scope the caller named the column, so an error is right,
  but under `scope = NULL` the package enumerates the candidates itself and a
  column it cannot model is its own choice to make better. Any data frame
  carrying a character or factor column -- which is most of them -- was
  otherwise unusable with the default scope.

* `hzr_bootstrap()` no longer returns a silent `n_success = 0` (and
  `n_failed = n_boot`, with no error and no warning) when the model was fitted
  inside a function. `hazard()` stored its call but not the environment that
  call was written in, so each replicate's refit resolved arguments passed by
  symbol -- `theta`, `phases`, `control` -- against the package namespace and
  `globalenv()` rather than the caller's locals. Fits built at the top level
  appeared to work by falling through to `globalenv()`; fits built inside a
  function failed on every replicate, and the per-replicate `tryCatch()`
  swallowed the error. `hazard()` now records the fitting environment, and each
  replicate is evaluated in a child of it that carries the resampled data and
  weights. Affects both `refit` and `select` modes.

* **`hzr_bootstrap(scope = ...)` selected nothing when the base fit's formula
  was passed by symbol.** `hazard()` records its call with `match.call()`, so a
  formula assigned to a variable first (`f <- Surv(t, d) ~ 1; hazard(f, ...)`)
  is stored as a *symbol* rather than a call. The scope-mutating refit
  recovered it with `as.formula(deparse(...))`, which turns that symbol into
  the string `"f"` and errors with `invalid formula "f": not a call`. Every
  post-entry refit therefore failed, no candidate ever entered, and the run
  reported `n_success = n_boot`, `n_failed = 0`, no error and no warning --
  with a summary holding only the base model's parameters. The stored formula
  is now evaluated in the fit's recorded calling environment, which handles
  the literal and by-symbol forms alike, and a stored formula that fails to
  resolve raises an error naming the problem instead of degrading to an empty
  screen. The same defect affected `hzr_stepwise()` directly. (#114)

* **A select-mode `hzr_bootstrap()` run that selects no covariate now warns.**
  The base model's own parameters appear in every replicate by construction,
  so they fill the summary at `pct = 100` and an empty screen reads as a set
  of perfectly reliable variables; nothing in the output prompted the reader
  to compare the parameter names against `names(coef(object))`. The warning
  names the likely causes: an entry criterion stricter than intended, a
  `scope` naming columns absent from the data, or a base fit whose stored call
  cannot be rewritten. Legitimate empty screens warn too -- an entry criterion
  no candidate can clear is also worth reporting. (#115)

* **A stepwise screen that could not score anything now says so, instead of
  looking like one that finished.** Under `criterion = "score"` a candidate
  whose Q statistic cannot be computed -- a degenerate or collinear column,
  or an information matrix that will not invert on this data -- yields `NA`
  and is dropped from consideration. When that happened to every remaining
  candidate the step returned exactly what a legitimate "no candidate met
  `slentry`" stop returns, so a screen that stopped because it was *unable
  to test* its candidates was indistinguishable from one that tested them
  and found nothing. The per-step diagnostic existed on the returned object
  the whole time and had no readers.

  `hzr_stepwise()` now warns when a run stops this way and reports
  `$criteria$n_uncomputable_scores`. Because `hzr_bootstrap()` runs each
  replicate under `suppressWarnings()` -- deliberately, so per-replicate
  numerical noise does not swamp the console -- that warning cannot surface
  in the mode where it matters most, so the count is aggregated instead:
  `hzr_bootstrap()` gains `$n_uncomputable_replicates` and warns once when it
  is non-zero. A replicate that scored nothing still counts toward
  `n_success` while contributing no selections, so it silently depresses
  every reported selection frequency -- which is the whole deliverable of a
  bootstrap screen.

  Found by a pre-release review pass, not by a failing test: the package's
  own `print.hzr_bootstrap` test runs a five-replicate screen in which four
  replicates cannot score a candidate and none selects anything, and it
  passed throughout because it only ever asserted the printed label.

* `hzr_bootstrap()` no longer floods the console with per-replicate numerical
  warnings (e.g. ill-conditioned-Hessian notes from unstable resamples), which
  are not individually actionable when the bootstrap aggregates over replicates.
  Structural problems (a mistyped `scope` column, an invalid scope) still
  surface once, up front.

* `hzr_bootstrap(scope = ..., trace = ...)` no longer errors with "formal
  argument matched by multiple actual arguments". Select-mode forwarded
  `...` to `hzr_stepwise()` alongside an explicit `trace = FALSE`, so any
  caller-supplied `trace=` collided with it.

* Multiphase models with a `"cdf"`/`"hazard"` phase whose shape sits exactly
  at the `m = 0` (Case 3L) or `nu = 0` (Case 2L) limiting-case boundary no
  longer lose their analytic Hessian. The finite-difference second
  derivative used to probe the *other* shape parameter's `-h` side, which
  can cross into the mathematically undefined `m < 0 && nu < 0` region and
  raise an error; this silently fell back to a numerical Hessian (or, if
  that also failed to invert, to `NA` standard errors) for every affected
  fit, not just `hzr_bootstrap()`'s Conservation-of-Events full-information
  recompute. The boundary direction now uses a one-sided finite difference
  instead.

* Multiphase fits with a single free parameter (a two-phase model with all
  shapes fixed, where Conservation of Events fixes one of the two `log_mu`)
  now use the analytic Hessian for standard errors instead of silently
  falling back to a numerical one. Restricting the Hessian to the lone free
  parameter dropped it from a 1x1 matrix to a scalar, which was rejected as
  non-conformant; it is now kept as a matrix (`drop = FALSE`).

# TemporalHazard 1.1.0

## New features

* `predict.hazard(type = "hazard")` now works for **multiphase** models,
  returning the instantaneous additive hazard
  `h(t|x) = sum_j mu_j(x) phi_j'(t)` (previously only single-distribution
  models supported `"hazard"`, via `exp(eta)`). Like `"survival"` /
  `"cumulative_hazard"` it is time-based (requires `newdata$time`), supports
  covariate `newdata`, and `se.fit = TRUE` (delta-method limits on the log
  scale via a numeric Jacobian of the hazard evaluator). `decompose = TRUE` is
  not supported for `"hazard"`. This gives the multiphase instantaneous hazard a
  public route (it was previously reachable only through internal functions).

* `predict.hazard(..., se.fit = TRUE, conf.type = "logit")` selects the survival
  confidence-limit transform. The default `"log-log"` builds limits on
  `log(-log S)` (the `survival::survfit` standard); `"logit"` builds them on
  `logit(1 - S)`, reproducing SAS HAZARD's `HAZPRED` survival limits. With the
  full-information vcov for CoE fits, `conf.type = "logit"` matches the SAS
  `hp.death.AVC` survival CLs to ~1e-5. Hazard / cumulative-hazard limits are
  unaffected (their log scale already matches HAZPRED).

* `predict.hazard(type = "cumulative_hazard", decompose = TRUE, se.fit = TRUE)`
  now returns per-phase **and** total delta-method confidence limits for
  multiphase models, as a long data frame
  (`time`, `component`, `fit`, `se.fit`, `lower`, `upper`). Each phase's CL uses
  only that phase's parameters, so per-phase limits do not sum to the total.
  Previously this combination raised an error.

## Changes

* **`hzr_deciles()` now matches the SAS `deciles.hazard` macro exactly.**
  Previously it excluded subjects censored before the horizon and defined the
  expected count as `sum(1 - S(horizon))`. It now follows the SAS method: **all**
  subjects are ranked into equal-sized risk groups by predicted survival at the
  horizon, and the expected count per group is the **sum of predicted cumulative
  hazard at each subject's own follow-up time** (so group totals sum to the total
  observed events under conservation of events). The `time` argument now only
  stratifies subjects into risk groups; it no longer restricts or excludes any
  subject, and the expected/observed totals are horizon-independent. Verified to
  reproduce the `hm.death.AVC.deciles` SAS decile table (CASES/EXPECTED/ACTUAL)
  to print precision. The output columns are unchanged; their definitions are
  updated in `?hzr_deciles`.

## Bug fixes

* **Conservation-of-Events fits now report the full-information variance.**
  CoE removes one phase's `log_mu` from the optimizer *search* (its score
  equation is the CoE constraint), but the previous code also dropped it from
  the *uncertainty* -- the conserved phase got an `NA` standard error, and
  anything depending on it (other SEs, `se(H)`, prediction confidence limits)
  was understated wherever that phase contributed. At the optimum the CoE
  solution is the unconstrained MLE, so `vcov()` is now recomputed from the
  unconstrained-objective Hessian over the full free set (including the
  conserved `log_mu`), matching an all-`mu`-free (`conserve = FALSE`) fit at the
  same point. On `hz.death.AVC` every parameter SE now matches the SAS HAZARD
  reference (e.g. the conserved early `log_mu`: 0.133 vs the previous ~0.059).
  The recomputation uses `numDeriv` (Suggests) and an invertible Hessian; if
  either is unavailable the fit emits a warning and the conserved `log_mu`
  retains an `NA` standard error (as before).

* **Conservation of Events ignored left-truncation (counting-process entry
  times).** For multiphase fits on `Surv(start, stop, event)` data, the CoE
  reparameterization conserved `Sum H(stop)` while the likelihood scores the
  intercepts on the entry-time scale, `Sum E = Sum [H(stop) - H(start)]`. The
  conserved phase therefore absorbed the spurious `Sum H(start)`, biasing its
  intercept and lowering the attained log-likelihood (the `hz.te123.OMC` fit-1
  parity offset, gap-list P1 #6). `.hzr_conserve_events()` and
  `.hzr_select_fixmu_phase()` now subtract the per-phase entry-time cumulative
  hazard, matching the likelihood and C HAZARD `setcoe` under `LCENSOR`/
  `STARTTME`. Plain right-censored fits (no `start` time) are unaffected.

* **`vcov()` was unusable for multiphase fits and returned an unnamed matrix.**
  `vcov.hazard()` collapsed the entire matrix to a scalar `NA` whenever any cell
  was `NA`. Multiphase fits legitimately have `NA` variance rows -- for
  parameters held fixed (e.g. early shapes) and for the
  Conservation-of-Events-conserved phase `log_mu` -- so the finite
  free-parameter block was discarded for almost every multiphase model. The
  method now returns the full matrix with `NA` rows preserved and labels rows
  and columns with the coefficient names (phase-prefixed for multiphase, e.g.
  `early.x` vs `constant.x`), so a covariate shared across phases resolves to
  distinct, name-addressable slots. A scalar `NA` is returned only when no
  covariance matrix is available.

* **Weibull analytic gradient produced `NaN` for right-censored `time = 0` rows.**
  `.hzr_gradient_weibull()` used an unguarded `log(time)` in the shape (`nu`)
  score; a legal right-censored row at `time = 0` made `0 * -Inf = NaN`, which
  poisoned the entire summed shape-gradient component (then silently zeroed by
  the optimizer, harming convergence). `log(time)` is now guarded with
  `log(pmax(time, .Machine$double.xmin))`, matching the analytic Hessian. The
  other families were audited: exponential (no `log(time)` in the score),
  log-normal (rejects `time = 0`), and multiphase (the decomposition clamps
  `time`) are unaffected.

* **Weibull event hazard was inconsistent with its cumulative hazard.**
  `.hzr_logl_weibull()` defined the event hazard as `mu*nu*t^(nu-1)*exp(eta)`
  while the cumulative hazard was `(mu*t)^nu*exp(eta)`; the former is missing a
  `mu^(nu-1)` factor (the exact derivative is `nu*mu^nu*t^(nu-1)*exp(eta) =
  (nu/t)*H`, Form A as in the C/SAS HAZARD reference). The natural-scale
  log-likelihood and its analytic gradient (`d/dmu`, `d/dnu` event terms) are
  corrected to match. Pure event/right-censored fits were already correct (they
  use the self-consistent internal reparameterization); the visible effect is on
  **mixed event + interval/left-censored Weibull fits**, which delegate to this
  likelihood and previously optimized a slightly mis-specified event term.

* **Weibull gradient attribute ignored observation weights.**
  `.hzr_logl_weibull(..., return_gradient = TRUE)` attached an unweighted
  gradient even when `weights` were supplied (the analytic gradient was off by
  the weight scale, e.g. halved under `weights = 2`). `weights` is now forwarded
  to the score computation. The model-fitting path was unaffected (it uses a
  separate internal weighted gradient); this only changes callers reading the
  `return_gradient = TRUE` attribute on weighted data.

* **`hzr_bootstrap()` was non-functional for weighted fits** (Phase 7c).
  The resample loop rewired only `data` in the refit call, leaving the
  original `weights` argument bound to a symbol in the *caller's* frame.
  The internal `eval()` could not resolve that symbol, so **every** replicate
  of a weighted model errored out (`n_success == 0`) regardless of `fraction`;
  even had it resolved, the un-resampled weights would have been misaligned
  with the bootstrapped rows.  `weights` is now evaluated once and resampled
  in lockstep with the data on each replicate (mirroring how `data` is
  handled).  Unweighted bootstraps are unaffected.  A regression test covers
  both the `fraction < 1` and full-size weighted paths in
  `test-diagnostics.R`.  Follow-up: `hzr_bootstrap()` now resamples the
  weights already stored on the fitted object (`object$data$weights`) rather
  than re-evaluating the call's `weights` expression in `parent.frame()`,
  which fails when the original symbol is no longer in scope (e.g. the fit
  was built inside a helper that has returned).  Caller-frame evaluation
  remains a fallback for objects fitted before weights were stored.  The same
  fragility applied to the call's `data` argument: `hazard()` now stores the
  evaluated `data` argument (the data frame passed to `hazard()`, not a
  `model.frame()` result) on the fitted object (`object$data$frame`), and
  `hzr_bootstrap()` resamples that stored frame instead of re-evaluating
  `cl$data` in `parent.frame()`, so bootstrap succeeds even when the original
  `data` symbol is out of scope.  Caller-frame evaluation remains a fallback
  for objects fitted before the frame was stored.

* **4-phase CoE fixmu-phase selection** (Phase 7d).
  `.hzr_select_fixmu_phase()` used `which.max()` over raw per-phase cumhaz
  at the starting theta.  G3 late phases with typical shape parameters have
  unnormalized cumhaz orders of magnitude larger than other phases, causing
  CoE to pin the G3 `log_mu` away from its true near-zero MLE.  Fixed by
  excluding phases whose cumhaz contribution exceeds 10× the median before
  selecting (falls back to `which.max` when all phases are outliers).  On the
  4-phase CABGKUL fit the CoE vs no-CoE LL gap closes from 6.9 to < 0.1
  units.  Six new tests cover the 4-phase code path in
  `test-conservation-of-events.R`.
* **`time_lower` dual-use bug in Weibull and multiphase likelihoods.**
  When `time_lower` was supplied for a mixed interval-censored + right-censored
  dataset, the Weibull LL interpreted `time_lower` as the counting-process
  *entry time* for right-censored rows, computing H(stop) − H(start) = 0 and
  silently zeroing those rows' likelihood contribution.  Fixed in
  `likelihood-weibull.R` (4 sites: LL, gradient, L-BFGS-B internal LL/gradient)
  and `likelihood-multiphase.R`: `start_vec` is now set from `time_lower` only
  for genuine epoch rows (`status %in% c(0L, 1L)` and `time_lower < time`).
  Two regression tests added to `test-interval-censoring-weibull.R`.

* **`hzr_decompos()` Case 3 corrected and `nu = 0, m >= 0` now fails loud**
  (Phase 7d).  Two issues in the early-phase (G1) sign dispatch:
    - **Case 3 (`m > 0, nu < 0`, "bounded cumulative") carried a spurious
      factor of `m`.**  Its `rho` used a bare `(2^m - 1)^nu` instead of the
      `((2^m - 1)/m)^nu` form used by Case 1, leaving an `m` factor on the
      `bt^(-1/nu)` term.  The CDF diverged from the C HAZARD G1 evaluator
      (`g1flag = 5`) by up to ~0.2 and was discontinuous with its `m -> 0`
      limit (Case 3L).  Adding the `/m` divisor makes the `m` factors cancel,
      reproducing the C evaluator exactly and restoring continuity (verified
      against `src/common/hzd_ln_G1_and_SG1.c`).  No shipped phase uses
      Case 3, so fitted models are unaffected; the synthetic 3-phase golden
      fixture was regenerated because its free-shape optimizer path crosses
      Case 3 territory.
    - **`nu = 0` with `m >= 0`** fell through every dispatch branch, leaving
      the CDF unassigned and raising the cryptic `object 'G' not found`.  The
      `nu -> 0` limit is defined only for `m < 0`; for `m >= 0` it is
      degenerate.  The function now raises a clear, explanatory error.
  New `test-decompos-boundary.R` locks in continuity of all limiting branches
  (Case 1 -> 1L, 2 -> 1L, 2 -> 2L, 3 -> 3L), Case 3 <-> C `g1flag=5` parity,
  `g = dG/dt` internal consistency, CDF sanity, and stability at extreme
  `t_half`.

## Improvements

* **Hardened Hessian inversion for standard errors (Phase 7c).**
  Post-fit variance-covariance estimation now symmetrizes the Hessian,
  checks its reciprocal condition number, inverts via Cholesky with a
  `solve()` fallback for non-positive-definite Hessians, and guards
  non-positive variances instead of silently emitting `NaN` standard
  errors. Ill-conditioned, non-positive-definite, and non-finite Hessians
  now raise specific, named warnings, and fits carry `rcond` / `pd`
  diagnostics that `summary()` surfaces as a note when a fit is flagged.
  This closes the "12+-parameter Hessian stability" hardening item for the
  inversion layer; analytic Hessians (more accurate standard errors) follow
  in subsequent releases.

* **Analytic Hessian for exponential standard errors (Phase 7c, Layer 2).**
  The exponential distribution now computes its post-fit Hessian in closed form
  (`X~' diag(wH) X~` over event + right-censored rows) rather than numerically,
  giving more accurate standard errors. The shared optimizer gained a
  `hessian_fn` hook that analytic Hessians for the remaining families will reuse;
  left/interval-censored exponential fits fall back to the numerical Hessian.
* **Analytic Hessian for Weibull standard errors (Phase 7c, Layer 2).**
  The Weibull distribution now computes its post-fit Hessian in closed form on
  the internal `(alpha, psi, beta)` optimization scale (then mapped to the
  natural scale by the existing delta method) rather than numerically, giving
  more accurate standard errors. Covers event + right-censored data (including
  counting-process start times); left/interval-censored fits fall back to the
  numerical Hessian.
* **Analytic Hessian for log-logistic standard errors (Phase 7c, Layer 2).**
  The log-logistic distribution now computes its post-fit Hessian in closed form
  on the internal `(log alpha, log beta, beta_coef)` scale rather than numerically,
  giving more accurate standard errors. Covers event + right-censored data;
  left/interval-censored fits fall back to the numerical Hessian.

* **Analytic Hessian for log-normal standard errors (Phase 7c, Layer 2).**
  The log-normal distribution now computes its post-fit Hessian in closed form
  on the internal `(mu, log_sigma, beta_coef)` scale rather than numerically,
  giving more accurate standard errors. Covers event + right-censored data;
  left/interval-censored fits fall back to the numerical Hessian.

* **Analytic Hessian for multiphase standard errors (Phase 7c, Layer 2 PR-6).**
  Post-fit standard errors for all multiphase fits now come from a closed-form
  Hessian of the negative log-likelihood rather than a numerical Richardson
  approximation. The Hessian is assembled from three terms: (A) a
  phase-block-diagonal curvature of Σᵢ wᵢ H(tᵢ), (B) a dense Fisher
  information outer product Σₑ (wᵢ/hᵢ²) ∇h ∇hᵀ capturing cross-phase
  parameter interactions, and (C) a phase-block-diagonal curvature of
  −Σₑ wᵢ log h(tᵢ). μ/β parameters use fully closed-form expressions;
  shape parameters (t_half, ν, m, and G3 parameters) use second-order
  central differences. The Conservation-of-Events full-information vcov
  path also switches to the analytic Hessian.
  Left/interval-censored fits fall back to the numerical Hessian.
  Completes the 6-PR analytic-Hessian rollout across all five families.

## Documentation

* `vignette("fitting-hazard-models")` gains an **Interval and left censoring**
  section covering: status coding reference (`-1`/`0`/`1`/`2`), a cardiac
  clinic-visit simulation with right- and interval-censored observations,
  the direct `time_lower`/`time_upper` API, and a comparison showing the
  interval-censored fit recovering `nu` close to 1.0 (true value) while the
  naive exact-at-upper fit incurs a shape bias of ~+0.45.  Includes a callout note on the correct
  use of `time_lower = 0` for right-censored rows.
* `vignette("fitting-hazard-models")` gains a **Convergence troubleshooting**
  section covering: reading the KM cumulative hazard for Weibull starting
  values (log-log plot), when to fix shape parameters vs. estimate freely,
  diagnosing overparameterization via near-zero phase scales and `NA` from
  `vcov()`, and `control` options (`n_starts`, `maxit`).
* Added a package-level overview help page (`?TemporalHazard`) giving the
  additive multiphase model, the phase-type vocabulary, the SAS/C HAZARD
  bridge, and a map of the main entry points.
* Expanded the mathematical content of the core help files in the style of
  `randomForestSRC`: explicit display equations for the generalized temporal
  decomposition `G(t)` (`?hzr_decompos`), the additive cumulative-hazard model
  on `?hzr_phase` and `?hazard`, and defining formulas plus the
  Mächler (2012) reference for the numerical primitives (`?hzr_log1pexp`,
  `?hzr_log1mexp`, `?hzr_clamp_prob`).
* Added methodological references to the nonparametric diagnostics
  (Kaplan-Meier/Greenwood, Nelson-Aalen, Aalen-Johansen) and filled in missing
  cross-references across the exported help pages.
* Explained the remaining enumerated options in the style of the `hzr_phase()`
  phase-type help. `?hazard` gains a **Baseline distributions** section
  describing each `dist` value (`"weibull"`, `"exponential"`, `"loglogistic"`,
  `"lognormal"`, `"multiphase"`) by its hazard shape and when to use it;
  `?hzr_stepwise` gains a **Selection direction and criterion** section
  explaining each `direction` (`"forward"`/`"backward"`/`"both"`) and
  `criterion` (`"wald"`/`"aic"`), including how Wald selection differs from
  C/SAS HAZARD's score-statistic path.

## Testing

* **Patient-specific HAZPRED prediction parity** (Group A fixtures
  `hp.death.AVC.hm1` / `hm2`).  New `test-sas-parity.R` blocks predict survival
  and instantaneous hazard -- with logit survival CLs and log hazard CLs at the
  SAS 1-SD level -- from the saved multivariable both-phase model
  (`hm.death.AVC` final fit, "HMDEATH") for two covariate profiles each
  (hm1: with/without an associated cardiac anomaly; hm2: complete vs partial
  canal by date of repair), matching SAS to ~5e-4 (survival) / ~8e-3 (hazard;
  the looser hazard tolerance reflects the near-singular 9-coefficient fit and
  the steep early-phase times).  Adds a header-driven
  `.hzr_parse_sas_nomogram_mv()` (parses the BY-group "digital nomogram" whose
  rows each carry their own covariate vector) and a shared
  `.hzr_fit_avc_hmdeath()` helper.

* **Stratified HAZPRED calibration parity** (Group A fixture
  `hs.death.AVC.hm1`).  New `test-sas-parity.R` blocks reproduce the
  population-averaged, stratified-by-`COM_IV` outputs from the same HMDEATH
  model: (1) the observed-vs-expected "predict number of deaths" table --
  per stratum, EXPECTED = sum of predicted cumulative hazard at each subject's
  own follow-up, PEXPECT = sum of predicted death probability, ACTUAL =
  observed deaths (totals conserve events, 14.76 + 55.24 = 70), to ~5e-3; and
  (2) the per-stratum mean survival curve (MSURVIV) at the digital time grid, to
  ~5e-4.  Adds `.hzr_parse_sas_calibration()` and
  `.hzr_parse_sas_strata_survival()`.

* **`hm.death.AVC` stepwise documented as a non-parity gap** (Group A).  The
  phase-aware forward `SELECTION SLE=0.2 SLS=0.1` fit's *final* selected model
  is the saved "HMDEATH" fit already verified by the `hm.death.AVC.deciles` /
  `hp.death.AVC.hm1` / `hm2` parity tests; its *selection path* cannot be
  reproduced (SAS uses approximate variances during selection while R's full
  Hessian is near-singular here; SAS's `/I` `/S` flags are phase-level but R's
  `force_in` is phase-blind; R oscillates at p ~ slstay and lands in a worse
  basin -- the same divergence already documented for `hm.deadp.VALVES`).
  `test-sas-parity.R` gains a regression-guard test that exercises the
  multiphase phase-aware stepwise path end-to-end on real data without
  asserting path parity; see `inst/dev/FIXTURE-GAP-LIST.md`.

* **`bs.death.AVC` bootstrap documented as a non-parity gap** (Group A).  SAS
  `%HAZBOOT` runs a fresh stepwise selection on each bootstrap resample and
  reports a variable-selection frequency; R's `hzr_bootstrap()` resamples and
  refits a *fixed* model (no embedded-selection mode), and reimplementing the
  SAS procedure would inherit the documented `hm.death.AVC` stepwise
  divergence.  `test-sas-parity.R` adds `.hzr_parse_sas_bootstrap()` and asserts
  the SAS reference selection frequencies in parseable form (so the parity test
  is half-written for a future bootstrap-with-selection capability), plus a
  regression guard that R's fixed-model bootstrap runs on the cohort; see
  `inst/dev/FIXTURE-GAP-LIST.md`.

* **Phase-specific covariate recovery tests** (Phase 7d).  New
  `test-phase-specific-covariates.R` confirms that `hzr_phase(formula = ~ ...)`
  is correct, not just runnable: simulation-based recovery tests verify that a
  covariate entered into one phase recovers its true coefficient, that the same
  covariate carries independent (here opposite-sign) effects across two phases,
  and that a covariate confined to one phase does not leak into another.  This
  is the honest substitute for a SAS parity fixture and guards against the
  "accepts the formal but never applies it" regression that has surfaced
  before with weights and counting-process times.
* Added fractional (non-integer) weight coverage to close the roadmap 7a gap.
  Prior weight tests verified weighting only via integer row duplication, which
  cannot express fractional (e.g. inverse-probability) weights. The new tests
  assert the two properties that define a correct per-row weighted
  log-likelihood: an **additive split** (a row of weight `a + b` equals two
  identical copies of weights `a` and `b`) and **linear scaling**
  (`L(theta; c*w) = c * L(theta; w)`, gradient likewise, MLE invariant), across
  the Weibull, exponential, and multiphase-with-covariates paths.
* Made the single-distribution weighted-fit tests exercise a real fit. They
  previously omitted `theta` start values, so `hazard(fit = TRUE)` took its
  unfitted branch and the assertions compared `NULL`/`NA` vacuously; they now
  supply starts and genuinely compare the weighted MLE to the duplicated-row
  MLE.
* Added interval-censoring coverage under the multiphase model (roadmap 7c).
  The multiphase likelihood's interval-/left-censored branch had a working code
  path but no isolated test. New R-only self-consistency invariants in
  `test-interval-censoring-multiphase.R` verify the interval contribution equals
  `log(S(lower) - S(upper))`, the left-censored term equals
  `log(1 - exp(-H(u)))`, right-censoring stays `-(H(stop) - H(start))`
  (including left truncation), invalid bounds (`lower > upper`) yield `-Inf`,
  integer weights match row duplication on interval rows, and an
  interval-censored multiphase fit converges.
* Added a SAS fractional-weight parity capture scaffold under
  `inst/extdata/weights-fixtures/` (roadmap 7a / FIXTURE-GAP-LIST B5): a
  `PROC HAZARD ... WEIGHT IPW` template, a deterministic non-integer weight
  dataset, a `.lst` parser, and `test-weights-sas-parity.R`. The parity test
  re-fits the SAS specification in R and compares covariate estimates and
  log-likelihood; it skips when the capture fixture is absent (as it is by
  default), so CI and installation are unaffected until a SAS run is dropped
  in. R-side fractional-weight correctness is already proven by the invariants
  above; this is the drop-in external SAS confirmation.

---

# TemporalHazard 1.0.3

## Bug fixes / CRAN compliance

* `hzr_bootstrap()` no longer touches `.GlobalEnv` directly. The 1.0.2
  `oldseed`/`on.exit()`/`assign(".Random.seed", ...)` save-restore wrapper
  added in 1.0.2 violated CRAN policy on writing to `.GlobalEnv` and has
  been removed. When `seed` is supplied the function simply calls
  `set.seed(seed)` (the documented R API for seeded reproducibility); the
  `@param seed` documentation now notes that the caller's RNG state is not
  restored on exit. With `seed = NULL` (the default) the function does
  not call `set.seed()` at entry, so it starts from the caller's current
  RNG state; the bootstrap still consumes random numbers and advances
  that state in the usual way.

# TemporalHazard 1.0.2

## Bug fixes / CRAN compliance

* The golden-fixture generators (`.hzr_create_*_golden_fixture()`,
  previously `R/golden_fixtures.R`) have been moved out of the package to
  `data-raw/golden_fixtures.R`. They are maintainer-only helpers for
  regenerating the bundled `inst/fixtures/*.rds` reference outputs and are
  not part of the installed package, so they are no longer shipped, checked,
  or user-reachable. This resolves the home-filespace concern at its root:
  the earlier fallback resolved to `system.file("fixtures", ...)` — i.e. the
  installed package directory — whenever the package was installed, so the
  1.0.1 "falls back to `tempdir()`" fix did not actually prevent writing to
  the user library. The bundled `.rds` fixtures still ship and the parity
  tests still read them via `system.file()`.
* `.hzr_generate_golden_fixture()` (the C-binary reference writer in
  `R/parity-helpers.R`, which shares a file with test-time helpers and so
  was kept in the package) now takes a required `output_dir` argument with
  no default path.
* Removed the remaining hardcoded `seed = 42` literals from the relocated
  generators; recorded fixture metadata reflects the actual `seed` argument
  passed (`NULL` by default, so no seed is set inside the function).
* `hzr_bootstrap()` no longer leaves the caller's random-number stream
  altered when `seed` is supplied: the global `.Random.seed` is saved before
  `set.seed()` and restored via `on.exit()`, matching the fixture generators.
  Bootstrap reproducibility under a given `seed` is unchanged.

# TemporalHazard 1.0.1

## Bug fixes / CRAN compliance

* Added `\value` documentation to all exported functions that were missing it:
  `hazard()`, `coef.hazard()`, `vcov.hazard()`, `print.hzr_calibrate()`,
  `print.hzr_deciles()`, `print.hzr_gof()`, and `print.hzr_kaplan()`.
* Internal fixture generators (`R/golden_fixtures.R`) no longer set a specific
  seed unconditionally. Generators now accept an optional `seed` argument;
  when provided, the global RNG state is saved and restored via `on.exit()`.
* Default `output_dir` for fixture generators falls back to `tempdir()` instead
  of the package source directory, keeping the home filespace unmodified.

# TemporalHazard 0.9.8

## New features

* **Delta-method confidence limits on `predict.hazard()`** — Phase 4g of
  the development plan lands. Two new arguments: `se.fit = FALSE` and
  `level = 0.95`. When `se.fit = TRUE`, the return value becomes a
  data frame with columns `fit`, `se.fit`, `lower`, `upper`.
  - **Weibull and multiphase use closed-form Jacobians**
    (`dH/dtheta`, `dexp(eta)/dtheta`, `deta/dtheta`); exponential /
    log-logistic / log-normal fall back to `numDeriv::jacobian` on a
    per-call cumhaz closure.
  - **Transforms match SAS HAZARD** (`hzp_calc_haz_CL.c` /
    `hzp_calc_srv_CL.c`): `hazard` and `cumulative_hazard` use
    log-scale CLs; `survival` uses log(-log S) CLs (equivalent to
    log-cumhaz) so 0 <= lower <= upper <= 1; `linear_predictor` is
    symmetric on the natural scale.
  - **Fixed-shape / CoE multiphase fits produce meaningful CLs** — the
    delta-method sandwich is restricted to the free-parameter submatrix
    of `vcov`, treating fixed parameters as known-with-zero-variance.
  - Backward compatible: `se.fit = FALSE` (default) preserves the
    pre-0.9.8 scalar-vector / decompose-data-frame return shape.

# TemporalHazard 0.9.7

## New features

* **Counting-process / repeating-events likelihood wired up** — Phase 4f
  of the development plan lands. `Surv(start, stop, event)` with any
  `start > 0` is now accepted. The Weibull and multiphase log-likelihoods
  apply `H(stop) - H(start)` to event and right-censored terms; the
  trivial `start = 0` case degenerates to `H(stop)` and recovers the
  plain-Surv fit exactly. Splitting each row into contiguous epochs
  preserves both the log-likelihood and the MLE to optimizer tolerance
  (split-invariance).
* **Weibull + multiphase analytic gradients handle H(start).** The
  closed-form Weibull score adds a `-d H(start)/d theta` term per row
  (guarded at `start = 0`). The multiphase analytic gradient computes
  per-phase `Phi_j(start)` and its shape derivatives, then adds
  `+w_H_start * mu_j * dPhi_j(start)` to each parameter's score; G3
  phase derivatives at `start` use the same finite-difference machinery
  as at `stop`.
* **0.9.5 narrowing removed.** The `hazard()` guard that rejected
  counting-process `Surv(start, stop, event)` with any `start > 0` is
  gone.

# TemporalHazard 0.9.6

## New features

* **`weights` now supported for all distributions** — Phase 4e of the
  development plan lands. The exponential, log-logistic, and
  log-normal likelihoods and their analytic gradients now apply row
  weights to every censoring term (event, right-censored,
  left-censored, interval-censored). The 0.9.5 guard in `hazard()`
  that rejected `weights` for `dist %in% c("exponential",
  "loglogistic", "lognormal")` has been removed. Fits with integer
  weights reproduce the row-duplicated fit to optimizer tolerance
  across all five distributions.
* **Conservation of Events now honours weights.**
  `.hzr_conserve_events()` and `.hzr_select_fixmu_phase()` take an
  optional `weights` argument; the multiphase optimizer threads it
  through so per-phase cumulative hazards are summed on the same
  scale as the (weighted) observed event count. CoE no longer
  auto-disables when weights are non-uniform — the dimension
  reduction stays on and the MLE matches the full-dim path.

## Bug fixes

* **Multiphase analytic gradient now applies `weights`.**
  `.hzr_gradient_multiphase()` accepted neither `weights` nor its
  downstream equivalents: the per-row score weights `w_H` / `inv_h`
  were set to ±1 and the interval-censored finite-difference
  correction summed an unweighted LL. Weighted multiphase fits
  therefore optimised a weighted objective with an unweighted score;
  BFGS line search still converged near the correct MLE but the
  final gradient norm did not go to zero. All three paths now honour
  row weights, and the optimizer's `gradient_fn` wrapper (including
  the all-zero numeric fallback and the CoE wrapper) forwards
  `weights` consistently. Regression test covers weighted analytic
  vs numerical gradient parity. Surfaced by Copilot review on PR #18.

# TemporalHazard 0.9.5

## New features

* **Stepwise covariate selection** — `hzr_stepwise()` runs forward,
  backward, or two-way stepwise selection on an existing `hazard` fit
  using Wald p-values or AIC deltas as the entry / retention criterion.
  Phase-specific entry is supported for multiphase models: a covariate
  can enter one phase and not another. Defaults match SAS `PROC HAZARD`
  (`SLENTRY = 0.30`, `SLSTAY = 0.20`); AIC mode uses `ΔAIC < 0`
  uniformly. SAS-style `MOVE` oscillation guard freezes variables that
  enter + exit more than `max_move` times. Returns an object of class
  `c("hzr_stepwise", "hazard")` with a `$steps` selection trace, scope
  record, and elapsed timer. Implements the core algorithm from C
  HAZARD `stepw.c` / `backw.c`.

## Bug fixes

* **Multiphase convergence after weights/repeating-events merge** —
  restored multiphase optimization that regressed in 0.9.4: three
  interacting defects in the new `weights` threading (dup-arg
  collision in the multiphase / Weibull closures, positional-arg
  corruption in every distribution's gradient call) made every
  optimizer iteration error silently inside `tryCatch`. Diagnosed and
  fixed via commit 73b4657.
* **Weibull analytic gradient now applies `weights`** — both
  `.hzr_gradient_weibull()` and the `grad_internal` closure inside
  `.hzr_optim_weibull()` accepted `weights` as a formal but did not
  apply it to the score vector. The optimizer still converged via
  line search on the (weighted) log-likelihood, but the gradient
  direction was wrong and the final gradient norm did not go to
  zero. Both gradient paths now weight the event indicator and
  cumulative hazard building blocks. Fits with integer weights
  reproduce the equivalent row-duplicated fit to optimizer tolerance.

## Scope change

* `weights` is now only accepted for `dist = "weibull"` and
  `dist = "multiphase"`. The 0.9.4 NEWS claimed weights were
  threaded through all distribution-specific likelihoods; in fact the
  exponential, log-logistic, and log-normal single-distribution paths
  accepted the formal but never applied it, so the fit was silently
  unweighted. `hazard()` now raises an explicit error when `weights`
  is supplied with one of those distributions rather than returning
  an unweighted fit. Full support for the remaining single-dist paths
  is tracked in `inst/dev/DEVELOPMENT-PLAN.md` Phase 4e.
* **Conservation of Events is auto-disabled when weights are not
  all 1.** `.hzr_conserve_events()` receives the weighted event count
  as its target but sums per-phase cumulative hazards across rows
  *without* applying weights, so Turner's adjustment comes out on a
  mismatched scale. The multiphase optimizer now detects non-unit
  weights and skips the CoE dimension reduction, falling through to
  the (correctly weighted) full-dimensional path. Fits are still
  correct; they just don't benefit from the one-parameter
  analytical closed-form solve. Weighted CoE wire-up is tracked
  alongside the other weights completion work in
  `inst/dev/DEVELOPMENT-PLAN.md` Phase 4e.
* **Repeating-events / counting-process notation narrowed.**
  `Surv(start, stop, event)` with `start > 0` is no longer accepted
  by `hazard()`. The 0.9.4 NEWS claimed each epoch contributed
  `H(stop) - H(start)` to the likelihood, but downstream likelihoods
  only read `time_lower` for interval-censored rows (`status == 2`);
  counting-process rows (`status` in `{0, 1}`) were silently scored
  with `H(stop)` alone, so any fit with nonzero entry times was
  silently wrong. `hazard()` now raises an explicit error. The
  trivial case `Surv(0, t, d)` -- equivalent to `Surv(t, d)` --
  continues to work. Full wire-up of `H(stop) - H(start)` for all
  distribution paths is tracked in
  `inst/dev/DEVELOPMENT-PLAN.md` Phase 4f.

# TemporalHazard 0.9.4

## New features

* **Observation weights** — `weights` argument in `hazard()` applies Fisher
  weighting to the log-likelihood for `dist = "weibull"` and
  `dist = "multiphase"`. Each observation's contribution is multiplied
  by its weight, enabling severity-weighted event analyses. Implements
  the SAS `WEIGHT` statement. _The original 0.9.4 entry claimed
  coverage of all distribution paths; the 0.9.5 patch corrected the
  claim and fixed a gradient wire-up bug in the Weibull path._
* **Repeating events** — `Surv(start, stop, event)` start-stop notation
  is parsed. _The original 0.9.4 entry claimed each epoch contributed
  `H(stop) - H(start)` to the likelihood, but the downstream
  likelihoods never applied the lower bound for counting-process rows;
  the 0.9.5 patch narrowed the feature to the trivial `start = 0`
  case and added an explicit error for nonzero starts._

# TemporalHazard 0.9.3

## New features

* `hzr_deciles()` — Decile-of-risk calibration function comparing observed
  vs. expected event counts across risk groups with chi-square GOF testing.
  Implements the SAS `deciles.hazard.sas` macro workflow.
* `hzr_gof()` — Goodness-of-fit function comparing parametric predictions
  against nonparametric (Kaplan-Meier) estimates with observed vs. expected
  event counting. Implements the SAS `hazplot.sas` macro workflow.
* `hzr_kaplan()` — Kaplan-Meier survival estimator with logit-transformed
  confidence limits that respect the [0, 1] boundary, interval hazard rate,
  density, and restricted mean survival time (life integral). Implements the
  SAS `kaplan.sas` macro output structure.
* `hzr_calibrate()` — Variable calibration function for assessing functional
  form before model entry. Groups a continuous covariate into quantile bins
  and applies logit, Gompertz, or Cox link transforms. Supports
  stratification via the `by` parameter. Implements the SAS `logit.sas` and
  `logitgr.sas` macros.
* `hzr_nelson()` — Wayne Nelson cumulative hazard estimator with lognormal
  confidence limits. Supports weighted events for severity-adjusted repeated
  event analyses. Implements the SAS `nelsonl.sas` macro.
* `hzr_bootstrap()` — Bootstrap resampling for hazard model coefficients with
  bagging support (fractional sampling). Returns per-replicate estimates and
  summary statistics (mean, SD, percentile CI). Implements the SAS
  `bootstrap.hazard.sas` macro workflow.
* `hzr_competing_risks()` — Competing risks cumulative incidence using the
  Aalen-Johansen estimator with Greenwood variance. Handles any number of
  competing event types. Implements the SAS `markov.sas` macro.
* **Conservation of Events (CoE)** — Turner's theorem is now integrated into
  the multiphase optimizer. One phase's log_mu scaling parameter is solved
  analytically at each iteration, reducing the optimization dimension by 1
  and improving numerical stability and convergence. Enabled by default;
  disable with `control = list(conserve = FALSE)`. Implements the core
  algorithm from C HAZARD `setcoe.c` / `consrv.c`.
* New vignette: "Complete Clinical Analysis Walkthrough" — end-to-end
  workflow from Kaplan-Meier baseline through validated multivariable model,
  mirroring the SAS HAZARD analytical sequence.

## Improvements

* Multi-start optimizer now respects user-set RNG seeds for reproducibility
  (removed `set.seed(NULL)` that was actively breaking determinism).
* Vignette metadata normalized to YAML `vignette:` key across all 8 files.
* `fit` parameter documentation corrected to state default is FALSE.
* README now includes key capabilities table and development plan link.

# TemporalHazard 0.9.1

## New features

* G3 late-phase decomposition (`hzr_phase("g3", ...)`) now fully integrated
  into the multiphase optimizer, Hessian, and prediction pipeline.
* `fixed = "shapes"` parameter in `hzr_phase()` allows fixing shape parameters
  during estimation (matching C/SAS HAZARD workflow of estimating only log-mu
  scale parameters).

## Bug fixes

* `summary.hazard()` now correctly reports standard errors when some
  parameters are fixed. Previously, `anyNA(vcov)` rejected the entire
  variance-covariance matrix when fixed parameters had NA entries.
* `print.summary.hazard()` coefficient table now shows the correct label
  for G3 phases (was printing empty parentheses).
* `print.summary.hazard()` phase listing now uses the phase name in
  CDF labels (e.g., "cdf (late risk)") instead of hardcoded "early risk".
* SAS missing value markers (`.`) in CSV datasets are now handled via
  `na.strings = c("NA", ".")` in `data-raw/make_data.R`, preventing
  numeric columns from being read as character.

## Documentation

* Seven Quarto vignettes: getting-started, fitting-hazard-models,
  prediction-visualization, inference-diagnostics, mathematical-foundations,
  package-architecture, and sas-to-r-migration.
* Roxygen examples now include both single-phase and multiphase models.
* README switched to self-contained CABGKUL examples with G3 late phase.
* Dataset axis labels corrected to "Months" (not "Years").

## Infrastructure

* CI workflows updated to use `roxygen2::load_pkgload` for lazy data
  compatibility.
* Added lintr CI workflow with `.lintr` configuration.
* pkgdown action bumped to `peaceiris/actions-gh-pages@v4`.
* Added `use-public-rspm: true` to all CI workflows.
* Added `lintr` to Suggests.

# TemporalHazard 0.9.0

## New features

* Multiphase engine: N-phase additive cumulative hazard models via
 `dist = "multiphase"` with `hzr_phase()` specification.
* `hzr_decompos()` parametric family implementing the three-parameter
  temporal decomposition of Blackstone, Naftel, and Turner (1986).
* Multi-start optimizer with Hessian-based variance-covariance estimation.
* C binary parity tests against the KUL CABG reference dataset.
* Five clinical reference datasets: `avc`, `cabgkul`, `omc`, `tga`, `valves`.

# TemporalHazard 0.1.0

## New features

* Single-phase engine: Weibull, exponential, log-logistic, and log-normal
  distributions with formula interface.
* `hazard()` API with `predict()`, `summary()`, `coef()`, `vcov()` S3 methods.
* Golden fixture regression testing system.
* Numerically stable helper primitives (`hzr_log1pexp`, `hzr_log1mexp`,
  `hzr_clamp_prob`).

# TemporalHazard 0.0.0.9000

* Initial package scaffold.
* Added numerically stable helper primitives.
* Added baseline unit tests and CI workflow.
