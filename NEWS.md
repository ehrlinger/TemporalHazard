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
  -159.15 and a perturbed start -158.30, and start 1 wins about one time in
  ten.

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
  step sequences. Pass `criterion = "wald"` to restore the previous behaviour
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

  One behaviour change comes with it: the guard on the adjusted variance is now
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
  from one that deliberately declines. Behaviour is unchanged -- the
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
