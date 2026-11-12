# A phase parameter running to its limit: restart, then refuse (#683, #729)

Approved by the maintainer 2026-11-10, against `main` `a12c4c32`. Design only; no code yet.

## The defect

A multiphase fit can walk a phase's shape parameter toward a limit of the
Blackstone-Naftel-Turner parameterisation and return the last iterate with
`converged = TRUE`. #729 measured it on `hm.death.patient` from the start
`t_half = 10`, `nu = 0.2`:

| Run | Where it stopped | Log-likelihood |
|---|---|---|
| `n_starts = 1` | early `nu` 0.067 | -2036.66 |
| `n_starts = 5` | early `log_mu` 113.6, `log_t_half` 275.8, `nu` -2.3 | -1817.64 |
| SAS / warm start | interior | -1792.07 |

Both bad fits report convergence, and the relative gradient passes at both. Neither is a
limit the data support: an interior maximum 25 log-likelihood units higher exists.
So refitting in the limiting form would return a cleaner wrong model. The response has to
find the interior maximum or say it could not.

## What SAS/C HAZARD does

Read in `ehrlinger/hazard` at `40bf657`. Paths are under `src/`.

- It optimises `log(t_half)`, `log|nu|` and `log|m|` (`common/hzd_early_p2t.c:35-40`).
  It never compares a transformed value with zero during the fit.
- The limiting form (the G1 case flag) is chosen once at setup, from the values the
  analyst fixed (`model/setg1.c:620-636`, `716-733`). Asking for a limit means fixing the
  parameter: FIXNU with nu = 0 and M < 0 selects the exponential-decay case.
- The only check during the fit is a fatal overflow error when a transformed parameter
  leaves double range (`common/hzd_early_t2p.c:37-63`).
- The only automatic retry is `RETRY_T`: zero the covariate coefficients, keep the shapes,
  run again (`hazard/hazrd3.c:343-361`).
- After the fit, a non-positive Hessian diagonal or a Cholesky perturbation stops the run
  with "Possible singularity" (`vcov/dcovar.c:207-234`, code 9585). Nothing is dropped,
  refixed or refitted.
- Side finding, not ours: `optim/umstop.c:215` reads `i = *consec++;`, which advances the
  pointer, not the count. The "five consecutive maximum steps" stop cannot fire.

So SAS never detects a limit and never switches a model. The analyst does both by hand.

## The limits

| Limit | What the phase becomes | Expressible today? |
|---|---|---|
| `nu -> 0`, `m < 0` | exponential decay (case 2L) | yes: fix `nu` at 0 |
| `nu -> 0`, `m >= 0` | a step at `t_half` | no; `hzr_decompos()` refuses it |
| `t_half -> 0` | `G = 1` at every observed time; the phase vanishes | no |
| `t_half -> Inf` | a power law, `G ~ t^(1 / (m nu))` | only by switching to a g3 phase with alpha fixed at 1 (#415) |
| `m -> +/-Inf` | truncated power law (the gradient fixture) | no; out of scope here |

## Decisions

1. **Restart, then refuse.** Never switch the model. Never return the drifted iterate.
2. **The trigger** is post-fit and fires on the data test or the weak flag.
3. **The restart grid** is data-anchored and deterministic.
4. **It is built as a post-fit wrapper in `hazard()`**, on by default, with one switch.
5. **Refusals reached through stepwise and bootstrap refits** are handled explicitly, in
   this change.

## Components

### Detection: `.hzr_limit_trigger(fit)`

This lives in a new file, `R/limit-restart.R`. It reads only what the fit object carries.
It returns a list of `(phase, parameter, limit)` records, or `NULL`. It fires on any of
these:

- `$weak` names a phase's shape or log-scale parameter. This covers #729's stall at
  `nu = 0.067`.
  - **Exception (amended 2026-11-10):** not for a phase that `$boundary` records as
    vanished (`coe_phase_vanished`). Its shapes are unidentified because its `mu` went to
    0. That is a different limit, already warned by #565.
  - Task 1 measured this on #444's `t_half = 3` fixture. Every grid start vanished the
    same way, at about -410.1.
- `$boundary` holds a `phase_discontinuity` record (#448: `nu -> 0`, a step).
- `$boundary` holds an `unbounded_phase` record (#444: `t_half` below the first observed
  time).
- **New:** `$boundary` holds a `phase_beyond_data` record. It fires when a G1 phase's `G`
  at the last observed time is still below `.hzr_step_resolution`, which means the rise
  begins after the data ends. That is `t_half -> Inf`, #729's `log_t_half` 275.8.
  - It is built beside the existing records in `hazard_api.R` (around lines 95-140),
    and reads every observed time, as #444 does.
  - Like #448, it thresholds no parameter and asks a question of the data.
- **Amended 2026-11-12 (#734): the span check, for the joint limit nothing records.**
  - The `m < 0` joint limit: `nu` and `t_half` run to infinity together at a fixed
    `c = ((1 - 2^m)^-nu - 1) / t_half`, and the phase becomes `K log(1 + ct)^(1/|m|)`,
    a power law over the data. `G(t_max)` can be far above `sqrt(eps)` (5.8e-7 to 0.125
    on the misses), so `phase_beyond_data` does not fire, and a usable Hessian leaves
    the weak check silent. The `weak_na` trigger above caught 10 of the 14 misses.
  - `.hzr_limit_span()` asks the likelihood. For each free `t_half` of a `cdf` or
    `hazard` phase, it moves `log t_half` by `L = log(t_max / t_min)` of the observed
    times, up and down. It refits the phase's own other free parameters, first to the
    phase's cumulative hazard over those times (in 10 steps, so a curved ridge is
    followed), then on the likelihood. It fires when either move changes the
    log-likelihood by less than `0.5 * qchisq(0.5, 1)` (0.227), or cannot be computed.
  - It is gated: it runs only on a fit nothing else triggered on, and only for a
    `t_half` outside the observed times or a phase with `m < 0`. It costs 30 to 2000
    likelihood evaluations, 0.4 to 2.5 s at n = 300.
  - It is **weak-only** (maintainer's decision, 2026-11-12). The margin is thin: the
    largest change on a drifted fit was 0.038, against 0.082 on one clean parity fit
    (the R-generated synthetic three-phase fixture's late phase). Weak-only, a false
    positive costs a grid and a `"kept"` record, never a refusal.
  - Under it, a grid start that does better wins only if its own span check is clean.
    On the seed-21 fixture three starts beat the original at the same limit.
  - The `nu * exp(-L)` move the measurement also tried is not used: it fires on fits
    with `nu` near 0 that the step check already handles, and on clean ones.
  - The KUL reference fit "not at an optimum" (open question 4) is the C binary's
    model, whose shapes the analyst fixed. Freeing the early phase's three shapes
    gains 11.27 (-3740.52 to -3729.26) at an interior point. A fixed `t_half` is never
    checked.
- **Amended 2026-11-11 (#734): `$weak` is `NA` because the Hessian was unusable.**
  - On 10 of the 14 drifted fits measured for #734, the fit converged into the
    `m < 0` joint limit (`nu` and `t_half` to infinity together, a power law). Its
    Hessian was not positive definite, or missing, so `$weak` was `NA`. The trigger read
    only a list, so nothing fired.
  - Such a fit is now a **weak-only** trigger, with source `"weak_na"`. `NA` names no
    parameter, so there is one record per free G1 (`cdf` or `hazard`) phase: `t_half`,
    limit "no usable Hessian (identification not checked)". The vanished-phase exclusion
    and the fixed mask apply as above. With two G1 phases, that is up to 18 refits.
  - `NA` triggers only on a **fitted** fit: converged, a finite objective, and no
    `fitting` entry on `$degraded`. That excludes `fit = FALSE` and imported fits.
  - It adds records only to a fit that nothing else triggers on. A hard trigger's grid,
    refusal text and rescue are unchanged.
  - Under this trigger alone, a grid start whose Hessian is also unusable is not
    clean, so it cannot be the rescue: it is in the same unchecked state. Every other
    trigger judges its starts as before. The first suite run judged all of them
    strictly, and turned 16 rescues into refusals. Those rescues are
    starts with an unusable Hessian of their own, which is a question for #734.
  - It does not trigger when the cause recorded in `$degraded_causes` is that
    numDeriv is not installed. A missing Suggests says nothing about the fit, and the
    estimates must not depend on whether it is installed. `$weak` itself is unchanged.
    The exclusion covers only the fit **without** numDeriv. A fit **with** numDeriv
    can still trigger here and be restarted, so its estimates can still differ from
    the same call without numDeriv.
  - The edge case: under a `weak_na` trigger, a grid start whose analytic Hessian failed
    with numDeriv absent records the numDeriv cause, so it does not trigger and reads as
    clean, although its identification was not checked either.
  - The rest of the misses (G21's `NULL` reading just above the rcond gate, and the
    `nu -> 0` exponential decay with `m < 0`) stay open in #734.

### Restart: `.hzr_limit_restart()`

**The grid.**
- Covariate coefficients are zeroed, as SAS's `RETRY_T` does.
- The triggered phase's `t_half` is set to the 25%, 50% and 75% quantiles of the event
  times, crossed with `nu` in {0.5, 1, 2}.
- `m` and every other parameter keep the user's start.

That is at most 9 refits per triggered phase. There is no RNG.

**Each start runs `hazard()`'s whole fit path.** That includes CoE, vcov, `$boundary`,
`$weak` and `$degraded`, so a rescued fit is an ordinary fit object. The check is switched
off inside a restart, so the grid cannot recurse.

**The winner** is the fit with the best objective **overall**, taken across the original
fit and every grid start.
- If that fit is clean, meaning its own `.hzr_limit_trigger()` is `NULL`, it is returned.
- If the best objective belongs to a fit at a limit, the data support the limit, and the
  fit refuses. It never returns a clean but worse fit.

Amended 2026-11-10 after Task 1 measured the #448 `avc` fixture. Its grid's best was
-21.66, at the step. The best fit that did not trigger was -208.31, 130 log-likelihood
units below the -78.15 the fit returns today. "Best that does not trigger" would have
returned it as a rescue.

### Hard and weak-only triggers (amended 2026-11-10, after Task 5)

The Task 5 suite census found the weak-flag source refusing 59 ordinary,
weakly identified test fits. The grid found nothing better for them, and each
fit paid 9 refits. So the two kinds of trigger now behave differently.

- **Hard trigger: any record from `$boundary`** (the step, the rise past the
  data, or below the first observed time). The rule is restart, then refuse, as
  described below. These are facts about the data.
- **Weak-only trigger: every record is the weak flag.** The grid is only tried.
  - **A grid start beats the original** by more than the relative tolerance: it
    is returned as a rescue.
  - **Otherwise the original fit is returned unchanged**, with its own weak
    warning, and `$limit` records the grid that was tried with
    `outcome = "kept"`. A weak-only trigger never refuses.
- **`control$limit_restart = FALSE`:** a hard trigger refuses at once. A
  weak-only trigger returns the original fit, with `$limit` `NULL`.

#729's two fits illustrate the split. The 1-start fit (nu 0.067) is weak-only,
and it is rescued because the grid beats it. The 5-start fit (log t_half 275.8)
is caught by `phase_beyond_data`, so its trigger is hard.

### Outcome, wired after the `$boundary` assembly (`hazard_api.R` around line 2010)

- **No trigger.** The object is unchanged and `$limit` is `NULL`. The parity fixtures
  must all take this path.
- **Trigger, and an interior restart wins.** The fit returned is the winning restart.
  - `$limit` records the trigger, each start with its objective and trigger, and the
    winner.
  - A warning names the original drift and the objective gained.
- **Trigger, and the best objective overall is at a limit** (the #448 step), or no clean
  restart, or `control$limit_restart = FALSE`. The fit stops
  with `stop()`, class `hzr_limit_refusal`.
  - The fields are `phase`, `parameter`, `limit` and `best_objective`.
  - The message names the phase, the parameter and the limit, and points to
    `hzr_phase(..., fixed = )`, which is SAS's FIXNU and FIXTHALF route.

The model never changes, so degrees of freedom, `logLik()` and `vcov()` keep their
ordinary meaning (#683 question 4). The only new state is `$limit`.

### Reporting

`$limit` is a tri-state, like `$weak` and `$boundary`:

- `NA`: not examined, because the fit did not run.
- `NULL`: examined, and nothing found.
- a list: a restart happened.

`print()` and `summary()` **always** show a line: `Limit check: none`,
`Limit check: restarted (<phase>.<param> -> <limit>)`, or `Limit check: not run`. A line
that appears only on bad news cannot be told from a line its author forgot to write
(`R/degraded.R`).

`$degraded` gets no new capability, because a rescued fit is not degraded.

## Refits that can now refuse

The plan phase censuses these by **executing** each route, not by grep. A static census
misses `do.call()` and generated code (#376).

- **`hzr_stepwise()`.** The warm and default refit pair is in `R/stepwise-refit.R`
  (lines 474-542). The candidate and drop refits are inside `tryCatch()` in
  `R/stepwise-step.R` (lines 273, 508, 638 and 1041).
- **`hzr_bootstrap()`.** The replicate refits are in `R/diagnostics.R`.
- **`hzr_evaluate()`** does not fit. It is listed so the census shows it was checked.

**Rules for each:**

- **Bootstrap.**
  - A refusing replicate is a failed replicate, counted in `n_failed` with reason
    `"limit"`.
  - A restarted replicate counts as a success and keeps its `$limit`.
  - The summary reports how many replicates were restarted.
- **Stepwise entry.**
  - A candidate whose refit refuses is recorded as unscored, reason `"limit"`, and the
    step warns, naming it.
  - "No candidate entered" must not hide "a candidate could not be scored".
- **Stepwise drop.**
  - Removals are tested on the current fit's Wald p-value, so only the accepted drop
    refits.
  - If that refit refuses, the step stops with the refusal. It does not keep a
    half-applied model.

Callers match on the class `hzr_limit_refusal`, never on message text. Mind the known
`tryCatch` traps: sibling handlers catch each other's re-raise, and `stop()` inside a
classed handler is swallowed by an outer `error =` handler.

## Testing and acceptance

**#729.**
- `hm.death.patient` from the blind start `t_half = 10`, `nu = 0.2`, with `n_starts = 1`
  and with `n_starts = 5`.
- Each must reach -1792.07 to a relative tolerance, with `$limit` a list naming the drift.
  Or it must refuse, naming the parameter.
- This replaces the warm start in `tests/testthat/test-sas-parity.R` (lines 205-300).

**One reproducer per trigger, each asserting the response, not just a warning:**
- `t_half -> Inf`: simulated power-law data fitted with a cdf phase. This tests the new
  `phase_beyond_data` record.
- `nu -> 0` as a step: the #448 fixture.
- `t_half -> 0`: the #444 fixture.
- A stall that only the weak flag catches.

**The restart must have varied.**
- Assert that the grid's objectives are not all identical.
- Assert that the winner differs from the drifted fit.

**Planted faults.** Disable each trigger branch in turn. Each time a named test must go
red, and the PR prints which test killed which mutant.

**Refusal.** Use a case where every grid start still triggers, and assert the class and its
fields.

**Callers.** Force a refusal through the classed error, then assert bootstrap's `n_failed`
with reason `"limit"`, and stepwise's unscored candidate with its warning.

**False-positive gate (#683 question 3).**
- Run `.hzr_limit_trigger()` over every parity-fixture fit and every multiphase fit in
  the suite.
- The parity fixtures must read **zero** triggers.
- Every other firing is listed and explained one by one in the PR, with the per-file
  pass-count delta stamped at its sha.

**Cost.**
- Refits run only on a trigger.
- A new test file that takes more than about a second with `NOT_CRAN=false` gets a
  file-level `skip_on_cran()`.
- The `--as-cran` check time is re-measured and stamped.

## Out of scope

- `m -> +/-Inf`, the gradient fixture's ridge. `m` is on the natural scale, so the weak
  reading does not name it. File it as its own issue.
- Refitting in a limiting form, which would be a model switch.
- Detection during the iterations. SAS has none, and it would need a hook into the
  optimiser.
- #589, the convergence criterion. #729 is evidence for it but does not need it.
- The `umstop.c:215` pointer defect. It is reported in `ehrlinger/hazard`, not fixed here.
