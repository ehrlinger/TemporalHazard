# Design: `objective = "sas"` -- PROC HAZARD's interval-censored likelihood

Status: approved 2026-08-20. Implements the form established by cross-study
parity work (temporal_hazard#142).

## Problem

`PROC HAZARD` does not use the interval-censored likelihood. For an
interval-censored row `(l, u]` carrying ICENSOR weight `d`, it accumulates

```
  d * log[ S(u) * (Lambda(u) - Lambda(l)) / (u - l) ]
= d * [ -Lambda(u) + log( mean hazard over (l, u] ) ]
```

-- the ordinary event-density term with the *instantaneous* hazard replaced by
the *interval-mean* hazard. The three row types are then one family:

| row type | contribution |
|---|---|
| right censored | `log S(u)` |
| exact event | `log S(u) + log h(u)` |
| interval censored | `log S(u) + log h_bar(l,u]` |

This package computes the statistically correct
`log P(interval) = log(S(l) - S(u))`. The two differ by
`log[dLambda/(exp(dLambda) - 1)] - log(u - l)` per row, so R cannot reproduce
either SAS's log-likelihood or -- where intervals are wide -- SAS's parameter
estimates.

### Evidence

Four references, four regimes, each evaluated at SAS's own printed final
estimates. `diff` is against SAS's printed log-likelihood; admissible error is
the 6-significant-figure half-ulp.

| fit | regime | SAS LL | **S(u)dL/w** | `log[P/w]` | `log f(u)` | `log P` |
|---|---|---|---|---|---|---|
| uslife table2023 `hz.icall` | w=1 yr, weighted, 124 all-interval rows | -410414 | **-0.002** | +4146 | +3271 | +4146 |
| esophagectomy `hz.phernia_int_censor_all` | 37 distinct widths, 0.19-12.3 yr | -213.794 | **+0.00046** | +0.377 | -9.77 | +22.2 |
| preserve_root `hz.dead_predim_avail` | w=1 day, l=0, eq 9-11 | -239.194 | **-0.00032** | +0.047 | +0.021 | -29.5 |
| preserve_root `..._3grp` g1 | w=1 day, l=0, eq 15-17 | -91.3489 | **-0.000037** | +0.021 | +0.401 | -11.8 |

4/4 inside print precision; every rival 0/4. Bit-identical under
TemporalHazard 1.2.0 and 1.2.1, so it is a property of SAS's likelihood rather
than of one package build.

## Scope

**In:** the multiphase interval-censored branch, its gradient, an `objective`
argument on `hazard()`, guards, a shipped public fixture, tests.

**Out:** conservation (`conserve` is a constraint, not the objective, and is
unaffected); the analytic Hessian (it already declines interval rows and falls
back to `numDeriv`); the exponential/Weibull/lognormal/loglogistic likelihoods
(none of them is a `PROC HAZARD` target); left-censoring (SAS has no statement
for it).

## Design

### 1. One helper, two branches

The interval contribution is currently written **twice** -- in
`.hzr_logl_multiphase()` and again in the finite-difference closure `logl_iv()`
inside `.hzr_gradient_multiphase()`. Those copies must agree or the optimiser
steps by the gradient of a different objective than it evaluates. Extract:

```r
.hzr_logl_interval <- function(cumhaz_lower, cumhaz_upper, lower, upper,
                               weights, objective = c("likelihood", "sas"))
```

Callers pass **only the interval rows** -- already subset by
`status == 2` -- so the helper never sees the mask and cannot disagree with a
caller about which rows are intervals. It returns the summed contribution.

| objective | per-row contribution | notes |
|---|---|---|
| `"likelihood"` (default) | `w * (-Lambda(l) + log1mexp(dLambda))` | today's behaviour, byte-identical |
| `"sas"` | `w * (-Lambda(u) + log(dLambda) - log(u - l))` | needs no `log1mexp` |

Both call sites delegate. This is the only place the formula appears.

### 2. Threading

`hazard(..., objective = c("likelihood", "sas"))`, resolved with `match.arg()`,
passed to `.hzr_fit_multiphase()` and on to both `.hzr_logl_multiphase()` and
`.hzr_gradient_multiphase()`. Exact-event, right-censored and left-censored
branches are untouched.

`objective` is a top-level argument, not a `control` list element: it changes
the estimand, and burying that among convergence tolerances makes it easy to
miss in review.

### 3. Guards

Two failure kinds, deliberately distinguished:

| condition | response | why |
|---|---|---|
| `u <= l` on an interval row, `objective="sas"` | `stop()`, naming row indices | data defect; the divisor is undefined |
| any `status == -1`, `objective="sas"` | `stop()` | `PROC HAZARD` has no left-censoring statement, so no SAS run corresponds to the result |
| `dLambda <= 0` | return `-Inf` | parameter infeasibility, not a data defect; matches how the existing code signals infeasible regions to the optimiser |

The `stop()`/`-Inf` split matters: `-Inf` is a value the optimiser is expected
to walk away from, whereas a data defect must not be survivable. A parity run
that silently substitutes a different formula for some rows produces a number
that looks like an answer.

### 4. Shipped fixture

`uslife2023` -- 124 rows x 3 columns (`age_l`, `age_u`, `d_all`), ~1 KB xz.
The published NCHS US Life Table 2023 on a synthetic 100,000 radix: aggregate
counts, no patient-level data, no PHI.

It is the ideal anchor. Every interval is exactly one year, so `log(u-l) = 0`
and the width term switches off, isolating the `S(u)*dLambda` core; and every
row is interval-censored, so nothing else dilutes the signal.

### 5. Tests

1. **End-to-end parity.** `objective="sas"` on `uslife2023` at SAS's final
   estimates reproduces `-410414` to print precision. `objective="likelihood"`
   gives `-406268` and must *not* match -- a test that only asserts the new
   path passes would also pass if the switch were ignored.
2. **Analytic identity.** On random parameter draws,
   `sas - likelihood == sum(w*log(dLambda/(exp(dLambda)-1))) - sum(w*log(u-l))`.
   Pins the algebra with no SAS dependency, so it runs on CRAN.
3. **Gradient/objective consistency.** `numDeriv::grad()` of the SAS objective
   matches `.hzr_gradient_multiphase()` under `objective="sas"`. This is the
   check that catches the duplicated-formula hazard the extraction removes.
4. **Guards.** Both `stop()` conditions error; `dLambda <= 0` yields `-Inf`.
5. **Default unchanged.** Existing fixtures produce identical values under the
   default, confirming the extraction is behaviour-preserving.
6. **Iteration-trace evaluation.** The SAS trace for `hz.icall` gives four
   off-optimum parameter vectors spanning 198 LL units (iterations 0/4/8/13,
   SAS values -410612 / -410419 / -410414 / -410414). Assert agreement at all
   four. At the optimum a wrong objective and a different optimum are
   indistinguishable; off-optimum iterates separate them, and this costs four
   extra objective evaluations with no refitting.

### 6. Documentation

`?hazard` gains an `@note`: the SAS objective exists to reproduce legacy
`PROC HAZARD` runs and **must not be used for new analyses**. It is a density,
not a probability -- inconsistent for wide intervals, where the two forms
differ by 22 log-likelihood units on the esophagectomy reference. Cross-ref
`inst/dev/SAS-INTERVAL-OBJECTIVE-DESIGN.md`.

### 7. Version

~~1.2.1 -> **1.2.2**~~ **Superseded 2026-08-24:** `main` already shipped 1.2.2
(vector-interface `data` scoping, `hzr_translate_sas()`) while this design sat
unimplemented. The target is now 1.2.2 -> **1.2.3** (patch), bumped in both
`DESCRIPTION` and `NEWS.md`, which a test greps for the exact `DESCRIPTION`
version. Still not a minor bump. **Bump when the feature lands, not while the
`"sas"` branch is a stub** -- a version claiming the feature that errors on use
is worse than no bump.

## 8. Implementation status (2026-08-24)

Scaffolded on `feat/sas-interval-objective-v2`:

- [x] §1 `.hzr_logl_interval()` extracted; both call sites delegate. Verified
      bit-identical to installed 1.2.2 (log-likelihood *and* gradient, three
      parameter vectors) on mixed-status weighted interval data.
- [x] §2 threading: `objective` on `hazard()`, `.hzr_optim_multiphase()`,
      `.hzr_logl_multiphase()`, `.hzr_gradient_multiphase()`.
- [x] §3 guards, all four verified firing.
- [x] **The `"sas"` per-row contribution itself.** Implemented in log space,
      not as `log(S(u)*dL/(u-l))`: `S(u)` underflows to 0 for `Lambda(u) > ~745`
      and would take the sum to `-Inf`. Agrees with the parity document's
      product form to 0 ULP, reproduces the closed-form identity exactly, and
      collapses correctly when every width is 1.
- [x] §4 `uslife2023` fixture. 124 rows x 3 columns, 1000 bytes xz -- the
      estimate was ~1 KB. Derived in `data-raw/make_data.R` and gated there on
      the `.lst`'s own printed figures. `inst/extdata/uslife2023.csv` is the
      durable source and round-trips bit-identically to the SAS read.
- [x] §5 tests 2-5, in `tests/testthat/test-sas-interval-objective.R`
      (14 blocks, 21 assertions). Built from deterministic sequences rather
      than `set.seed()`, so an R RNG change cannot move an assertion, and
      asserted against algebraic identities rather than stored numbers, so
      they pin the formula rather than one build's output. Two additions
      beyond the design: the underflow case that motivates the log-space
      form, and an explicit "the two objectives disagree" assertion without
      which every `"sas"` test would still pass if the switch were ignored.
- [x] §5 tests 1 and 6. Reproduce SAS's -410414 at the printed natural
      estimates (achieved -410414.0025, against a +-0.5 print precision) and
      at iterations 0/4/8/13 of SAS's own trace, spanning 198 LL units. The
      conservation check `sum(d*Lambda(u))` lands at 100000.030 against SAS's
      reported 100000, validating the whole cumulative hazard rather than one
      point of it.
- [x] §6 `@note` on `?hazard`, plus a `NEWS.md` entry stating that this is a
      different estimator and not for new analyses.
- [x] §7 version bumped 1.2.2 -> 1.2.3 in `DESCRIPTION` and `NEWS.md`.

**The design is now fully implemented.** What remains is release hygiene, not
design: `R CMD check --as-cran` with the manual, and a PR.

### Note on the `objective` name

`hazard()` already returns `fit$objective`, a *number* (the log-likelihood at
the optimum). The new argument is a *choice of estimand* with the same name.
They never collide in code, but a reader of `summary()` output may reasonably
read one as the other. Flagged rather than renamed: the design was approved
with this name, and renaming is the author's call.

### Guard ordering is load-bearing

`u <= l` must be checked **before** `delta_h <= 0`. A row with `u <= l` almost
always also produces `delta_h <= 0`, so the feasibility guard would return
`-Inf` first and silently swallow the data defect -- the optimiser would walk
away from a corrupt row instead of stopping on it. This is not obvious from
the §3 table, which lists the conditions but not their order.

## Risks

- **Misuse.** The switch makes an inconsistent estimator fittable. Mitigated by
  the explicit argument, the `@note`, and the guards -- not eliminated.
- **Gradient drift.** The FD closure is the historical duplication point. The
  extraction plus test 3 close it.
- **Fixture provenance.** `uslife2023` must be reproducible from
  `/studies/general/uslife/table2023/datasets/built.sas7bdat` after
  `IF D_ALL=0 THEN DELETE`; record the derivation in the data documentation so
  it survives the SAS licence lapse.
