# TemporalHazard

CRAN R package. A native R implementation of the multiphase parametric hazard model of
Blackstone, Naftel and Turner (1986) — the model the SAS/C `HAZARD` program implements. The
public API is `hazard()` plus the `hzr_*` family (20 exports, 18 S3 methods, as declared in
`NAMESPACE`).

The package exists to **reproduce a reference implementation**. That shapes almost every rule
below: where R and SAS disagree, the disagreement is a finding to be decided on first
principles and documented, not a bug to be papered over in either direction.

This file is the operational contract for any agent working in this repo.

## Definition of done

A change is not done until these pass, in this order:

```bash
Rscript -e 'devtools::document()'
Rscript -e 'lintr::lint_package()'      # must be 0 lints
Rscript -e 'devtools::test()'           # 0 failures
```

Then once per PR, not once per edit:

```bash
git archive HEAD | tar -x -C "$TMPDIR/tree"     # commit first — see below
R CMD build "$TMPDIR/tree"
R CMD check --as-cran TemporalHazard_X.Y.Z.tar.gz
```

Run the commands. Reading the code is not evidence, and neither is a subagent's report.

Four details there are load bearing:

- **`document()` runs first, every time.** `man/` and `NAMESPACE` are generated. A stale
  `NAMESPACE` makes the test run answer a question about the previous commit.
- **Lint runs before tests** because it is seconds against about a minute for the suite. Cheap
  failures first.
- **`R CMD check` does not run the whole suite.** It runs with `NOT_CRAN` false, so
  `skip_on_cran()` tests are skipped: 111 skips and 1500 passes under check, against 96 skips
  and **1538** passes locally (both measured 2026-08-17 on `dev`). A green check is **not**
  evidence that those tests pass; only the local `devtools::test()` line is.
- **Commit before you `git archive`.** It exports the committed tree, so an uncommitted fix is
  silently absent and the check answers a question about the wrong code. This has already
  cost one wrong conclusion. Nothing in the output tells you.

Check with the manual — never `--no-manual`. The PDF step is the only thing that catches
raw Unicode in `Rd`.

## The one thing that destroys work

**This package's signature defect is an output that looks like a result and is not.** Not a
crash, not a failing test: a populated summary, a full progress bar, a green badge, over
nothing. Every instance below shipped and was caught later, usually by a production run:

| What it produced | What was actually there |
|---|---|
| 500 bootstrap replicates, `n_success = 500`, `n_failed = 0`, no warning | 500 **identical** fits; `sd` exactly 0 on every parameter |
| A parsed life table, badge PASS | **97 of 844 rows** — page one of a paginated listing |
| A fitted object with `vcov()` present | a bare `logical`; no standard errors, `rcond = NA` |
| `expect_true(all(status %in% c(0,1,2,3)))` passing for years | the full set of codes `Surv()` can emit — the assertion could not fail |
| A stepwise screen that stopped cleanly | it could not *score* any candidate, which is indistinguishable from finding none good enough |

The defence that works is **asserting a computation varied, or covered what it claims** — not
that it returned without error. Concretely, and these have all earned their place:

- assert row coverage before comparing values;
- assert a free parameter actually varies across resamples;
- warn when a maximum discrepancy is exactly zero, because that means nothing was compared;
- compare replicates, never a summary statistic of them;
- prefer an assertion that can fail: `expect_equal(status, c(1, -1))`, not
  `expect_true(all(status %in% <every possible value>))`.

Three shapes of "absent", in ascending danger: `NULL` (caught by `is.null`), an empty
container (defeats `is.null`), and a **hollow object** — right shape, empty inside — which
defeats both `is.null` and a length check.

## The automated gates

| Gate | When | What it runs |
|---|---|---|
| `.github/workflows/lint.yaml` | PR | `lintr::lint_package()`, `LINTR_ERROR_ON_LINT=true` |
| `R-CMD-check.yaml` | PR | 5 platforms: ubuntu devel/release/oldrel-1, macOS, Windows |
| `test-coverage.yaml`, `pkgdown.yaml` | PR | coverage; docs site |

**There are no git hooks and no Claude Code hooks in this repo.** Nothing runs locally on
your behalf. The definition of done above is entirely manual — run it yourself.

CI is not a substitute for the local suite, for the reason in the previous section: CI's
`R CMD check` skips 111 tests that the local run exercises.

## Before you touch code

Orient on the public API surface and where things live **before** editing; do not infer
structure from a partial file read. In Claude Code the codemap is the fast path — see
`CLAUDE.md`.

The multiphase path is where the bugs are. A fit reaches the optimizer through
`hazard()` → `.hzr_optim_multiphase()` → the likelihood in `R/likelihood-multiphase.R`, with
shape handling split across `R/phase-spec.R` and the decomposition helpers, so the function
you found is often not the one that runs.

## Generated files: never hand-edit

| Path | Generated by |
|---|---|
| `man/`, `NAMESPACE` | roxygen2, via `devtools::document()` |
| `.claude/house-style.md` | the `ehrlinger/house-style` composer |

`.claude/house-style.md` carries a `DO NOT EDIT` banner and the `house-style` CI job **fails
the build** when it drifts from its vault sources. Editing it reddens CI and the next
recompose reverts you. Edit the sources under `~/Documents/ObsidianVault/memory/`, then:

```bash
Rscript ~/Documents/GitHub/house-style/compose-house-style.R --repo TemporalHazard
```

Note the registry name is `TemporalHazard` while the directory is `temporal_hazard`.

## Rules for this repo

- **Censoring status is coded `-1` left, `0` right, `1` event, `2` interval.** `survival::Surv()`
  uses *different* integers for the same meanings — under `type = "interval"` it is `0`/`1`/`2`/`3`
  for right/event/left/interval. The formula path translates in `.hzr_parse_formula()`; the
  vector path does not. Never carry one coding into the other. Passing them through unchanged
  is a real bug this package shipped: `Surv(type = "left")` read left-censored rows as
  *right*-censored, a wrong answer with no error.
- **Two interfaces, not interchangeable in the details.** `hazard(formula, data)` stores an
  unevaluated call; `hazard(time =, status =)` stores evaluated vectors. Anything that
  rewrites or resamples a stored call has to handle both — this asymmetry has produced three
  separate defects.
- **`stats::` prefixing is the house style.** `R/diagnostics.R` alone uses it ten times.
  An unqualified `stats` generic is an `R CMD check` NOTE waiting to happen; it will not fail
  a test, because it resolves fine at runtime.
- **`numDeriv` is a `Suggests`.** The analytic Hessian declines for left- and interval-censored
  rows by design and falls back to it, so on an install without Suggests those fits silently
  produce no standard errors. Add it to every environment audit.
- **Conservation of Events auto-disables** when any status is outside `{0, 1}`
  (`coe_supported_data`). That is deliberate — the CoE identity counts only exact events.
- **The score criterion is an *entry* criterion.** The drop path never refits per candidate
  under either criterion: removals are tested on the current model's Wald p-value against
  `slstay`, as SAS does, and only the accepted drop is refitted.
- Anything slow gets `skip_on_cran()`. There are 92 calls today.
- No `browser()`, no bare `print()`, no `library()` inside `R/`. `cat()` belongs only in a
  `print.*` method.

### Where the check's time actually goes

Overall `R CMD check --as-cran` with the manual: **3m 44s**, tarball **2.7 MB** (measured
2026-08-17 on `dev`). CRAN's ceiling is about 10 minutes and it rejects on that even at 0/0/0
— that is what got ggRandomForests archived in June 2026.

| Component | Time |
|---|---|
| Tests | 67s / 77s |
| Vignette rebuild | 50s / 58s |
| Everything else | the remainder |

Both dominant costs are the ones that grow silently. Watch the budget when adding a vignette
chunk or an unskipped slow test.

## Change discipline

1. **Think before coding.** Do not assume, ask. If the request is ambiguous or a name, path or
   signature is uncertain, surface the confusion instead of running with a guess. One good
   clarifying question beats a confident wrong edit.
2. **Simplicity first.** Write the minimum code that solves the stated problem. No speculative
   abstractions, no "while I'm here" generalising. For this scientific code, prefer the plain
   readable form a future reader can follow over the clever one.
3. **Surgical changes.** Touch only what the task requires. Do not refactor, reformat or
   re-style adjacent code, and do not reorganise imports or rename things that were not asked
   for. If you spot something worth changing nearby, note it separately rather than folding it
   in.
4. **Define "done" as a passing test.** State what done looks like before you start. If no test
   covers the change, add or propose one rather than declaring success from inspection.

A new dependency is a CRAN cost. Ask first.

## Git and versioning

- **Never push to `main` or `dev` directly.** Branch, commit, push the branch, open a PR, then
  stop. The maintainer merges.
- **Never roll the MINOR or MAJOR digit.** That is the maintainer's call, made when a feature
  set is consolidated into a release. Patch bumps are fine for incremental work; say so when
  you make one.
- **Always a plain three-digit version. No `.9000`, no fourth digit, ever** — not for dev
  snapshots, not to satisfy a version grep. `RELEASE.md` carries the full convention.
- ⚠️ **Unlike ggRandomForests, there is no test here that greps `NEWS.md` for the `DESCRIPTION`
  version.** Nothing will catch a mismatch. Keep `DESCRIPTION` `Version:` and the `NEWS.md`
  top heading in sync **by hand** on every bump.
- ⚠️ **`Closes #123` does not fire on this repo's normal flow.** GitHub auto-closes only on
  merge into the *default* branch, and routine work merges to `dev`. Close the issue manually
  after the merge, or it stays open silently.

## Prose

Documentation prose (vignettes, README, roxygen `@description` and `@details`, release copy)
follows the house style in `.claude/house-style.md`: a specific voice, reader persona and
project context. Read it before writing user-facing text.

The default persona for this package is **(d), the public CRAN user** — no institutional
context, reading `?fn` or a vignette. The registered secondary is **(c), the bilingual
SAS-to-R migrant**, which is the persona `sas-to-r-migration.qmd` is written for. Write for
one persona at a time, not a blend.

## Gotchas

- **Build the release tarball from a clean `git archive` export, not the working tree.** An
  empty `inst/doc` fabricates two vignette WARNINGs, and in a git worktree `.git` is a *file*,
  so `R CMD build`'s VCS exclusion misses it and it lands in the tarball as a spurious
  hidden-files NOTE. Both look like package defects and are not. Commit first.
- **The SAS `.lst` parsers encode what one study's listings happen to print.** Expect roughly
  one gap per parser per new study: a hardcoded time-variable name, a missing column, a
  colon-terminated count line, a paginated table read only to page one. Read the header, do
  not assert it.
- **`.Rbuildignore` must exclude `^\.claude$`, `^CLAUDE\.md$` and `^AGENTS\.md$`.** These are
  developer files and have no business in a CRAN tarball. Confirm with
  `tar tzf <tarball> | grep -iE 'CLAUDE|AGENTS'` after building.
- The CRAN incoming aspell check is **not** run by a local `--as-cran` unless `aspell` is
  installed, so local checks under-report it. `devtools::check_win_devel()` is the source of
  truth for that NOTE, and `inst/WORDLIST` feeds only `devtools::spell_check()` — it has no
  effect on CRAN's check.
- `n_starts` explores a neighbourhood, not the space. `n_starts = 5` can be bit-identical to
  `n_starts = 1` from the same start, while a genuinely different start finds a worse optimum.
  The multiphase likelihood is multimodal.
