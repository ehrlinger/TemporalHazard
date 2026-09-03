# TemporalHazard

CRAN R package. A native R implementation of the multiphase parametric
hazard model of Blackstone, Naftel and Turner (1986) — the model the
SAS/C `HAZARD` program implements. The public API is
[`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md)
plus the `hzr_*` family (23 exports, 20 S3 methods, as declared in
`NAMESPACE`).

The package exists to **reproduce a reference implementation**. That
shapes almost every rule below: where R and SAS disagree, the
disagreement is a finding to be decided on first principles and
documented, not a bug to be papered over in either direction.

This file is the operational contract for any agent working in this
repo.

## Definition of done

A change is not done until these pass, in this order:

``` bash
Rscript -e 'devtools::document()'
Rscript -e 'lintr::lint_package()'      # must be 0 lints
Rscript -e 'devtools::test()'           # 0 failures
```

Then once per PR, not once per edit:

``` bash
mkdir -p "$TMPDIR/tree"                          # tar -C will not create it
git archive HEAD | tar -x -C "$TMPDIR/tree"      # commit first — see below
R CMD build "$TMPDIR/tree"
R CMD check --as-cran TemporalHazard_X.Y.Z.tar.gz
```

Run the commands. Reading the code is not evidence, and neither is a
subagent’s report.

Before opening the PR, and whenever the change touches an exported
function or the likelihood, run the `r-reviewer` agent
(`.claude/agents/r-reviewer.md`) over the diff. It is advisory: verify
each finding against the code rather than acting on it, and do not let a
clean report stand in for the commands above. `RELEASE.md` step 5 runs
the same agent over the accumulated diff from the previous release tag
to `main` (`git diff v1.2.2...main` at the time of writing) at release
time, because per-PR review never sees the release as a whole.

Four details there are load bearing:

- **`document()` runs first, every time.** `man/` and `NAMESPACE` are
  generated. A stale `NAMESPACE` makes the test run answer a question
  about the previous commit.
- **Lint runs before tests** because it is seconds against about a
  minute for the suite. Cheap failures first.
- **`R CMD check` does not run the whole suite.** It runs with
  `NOT_CRAN` false, so `skip_on_cran()` tests are skipped: 158 skips and
  2473 passes under check, against 6 skips and **3373** passes locally
  (both measured 2026-08-31 at `3f4a506`, with no environment variables
  set; the local figure needs the SAS fixture checkouts present under
  `~/Documents/GitHub/hazard`). A green check is **not** evidence that
  those tests pass; only the local `devtools::test()` line is.
- **Commit before you `git archive`.** It exports the committed tree, so
  an uncommitted fix is silently absent and the check answers a question
  about the wrong code. This has already cost one wrong conclusion.
  Nothing in the output tells you.

Check with the manual — never `--no-manual`. The PDF step is the only
thing that catches raw Unicode in `Rd`.

## The one thing that destroys work

**This package’s signature defect is an output that looks like a result
and is not.** Not a crash, not a failing test: a populated summary, a
full progress bar, a green badge, over nothing. Every instance below
shipped and was caught later, usually by a production run:

| What it produced | What was actually there |
|----|----|
| 500 bootstrap replicates, `n_success = 500`, `n_failed = 0`, no warning | 500 **identical** fits; `sd` exactly 0 on every parameter |
| A parsed life table, badge PASS | **97 of 844 rows** — page one of a paginated listing |
| A fitted object with [`vcov()`](https://rdrr.io/r/stats/vcov.html) present | a bare `logical`; no standard errors, `rcond = NA` |
| `expect_true(all(status %in% c(0,1,2,3)))` passing for years | the full set of codes `Surv()` can emit — the assertion could not fail |
| A stepwise screen that stopped cleanly | it could not *score* any candidate, which is indistinguishable from finding none good enough |

The defence that works is **asserting a computation varied, or covered
what it claims** — not that it returned without error. Concretely, and
these have all earned their place:

- assert row coverage before comparing values;
- assert a free parameter actually varies across resamples;
- warn when a maximum discrepancy is exactly zero, because that means
  nothing was compared;
- compare replicates, never a summary statistic of them;
- prefer an assertion that can fail: `expect_equal(status, c(1, -1))`,
  not `expect_true(all(status %in% <every possible value>))`.

Three shapes of “absent”, in ascending danger: `NULL` (caught by
`is.null`), an empty container (defeats `is.null`), and a **hollow
object** — right shape, empty inside — which defeats both `is.null` and
a length check.

## The automated gates

Eight of these are **required status checks** on `main`: a PR cannot
merge until they pass. The rest run and report but do not block.

| Gate | When | Blocks merge | What it runs |
|----|----|----|----|
| `lint.yaml` → `lint` | PR, push | **yes** | [`lintr::lint_package()`](https://lintr.r-lib.org/reference/lint.html), `LINTR_ERROR_ON_LINT=true` |
| `lint.yaml` → `house-style` | PR, push | **yes** | fails when `.claude/house-style.md` has drifted from its vault sources |
| `lint.yaml` → `docs-current` | PR, push | no | `git diff --exit-code man/ NAMESPACE DESCRIPTION` after `document()` |
| `spelling.yaml` | PR, push | **yes** | `spelling::spell_check_package(use_wordlist = TRUE)` |
| `R-CMD-check.yaml` | PR, push, release | **yes**, all five | ubuntu devel/release/oldrel-1, macOS, Windows |
| `test-coverage.yaml` | PR, push, release | no | coverage upload |
| `pkgdown.yaml` → `build-and-deploy` | PR, push, release | no | docs site |
| `check-manual.yaml` | push to `main`, release | **cannot** | the PDF manual — the only thing that catches raw Unicode in `Rd` |
| `check-release.yaml` | release published | no | `R CMD check --as-cran` |

`check-manual` says *cannot* rather than *no*: it deliberately does not
run on pull requests, because building the manual is slow and a check
that makes every PR wait is one people learn to route around. A check
that never reports on a PR can never be required — making it one would
block every merge permanently. It runs after the merge instead, so a
raw-Unicode `Rd` is caught on `main`, not before it lands.

The rules live in the repository **ruleset** `protect main`, not in the
legacy branch-protection settings — the two are separate systems, and
the `branches/main/protection` API returns 404 here even though `main`
is protected. Alongside the required checks the ruleset blocks deletion
and force-push, requires a pull request, auto-requests Copilot review,
and requires **one approving review**.

That last one is what actually blocks a merge, and it is easy to miss: a
PR with all eight required checks green still sits at
`mergeStateStatus: BLOCKED` and `reviewDecision: REVIEW_REQUIRED` until
someone approves it. Copilot does not satisfy it — its reviews come back
`COMMENTED`, never `APPROVED`.

⚠️ **The maintainer can bypass all of it.** An earlier version of this
paragraph said there were no bypass actors and the rules therefore
applied to the maintainer too. That was wrong on both halves, and it
mattered, because “nobody can bypass this” is the sentence that makes
*branch, PR, stop* feel non-negotiable. Verified 2026-09-03:

``` bash
gh api repos/ehrlinger/TemporalHazard/rulesets/15010037 \
  --jq '{bypass_actors, current_user_can_bypass}'
# {"bypass_actors":[{"actor_id":5,"actor_type":"RepositoryRole","bypass_mode":"always"}],
#  "current_user_can_bypass":"always"}
```

Read `current_user_can_bypass`, not `actor_id`. GitHub does not publish
what its `RepositoryRole` ids mean, and no endpoint resolves them on a
personal repo — `repos/:o/:r/roles` and
`orgs/:o/custom-repository-roles` both 404 here. The documentation says
only that a `RepositoryRole` bypass actor is an admin, a maintain or
write role, or a custom role built on write, so `5` is one of those four
and which one is not checkable from here. Naming it would be a guess
wearing the costume of a fact, which is the house failure mode.
`current_user_can_bypass` answers the question that actually matters and
is evaluated for whoever holds the token.

It is not theoretical: PR \#222 merged with **zero** approving reviews.
So the rule is a matter of practice, not of enforcement — which is the
stronger reason to follow it, not a licence to skip it. **An agent still
never merges and never bypasses**; the maintainer decides when to.

Required checks are *not* strict: a PR is not forced to re-run the
matrix every time `main` moves. That is a deliberate trade against a
roughly 45-minute Windows job, and it means a branch can merge green
having never been tested against the current `main`. When a change
depends on something `main` gained after the branch was cut, merge
`main` in and re-run rather than trusting the badge.

**There are no git hooks and no Claude Code hooks in this repo.**
Nothing runs locally on your behalf. The definition of done above is
entirely manual — run it yourself.

**CI does not skip `skip_on_cran()` tests, whatever a local `--as-cran`
run does.** `r-lib/actions/check-r-package` sets `NOT_CRAN: true`
itself, so every one of them runs on all five platforms. This was
verified the hard way on 2026-09-01: a `skip_on_cran()` test added in PR
\#200 passed locally and *failed* on Linux and Windows.

That matters mostly for how you read a red CI. The previous text here
claimed CI skipped those tests, which would lead you to dismiss a
CI-only test failure as impossible — and the CI-only failures are the
valuable ones, because they are the platform differences a single
machine cannot show you. CI does skip more than the local run (42
against 6 on that commit), but those are the SAS fixture-availability
skips: the local machine has the checkouts under
`~/Documents/GitHub/hazard` and the runners do not.

CI is still not a substitute for the local suite, for the opposite
reason to the one previously given: it exercises *more* of the suite
than a local `--as-cran` check does, and runs it on platforms you do not
have — so a green local check and a green CI answer different questions,
and you want both.

## Before you touch code

Orient on the public API surface and where things live **before**
editing; do not infer structure from a partial file read. In Claude Code
the codemap is the fast path — see `CLAUDE.md`.

The multiphase path is where the bugs are. A fit reaches the optimizer
through
[`hazard()`](https://ehrlinger.github.io/TemporalHazard/reference/hazard.md)
→
[`.hzr_optim_multiphase()`](https://ehrlinger.github.io/TemporalHazard/reference/dot-hzr_optim_multiphase.md)
→ the likelihood in `R/likelihood-multiphase.R`, with shape handling
split across `R/phase-spec.R` and the decomposition helpers, so the
function you found is often not the one that runs.

## Generated files: never hand-edit

| Path                     | Generated by                         |
|--------------------------|--------------------------------------|
| `man/`, `NAMESPACE`      | roxygen2, via `devtools::document()` |
| `.claude/house-style.md` | the `ehrlinger/house-style` composer |

`.claude/house-style.md` carries a `DO NOT EDIT` banner and the
`house-style` CI job **fails the build** when it drifts from its vault
sources. Editing it reddens CI and the next recompose reverts you. Edit
the sources under `~/Documents/ObsidianVault/memory/`, then:

``` bash
Rscript ~/Documents/GitHub/house-style/compose-house-style.R --repo TemporalHazard
```

## Rules for this repo

- **Censoring status is coded `-1` left, `0` right, `1` event, `2`
  interval.**
  [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) uses
  *different* integers for the same meanings — under `type = "interval"`
  it is `0`/`1`/`2`/`3` for right/event/left/interval. The formula path
  translates in
  [`.hzr_parse_formula()`](https://ehrlinger.github.io/TemporalHazard/reference/dot-hzr_parse_formula.md);
  the vector path does not. Never carry one coding into the other.
  Passing them through unchanged is a real bug this package shipped:
  `Surv(type = "left")` read left-censored rows as *right*-censored, a
  wrong answer with no error.
- **Two interfaces, not interchangeable in the details.**
  `hazard(formula, data)` stores an unevaluated call;
  `hazard(time =, status =)` stores evaluated vectors. Anything that
  rewrites or resamples a stored call has to handle both — this
  asymmetry has produced three separate defects.
- **`stats::` prefixing is the house style.** `R/diagnostics.R` alone
  uses it twelve times. An unqualified `stats` generic is an
  `R CMD check` NOTE waiting to happen; it will not fail a test, because
  it resolves fine at runtime.
- **`numDeriv` is a `Suggests`.** The analytic Hessian declines for
  left- and interval-censored rows by design and falls back to it, so on
  an install without Suggests those fits silently produce no standard
  errors. Add it to every environment audit.
- **Conservation of Events auto-disables** when any status is outside
  `{0, 1}` (`coe_supported_data`). That is deliberate — the CoE identity
  counts only exact events.
- **The score criterion is an *entry* criterion.** The drop path never
  refits per candidate under either criterion: removals are tested on
  the current model’s Wald p-value against `slstay`, as SAS does, and
  only the accepted drop is refitted.
- Anything slow gets `skip_on_cran()`. There are 127 calls today.
- No [`browser()`](https://rdrr.io/r/base/browser.html), no bare
  [`print()`](https://rdrr.io/r/base/print.html), no
  [`library()`](https://rdrr.io/r/base/library.html) inside `R/`.
  [`cat()`](https://rdrr.io/r/base/cat.html) belongs only in a `print.*`
  method.

### Where the check’s time actually goes

Overall `R CMD check --as-cran` with the manual: **3m 33s** (213s),
tarball **2.9 MB** (measured 2026-09-02 at `f790c73`, the `v1.2.8`
release artifact). CRAN’s ceiling is about 10 minutes and it rejects on
that even at 0/0/0, which is what got ggRandomForests archived in June
2026.

| Component        | Time          |
|------------------|---------------|
| Tests            | 91s / 97s     |
| Vignette rebuild | 44s / 46s     |
| Everything else  | the remainder |

Re-measure these when you change them, and stamp the commit.

⚠️ **One pair of measurements is not a trend.** Five `--as-cran` runs on
2026-09-02, on the same machine, gave **199s, 206s, 212s, 213s and
245s**: an observed range of 199s to 245s, a **46-second spread**,
around a median of 212s. Two of those runs, `d5212d9` at 245s and
`8f68f2a` at 206s, differ by 39s while their test content differs by six
assertions, so the spread is the machine and not the code. Quote the
range rather than a mean and a band; five runs do not support a
distribution. The practical consequence is that any two runs can differ
by more than a release cycle’s real growth, so read the dated series,
not the last delta:

| Measured   | Commit        | Overall | Tests |
|------------|---------------|---------|-------|
| 2026-08-22 | (not stamped) | 2m 52s  | 53s   |
| 2026-08-23 | `c9f6c28`     | 3m 29s  | 84s   |
| 2026-09-02 | `f790c73`     | 3m 33s  | 91s   |

Across those three the total moved about 41 seconds and the tests about
38, all of it in tests, over two release cycles. Real, slow, and still a
wide margin against CRAN. The 2026-08-23 entry read the 2026-08-22 delta
as a 37-second jump that “nothing flagged”; against a 46-second spread
on unchanged code, a delta that size can be entirely the machine. The
direction is what is worth watching, and it only shows against several
dated points.

Both dominant costs are the ones that grow silently. Watch the budget
when adding a vignette chunk or an unskipped slow test.

## Change discipline

1.  **Think before coding.** Do not assume, ask. If the request is
    ambiguous or a name, path or signature is uncertain, surface the
    confusion instead of running with a guess. One good clarifying
    question beats a confident wrong edit.
2.  **Simplicity first.** Write the minimum code that solves the stated
    problem. No speculative abstractions, no “while I’m here”
    generalising. For this scientific code, prefer the plain readable
    form a future reader can follow over the clever one.
3.  **Surgical changes.** Touch only what the task requires. Do not
    refactor, reformat or re-style adjacent code, and do not reorganise
    imports or rename things that were not asked for. If you spot
    something worth changing nearby, note it separately rather than
    folding it in.
4.  **Define “done” as a passing test.** State what done looks like
    before you start. If no test covers the change, add or propose one
    rather than declaring success from inspection.

A new dependency is a CRAN cost. Ask first.

## Git and versioning

- **Never push to `main` directly.** Branch, commit, push the branch,
  open a PR, then stop. The maintainer merges. `main` is the only
  long-lived branch: there is no `dev`, and every kind of change takes
  the same route (see `RELEASE.md`, “Branch model”).
- **Never roll the MINOR or MAJOR digit.** That is the maintainer’s
  call, made when a feature set is consolidated into a release. Patch
  bumps are fine for incremental work; say so when you make one.
- **Always a plain three-digit version. No `.9000`, no fourth digit,
  ever** — not for dev snapshots, not to satisfy a version grep.
  `RELEASE.md` carries the full convention.
- ⚠️ **Unlike ggRandomForests, there is no test here that greps
  `NEWS.md` for the `DESCRIPTION` version.** Nothing will catch a
  mismatch. Keep `DESCRIPTION` `Version:` and the `NEWS.md` top heading
  in sync **by hand** on every bump.
- **`Closes #123` fires on the normal flow.** GitHub auto-closes only on
  merge into the *default* branch, and routine work now merges to
  `main`, which is it. This inverted on 2026-08-24 when the `dev` branch
  was dropped; before that the keyword silently did nothing and issues
  had to be closed by hand.

## Prose

Documentation prose (vignettes, README, roxygen `@description` and
`@details`, release copy) follows the house style in
`.claude/house-style.md`: a specific voice, reader persona and project
context. Read it before writing user-facing text.

The default persona for this package is **(d), the public CRAN user** —
no institutional context, reading `?fn` or a vignette. The registered
secondary is **(c), the bilingual SAS-to-R migrant**, which is the
persona `sas-to-r-migration.qmd` is written for. Write for one persona
at a time, not a blend.

## Gotchas

- **Build the release tarball from a clean `git archive` export, not the
  working tree.** An empty `inst/doc` fabricates two vignette WARNINGs,
  and in a git worktree `.git` is a *file*, so `R CMD build`’s VCS
  exclusion misses it and it lands in the tarball as a spurious
  hidden-files NOTE. Both look like package defects and are not. Commit
  first.
- **The SAS `.lst` parsers encode what one study’s listings happen to
  print.** Expect roughly one gap per parser per new study: a hardcoded
  time-variable name, a missing column, a colon-terminated count line, a
  paginated table read only to page one. Read the header, do not assert
  it.
- **`.Rbuildignore` must exclude `^\.claude$`, `^CLAUDE\.md$` and
  `^AGENTS\.md$`.** These are developer files and have no business in a
  CRAN tarball. Confirm with
  `tar tzf TemporalHazard_X.Y.Z.tar.gz | grep -iE 'CLAUDE|AGENTS'` after
  building.
- The CRAN incoming aspell check is **not** run by a local `--as-cran`
  unless `aspell` is installed, so local checks under-report it.
  `devtools::check_win_devel()` is the source of truth for that NOTE,
  and `inst/WORDLIST` feeds only `devtools::spell_check()` — it has no
  effect on CRAN’s check.
- `n_starts` explores a neighbourhood, not the space. `n_starts = 5` can
  be bit-identical to `n_starts = 1` from the same start, while a
  genuinely different start finds a worse optimum. The multiphase
  likelihood is multimodal.
