# Release process

Maintainer checklist for TemporalHazard. Not shipped (`.Rbuildignore`d).

## Versioning convention (three digits, decide the bump at release)

- **Every version is a straight three-part number**: `1.1.0`, `1.2.0`, `1.2.1`.
  Never a fourth component, and never a `.9000` dev suffix — not for dev
  snapshots, not to satisfy a version-grep test, not anywhere.
- **Acceptance itself never bumps the version.** `dev` stays at the number
  that shipped; there is no placeholder to strip later and no bump to
  remember. The version moves only when work moves it, per the next two
  bullets.
- **The patch digit moves when a fix actually lands**, not before: `1.2.0` →
  `1.2.1` once there is work behind it.
- **The minor and major digits are the maintainer's call, taken at release
  time** from what is actually in NEWS — breaking API → major, back-compatible
  features → minor, fixes only → patch. Accumulate under the current minor and
  let one consolidated minor carry a whole feature set rather than spending a
  new one per batch.
- **The major digit is reserved for a deliberate product milestone**, and this
  overrides the mechanical rule above. A breaking change does **not** by itself
  force a major bump while the package is still working toward its first
  production release — the whole `1.x` line is the run-up to that milestone,
  and spending the major digit on a single default flip mid-run-up leaves
  nothing to mark the milestone with.

  This is not hypothetical. On 2026-08-04 the mechanical reading took `dev` to
  `2.0.0` for one changed default in `hzr_stepwise()` (`criterion = "wald"` →
  `"score"`). That work was folded back into `1.2.0`, where it belonged: the
  old default deviated from the SAS/C reference this package exists to
  reproduce, so the change is closer to a correction than to a redesign, and
  the `1.x` line had not yet reached production.

  When a breaking change ships in a minor for this reason, say so plainly in
  NEWS under a **Breaking changes** heading. The version number is a release
  decision; the warning to users is not negotiable and is not softened by it.

So the per-release cycle is:

```
1.1.0 (released)  ->  dev stays 1.1.0  ->  accumulate work, patch-bumping
                                            as fixes land (1.1.1, 1.1.2 ...)
  ->  at release, decide from NEWS: 1.1.3 | 1.2.0 | 2.0.0  ->  submit
  ->  accepted  ->  dev stays at that number until work moves it ...
```

**Do not pre-stamp a future release number into the version.** Labelling `dev`
with the number you expect to ship silently commits a decision that should be
made last, from the evidence. That is the mistake behind the 2026-06 version
confusion, when `dev` was labelled for a minor release before `1.1.0` had even
shipped and before the work existed to justify it.

## Branch model (default: single line)

| Branch | Version | Role |
|--------|---------|------|
| `main` | the clean CRAN release (e.g. `1.1.0`) | the released/stable line; always == the latest version on CRAN. Tagged at each release. |
| `dev`  | the same three-part number as `main`, patch-bumped as fixes land | the working line; accumulates the next release. Carries no dev marker. |

Both branches carry a clean three-part version. A README/badge tweak on `main`
at the same version number is normal and fine — it does not need a version
marker. (Historically `main` carried a fourth `.9000` component; that was
dropped after `1.1.0`, and the convention was dropped entirely in 2026-08.)

**Routing:**

- **Routine work** (R code, `man/`, vignettes, tests, docs) — branch off `dev`,
  PR to `dev`. It reaches `main` only via a release.
- **Urgent patch / hotfix** to the released line — branch off `main`, fix, set
  `DESCRIPTION` to the next *patch* (`X.Y.(Z+1)`), PR to `main`, release it,
  then merge `main -> dev` so the fix flows into the dev line. (This is the
  "bug fixes go to `main`" path: ship a patch without dragging in `dev`'s
  unfinished features.)
- **Cosmetic / GitHub-facing fixes** (README prose, badges, figures, pkgdown
  config) — branch off `main`, PR to `main`. No version change; `main` stays at
  the released number.
- **Release ops** (tag, GitHub Release) — directly on `main` (the documented
  exception to the branch+PR rule).

**Keep them in sync:** after any `main`-side fix, merge `main -> dev` so it
flows into the dev line. Skipping this is what produces the "same fix, two
parallel PRs" situation.

### When to diverge into a two-branch model

Reserve the heavier git-flow split — `main` keeps shipping `1.1.x` patches
while `dev` becomes a diverged **next-major** line (e.g. `2.0.0`) — for a
**genuine breaking rewrite** running in parallel with patch maintenance (the
ggRandomForests-v4 situation). It carries a real cost: every `main` hotfix must
be back-merged into `dev` or `dev` drifts from the patched line. Do **not** use
it as the default — it is not worth that tax for ordinary incremental cycles.

**Release flow (default single line):**

```
dev (X.Y.Z) accumulates, patch-bumping as fixes land
  -> at release, decide the next number N from NEWS scope
  -> set DESCRIPTION + NEWS top heading to N on dev
  -> run the pre-submission gate (below), submit to CRAN
  -> on acceptance: merge dev -> main (main = N), tag vN, GitHub Release
  -> dev stays at N; it moves again only when work lands (N patch-bumped)
```

## Pre-submission checklist

1. `DESCRIPTION` `Version:` is a three-part `X.Y.Z` (no fourth component).
2. `NEWS.md` top heading is `# TemporalHazard X.Y.Z` and describes the
   user-visible changes since the last CRAN release.
3. `cran-comments.md` version heading matches `DESCRIPTION`; if a
   resubmission, each reviewer point is itemised with how it was addressed.
4. `devtools::document()` is clean and `man/` is in sync.
5. **Adversarial review of the accumulated release diff.** Every PR is reviewed
   on its own; nobody reviews the release as a whole, and by submission time
   `dev` is typically dozens of commits ahead of `main`. Run the `r-reviewer`
   agent (`.claude/agents/r-reviewer.md`) over `git diff main...dev`.

   This step is **advisory**, unlike the rest of this checklist. It returns
   findings to triage, not a pass/fail. Verify each one against the code before
   acting -- and equally, do not let a clean report substitute for step 6. The
   two catch different things: the check finds what static analysis sees in a
   built tarball, the review finds what is wrong but passing, including prose
   that contradicts the code. Run it *before* the check so fixes land first.
6. `R CMD check --as-cran` **with the manual built** → 0 errors, 0 warnings,
   0 notes. Build the real tarball and check it, *not* `--no-manual`:

   ```sh
   R CMD build .
   R CMD check --as-cran TemporalHazard_X.Y.Z.tar.gz
   ```

   The PDF-manual step catches raw-Unicode-in-Rd that `--no-manual` skips. Also
   confirm overall check time is well under the CRAN ~10-min budget and the
   tarball is < 5 MB. The built tarball must contain `build/vignette.rds` — if
   it is missing you get a "no prebuilt vignette index" NOTE (do **not** re-add
   `^build$` to `.Rbuildignore`; that strips the index).
7. `devtools::check_win_devel()` (and optionally `rhub::rhub_check()`).
   **This is the source of truth for the aspell NOTE** — the local `--as-cran`
   does *not* run the CRAN incoming aspell step unless `aspell` is installed,
   so it under-reports. Reconcile the `## NOTE disposition` section of
   `cran-comments.md` against the *win-builder* `00check.log`. See "Known
   benign NOTE" below.
8. `urlchecker::url_check()` and `tools::package_dependencies(reverse = TRUE)`
   (revdeps must be handled; currently 0). doi.org links may 403 to automated
   checkers but resolve in browsers — note, don't chase.
9. Submit: `devtools::submit_cran()` (writes `CRAN-SUBMISSION`). Submit from a
   checkout of the release tree (normally `dev` at the release number, or
   `main` once merged) — the shipped tarball content is what matters, not which
   branch the working tree was on.

### Known benign NOTE

The CRAN incoming check reports one expected NOTE the local `--as-cran` does
not show (without `aspell` installed):

- **Possibly misspelled words in DESCRIPTION** — proper nouns (`Naftel`,
  `Rajeswaran`), the acronym `UAB`, the `et al.` citation, and the domain
  term `multiphase`. All intentional; not misspellings. `inst/WORDLIST`
  feeds only the `spelling` test, **not** the CRAN aspell check, so this
  recurs by design — document it, don't try to suppress it.

The canonical wording lives in `cran-comments.md` `## NOTE disposition`;
keep the two in agreement.

## On CRAN acceptance

1. Tag the release and push the tag (annotated, mirroring existing tags), on
   `main` (the released line):

   ```sh
   git tag -a vX.Y.Z -m "TemporalHazard vX.Y.Z — accepted & published on CRAN <date>" <main-HEAD>
   git push origin vX.Y.Z
   ```

   `CRAN-SUBMISSION`'s `SHA:` records the tree `submit_cran()` ran from; if you
   submitted from `dev`, the shipped content still matches `main` HEAD (they
   differ only in `.Rbuildignore`'d files), so tag `main`.

2. Publish a GitHub Release from the tag, using that version's `NEWS.md`
   section as the body, marked as the latest release:

   ```sh
   gh release create vX.Y.Z --verify-tag --latest \
     --title "TemporalHazard X.Y.Z" --notes-file <news-section.md>
   ```

3. **Nothing to bump at this step.** `dev` stays at `X.Y.Z`, the number that
   just shipped. Do not run `usethis::use_dev_version()` — it writes a fourth
   `.9000` component, which this package does not use. The version moves later,
   when work moves it: the patch digit once a fix actually lands.

   Steps 1 and 2 (tag push, GitHub Release) are the standard release-op
   exception to the branch + PR rule.

The badge block (CRAN status + cranlogs + R-CMD-check / codecov / lint /
pkgdown) was aligned at the first CRAN release (`1.0.3`) and the manual
version-badge CI (`update-version-badge.yaml`) retired then — no badge work on
subsequent releases.

All subsequent development happens at `X.Y.Z` on `dev`, patch-bumping as
fixes land, until the next release decides the number to submit.

## Historical note

- Releases through `1.0.2` used clean three-part versions only.
- A `.9000` dev convention was adopted after `1.0.2`; `1.0.3` (accepted
  2026-05-29) was the first to use it, under a **two-branch** model where
  `main` carried the released `.9000` and `dev` the next-version `.9000`.
  **That convention was dropped entirely on 2026-08-17** — see the versioning
  section above. Versions are three digits, always.
- **2026-06-12:** moved to the **single-line, decide-the-bump-at-release**
  model documented above, after the `dev = 1.2.0.9000` mislabel (a minor bump
  pre-stamped before the work existed, while `1.1.0` was frozen on `main` but
  never shipped). `1.1.0` folded that work in and was accepted & published on
  CRAN 2026-06-12; `main` stays at the released number, and since 2026-08-17
  so does `dev`.
- Existing tags: `v0.1.0`, `v0.9.3`, `v1.0.0`, `v1.0.1`, `v1.0.3`, `v1.1.0`.
