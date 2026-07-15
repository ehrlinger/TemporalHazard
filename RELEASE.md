# Release process

Maintainer checklist for TemporalHazard. Not shipped (`.Rbuildignore`d).

## Versioning convention (single-line, decide the bump at release)

- **CRAN releases** carry a clean three-part version: `1.1.0`.
- **Development versions** carry a fourth `.9000` component (`1.1.0.9000`), so
  it is unambiguous whether a given install is the CRAN release or a dev build.
- After a release `X.Y.Z` is accepted, the dev version becomes **`X.Y.Z.9000`**
  — the *just-released* version plus `.9000`. It deliberately does **not**
  encode the next release number.
- **Decide the next number at release time, from what is actually in NEWS** —
  breaking API → major, back-compatible features → minor, fixes only → patch.
  Then strip `.9000` and set that number. You never submit a `.9000` version
  to CRAN.

So the per-release cycle is:

```
1.1.0 (released)  ->  1.1.0.9000 (dev, number TBD)  ->  accumulate work
  ->  at release, decide from NEWS: 1.1.1 | 1.2.0 | 2.0.0  ->  submit
  ->  accepted  ->  <new>.9000 (dev) ...
```

**Do not pre-stamp the next number into the dev label.** A dev version of
`1.2.0.9000` (or `1.1.1.9000`) silently commits a decision that should be made
last, from the evidence. This is the mistake that produced the 2026-06 version
confusion (`dev` was labelled `1.2.0.9000` for what was really incremental
work, before `1.1.0` had even shipped). `1.1.0.9000` makes no such claim.

## Branch model (default: single line)

| Branch | Version | Role |
|--------|---------|------|
| `main` | the clean CRAN release (e.g. `1.1.0`) | the released/stable line; always == the latest version on CRAN. Tagged at each release. |
| `dev`  | `<released>.9000` (e.g. `1.1.0.9000`) | the working line; accumulates the next release. The `.9000` just means "ahead of the last release, number TBD". |

`main` carries the **clean** released version — no `.9000`. (Earlier releases
put `.9000` on `main`; that was dropped after `1.1.0`. A README/badge tweak on
`main` at the same version number is normal and fine — it does not need a
version marker.)

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
while `dev` becomes a diverged **next-major** line (e.g. `2.0.0.9000`) — for a
**genuine breaking rewrite** running in parallel with patch maintenance (the
ggRandomForests-v4 situation). It carries a real cost: every `main` hotfix must
be back-merged into `dev` or `dev` drifts from the patched line. Do **not** use
it as the default — it is not worth that tax for ordinary incremental cycles.

**Release flow (default single line):**

```
dev (X.Y.Z.9000) accumulates
  -> at release, decide the next number N from NEWS scope
  -> set DESCRIPTION + NEWS top heading to N (drop .9000) on dev
  -> run the pre-submission gate (below), submit to CRAN
  -> on acceptance: merge dev -> main (main = N, clean), tag vN, GitHub Release
  -> bump dev -> N.9000, open the next cycle
```

## Pre-submission checklist

1. `DESCRIPTION` `Version:` is a clean `X.Y.Z` (no `.9000`).
2. `NEWS.md` top heading is `# TemporalHazard X.Y.Z` and describes the
   user-visible changes since the last CRAN release.
3. `cran-comments.md` version heading matches `DESCRIPTION`; if a
   resubmission, each reviewer point is itemised with how it was addressed.
4. `devtools::document()` is clean and `man/` is in sync.
5. `R CMD check --as-cran` **with the manual built** → 0 errors, 0 warnings,
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
6. `devtools::check_win_devel()` (and optionally `rhub::rhub_check()`).
   **This is the source of truth for the aspell NOTE** — the local `--as-cran`
   does *not* run the CRAN incoming aspell step unless `aspell` is installed,
   so it under-reports. Reconcile the `## NOTE disposition` section of
   `cran-comments.md` against the *win-builder* `00check.log`. See "Known
   benign NOTE" below.
7. `urlchecker::url_check()` and `tools::package_dependencies(reverse = TRUE)`
   (revdeps must be handled; currently 0). doi.org links may 403 to automated
   checkers but resolve in browsers — note, don't chase.
8. Submit: `devtools::submit_cran()` (writes `CRAN-SUBMISSION`). Submit from a
   checkout of the release tree (normally `dev` with the suffix stripped, or
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

3. Bump **`dev`** (not `main`) to the next development version — `main` stays
   clean at `X.Y.Z`:

   ```r
   # on a branch off dev:
   usethis::use_dev_version()   # DESCRIPTION -> X.Y.Z.9000
   ```

   or manually: set `DESCRIPTION` `Version:` to `X.Y.Z.9000` and add a NEWS
   heading `# TemporalHazard X.Y.Z.9000 (development version)`. PR to `dev`.

4. Step 3 changes tracked files, so per the repo's mandatory git workflow it
   goes through a branch + PR. Steps 1 and 2 (tag push, GitHub Release) are the
   standard release-op exception.

The badge block (CRAN status + cranlogs + R-CMD-check / codecov / lint /
pkgdown) was aligned at the first CRAN release (`1.0.3`) and the manual
version-badge CI (`update-version-badge.yaml`) retired then — no badge work on
subsequent releases.

All subsequent development happens at `X.Y.Z.9000` on `dev` until the next
release decides and strips the number for submission.

## Historical note

- Releases through `1.0.2` used clean three-part versions only.
- The `.9000` dev convention was adopted after `1.0.2`; `1.0.3` (accepted
  2026-05-29) was the first to use it, under a **two-branch** model where
  `main` carried the released `.9000` and `dev` the next-version `.9000`.
- **2026-06-12:** moved to the **single-line, decide-the-bump-at-release**
  model documented above, after the `dev = 1.2.0.9000` mislabel (a minor bump
  pre-stamped before the work existed, while `1.1.0` was frozen on `main` but
  never shipped). `1.1.0` folded that work in and was accepted & published on
  CRAN 2026-06-12; `main` now stays clean at the released number and `dev` is
  `<released>.9000`.
- Existing tags: `v0.1.0`, `v0.9.3`, `v1.0.0`, `v1.0.1`, `v1.0.3`, `v1.1.0`.
