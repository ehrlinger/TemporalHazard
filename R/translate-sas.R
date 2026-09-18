# translate-sas.R -- the public entry point that wires the SAS lexer, block
# extraction, PROC HAZARD/HAZPRED parsers and the Quarto renderer together.
#
# This is the feature's only export. Everything upstream (sas-lex.R,
# sas-parse-job.R, sas-parse-parms.R, sas-render-qmd.R, sas-job.R) is
# internal and reachable only through here.

#' Reserve a stable, unique name for a chunk/call.
#'
#' The first block of a kind keeps the bare name (`fit`, `pred`, ...), so
#' existing chunk labels, `job$calls$fit` accesses and tests are unaffected;
#' a second block of the same kind becomes `fit_2`, a third `fit_3`, and so
#' on. `sas-render-qmd.R` uses `names(job$calls)` directly as chunk labels,
#' and Quarto rejects duplicate labels, so this is the single source of
#' truth for both the `calls` list key and the label.
#' @noRd
.hzr_next_call_name <- function(calls, base) {
  if (is.null(calls[[base]])) return(base)
  i <- 2L
  while (!is.null(calls[[paste0(base, "_", i)]])) i <- i + 1L
  paste0(base, "_", i)
}

#' Point a `predict()` call at a different fitted-model variable.
#'
#' `.hzr_parse_hazpred()` always builds `predict(fit, ...)`; the literal
#' name `fit` is a placeholder, first positional argument. With more than
#' one fit in a job, a given `PROC HAZPRED` block's `INHAZ=` may resolve to
#' `fit_2` instead; this swaps the placeholder for the resolved name without
#' having to thread the name through `.hzr_parse_hazpred()` itself.
#' @noRd
.hzr_retarget_fit <- function(call, fit_ref) {
  if (is.null(call)) return(NULL)
  # A refused grid makes this a stop() chunk, whose first argument is the
  # message, not the fit. Overwriting it would replace the explanation with
  # a bare symbol -- a refusal that no longer says what it refused.
  if (!identical(call[[1L]], as.name("predict"))) return(call)
  call[[2L]] <- as.name(fit_ref)
  call
}

#' Translate a SAS HAZARD job into a Quarto document
#'
#' Reads a SAS program containing `PROC HAZARD` and/or `PROC HAZPRED` blocks and
#' emits a Quarto document of the equivalent [hazard()] and [predict.hazard()]
#' calls.
#'
#' **Experimental:** a job that translates does render (the emitted
#' `hazard()` chunk binds its fit and asks for an actual fit), but this is a
#' translation aid, not a turnkey reproduction, and some SAS constructs are
#' refused rather than translated. See the Experimental section below.
#'
#' Constructs the translator does not cover are recorded on the returned object
#' and rendered as visible callouts, never dropped. A `PROC HAZPRED` job whose
#' `INHAZ=` fitted model cannot be located emits a `stop()`, so the document
#' fails to render rather than reporting over a model it did not load.
#'
#' A job may contain more than one `PROC HAZARD` and/or `PROC HAZPRED` block.
#' Every block is preserved: the first of a kind keeps the bare chunk name
#' (`fit`, `pred`, `pred_haz`), later ones get `fit_2`, `fit_3`, `pred_2`, and
#' so on, so the emitted document contains every model, not just the last one
#' seen. `job$outhaz` is therefore a character vector, one element per `PROC
#' HAZARD` block that set `OUTHAZ=` (in the order the blocks appear). Each
#' `PROC HAZPRED` block's `INHAZ=` is resolved independently against that
#' vector, matching the most recently written `OUTHAZ=` at that point in the
#' file, mirroring SAS itself, where a later `OUTHAZ=` write overwrites the
#' dataset an earlier one wrote under the same name. If a job's own `OUTHAZ=`
#' values don't cover it, `librefs` is tried next; distinct external
#' `INHAZ=` values each get their own loaded-fit chunk. When neither
#' resolves and the job holds more than one local fit, which fit a
#' `predict()` call belongs to is genuinely unknown; the emitted call falls
#' back to referencing `fit` and the ambiguity is recorded in
#' `untranslated`, never guessed at silently.
#'
#' @param path Path to a `.sas` file.
#' @param out_dir Directory to write the `.qmd` into. `NULL` (default) parses
#'   without writing.
#' @param librefs Optional named character vector mapping SAS librefs to
#'   directories, e.g. `c(EX = "estimates")`, used to resolve `INHAZ=`. The
#'   resolved member is read as `<member>.sas7bdat`; give a librefs value
#'   that already ends in `.rds` or `.sas7bdat` (a specific file, not a
#'   directory) to name an already-converted fit directly, e.g.
#'   `c(EX = "estimates/hzdeath.rds")`.
#' @return An `hzr_sas_job` object, invisibly. `$calls` holds the emitted
#'   calls keyed by chunk label, `$grid` the last prediction grid seen and
#'   `$inhaz` the first unresolved `INHAZ=` (not all of each, when a job has
#'   several), `$untranslated` the recorded gaps, `$coverage` the token
#'   counts, and `$notes` any callout attached to a chunk by label, emitted
#'   immediately above that chunk in the rendered document.
#'
#' @section What the emitted document needs to run:
#' Base R at the version this package declares in `DESCRIPTION` (`Depends`),
#' plus TemporalHazard itself. The emitted chunks never reach into this
#' package's internals, so the document renders in your session, and in a
#' colleague's, without anything further installed, with one exception. A
#' translated `SELECTION` screen on a multiphase job with interval-censored
#' rows (an `ICENSOR` job; `LCENSOR` is left truncation and does not need
#' it) needs the suggested package \pkg{numDeriv}. Without it a screen that
#' tests an entry stops and says so. A screen that completes but cannot test
#' a removal (a `BACKWARD` screen, whose base has no usable variance) is
#' caught by the emitted check, which warns and names the variables.
#'
#' @section Experimental:
#' The emitted document renders: the `hazard()` chunk binds its fit to a name
#' and passes `fit = TRUE`, so the `predict()` chunks have something to
#' predict from. A `SELECTION` statement is translated into an
#' `hzr_stepwise()` call carrying the job's own candidates, per-variable
#' flags and thresholds (#160); the screen is real, and may select a
#' different model than `PROC HAZARD` did, which the emitted document says
#' in a callout above the chunk. `SELECTION` options with no faithful
#' translation are refused outright rather than mistranslated, emitting a
#' `stop()` in place of the fit: `FAST`, `MAXVARS`, `RESTRICT`, a
#' per-variable `MOVE=` or `ORDER=`, and a variable held by
#' `/I` in one phase but movable in another. So is `LCENSOR`
#' combined with `ICENSOR`, which one `time_lower` argument cannot express
#' (#155). Prediction grids the parser cannot resolve are refused whole, and
#' the `predict()` chunks that would have read such a grid become a `stop()`
#' naming it, rather than a `predict(newdata = )` over a name no chunk
#' builds. An unresolved `INHAZ=` stops the render on purpose.
#'
#' On a fit loaded from an external `INHAZ=` dataset, point predictions work
#' but `se.fit = TRUE` is refused when `PROC HAZARD` estimated a late shape
#' parameter on a composite scale (the generic unconstrained three-phase
#' case, not an exotic one). A translated `PROC HAZPRED` block asks for
#' confidence limits unless the SAS job says `NOCL`, so such a job stops at
#' its `predict()` chunks.
#'
#' @section Repeated events:
#' A job that builds its fit input with the SAS macro `%repeat` has that call
#' translated to [hzr_repeated_events()]. The emitted chunk renames the
#' function's output columns to the names the job gives them, upper-cased like
#' every name the translator emits, so the fit reads them. An input column named
#' like one of those outputs is dropped first, with a warning, because the macro
#' overwrites it. The macro's input is built by the job's own DATA steps, which
#' are not translated, so the document stops until that input is assigned.
#' Any step between the macro and the fit that names the macro's output, or
#' uses a macro variable that might, is not translated either; a plain
#' `PROC SORT` that only reorders rows is the one step let through. Its chunk
#' stops and quotes the step, for the reader to replace with R code or to delete
#' if the step leaves the output unchanged. The same holds between two
#' `%repeat` calls when the second reads the first's output. A step that changes
#' the output without naming it, such as a macro that writes it internally, is
#' not detected.
#'
#' @section Comparing a translated fit against a SAS listing:
#' The emitted `hzr_phase("g3", ...)` call carries the `TAU`, `GAMMA`,
#' `ALPHA` and `ETA` that `PARMS` specified. `PROC HAZARD` does not always
#' report the values it was given: when `ALPHA` is fixed at 1 it rewrites the
#' late phase before fitting, pinning `TAU` at 1 and folding the two shape
#' parameters into one exponent, because at `alpha = 1, tau = 1` the G3 form
#' collapses to \eqn{t^{\gamma\eta}} and only that product is identified.
#' With `ETA` fixed it reports `GAMMA` as \eqn{\gamma\eta} and `ETA` as 1;
#' otherwise it reports `ETA` as \eqn{\gamma\eta} and `GAMMA` as 1.
#'
#' So a listing can print a `GAMMA`/`ETA` pair that differs from the job's own
#' `PARMS` text, and from the emitted call, while describing the same fitted
#' model. Compare the **product** \eqn{\gamma\eta}, not the two parameters
#' separately, and expect `TAU` to read 1 and fixed. The log-likelihood,
#' `MUL` and every prediction agree regardless. Most jobs are unaffected in
#' practice: `GAMMA = 1 FIXGAMMA` is the usual companion to `ALPHA = 1
#' FIXALPHA`, and it makes the rewrite an identity.
#'
#' The one case that is a genuine model difference (`ALPHA` fixed at 1 with
#' `GAMMA` and `ETA` *both* estimated, where `PROC HAZARD` fixes `ETA` and
#' fits one parameter fewer than [hazard()] would) is recorded in
#' `$untranslated` rather than left to be discovered in the comparison.
#'
#' `$coverage` counts tokens the parser recognised; it is not evidence that
#' the emitted calls execute, and a job can report full coverage with an
#' empty `$untranslated` while its document still errors on render. See the
#' 1.2.2 `NEWS.md` entry. The function's API, the `hzr_sas_job` field layout
#' and the emitted document format are all expected to change.
#' @examples
#' \donttest{
#' job <- hzr_translate_sas(
#'   system.file("extdata", "hz-example.sas", package = "TemporalHazard")
#' )
#' }
#' @export
hzr_translate_sas <- function(path, out_dir = NULL, librefs = NULL) {
  stopifnot(is.character(path), length(path) == 1L)
  if (!file.exists(path)) stop("no such file: ", path, call. = FALSE)
  # anyNA() is a separate term from the nzchar() check on purpose: nzchar()
  # defaults to keepNA = FALSE, so nzchar(NA_character_) is TRUE and an NA
  # name would otherwise pass a check whose message promises a named vector.
  if (!is.null(librefs) &&
      (!is.character(librefs) || is.null(names(librefs)) ||
       anyNA(names(librefs)) || any(!nzchar(names(librefs))))) {
    stop('librefs must be a named character vector, e.g. c(EX = "estimates").',
         call. = FALSE)
  }

  txt <- .hzr_sas_normalise(readLines(path, warn = FALSE))
  blocks <- .hzr_sas_blocks(txt)
  # A %repeat call alone is a dataset-building job (the tp.bd.* templates),
  # not a HAZARD job, and there is nothing to fit.
  if (!any(vapply(blocks, function(b) b$proc != "REPEAT", logical(1L)))) {
    stop("no HAZARD or HAZPRED block found in ", path, call. = FALSE)
  }

  calls <- list()
  untr <- .hzr_untranslated_frame()
  # Callouts attached to a named chunk, rendered immediately ABOVE it.
  notes <- list()
  seen <- 0L
  mapped <- 0L
  grid <- NULL

  fits <- list()             # list(list(slot = <chr>, outhaz = <chr|NULL>))
  guarded_data <- character(0L)
  loaded_ext <- list()       # raw INHAZ string -> loaded chunk's slot name
  first_unresolved_inhaz <- NULL
  n_unresolved_inhaz <- 0L
  repeat_scan <- list()      # %repeat OUT= -> offset its rewrite scan resumes from

  for (b in blocks) {
    if (identical(b$proc, "REPEAT")) {
      r <- .hzr_parse_repeat(b)
      # A %repeat whose IN= is an earlier %repeat's OUT= reads that dataset as
      # the job left it, so a step between the two that changes it is the same
      # untranslated rewrite a fit would read; stop on it here, ahead of the
      # second macro's chunk.
      if (!is.null(r$in_name) && !is.null(repeat_scan[[r$in_name]])) {
        rs <- .hzr_rewrite_stops(substring(txt, repeat_scan[[r$in_name]] + 1L, b$start - 1L), r$in_name)
        for (cl in rs$calls) calls[[.hzr_next_call_name(calls, "rewrite")]] <- cl
        untr <- rbind(untr, rs$untranslated)
        repeat_scan[[r$in_name]] <- b$end
      }
      # The macro's input is built by the job's own DATA steps, which this
      # translator does not translate: the same loud guard a PROC HAZARD
      # DATA= gets, once per name.
      if (!is.null(r$in_name) && !(r$in_name %in% guarded_data)) {
        calls[[.hzr_next_call_name(calls, "data")]] <- bquote(
          if (!exists(.(r$in_name))) {
            stop("This job built ", .(r$in_name), " in SAS DATA steps, which ",
                 "hzr_translate_sas() does not translate. Assign ", .(r$in_name),
                 " as it stood at the job's %repeat call, with its columns named ",
                 "as the job spells them, in upper case, before rendering.")
          }
        )
        guarded_data <- c(guarded_data, r$in_name)
      }
      calls[[.hzr_next_call_name(calls, "repeated")]] <- r$call
      # OUT= is now built by a chunk, so a fit reading it needs no guard.
      if (!is.null(r$out_name)) {
        guarded_data <- c(guarded_data, r$out_name)
        repeat_scan[[r$out_name]] <- b$end
      }
    } else if (identical(b$proc, "HAZARD")) {
      r <- tryCatch(.hzr_parse_hazard(b), error = function(e) {
        stop("failed to parse PROC HAZARD block in ", basename(path), ": ",
             conditionMessage(e), call. = FALSE)
      })

      # A step between %repeat and this fit that may change the macro's
      # OUT= is job code this translator does not fold in. Fitting without it
      # fits data SAS did not fit, so stop here -- ahead of the status chunk,
      # which writes into the same data frame. The scan resumes where the
      # last one ended, so one rewrite stops once however many fits follow.
      dname <- if (is.null(r$call[["data"]])) NULL else as.character(r$call[["data"]])
      if (!is.null(dname) && !is.null(repeat_scan[[dname]])) {
        rs <- .hzr_rewrite_stops(substring(txt, repeat_scan[[dname]] + 1L, b$start - 1L), dname)
        for (cl in rs$calls) calls[[.hzr_next_call_name(calls, "rewrite")]] <- cl
        untr <- rbind(untr, rs$untranslated)
        repeat_scan[[dname]] <- b$end
      }

      # A hazard() call that reads a SAS DATA= dataset by name cannot
      # actually run: the DATA step that built it is out of scope for this
      # translator. Rather than let the document fail later with an obscure
      # "object 'AVCS' not found", insert a chunk that fails loudly and
      # explains what the reader still has to supply, right before this
      # fit -- once per distinct dataset name, however many fits reference it.
      if (!is.null(dname)) {
        if (!(dname %in% guarded_data)) {
          guard_slot <- .hzr_next_call_name(calls, "data")
          calls[[guard_slot]] <- bquote(
            if (!exists(.(dname))) {
              stop("This job read ", .(dname), " from a SAS DATA step, which ",
                   "hzr_translate_sas() does not translate. Assign ", .(dname),
                   " before rendering.")
            }
          )
          guarded_data <- c(guarded_data, dname)
        }
      }

      if (!is.null(r$status_call)) {
        status_slot <- .hzr_next_call_name(calls, "status")
        calls[[status_slot]] <- r$status_call
      }
      fit_slot <- .hzr_next_call_name(calls, "fit")
      # Bind the fit: predict() chunks reference the fit by its slot name, and
      # a bare hazard(...) call binds nothing, so those chunks failed with
      # "object 'fit' not found" -- or worse, silently used an unrelated
      # object of that name already in the rendering session (#151).
      if (is.null(r$stepwise_call)) {
        calls[[fit_slot]] <- call("<-", as.name(fit_slot), r$call)
      } else {
        # A SELECTION job is two chunks: the shape-fixed base fit, then the
        # screen. `fit <- hzr_stepwise(fit, ...)` would be self-referential,
        # so the base is bound under <slot>_base and the screen's RESULT
        # takes the slot name every predict() chunk already references.
        base_slot <- paste0(fit_slot, "_base")
        calls[[base_slot]] <- call("<-", as.name(base_slot), r$call)
        sw <- r$stepwise_call
        sw[[2L]] <- as.name(base_slot)
        calls[[fit_slot]] <- call("<-", as.name(fit_slot), sw)
        notes[[fit_slot]] <- .hzr_selection_divergence_note()
        if (!is.null(r$screen_check_call)) {
          chk <- do.call(substitute,
                         list(r$screen_check_call,
                              list(fit = as.name(fit_slot),
                                   fit_label = fit_slot)))
          calls[[.hzr_next_call_name(calls, "screen_check")]] <- chk
        }
      }
      fits[[length(fits) + 1L]] <- list(slot = fit_slot, outhaz = r$outhaz)
    } else {
      r <- tryCatch(.hzr_parse_hazpred(b, txt), error = function(e) {
        stop("failed to parse PROC HAZPRED block in ", basename(path), ": ",
             conditionMessage(e), call. = FALSE)
      })

      # --- which fit does this predict() belong to? --------------------
      # 1. this job's own OUTHAZ, most-recently-written match wins (mirrors
      #    SAS: a later OUTHAZ= write overwrites the dataset an earlier one
      #    wrote under the same name, and blocks are visited in file order).
      # 2. else librefs, one loaded-fit chunk per distinct external INHAZ.
      # 3. else: unresolved. With more than one local fit this is genuinely
      #    ambiguous -- record it rather than guess. With 0 or 1 local fits
      #    it is the pre-existing "unresolved INHAZ" case, reported via
      #    job$inhaz/job$inhaz_resolved and the stop() chunk in the render.
      fit_ref <- "fit"
      this_inhaz <- r$inhaz
      if (!is.null(this_inhaz)) {
        match_slot <- NULL
        if (length(fits)) {
          for (fi in rev(seq_along(fits))) {
            if (identical(fits[[fi]]$outhaz, this_inhaz)) {
              match_slot <- fits[[fi]]$slot
              break
            }
          }
        }
        if (!is.null(match_slot)) {
          fit_ref <- match_slot
        } else if (!is.null(librefs) &&
                   sub("[.].*$", "", this_inhaz) %in% names(librefs)) {
          if (!is.null(loaded_ext[[this_inhaz]])) {
            fit_ref <- loaded_ext[[this_inhaz]]
          } else {
            lib <- sub("[.].*$", "", this_inhaz)
            mem <- tolower(sub("^[^.]*[.]", "", this_inhaz))
            dir <- unname(librefs[[lib]])
            # A libref names a directory of SAS datasets, so the member's
            # on-disk form is <member>.sas7bdat -- a bare member name
            # matches neither branch hzr_read_outhaz() dispatches on and
            # cannot run. If the librefs value already carries a recognised
            # extension, an analyst who has already converted the fit (e.g.
            # to .rds) is naming that file directly, verbatim.
            read_expr <- if (grepl("[.](rds|sas7bdat)$", dir,
                                    ignore.case = TRUE)) {
              bquote(hzr_read_outhaz(.(dir)))
            } else {
              bquote(hzr_read_outhaz(
                file.path(.(dir), .(paste0(mem, ".sas7bdat")))))
            }
            fit_ref <- .hzr_next_call_name(calls, "fit")
            calls[[fit_ref]] <- call("<-", as.name(fit_ref), read_expr)
            loaded_ext[[this_inhaz]] <- fit_ref
          }
        } else {
          n_unresolved_inhaz <- n_unresolved_inhaz + 1L
          if (is.null(first_unresolved_inhaz)) {
            first_unresolved_inhaz <- this_inhaz
          }
          if (length(fits) > 1L) {
            untr <- rbind(untr, .hzr_untranslated_frame(
              NA_integer_, paste0("INHAZ=", this_inhaz),
              sprintf(paste("does not match any of this job's %d OUTHAZ",
                             "values; ambiguous which fit() this predict()",
                             "belongs to, defaulting to `fit`"),
                      length(fits))
            ))
          }
        }
      }

      pred_call <- .hzr_retarget_fit(r$call, fit_ref)
      pred_haz_call <- .hzr_retarget_fit(r$call_haz, fit_ref)

      # The grid is an assignment, placed right before the predict() chunk
      # that uses it.
      if (!is.null(r$grid) && !is.null(pred_call[["newdata"]])) {
        nm <- as.character(pred_call[["newdata"]])
        grid_slot <- .hzr_next_call_name(calls, "grid")
        calls[[grid_slot]] <- call("<-", as.name(nm), r$grid)
        grid <- r$grid
      }

      calls[[.hzr_next_call_name(calls, "pred")]] <- pred_call
      if (!is.null(pred_haz_call)) {
        calls[[.hzr_next_call_name(calls, "pred_haz")]] <- pred_haz_call
      }
    }
    untr <- rbind(untr, r$untranslated)
    seen <- seen + r$tokens_seen
    mapped <- mapped + r$tokens_mapped
  }

  outhaz_vec <- unlist(lapply(fits, function(x) x$outhaz))
  if (!length(outhaz_vec)) outhaz_vec <- NULL

  job <- .hzr_sas_job(
    source = list(path = path,
                  checksum = unname(tools::md5sum(path))),
    calls = calls, grid = grid, inhaz = first_unresolved_inhaz,
    outhaz = outhaz_vec, untranslated = untr,
    coverage = list(tokens_seen = seen, tokens_mapped = mapped),
    notes = notes
  )
  job$inhaz_resolved <- is.null(first_unresolved_inhaz)
  .hzr_validate_sas_job(job)

  if (!is.null(first_unresolved_inhaz)) {
    msg <- if (n_unresolved_inhaz > 1L) {
      sprintf(paste("%d unresolved INHAZ values (first: %s) in %s;",
                     "the emitted document will stop() rather than render."),
              n_unresolved_inhaz, first_unresolved_inhaz, basename(path))
    } else {
      sprintf(paste("unresolved INHAZ=%s in %s;",
                     "the emitted document will stop() rather than render."),
              first_unresolved_inhaz, basename(path))
    }
    warning(msg, call. = FALSE)
  }
  if (nrow(untr)) {
    warning(nrow(untr), " untranslated construct(s) in ", basename(path), ": ",
            paste(unique(untr$construct), collapse = ", "), call. = FALSE)
  }

  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    writeLines(
      .hzr_render_qmd(job),
      file.path(out_dir, sub("[.]sas$", ".qmd", basename(path)))
    )
  }
  invisible(job)
}


#' The callout a translated SELECTION job carries, above its screen.
#'
#' `hzr_stepwise()` runs the job's own candidates, flags and thresholds, but
#' it cannot be expected to reach `PROC HAZARD`'s selected set: SAS uses
#' approximate variances during selection, which the entry statistic here
#' reproduces (except for a candidate refitted because its information is
#' indefinite) but the Wald removal tests do not, and `force_in` is not
#' phase-keyed where SAS's `/I` is. The
#' divergence is recorded against a real fixture in
#' `tests/testthat/test-sas-parity.R` (hm.death.AVC). The reader meets this
#' before the code, which is why it is a note on the chunk rather than a row
#' in `$untranslated` (#160).
#' @noRd
.hzr_selection_divergence_note <- function() {
  list(
    title = paste("SELECTION: this screen may select a different model than",
                  "PROC HAZARD did"),
    body = paste(
      "This job's SELECTION statement is translated into hzr_stepwise() with",
      "the job's own candidates, per-variable flags and SLENTRY/SLSTAY",
      "thresholds. The screen is real, and the selected model may still",
      "differ from the one PROC HAZARD chose, for reasons that cannot be",
      "tuned away. PROC HAZARD uses approximate variances during selection",
      "(it ignores the shaping-parameter covariances). The entry statistic",
      "here reproduces that approximation, except for a candidate whose",
      "information is indefinite, which is refitted and Wald-tested instead;",
      "the Wald tests behind the drop decisions use the full Hessian, so",
      "removals can differ too; and SAS's /I",
      "holds a variable in ONE phase,",
      "while hzr_stepwise()'s force_in is keyed by variable name across every",
      "phase, which is why a job whose /I variable is movable in another",
      "phase is refused rather than screened here. And PROC HAZARD's MOVE",
      "limit counts a variable's DELETIONS, separately for each phase: at",
      "its default of 1 a variable removed from a phase can never return to",
      "it, while this package's oscillation guard counts entries and exits",
      "together across every phase and lets a removed variable re-enter. So",
      "the screen here can re-enter variables PROC HAZARD would have kept",
      "out, and occasionally freeze one PROC HAZARD would still move. The",
      "divergence is recorded against the hm.death.AVC",
      "fixture in",
      "tests/testthat/test-sas-parity.R. Read the selected model as this",
      "package's screen of this job's candidates, not as a reproduction of",
      "the SAS run, and compare it against the SAS listing before relying on",
      "it. Any candidate the screen could not score is reported by the",
      "uncomputable-score check below it."
    )
  )
}
