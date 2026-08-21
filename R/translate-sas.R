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
#' `.hzr_parse_hazpred()` always builds `predict(fit, ...)` -- the literal
#' name `fit` is a placeholder, first positional argument. With more than
#' one fit in a job, a given `PROC HAZPRED` block's `INHAZ=` may resolve to
#' `fit_2` instead; this swaps the placeholder for the resolved name without
#' having to thread the name through `.hzr_parse_hazpred()` itself.
#' @noRd
.hzr_retarget_fit <- function(call, fit_ref) {
  if (is.null(call)) return(NULL)
  call[[2L]] <- as.name(fit_ref)
  call
}

#' Translate a SAS HAZARD job into a Quarto document
#'
#' Reads a SAS program containing `PROC HAZARD` and/or `PROC HAZPRED` blocks and
#' emits a Quarto document that reproduces the analysis with [hazard()] and
#' [predict.hazard()].
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
#' file -- mirroring SAS itself, where a later `OUTHAZ=` write overwrites the
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
#' @return An `hzr_sas_job` object, invisibly.
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
  if (!length(blocks)) {
    stop("no HAZARD or HAZPRED block found in ", path, call. = FALSE)
  }

  calls <- list()
  untr <- .hzr_untranslated_frame()
  seen <- 0L
  mapped <- 0L
  grid <- NULL

  fits <- list()             # list(list(slot = <chr>, outhaz = <chr|NULL>))
  guarded_data <- character(0L)
  loaded_ext <- list()       # raw INHAZ string -> loaded chunk's slot name
  first_unresolved_inhaz <- NULL
  n_unresolved_inhaz <- 0L

  for (b in blocks) {
    if (identical(b$proc, "HAZARD")) {
      r <- tryCatch(.hzr_parse_hazard(b), error = function(e) {
        stop("failed to parse PROC HAZARD block in ", basename(path), ": ",
             conditionMessage(e), call. = FALSE)
      })

      # A hazard() call that reads a SAS DATA= dataset by name cannot
      # actually run: the DATA step that built it is out of scope for this
      # translator. Rather than let the document fail later with an obscure
      # "object 'AVCS' not found", insert a chunk that fails loudly and
      # explains what the reader still has to supply, right before this
      # fit -- once per distinct dataset name, however many fits reference it.
      if (!is.null(r$call[["data"]])) {
        dname <- as.character(r$call[["data"]])
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
      calls[[fit_slot]] <- r$call
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
    coverage = list(tokens_seen = seen, tokens_mapped = mapped)
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
