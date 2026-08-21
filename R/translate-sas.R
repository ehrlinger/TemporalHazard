# translate-sas.R -- the public entry point that wires the SAS lexer, block
# extraction, PROC HAZARD/HAZPRED parsers and the Quarto renderer together.
#
# This is the feature's only export. Everything upstream (sas-lex.R,
# sas-parse-job.R, sas-parse-parms.R, sas-render-qmd.R, sas-job.R) is
# internal and reachable only through here.

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
  if (!is.null(librefs) &&
      (!is.character(librefs) || is.null(names(librefs)) ||
       any(!nzchar(names(librefs))))) {
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
  inhaz <- NULL
  outhaz <- NULL
  grid <- NULL

  for (b in blocks) {
    if (identical(b$proc, "HAZARD")) {
      r <- tryCatch(.hzr_parse_hazard(b), error = function(e) {
        stop("failed to parse PROC HAZARD block in ", basename(path), ": ",
             conditionMessage(e), call. = FALSE)
      })
      if (!is.null(r$status_call)) calls$status <- r$status_call
      calls$fit <- r$call
      outhaz <- r$outhaz
    } else {
      r <- tryCatch(.hzr_parse_hazpred(b, txt), error = function(e) {
        stop("failed to parse PROC HAZPRED block in ", basename(path), ": ",
             conditionMessage(e), call. = FALSE)
      })
      inhaz <- r$inhaz
      grid <- r$grid
      calls$pred <- r$call
      if (!is.null(r$call_haz)) calls$pred_haz <- r$call_haz
    }
    untr <- rbind(untr, r$untranslated)
    seen <- seen + r$tokens_seen
    mapped <- mapped + r$tokens_mapped
  }

  # The grid is an assignment, placed before the predict() chunks that use it.
  if (!is.null(grid) && !is.null(calls$pred[["newdata"]])) {
    nm <- as.character(calls$pred[["newdata"]])
    calls <- c(list(grid = call("<-", as.name(nm), grid)), calls)
  }

  # --- INHAZ resolution: this job's own OUTHAZ, then librefs, then fail ----
  resolved <- FALSE
  if (!is.null(inhaz)) {
    if (!is.null(outhaz) && identical(inhaz, outhaz)) {
      resolved <- TRUE
    } else if (!is.null(librefs)) {
      lib <- sub("[.].*$", "", inhaz)
      mem <- tolower(sub("^[^.]*[.]", "", inhaz))
      if (lib %in% names(librefs)) {
        dir <- unname(librefs[[lib]])
        # A libref names a directory of SAS datasets, so the member's
        # on-disk form is <member>.sas7bdat -- a bare member name matches
        # neither branch hzr_read_outhaz() dispatches on and cannot run.
        # If the librefs value already carries a recognised extension, an
        # analyst who has already converted the fit (e.g. to .rds) is
        # naming that file directly, verbatim.
        read_expr <- if (grepl("[.](rds|sas7bdat)$", dir, ignore.case = TRUE)) {
          bquote(hzr_read_outhaz(.(dir)))
        } else {
          bquote(hzr_read_outhaz(
            file.path(.(dir), .(paste0(mem, ".sas7bdat")))))
        }
        calls <- c(
          list(fit = call("<-", as.name("fit"), read_expr)),
          calls
        )
        resolved <- TRUE
      }
    }
  }

  # A hazard() call that reads a SAS DATA= dataset by name cannot actually
  # run: the DATA step that built it is out of scope for this translator.
  # Rather than let the document fail later with an obscure
  # "object 'AVCS' not found", prepend a chunk that fails loudly and explains
  # what the reader still has to supply. This must come before every other
  # chunk, including "fit" and the grid assignment above.
  if (!is.null(calls$fit) && !is.null(calls$fit[["data"]])) {
    dname <- as.character(calls$fit[["data"]])
    guard <- bquote(
      if (!exists(.(dname))) {
        stop("This job read ", .(dname), " from a SAS DATA step, which ",
             "hzr_translate_sas() does not translate. Assign ", .(dname),
             " before rendering.")
      }
    )
    calls <- c(list(data = guard), calls)
  }

  job <- .hzr_sas_job(
    source = list(path = path,
                  checksum = unname(tools::md5sum(path))),
    calls = calls, grid = grid, inhaz = inhaz, outhaz = outhaz,
    untranslated = untr,
    coverage = list(tokens_seen = seen, tokens_mapped = mapped)
  )
  job$inhaz_resolved <- resolved
  .hzr_validate_sas_job(job)

  if (!is.null(inhaz) && !resolved) {
    warning("unresolved INHAZ=", inhaz, " in ", basename(path),
            "; the emitted document will stop() rather than render.",
            call. = FALSE)
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
