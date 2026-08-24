# sas-render-qmd.R -- turn a parsed hzr_sas_job into Quarto source.
#
# Calls are stored unevaluated on the job, so rendering is deparse(), not
# string templating: the parser already produced real R calls, and this file
# only lays them out as .qmd chunks. The one behaviour that must not soften:
# an unresolved INHAZ= emits stop(), not a comment. A PROC HAZPRED job whose
# fitted model cannot be located must produce a document that fails to
# render, rather than a beautiful report computed over a model it never
# loaded.

#' Render an hzr_sas_job as Quarto source.
#'
#' Calls are stored unevaluated, so rendering is deparse() rather than string
#' templating. An unresolved INHAZ emits stop(), not a comment: the document
#' must fail to render rather than produce a report over a model it never
#' loaded.
#' @noRd
.hzr_render_qmd <- function(job) {
  out <- character(0)
  add <- function(...) out <<- c(out, ...)
  base <- basename(job$source$path)

  add("---",
      sprintf('title: "%s"', sub("[.]sas$", "", base)),
      "format: html",
      "---",
      "")
  add(sprintf("<!-- Translated from %s (md5 %s) by hzr_translate_sas(). -->",
              base, substr(job$source$checksum, 1L, 12L)),
      "")
  add("```{r}", "#| label: setup", "library(TemporalHazard)", "```", "")

  if (!is.null(job$inhaz) && !isTRUE(job$inhaz_resolved)) {
    lib <- sub("[.].*$", "", job$inhaz)
    # Quote the libref. A SAS libref is nearly always a syntactic R name, but
    # one that collides with a reserved word (TRUE, if) would make the remedy
    # invalid R, and the whole point of the line is that it can be pasted.
    msg <- paste0(
      "unresolved INHAZ=", job$inhaz, ' -- pass librefs = c("', lib,
      '" = "<path>") to hzr_translate_sas()'
    )
    stop_call <- bquote(stop(.(msg)))
    add("```{r}",
        "#| label: inhaz-unresolved",
        deparse(stop_call, width.cutoff = 500L),
        "```", "")
  }

  for (nm in names(job$calls)) {
    add("```{r}",
        sprintf("#| label: %s", nm),
        deparse(job$calls[[nm]], width.cutoff = 60L),
        "```", "")
  }

  for (i in seq_len(nrow(job$untranslated))) {
    line <- job$untranslated$line[i]
    where <- if (is.na(line)) "" else sprintf(" (line %d)", line)
    add("::: {.callout-warning}",
        sprintf("## UNTRANSLATED: %s%s", job$untranslated$construct[i], where),
        job$untranslated$reason[i],
        ":::", "")
  }

  add(sprintf("<!-- coverage: %d/%d tokens mapped -->",
              job$coverage$tokens_mapped, job$coverage$tokens_seen))
  out
}
