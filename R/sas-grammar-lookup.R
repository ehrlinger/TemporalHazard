#' Resolve a SAS keyword to its parser token within a proc and lexer context.
#'
#' Returns NA_character_ when the keyword is unknown in that context. Callers
#' must treat NA as untranslated rather than falling back to a context-free
#' match: `M` is the early-phase shape parameter in PARM context and MOVE in
#' PHOP, and guessing between them silently changes the model.
#' @noRd
.hzr_sas_token <- function(keyword, proc, context) {
  g <- .hzr_sas_grammar
  hit <- g$keyword == keyword & g$proc == proc & g$context == context
  if (!any(hit)) return(NA_character_)
  unique(g$token[hit])[[1L]]
}
