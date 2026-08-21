# Evaluate a translated job's emitted chunks the way rendering the .qmd
# would: in a clean environment holding only the input data. Asserting the
# SHAPE of emitted calls cannot distinguish a document that runs from one
# that errors -- every defect in #151/#152 was invisible to shape assertions
# and visible to this in one line.

render_sim <- function(job, data = list()) {
  env <- new.env(parent = globalenv())
  for (nm in names(data)) assign(nm, data[[nm]], envir = env)
  res <- character(0)
  for (nm in names(job$calls)) {
    # job$calls holds quoted calls emitted by this package's own
    # hzr_translate_sas(), not arbitrary external input -- evaluating them
    # is exactly what rendering the translated .qmd does.
    out <- tryCatch({
      eval(job$calls[[nm]], env)
      "ok"
    }, error = function(e) paste0("ERROR: ", conditionMessage(e)))
    res[[nm]] <- out
  }
  # all(character(0) == "ok") is TRUE, so an empty `res` would report ok
  # having evaluated nothing -- require at least one call ran.
  list(ok = length(res) > 0 && all(res == "ok"), results = res, env = env)
}
