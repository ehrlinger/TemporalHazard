# Runs a Quarto vignette's code up to and including the chunk labelled
# `label`, and returns the environment it ran in.
#
# Used by test-vignette-fit-step.R and by data-raw/vignette-fit-step.R, so
# both refit exactly the code the vignette shows. The labelled chunk may be
# `eval: false`, which purl() would comment out, so its body is read from the
# .qmd directly. eval() runs this package's own vignette code, nothing
# external.
run_vignette_through <- function(qmd, label) {
  src <- readLines(qmd, warn = FALSE)
  at <- grep(paste0("^#\\| label: ", label, "$"), src)
  stopifnot(length(at) == 1L)
  open <- max(grep("^```\\{r", src[seq_len(at)]))
  close <- at + min(grep("^```\\s*$", src[-seq_len(at)]))
  body <- src[(at + 1L):(close - 1L)]
  body <- body[!grepl("^#\\|", body)]

  before <- knitr::purl(
    text = paste(src[seq_len(open - 1L)], collapse = "\n"),
    quiet = TRUE, documentation = 0
  )
  # Parent on the global environment, as a vignette's code has. testthat's
  # helper environment hangs below the search path entries that existed when
  # it was made, so a package the vignette attaches with library() (ggplot2)
  # would not be visible from it.
  env <- new.env(parent = globalenv())
  eval(parse(text = before), envir = env)
  eval(parse(text = body), envir = env)
  env
}
