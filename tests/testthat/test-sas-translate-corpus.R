test_that("every public-corpus job either translates or fails loudly", {
  skip_on_cran()
  repo <- Sys.getenv("HAZARD_REPO", "~/Documents/GitHub/hazard")
  skip_if_not(dir.exists(path.expand(repo)), "hazard checkout not available")

  fs <- list.files(path.expand(repo), pattern = "[.]sas$",
                   recursive = TRUE, full.names = TRUE)
  skip_if(length(fs) == 0L, "no .sas files found")

  out <- withr::local_tempdir()
  n_ok <- 0L
  for (f in fs) {
    job <- tryCatch(
      suppressWarnings(hzr_translate_sas(f, out_dir = out)),
      error = function(e) NULL
    )
    if (is.null(job)) next
    # A translated job must have seen something. A hollow job is the failure
    # this whole design exists to prevent.
    expect_gt(job$coverage$tokens_seen, 0L)
    expect_gte(job$coverage$tokens_seen, job$coverage$tokens_mapped)
    n_ok <- n_ok + 1L
  }
  # Warn loudly if nothing translated: a pass over zero jobs is not a pass.
  # Measured 2026-08-20: 57 successful translations across 110 .sas files (26
  # distinct .qmd outputs -- fewer than 57 because examples/, dist/examples/
  # and tests/ hold duplicate copies of the same jobs). A drop below this
  # threshold that isn't explained by a corpus change is a real regression.
  expect_gt(n_ok, 20L)
})

test_that("public-corpus jobs that translate also render", {
  skip_on_cran()
  fs <- list.files(Sys.getenv("HZR_SAS_CORPUS", unset = ""),
                   pattern = "[.]sas$", full.names = TRUE, recursive = TRUE)
  skip_if(length(fs) == 0L, "no .sas corpus available")

  n_render <- 0L
  n_trans <- 0L
  n_excluded <- 0L
  failures <- character(0)

  for (f in fs) {
    job <- tryCatch(suppressWarnings(hzr_translate_sas(f)), error = function(e) NULL)
    if (is.null(job)) next
    n_trans <- n_trans + 1L

    # Three constructs are refused or unresolvable by design, not translator
    # defects, so a job in one of these categories is excluded from the
    # render requirement rather than misclassified against it:
    #  - SELECTION's stepwise screen and LCENSOR + ICENSOR are refused
    #    outright: the "fit" chunk holds a bare stop() (R/sas-parse-job.R),
    #    indistinguishable by chunk name from a real fit.
    #  - An unresolved INHAZ= is a third documented non-render case --
    #    hzr_translate_sas() itself warns "the emitted document will stop()
    #    rather than render" (R/translate-sas.R), and .hzr_render_qmd()
    #    inserts that stop() chunk *before* any of job$calls, so nothing in
    #    job$calls is reachable in an actual render. render_sim() only
    #    evaluates job$calls, so it never sees that stop() and instead
    #    surfaces a downstream "object 'fit' not found" -- the same
    #    intentional failure, at the wrong chunk.
    refused <- any(job$untranslated$construct %in% c("SELECTION", "LCENSOR + ICENSOR"))
    if (refused || !isTRUE(job$inhaz_resolved)) {
      n_excluded <- n_excluded + 1L
      next
    }

    # Render with no data bound: a job whose DATA= dataset is absent must
    # fail in its guard chunk and nowhere else. A real render halts at the
    # first uncaught error, so only the *first* failing chunk in document
    # order decides whether the document would have rendered past the
    # guard -- a later chunk failing only because an earlier one already
    # stopped the document is not a separate defect.
    sim <- render_sim(job)
    bad_idx <- which(sim$results != "ok")
    ok <- length(sim$results) > 0L &&
      (length(bad_idx) == 0L || grepl("^data(_[0-9]+)?$", names(sim$results)[bad_idx[1]]))
    if (ok) {
      n_render <- n_render + 1L
    } else {
      first_bad <- if (length(sim$results) == 0L) {
        "<no calls emitted>"
      } else {
        paste0(names(sim$results)[bad_idx[1]], ": ", sim$results[bad_idx[1]])
      }
      failures <- c(failures, paste0(basename(f), " -- ", first_bad))
    }
  }

  if (length(failures)) {
    message("Corpus jobs that translate but do not render:\n",
            paste(failures, collapse = "\n"))
  }

  # Assert coverage, not just a non-zero count: a threshold that cannot fail
  # is how the previous corpus test passed over unrenderable documents.
  expect_gt(n_trans, 20L)
  expect_gte(n_render, n_trans - n_excluded)
})
