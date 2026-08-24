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
  corpus <- Sys.getenv(
    "HZR_SAS_CORPUS",
    Sys.getenv("HAZARD_REPO", "~/Documents/GitHub/hazard")
  )
  skip_if_not(dir.exists(path.expand(corpus)), "hazard checkout not available")
  fs <- list.files(path.expand(corpus),
                   pattern = "[.]sas$", full.names = TRUE, recursive = TRUE)
  skip_if(length(fs) == 0L, "no .sas corpus available")

  n_trans <- 0L        # .sas files that translated at all
  n_eligible <- 0L     # distinct jobs synthetic data can drive end to end
  n_rendered <- 0L     # of those, jobs whose every chunk ran
  n_partial <- 0L      # distinct jobs checked only up to their fit chunks
  failures <- character(0)
  ineligible <- character(0)
  seen <- new.env(parent = emptyenv())

  for (f in fs) {
    job <- tryCatch(suppressWarnings(hzr_translate_sas(f)), error = function(e) NULL)
    if (is.null(job)) next
    n_trans <- n_trans + 1L

    # The corpus holds the same jobs under examples/, dist/examples/ and
    # tests/, so 57 translations are 22 distinct documents. Fitting each
    # one once keeps the file inside its time budget; the duplicates would
    # exercise byte-identical calls.
    key <- paste(vapply(job$calls, function(x) paste(deparse(x), collapse = " "), ""),
                 collapse = " ;; ")
    if (!is.null(seen[[key]])) next
    assign(key, TRUE, envir = seen)

    shape <- sas_job_shape(job)
    # A meaningless fit is fine here -- the assertion is that the chunks run
    # -- but a nonsense fit warns freely (no convergence, no standard
    # errors). Errors are the signal; warnings are noise, so drop them.
    sim <- suppressWarnings(render_sim(job, sas_synth_data(job)))

    # Chunks the translator emits as a deliberate stop() must actually stop:
    # a SELECTION or LCENSOR + ICENSOR refusal that quietly succeeded would
    # be a fit over a model this package refuses to specify.
    for (nm in shape$refusals) {
      expect_match(sim$results[[nm]], "^ERROR: ", info = paste(basename(f), nm))
    }

    if (!shape$eligible) {
      # No local fit to bind (an unresolved INHAZ= job is the usual case),
      # so predict() chunks cannot run here. Everything else still must:
      # excluding the whole job would hide a defect in its grid or guard
      # chunks, which is what the previous version of this exclusion did.
      n_partial <- n_partial + 1L
      ineligible <- c(ineligible, basename(f))
      other <- setdiff(names(sim$results), c(shape$preds, shape$refusals))
      bad <- other[sim$results[other] != "ok"]
      if (length(bad)) {
        failures <- c(failures, paste0(basename(f), " [", bad[1], "] ", sim$results[[bad[1]]]))
      }
      expect_true(length(bad) == 0L,
                  info = paste(basename(f), paste(bad, collapse = ", ")))
      next
    }

    n_eligible <- n_eligible + 1L
    # Assert the fit was BOUND, not merely evaluated. A chunk that calls
    # hazard() without assigning it returns "ok" from render_sim and leaves
    # every downstream predict() chunk referencing an object that is not
    # there -- the defect (#151, R/translate-sas.R) this test exists for.
    for (nm in shape$hazard_chunks) {
      expect_true(exists(nm, envir = sim$env, inherits = FALSE),
                  info = paste(basename(f), nm, "not bound"))
      if (exists(nm, envir = sim$env, inherits = FALSE)) {
        expect_s3_class(get(nm, envir = sim$env, inherits = FALSE), "hazard")
      }
    }
    bad <- setdiff(names(sim$results)[sim$results != "ok"], shape$refusals)
    if (length(bad)) {
      failures <- c(failures, paste0(basename(f), " [", bad[1], "] ", sim$results[[bad[1]]]))
    } else {
      n_rendered <- n_rendered + 1L
    }
  }

  if (length(failures)) {
    message("Corpus jobs that translate but do not render:\n",
            paste(failures, collapse = "\n"))
  }
  if (length(ineligible)) {
    message(length(ineligible), " corpus job(s) checked only up to their fit chunks: ",
            paste(ineligible, collapse = "; "))
  }

  # Measured 2026-08-22 against ~/Documents/GitHub/hazard: 57 translations,
  # 22 distinct documents, 11 of them eligible and all 11 rendering. The
  # other 11 are PROC HAZPRED-only jobs -- an unresolved INHAZ= or a
  # refused SELECTION leaves no local fit to bind -- and are checked up to
  # their fit chunks only. That is the honest measurement: half the corpus
  # gets its predict() chunks exercised, half does not, and the eligible
  # count is asserted so a shrinking eligible set fails rather than
  # quietly reducing what is tested. Floors, not equalities, so a corpus
  # that grows does not fail.
  expect_gt(n_trans, 20L)
  expect_gte(n_eligible, 11L)
  expect_gte(n_partial, 11L)
  expect_equal(n_rendered, n_eligible)
})
