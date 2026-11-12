test_that("every public-corpus job either translates or fails loudly", {
  skip_on_cran()
  repo <- Sys.getenv("HAZARD_REPO", "~/Documents/GitHub/hazard")
  skip_if_not(dir.exists(path.expand(repo)), "hazard checkout not available")

  fs <- list.files(path.expand(repo), pattern = "[.]sas$",
                   recursive = TRUE, full.names = TRUE)
  # dist/ is a local `make install --prefix=./dist` tree, gitignored by the
  # hazard repo: not public corpus, and absent from a clean checkout such as
  # the one CI uses. Counting it makes a local run measure a different corpus.
  fs <- fs[!startsWith(fs, file.path(path.expand(repo), "dist", ""))]
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
  # Measured 2026-09-10 on a clean checkout of the hazard commit CI pins
  # (dad7978): 37 successful translations across 66 .sas files. examples/
  # and tests/ hold duplicate copies of the same jobs. A drop below this
  # threshold that isn't explained by a corpus change is a real regression.
  # (The 2026-08-20 figure, 57 across 110, also counted a local dist/ tree.)
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
  # Same exclusion as above: dist/ is a local install tree, not corpus.
  fs <- fs[!startsWith(fs, file.path(path.expand(corpus), "dist", ""))]
  skip_if(length(fs) == 0L, "no .sas corpus available")

  n_trans <- 0L        # .sas files that translated at all
  n_eligible <- 0L     # distinct jobs synthetic data can drive end to end
  n_rendered <- 0L     # of those, jobs whose every chunk ran
  n_partial <- 0L      # distinct jobs checked only up to their fit chunks
  failures <- character(0)
  ineligible <- character(0)
  # Keys are whole deparsed jobs; a character vector, because assign() caps
  # a name at 10000 bytes and the longest corpus job passed that (#340).
  seen <- character(0)
  # Jobs refused because PROC HAZARD itself rejects their phase statements
  # (#340). The corpus carries only single-flag options (/E, /I, /S), all
  # valid SAS, so none may be refused: this is the false-refusal check the
  # corpus can give that the unit tests cannot.
  n_parse_rejected <- 0L
  # Jobs whose fit chunk hazard() refuses on THIS synthetic data (#683), by
  # job and chunk. sas_synth_data() spaces the times evenly to 60 with no
  # early events, so a job's early phase has nothing to fit: its t_half runs
  # past the data, no restart reaches a fit off that limit, and hazard()
  # stops with class hzr_limit_refusal rather than return the drifted fit.
  # That is the test data, not the job: SAS's own hm.death.AVC fit, on the
  # study data, has THALF 0.151, inside its data. The refusal is asserted
  # for exactly these, with its message, so a new refusal, or one of these
  # starting to render, fails here. Measured at 1906f75b against the hazard
  # commit CI pins (dad7978), whose .sas files match the local checkout.
  limit_refused <- c("hm.death.AVC.deciles.sas" = "fit",
                     "hm.death.AVC.sas" = "fit_2",
                     "hzr.J24407.X32.sas" = "fit",
                     "hzr.J24407.X38.sas" = "fit")
  # A job on the edge of that limit: its early phase has nothing to fit here
  # either, and whether the optimizer goes far enough for the check to call
  # it a limit depends on the platform. On macOS fit_3 stops at log t_half
  # 7.30, where G at the last time is 3.9e-8, just above the 1.5e-8 the check
  # reads as 0, and renders; on Linux it runs to Inf and is refused (#737 CI).
  # Either is accepted, the refusal only with its message.
  limit_may_refuse <- c("hm.dead.n1surg.sas" = "fit_3")
  refused_seen <- character(0)

  for (f in fs) {
    job <- tryCatch(suppressWarnings(hzr_translate_sas(f)), error = function(e) NULL)
    if (is.null(job)) next
    n_trans <- n_trans + 1L
    src <- paste(unlist(lapply(job$calls, deparse)), collapse = " ")
    if (grepl("PROC HAZARD does not run this job", src, fixed = TRUE)) {
      n_parse_rejected <- n_parse_rejected + 1L
    }

    # The corpus holds the same jobs under examples/ and tests/, so 37
    # translations are 21 distinct documents. Fitting each one once keeps
    # the file inside its time budget; the duplicates would exercise
    # byte-identical calls.
    key <- paste(vapply(job$calls, function(x) paste(deparse(x), collapse = " "), ""),
                 collapse = " ;; ")
    if (key %in% seen) next
    seen <- c(seen, key)

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
    refused <- limit_refused[basename(f)]
    edge <- limit_may_refuse[basename(f)]
    if (is.na(refused) && !is.na(edge) &&
          !identical(sim$results[[edge]], "ok")) {
      refused <- edge
    }
    if (!is.na(refused)) {
      res <- sim$results[[refused]]
      expect_match(res, "^ERROR: hazard\\(\\): phase_1\\.t_half ran to Inf ",
                   info = basename(f))
      expect_match(res, "The drifted estimate is not returned.", fixed = TRUE,
                   info = basename(f))
      expect_identical(setdiff(names(sim$results)[sim$results != "ok"],
                               shape$refusals), unname(refused),
                       info = basename(f))
      refused_seen <- c(refused_seen, basename(f))
      next
    }
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
  #
  # Re-measured 2026-09-10 on a clean checkout of the hazard commit CI pins
  # (dad7978), dist/ excluded above: 37 translations, 21 distinct
  # documents, 10 eligible and all 10 rendering, 11 partial. The 2026-08-22
  # figures also counted a local dist/ install tree. That tree is the only
  # place hm.dthar.TGA.sas exists, the job whose `?` PARMS template now emits
  # a stop() instead of an unfitted fit; that shape is tested without the
  # corpus in test-sas-translate-fits.R.
  expect_gt(n_trans, 20L)
  expect_identical(n_parse_rejected, 0L)
  expect_gte(n_eligible, 10L)
  expect_gte(n_partial, 11L)
  # Every eligible job renders, except the refusals named above, each of
  # which was seen (#683).
  expect_true(all(names(limit_refused) %in% refused_seen))
  expect_true(all(refused_seen %in% c(names(limit_refused),
                                      names(limit_may_refuse))))
  expect_equal(n_rendered + length(refused_seen), n_eligible)
})
