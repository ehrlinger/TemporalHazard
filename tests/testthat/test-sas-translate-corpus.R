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
