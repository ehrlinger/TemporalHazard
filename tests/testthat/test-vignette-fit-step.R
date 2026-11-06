# The clinical-analysis-walkthrough vignette loads a stored stepwise fit
# instead of refitting it (CRAN check time). This refits the vignette's own
# code and compares, so a fit the code no longer produces fails here rather
# than shipping in a vignette. Regenerate with data-raw/vignette-fit-step.R.
skip_on_cran()

test_that("the vignette's stored stepwise fit matches a refit", {
  skip_if_not_installed("knitr")
  name <- "clinical-analysis-walkthrough.qmd"
  qmd <- c(
    system.file("doc", name, package = "TemporalHazard"),
    test_path("..", "..", "vignettes", name)
  )
  qmd <- qmd[nzchar(qmd) & file.exists(qmd)]
  expect_gt(length(qmd), 0L)
  stored_path <- system.file("extdata", "vignette-fit-step.rds",
                             package = "TemporalHazard")
  expect_true(nzchar(stored_path))
  stored <- readRDS(stored_path)

  # The vignette's chunks print (str() of the data and so on); keep that out
  # of the test log.
  utils::capture.output(refit <- suppressWarnings(
    run_vignette_through(qmd[1], "stepwise-forward")$fit_step
  ))

  # The selection the vignette's prose describes, step for step.
  cols <- c("step_num", "action", "variable", "phase")
  expect_gt(nrow(refit$steps), 0L)
  expect_identical(refit$steps[, cols], stored$steps[, cols])
  # Every column of the steps table, since the vignette prints it: AIC,
  # statistics, logLik and counts as well as the identifiers.
  expect_identical(names(refit$steps), names(stored$steps))
  expect_equal(refit$steps, stored$steps, tolerance = 1e-4)
  # What print(fit_step) shows, and what the AIC comparison chunk reads.
  expect_identical(utils::capture.output(print(refit)),
                   utils::capture.output(print(stored)))
  expect_identical(refit$fit$fixed_mask, stored$fit$fixed_mask)
  expect_equal(refit$fit$objective, stored$fit$objective, tolerance = 1e-6)
  expect_equal(unname(refit$fit$theta), unname(stored$fit$theta),
               tolerance = 1e-4)
})
