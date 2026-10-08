# The unfitted-prediction warning is switched off only inside the test files
# that predict from unfitted models on purpose, each for its own scope
# (#398). Nothing switches it off for the whole suite: a helper that did so
# would also silence every developer session that calls load_all(), which
# sources the test helpers, and would hide a new unfitted-prediction path in
# any other file.

test_that("no helper switches the unfitted-prediction warning off suite-wide (#398)", {
  expect_null(getOption("TemporalHazard.warn_unfitted_prediction"))
  helpers <- list.files(test_path(), pattern = "^helper.*[.][Rr]$", full.names = TRUE)
  expect_gt(length(helpers), 0L)
  uses <- vapply(helpers, function(f) {
    any(grepl("warn_unfitted_prediction", readLines(f, warn = FALSE), fixed = TRUE))
  }, logical(1))
  expect_false(any(uses), info = paste(basename(helpers[uses]), collapse = ", "))
})
