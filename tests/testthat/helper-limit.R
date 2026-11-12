# What the limit layer did to a fixture (#735).
#
# Many fixtures outside the test-limit-* files are restarted, or kept, by the
# limit layer without the test saying so. Pinning the outcome and the trigger
# source makes a change in what the layer does to that fixture fail loudly,
# instead of shifting the fit the rest of the test reads.
expect_limit <- function(fit, outcome, source) {
  lim <- fit$fit$limit
  testthat::expect_true(is.list(lim), label = "fit$fit$limit is a list")
  testthat::expect_identical(lim$outcome, outcome)
  testthat::expect_identical(
    unique(vapply(lim$trigger, `[[`, "", "source")), source
  )
}
