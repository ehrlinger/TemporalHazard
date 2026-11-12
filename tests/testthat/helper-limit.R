# What the limit layer did to a fixture (#735).
#
# Many fixtures outside the test-limit-* files are restarted, or kept, by the
# limit layer without the test saying so. Pinning the outcome and the trigger
# source makes a change in what the layer does to that fixture fail loudly,
# instead of shifting the fit the rest of the test reads.
#
# A fixture whose verdict differs by BLAS build (whether the Hessian is usable,
# whether a grid start wins by a relative 1e-6) takes every verdict seen on
# some machine: `outcome` and `source` then list the allowed values.
expect_limit <- function(fit, outcome, source) {
  lim <- fit$fit$limit
  testthat::expect_true(is.list(lim), label = "fit$fit$limit is a list")
  src <- unique(vapply(lim$trigger, `[[`, "", "source"))
  if (length(outcome) == 1L) {
    testthat::expect_identical(lim$outcome, outcome)
  } else {
    testthat::expect_true(isTRUE(lim$outcome %in% outcome),
                          label = paste("outcome", lim$outcome))
  }
  if (length(source) == 1L) {
    testthat::expect_identical(src, source)
  } else {
    testthat::expect_true(length(src) == 1L && src %in% source,
                          label = paste("source", paste(src, collapse = ",")))
  }
}
