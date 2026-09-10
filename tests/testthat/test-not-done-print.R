# The printed block is checked against the record by parsing the printed text
# back into entries (#242). The check never calls the printer's formatter, so
# a printer that drops or rewrites an entry cannot also rewrite what it is
# checked against. The last three tests make the printer wrong on purpose and
# assert that the check says so: a check that cannot fail proves nothing.

clean_fit <- function() {
  set.seed(7)
  dat <- data.frame(t = stats::rweibull(80, 1.4, 3), d = rep(1L, 80))
  hazard(survival::Surv(t, d) ~ 1, data = dat, dist = "weibull",
         fit = TRUE, theta = c(mu = 1, nu = 1))
}

unfitted <- function() {
  set.seed(7)
  hazard(time = stats::rweibull(80, 1.4, 3), status = rep(1L, 80),
         dist = "weibull", fit = FALSE)
}

printed <- function(x) {
  list(print = capture.output(print(x)),
       summary = capture.output(print(summary(x))))
}

test_that("print() and summary() say 'none' for a clean fit", {
  fit <- clean_fit()
  expect_identical(fit$degraded, character(0))     # guard
  for (out in printed(fit)) {
    expect_true("  Not done in this run: none" %in% out)
    expect_identical(.hzr_check_not_done_output(out, fit), character(0))
  }
})

test_that("print() and summary() show every entry of a degraded record", {
  fit0 <- unfitted()
  expect_length(fit0$degraded, 3L)                 # guard: something to print
  for (out in printed(fit0)) {
    expect_identical(.hzr_check_not_done_output(out, fit0), character(0))
  }
})

test_that("an object from before the record existed says so, not 'none'", {
  fit <- clean_fit()
  fit$degraded <- NULL
  fit$degraded_causes <- NULL
  for (out in printed(fit)) {
    expect_true(paste0("  Not done in this run: not recorded ",
                       "(object built before this record existed)") %in% out)
    expect_identical(.hzr_check_not_done_output(out, fit), character(0))
  }
})

test_that("the check fails when the printer drops an entry", {
  fit0 <- unfitted()
  local_mocked_bindings(.hzr_format_not_done = function(degraded, causes) {
    keep <- degraded[-1]
    c("  Not done in this run:", paste0("    ", keep, ": ", unname(causes[keep])))
  })
  out <- capture.output(print(fit0))
  expect_true("entry missing: fitting" %in% .hzr_check_not_done_output(out, fit0))
})

test_that("the check fails when the printer always says 'none'", {
  fit0 <- unfitted()
  local_mocked_bindings(.hzr_format_not_done = function(degraded, causes) {
    "  Not done in this run: none"
  })
  out <- capture.output(print(summary(fit0)))
  expect_true("'none' printed over a non-empty record" %in%
                .hzr_check_not_done_output(out, fit0))
})

test_that("the check fails when the printer leaves the block out", {
  fit <- clean_fit()
  local_mocked_bindings(.hzr_format_not_done = function(degraded, causes) {
    character(0)
  })
  expect_identical(.hzr_check_not_done_output(capture.output(print(fit)), fit),
                   "heading missing")
})
