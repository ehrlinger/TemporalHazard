test_that("a job with no statements seen is rejected, not returned empty", {
  # Three shapes of absent, ascending danger: NULL, empty container, and a
  # hollow object -- right shape, empty inside -- which defeats is.null() and
  # a length check both. This is the hollow case.
  j <- .hzr_sas_job(
    source = list(path = "x.sas", checksum = "abc"),
    calls = list(), grid = NULL, inhaz = NULL, outhaz = NULL,
    untranslated = .hzr_untranslated_frame(),
    coverage = list(tokens_seen = 0L, tokens_mapped = 0L)
  )
  expect_error(.hzr_validate_sas_job(j), "no SAS statements")
})

test_that("a populated job validates and prints its coverage", {
  j <- .hzr_sas_job(
    source = list(path = "x.sas", checksum = "abc"),
    calls = list(fit = quote(hazard(time = t, status = s))),
    grid = NULL, inhaz = NULL, outhaz = "LIB.D",
    untranslated = .hzr_untranslated_frame(),
    coverage = list(tokens_seen = 4L, tokens_mapped = 4L)
  )
  expect_silent(.hzr_validate_sas_job(j))
  expect_output(print(j), "4/4")
})

test_that("the untranslated frame keeps construct and line, not a count", {
  f <- .hzr_untranslated_frame(line = 12L, construct = "RESTRICT",
                               reason = "no R equivalent")
  expect_equal(f$construct, "RESTRICT")
  expect_equal(f$line, 12L)
})

test_that("the hollow-object guard rejects every malformed tokens_seen", {
  mk <- function(seen) {
    j <- .hzr_sas_job(
      source = list(path = "x.sas", checksum = "abc"),
      calls = list(fit = quote(hazard(time = t, status = s))),
      grid = NULL, inhaz = NULL, outhaz = NULL,
      untranslated = .hzr_untranslated_frame(),
      coverage = list(tokens_seen = seen, tokens_mapped = 0L)
    )
    j
  }
  for (bad in list(0L, NA, NA_integer_, NA_real_, "4", c(1L, 2L), Inf)) {
    expect_error(.hzr_validate_sas_job(mk(bad)), "no SAS statements",
                 info = paste("tokens_seen =", deparse(bad)))
  }
  # Test NULL separately since it is dropped by list(...)
  j_null <- .hzr_sas_job(
    source = list(path = "x.sas", checksum = "abc"),
    calls = list(fit = quote(hazard(time = t, status = s))),
    grid = NULL, inhaz = NULL, outhaz = NULL,
    untranslated = .hzr_untranslated_frame(),
    coverage = list(tokens_mapped = 0L)
  )
  expect_error(.hzr_validate_sas_job(j_null), "no SAS statements",
               info = "tokens_seen = NULL")
  # And a well-formed one still passes.
  expect_silent(.hzr_validate_sas_job(mk(4L)))
})

test_that("print.hzr_sas_job formats coverage with mapped/seen tokens", {
  j <- .hzr_sas_job(
    source = list(path = "x.sas", checksum = "abc"),
    calls = list(fit = quote(hazard(time = t, status = s))),
    grid = NULL, inhaz = NULL, outhaz = "LIB.D",
    untranslated = .hzr_untranslated_frame(),
    coverage = list(tokens_seen = 7L, tokens_mapped = 6L)
  )
  expect_output(print(j), "6/7")
})
