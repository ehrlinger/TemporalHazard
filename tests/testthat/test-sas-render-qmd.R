test_that("calls render as R chunks", {
  j <- .hzr_sas_job(
    source = list(path = "hz.death.AVC.sas", checksum = "abc"),
    calls = list(fit = quote(hazard(time = Time, status = D))),
    grid = NULL, inhaz = NULL, outhaz = "EX.HZD",
    untranslated = .hzr_untranslated_frame(),
    coverage = list(tokens_seen = 6L, tokens_mapped = 6L)
  )
  out <- .hzr_render_qmd(j)
  expect_true(any(grepl("^```\\{r\\}", out)))
  expect_true(any(grepl("hazard(time = Time, status = D)", out, fixed = TRUE)))
})

test_that("an untranslated construct becomes a visible callout", {
  j <- .hzr_sas_job(
    source = list(path = "x.sas", checksum = "abc"),
    calls = list(fit = quote(hazard(time = Time, status = D))),
    grid = NULL, inhaz = NULL, outhaz = NULL,
    untranslated = .hzr_untranslated_frame(12L, "STEEPEST",
                                           "no R equivalent (see #145)"),
    coverage = list(tokens_seen = 7L, tokens_mapped = 6L)
  )
  out <- .hzr_render_qmd(j)
  expect_true(any(grepl("UNTRANSLATED", out)))
  expect_true(any(grepl("STEEPEST", out)))
})

test_that("an unresolved INHAZ renders a stop(), not a comment", {
  # The document must FAIL to render rather than produce a report over a model
  # it never loaded.
  j <- .hzr_sas_job(
    source = list(path = "hp.sas", checksum = "abc"),
    calls = list(pred = quote(predict(fit, newdata = PREDICT))),
    grid = NULL, inhaz = "CABGKUL.HMDEADP", outhaz = NULL,
    untranslated = .hzr_untranslated_frame(),
    coverage = list(tokens_seen = 4L, tokens_mapped = 4L)
  )
  out <- .hzr_render_qmd(j)
  expect_true(any(grepl("stop(", out, fixed = TRUE)))
  expect_true(any(grepl("CABGKUL.HMDEADP", out, fixed = TRUE)))
})

test_that("a resolved INHAZ does not emit a stop()", {
  j <- .hzr_sas_job(
    source = list(path = "hp.sas", checksum = "abc"),
    calls = list(pred = quote(predict(fit, newdata = PREDICT))),
    grid = NULL, inhaz = "CABGKUL.HMDEADP", outhaz = NULL,
    untranslated = .hzr_untranslated_frame(),
    coverage = list(tokens_seen = 4L, tokens_mapped = 4L)
  )
  j$inhaz_resolved <- TRUE
  out <- .hzr_render_qmd(j)
  expect_false(any(grepl("stop(", out, fixed = TRUE)))
})

test_that(".hzr_render_qmd returns a character vector, not a printed side effect", {
  j <- .hzr_sas_job(
    source = list(path = "x.sas", checksum = "abc"),
    calls = list(fit = quote(hazard(time = Time, status = D))),
    grid = NULL, inhaz = NULL, outhaz = NULL,
    untranslated = .hzr_untranslated_frame(),
    coverage = list(tokens_seen = 1L, tokens_mapped = 1L)
  )
  out <- .hzr_render_qmd(j)
  expect_type(out, "character")
  expect_true(length(out) > 1L)
})
