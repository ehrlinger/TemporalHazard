test_that("render_sim reports per-chunk success and failure", {
  job <- list(calls = list(
    good = quote(x <- 1 + 1),
    bad  = quote(stop("boom")),
    uses = quote(y <- x * 2)
  ))
  got <- render_sim(job)
  expect_false(got$ok)
  expect_equal(unname(got$results[["good"]]), "ok")
  expect_match(got$results[["bad"]], "^ERROR: boom")
  expect_equal(got$env$y, 4)
})

test_that("render_sim binds supplied data and nothing else", {
  job <- list(calls = list(a = quote(n <- nrow(D)), b = quote(z <- MISSING_VAR)))
  got <- render_sim(job, data = list(D = data.frame(x = 1:3)))
  expect_equal(got$env$n, 3L)
  expect_match(got$results[["b"]], "object 'MISSING_VAR' not found")
})

test_that("render_sim does not report ok when nothing was evaluated", {
  # all(character(0) == "ok") is TRUE -- a job with zero emitted calls must
  # not be indistinguishable from one where every call succeeded.
  got <- render_sim(list(calls = list()))
  expect_false(got$ok)
})

test_that("a minimal translated job renders end to end", {
  skip_on_cran()
  set.seed(1)
  AVCS <- data.frame(INT_DEAD = stats::rexp(200, 0.2),
                     DEAD = rep(c(1, 0), length.out = 200))
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=AVCS CONDITION=14;",
    "EVENT DEAD; TIME INT_DEAD;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  got <- render_sim(job, data = list(AVCS = AVCS))
  expect_true(got$ok, info = paste(names(got$results), got$results,
                                   collapse = " | "))
  expect_s3_class(got$env$fit, "hazard")
  expect_true(isTRUE(got$env$fit$fit$converged))
})

test_that("a censoring job's status chunk renders", {
  skip_on_cran()
  set.seed(2)
  n <- 120
  AVCS <- data.frame(INT_DEAD = stats::rexp(n, 0.2),
                     DEAD = rep(c(1, 0, 0), length.out = n))
  AVCS$C3FLAG <- as.numeric(seq_len(n) %% 5 == 0)
  AVCS$ICTIME <- ifelse(AVCS$C3FLAG > 0, AVCS$INT_DEAD * 0.5, NA)
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=AVCS CONDITION=14;",
    "EVENT DEAD; TIME INT_DEAD; ICENSOR C3FLAG = ICTIME;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  got <- render_sim(job, data = list(AVCS = AVCS))
  expect_true(got$ok, info = paste(names(got$results), got$results,
                                   collapse = " | "))
  expect_true(any(got$env$.hzr_status == 2))
})
