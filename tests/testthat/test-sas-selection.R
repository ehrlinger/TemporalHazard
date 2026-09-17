test_that("SELECTION FORWARD maps to direction = both, not forward", {
  # The obvious mapping is wrong and fails silently on every stepwise job.
  got <- .hzr_selection_spec(c("FORWARD", "SLENTRY=0.05", "SLSTAY=0.10"))
  expect_equal(got$direction, "both")
  expect_equal(got$slentry, 0.05)
  expect_equal(got$slstay, 0.10)
})

test_that("every forward/stepwise spelling reaches the same direction", {
  for (kw in c("FORWARD", "FW", "SW", "SELECT", "STEPWISE")) {
    expect_equal(.hzr_selection_spec(kw)$direction, "both",
                 info = paste("spelling:", kw))
  }
})

test_that("BACKWARD and ONEWAY are distinct from stepwise", {
  expect_equal(.hzr_selection_spec("BACKWARD")$direction, "backward")
  expect_equal(.hzr_selection_spec("BW")$direction, "backward")
  # ONEWAY/NOSTEPWISE is still a screen, forward only: the SELECTION statement
  # itself sets sw = 1 (hazard_y.y setopt(33); stpwprc.c), and NOSTEPWISE only
  # sets nosw, which caps each variable at one move (#342 review). Reading it
  # as "no screen" put every candidate into a plain fit.
  for (kw in c("NOSTEPWISE", "NOSW")) {
    got <- .hzr_selection_spec(kw)
    expect_true(got$stepwise, info = kw)
    expect_equal(got$direction, "forward", info = kw)
  }
})

test_that("a SELECTION block emits a screen, not a refusal (#160)", {
  # This test used to assert a stop() chunk. The refusal existed because the
  # refit path required a formula-interface base fit, so every candidate
  # refit errored and the screen reported zero steps -- indistinguishable
  # from "nothing met slentry" (#159). .hzr_refit_blocker() is phase-aware
  # now, so a vector-interface multiphase fit refits, and #160 translates
  # the statement instead. What must NOT come back is a screen that cannot
  # screen, which the executed tests below pin.
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=14;",
    "EVENT D; TIME T; EARLY X1=1, X2=1, X3=1;",
    "PARMS MUE=1 THALF=1 NU=1;",
    "SELECTION FORWARD SLENTRY=0.05 SLSTAY=0.1; );"
  ))
  got <- .hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]])
  expect_identical(got$call[[1L]], as.name("hazard"))
  expect_identical(got$stepwise_call[[1L]], as.name("hzr_stepwise"))
  expect_equal(got$stepwise_call[["slentry"]], 0.05)
  expect_equal(got$stepwise_call[["slstay"]], 0.1)
  # FORWARD is SAS option 21, which is two-way (see .hzr_selection_spec()).
  expect_equal(got$stepwise_call[["direction"]], "both")
  expect_false(any(grepl("SELECTION", got$untranslated$construct, fixed = TRUE)))
})

test_that("a non-numeric SLENTRY value is untranslated, not coerced to NA silently", {
  got <- .hzr_selection_spec(c("STEPWISE", "SLENTRY=oops"))
  expect_true(any(grepl("SLENTRY", got$untranslated$construct)))
  # The value is discarded and PROC HAZARD's own default stands in its
  # place, so the emitted call never inherits hzr_stepwise()'s different
  # default (stpwprc.c: SLE 0.3).
  expect_equal(got$slentry, 0.3)
})

test_that("an unknown SELECTION option is recorded, not dropped", {
  got <- .hzr_selection_spec("BOGUS")
  expect_true(any(grepl("BOGUS", got$untranslated$construct)))
})

test_that("a SELECTION statement with no direction still enables stepwise", {
  # stepwisestmt : STEPWISE stepwiseopts { setopt(33); } -- stepwiseopts may be
  # empty, so the statement itself turns stepwise on.
  got <- .hzr_selection_spec(c("SLENTRY=0.05", "SLSTAY=0.1"))
  expect_true(got$stepwise)
  expect_equal(got$direction, "both")
  expect_equal(got$slentry, 0.05)
  expect_equal(got$slstay, 0.1)
})

test_that("a bare SELECTION with no options at all enables stepwise", {
  got <- .hzr_selection_spec(character(0))
  expect_true(got$stepwise)
  expect_equal(got$direction, "both")
})

test_that("NOSTEPWISE keeps its entry threshold, because it still screens", {
  got <- .hzr_selection_spec(c("NOSTEPWISE", "SLENTRY=0.05"))
  expect_true(got$stepwise)
  expect_equal(got$slentry, 0.05)
  expect_false(any(grepl("SLENTRY", got$untranslated$construct)))
})

test_that("a NOSTEPWISE job screens forward, with candidates out of the base (#342, #160)", {
  # Under SELECTION a bare phase variable starts OUT of the model
  # (setstat.c), so emitting ~AGE + X + Z as a plain fit was a wrong model.
  # #342 refused the job; #160 translates it: NOSTEPWISE is a forward-only
  # screen, because the SELECTION statement sets sw = 1 and nosw only caps
  # moves. The assertion that matters is unchanged -- the candidates are not
  # forced into the base model.
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=D CONDITION=14; EVENT DEAD; TIME TT;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; SELECTION NOSW;",
    "EARLY AGE, X/I, Z/S; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  cl <- job$calls$fit[[3L]]
  expect_identical(cl[[1L]], as.name("hzr_stepwise"))
  expect_equal(cl[["direction"]], "forward")
  # AGE is a bare candidate: offered as scope, absent from the base.
  expect_equal(deparse(as.list(cl[["scope"]])$phase_1), "~AGE")
  base_formula <- job$calls$fit_base[[3L]][["phases"]][[2L]][["formula"]]
  expect_equal(deparse(base_formula), "~X + Z")
  expect_equal(eval(cl[["force_in"]]), "X")
})

test_that("a directionless SELECTION block emits a screen too (#160)", {
  # A bare SELECTION (stepwisestmt with empty stepwiseopts) still enables
  # stepwise, so it translates like any other; it used to be refused.
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=14;",
    "EVENT D; TIME T; EARLY X1=0.5, X2=0.25;",
    "PARMS MUE=1 THALF=1 NU=1;",
    "SELECTION SLENTRY=0.05 SLSTAY=0.1; );"
  ))
  got <- .hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]])
  expect_identical(got$stepwise_call[[1L]], as.name("hzr_stepwise"))
  expect_equal(got$stepwise_call[["slentry"]], 0.05)
  expect_equal(got$stepwise_call[["slstay"]], 0.1)
})

test_that("a SELECTION job emits a screen that is not a no-op (#160)", {
  # The refusal this replaces existed to prevent a screen that runs and
  # selects nothing, which reads exactly like an honest "nothing met
  # slentry" (#159). So assert the emitted call carries what a screen needs:
  # a scope with candidates in it, and the thresholds the job asked for.
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=14;",
    "EVENT DEAD; TIME T; EARLY X1, X2, X3;",
    "SELECTION STEPWISE SLENTRY=0.05 SLSTAY=0.1;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  expect_false(any(grepl("SELECTION", job$untranslated$construct, fixed = TRUE)))
  cl <- job$calls$fit[[3L]]
  expect_identical(cl[[1L]], as.name("hzr_stepwise"))
  expect_equal(deparse(as.list(cl[["scope"]])$phase_1), "~X1 + X2 + X3")
  # The candidates are NOT also baked into the base model, which would
  # invert the statement's meaning: the base is intercept-only here.
  expect_null(job$calls$fit_base[[3L]][["phases"]][[2L]][["formula"]])
})

test_that("a job with no SELECTION statement is unaffected", {
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=14;",
    "EVENT DEAD; TIME T; EARLY X1, X2, X3;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  expect_false(any(grepl("SELECTION", job$untranslated$construct, fixed = TRUE)))
  expect_equal(job$calls$fit[[3L]][[1L]], as.name("hazard"))
})

# --- #160: SELECTION translates to hzr_stepwise() ---------------------------
#
# The acceptance condition from the issue is that a step is actually TAKEN:
# a screen that runs and selects nothing is indistinguishable from an honest
# "nothing met slentry", which is the defect the refusal existed to avoid.

.sel_job <- function(stmts, env = parent.frame()) {
  f <- withr::local_tempfile(fileext = ".sas", .local_envir = env)
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=D CONDITION=14; EVENT DEAD; TIME TT;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005;", stmts, ");"
  ), f)
  suppressWarnings(hzr_translate_sas(f))
}

.sel_data <- function(n = 400, seed = 31) {
  set.seed(seed)
  d <- data.frame(STRONG = stats::rnorm(n), NOISE = stats::rnorm(n),
                  MAL = stats::rbinom(n, 1, 0.4), AGE = stats::rnorm(n))
  d$TT <- stats::rexp(n, rate = 0.15 * exp(0.9 * d$STRONG))
  d$DEAD <- stats::rbinom(n, 1, 0.85)
  d
}

test_that("a SELECTION job translates to a screen that TAKES a step (#160)", {
  skip_on_cran()
  job <- .sel_job("SELECTION SLE=0.2 SLS=0.1; EARLY STRONG, NOISE, MAL/I;")
  expect_false(any(job$untranslated$construct == "SELECTION"))
  # Two chunks: the base fit, then the screen, so nothing is self-referential.
  expect_true(all(c("fit_base", "fit") %in% names(job$calls)))
  cl <- job$calls$fit[[3L]]
  expect_identical(cl[[1L]], as.name("hzr_stepwise"))
  expect_identical(cl[[2L]], as.name("fit_base"))
  # Candidates are WITHHELD from the phase formula and offered as scope,
  # keyed positionally, because the emitted phases list is unnamed and the
  # fit therefore carries phase_1/phase_2 (keying on "early" errors).
  base_call <- job$calls$fit_base[[3L]]
  ph1 <- base_call[["phases"]][[2L]]
  expect_equal(deparse(ph1[["formula"]]), "~MAL")
  expect_equal(names(as.list(cl[["scope"]])[-1L]), "phase_1")
  expect_equal(deparse(as.list(cl[["scope"]])$phase_1), "~STRONG + NOISE")
  expect_equal(eval(cl[["force_in"]]), "MAL")
  expect_equal(cl[["slentry"]], 0.2)
  expect_equal(cl[["slstay"]], 0.1)
  expect_equal(cl[["criterion"]], "score")
  expect_equal(cl[["direction"]], "both")
  expect_equal(cl[["max_move"]], 1)

  res <- suppressWarnings(render_sim(job, list(D = .sel_data())))
  expect_true(res$ok, info = paste(res$results, collapse = "; "))
  sw <- res$env$fit
  expect_s3_class(sw, "hzr_stepwise")
  steps <- as.data.frame(sw)
  expect_gt(nrow(steps), 0L)
  entered <- steps$variable[steps$action == "enter"]
  expect_true("STRONG" %in% entered, info = paste(entered, collapse = ","))
  expect_true(any(grepl("STRONG", names(stats::coef(sw)))))
  # MAL was /I: forced in, never dropped.
  expect_false("MAL" %in% steps$variable[steps$action == "drop"])
})

test_that("the divergence callout sits above the stepwise chunk (#160)", {
  job <- .sel_job("SELECTION SLE=0.2; EARLY STRONG, NOISE;")
  doc <- TemporalHazard:::.hzr_render_qmd(job)
  hit <- grep("may select a different", doc)
  expect_length(hit, 1L)
  # Above the chunk, not after it: the reader meets it before the code.
  chunk <- grep("^#\\| label: fit$", doc)
  expect_length(chunk, 1L)
  expect_lt(hit, chunk)
  # It must say WHAT differs and WHY, and cite the recorded gap.
  txt <- paste(doc, collapse = " ")
  expect_match(txt, "approximate variances")
  expect_match(txt, "full Hessian")
  expect_match(txt, "force_in")
  expect_match(txt, "hm[.]death[.]AVC")
  # It must not claim reproduction.
  expect_no_match(txt, "reproduces PROC HAZARD")

  # Absent for a job with no SELECTION statement.
  plain <- .sel_job("EARLY STRONG, NOISE;")
  expect_length(grep("may select a different",
                     TemporalHazard:::.hzr_render_qmd(plain)), 0L)
})

test_that("a BACKWARD job emits no scope and drops from a full base (#160, #343)", {
  skip_on_cran()
  job <- .sel_job("SELECTION BACKWARD; EARLY STRONG, NOISE, MAL/I;")
  cl <- job$calls$fit[[3L]]
  expect_identical(cl[[1L]], as.name("hzr_stepwise"))
  expect_null(cl[["scope"]])
  expect_equal(cl[["direction"]], "backward")
  # SAS's SLS default under BACKWARD is 0.05, not the 0.2 of a stepwise run.
  expect_equal(cl[["slstay"]], 0.05)
  # Under BACKWARD a bare variable starts IN the model (setstat.c).
  base_call <- job$calls$fit_base[[3L]]
  expect_equal(deparse(base_call[["phases"]][[2L]][["formula"]]),
               "~STRONG + NOISE + MAL")
  res <- suppressWarnings(render_sim(job, list(D = .sel_data())))
  expect_true(res$ok, info = paste(res$results, collapse = "; "))
  steps <- as.data.frame(res$env$fit)
  expect_true(any(steps$action == "drop"), info = paste(steps$action, collapse = ","))
  expect_false("MAL" %in% steps$variable[steps$action == "drop"])
})

test_that("SELECTION constructs with no faithful translation are refused (#160)", {
  refused <- function(stmts, what) {
    job <- .sel_job(stmts)
    expect_identical(job$calls$fit[[3L]][[1L]], as.name("stop"), info = what)
    expect_error(eval(job$calls$fit), what, info = what)
    expect_true(any(grepl(what, job$untranslated$reason)), info = what)
  }
  refused("SELECTION FAST; EARLY STRONG, NOISE;", "FAST")
  refused("SELECTION MAXVARS=2; EARLY STRONG, NOISE;", "MAXVARS")
  refused("SELECTION ROBUST; EARLY STRONG, NOISE;", "ROBUST")
  refused("SELECTION SEMIROBUST; EARLY STRONG, NOISE;", "SEMIROBUST")
  refused("SELECTION; EARLY STRONG/MOVE=2, NOISE;", "MOVE=")
  refused("SELECTION; EARLY STRONG/ORDER=1, NOISE;", "ORDER=")
  # /I in one phase and movable in another: force_in has no phase, so it
  # would be pinned in both, which is a wrong model, not a path difference.
  refused("SELECTION; EARLY MAL/I, STRONG; CONSTANT MAL, NOISE;", "MAL")
  # RESTRICT constrains which variables the screen may select (hazrd4.c
  # rsttbl). A screen that ignores it selects by a different rule than the
  # job asked for, so it is refused rather than recorded.
  refused("SELECTION; EARLY STRONG, NOISE; RESTRICT STRONG;", "RESTRICT")
})

test_that("printing options are recorded, not refused (#160)", {
  job <- .sel_job("SELECTION NOPRINTS NOPRINTQ; EARLY STRONG, NOISE;")
  expect_identical(job$calls$fit[[3L]][[1L]], as.name("hzr_stepwise"))
  expect_true(any(grepl("NOPRINTS", job$untranslated$construct)))
})

test_that("the emitted screen reports candidates it could not score (#160)", {
  skip_on_cran()
  job <- .sel_job("SELECTION SLE=0.2; EARLY STRONG, NOISE;")
  # The uncomputable-score tally is surfaced in the document, so a screen
  # that stopped without testing anything is not read as "nothing qualified".
  doc <- TemporalHazard:::.hzr_render_qmd(job)
  expect_true(any(grepl("uncomputable", doc, ignore.case = TRUE)))
})
