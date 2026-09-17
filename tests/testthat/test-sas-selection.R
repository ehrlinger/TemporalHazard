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
  # scope is always explicit, with a NULL for a phase offering nothing.
  expect_equal(names(as.list(cl[["scope"]])[-1L]), c("phase_1", "phase_2"))
  expect_equal(deparse(as.list(cl[["scope"]])$phase_1), "~STRONG + NOISE")
  expect_null(as.list(cl[["scope"]])$phase_2)
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

test_that("every mapped SELECTION option reaches the emitted call (#160)", {
  # A line-by-line check of the approved mapping against what is emitted.
  # RESTRICT was implemented as "record" where the design said "refuse",
  # which no behavioural test could catch, so each item is pinned here.
  cl <- function(j) j$calls$fit[[3L]]
  expect_equal(cl(.sel_job("SELECTION MAXSTEPS=7; EARLY A, B;"))[["max_steps"]], 7)
  expect_equal(cl(.sel_job("SELECTION MOVE=3; EARLY A, B;"))[["max_move"]], 3)
  d <- .sel_job("SELECTION; EARLY A, B;")
  expect_equal(cl(d)[["criterion"]], "score")
  expect_equal(cl(d)[["direction"]], "both")
  expect_equal(cl(d)[["slentry"]], 0.3)
  expect_equal(cl(d)[["slstay"]], 0.2)
  expect_equal(cl(d)[["max_move"]], 1)

  # /E under SELECTION: out of the base AND out of scope, but still in the
  # listwise guard, because PROC HAZARD deletes rows where it is missing.
  j <- .sel_job("SELECTION; EARLY A, B/E, C/S;")
  expect_equal(deparse(j$calls$fit_base[[3L]][["phases"]][[2L]][["formula"]]), "~C")
  expect_equal(deparse(cl(j)[["scope"]]), "list(phase_1 = ~A, phase_2 = NULL)")
  expect_true(any(grepl("\\bB\\b", deparse(j$calls$status))))
})

test_that("ROBUST translates with a loud row, it does not refuse (#160)", {
  skip_on_cran()
  # This INVERTS an earlier test. ROBUST used to refuse, because it changes
  # the variance the removal test is computed from. It is the same CLASS of
  # divergence the translation already ships and documents (PROC HAZARD uses
  # approximate variances while selecting; this package uses the full
  # Hessian), and it is on 90.5% of the SELECTION jobs in the production
  # corpus, so it is said loudly instead of refused. Executed, not just
  # asserted on the row text.
  for (kw in c("ROBUST", "SEMIROBUST")) {
    job <- .sel_job(paste0("SELECTION ", kw, " SLE=0.2; EARLY STRONG, NOISE;"))
    expect_identical(job$calls$fit[[3L]][[1L]], as.name("hzr_stepwise"),
                     info = kw)
    u <- job$untranslated
    expect_true(any(u$construct == kw), info = kw)
    expect_match(u$reason[u$construct == kw][1L], "removal tests", info = kw)
    res <- suppressWarnings(render_sim(job, list(D = .sel_data())))
    expect_true(res$ok, info = paste(kw, paste(res$results, collapse = "; ")))
    expect_s3_class(res$env$fit, "hzr_stepwise")
  }
})

test_that("the callout names ROBUST when the job asks for it (#160)", {
  # The reason most likely to apply, so it is asserted like the others: this
  # is the acceptance condition for translating ROBUST rather than refusing.
  doc <- TemporalHazard:::.hzr_render_qmd(
    .sel_job("SELECTION ROBUST SLE=0.2; EARLY STRONG, NOISE;"))
  hit <- grep("asks for a ROBUST", doc)
  expect_length(hit, 1L)
  chunk <- grep("^#\\| label: fit$", doc)
  expect_lt(hit, chunk)
  expect_match(paste(doc, collapse = " "), "removed, and at which step")
  expect_no_match(paste(doc, collapse = " "), "reproduces PROC HAZARD")
  # Absent for a job that does not ask for it.
  plain <- TemporalHazard:::.hzr_render_qmd(
    .sel_job("SELECTION SLE=0.2; EARLY STRONG, NOISE;"))
  expect_length(grep("asks for a ROBUST", plain), 0L)
  expect_length(grep("asks for a SEMIROBUST", plain), 0L)
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

# --- #160 r-reviewer findings ------------------------------------------------

test_that("scope is always emitted, so the screen never ranges over the data (#160)", {
  skip_on_cran()
  # hzr_stepwise(scope = NULL) enumerates EVERY data-frame column not already
  # in the model. Omitting the argument when no phase offers a candidate
  # therefore handed the screen the whole SAS dataset: the reviewer's run
  # entered DEAD, the EVENT count, as a covariate, under a callout saying
  # the screen used this job's candidates.
  job <- .sel_job("SELECTION SLE=0.5; EARLY STRONG/S, MAL/I;")
  cl <- job$calls$fit[[3L]]
  expect_false(is.null(cl[["scope"]]))
  sc <- as.list(cl[["scope"]])[-1L]
  expect_true(all(names(sc) %in% c("phase_1", "phase_2")))
  res <- suppressWarnings(render_sim(job, list(D = .sel_data())))
  expect_true(res$ok, info = paste(res$results, collapse = "; "))
  steps <- as.data.frame(res$env$fit)
  entered <- steps$variable[steps$action == "enter"]
  expect_false(any(c("DEAD", "TT", ".hzr_status") %in% entered),
               info = paste(entered, collapse = ","))
  expect_length(entered, 0L)
})

test_that("BACKWARD wins over any other direction keyword, in any order (#160)", {
  # stpwprc.c:16-22 is order-independent: STEPWISE sets sw = 1, and BACKWARD
  # then sets sw = 0 and bw = 1 whatever the order. A switch() over the
  # operands made it last-wins, so "SELECTION BACKWARD STEPWISE;" ran a
  # two-way screen from a base SAS never uses, at the wrong slstay.
  for (ops in list(c("BACKWARD", "STEPWISE"), c("STEPWISE", "BACKWARD"),
                   c("BACKWARD", "NOSTEPWISE"), c("NOSTEPWISE", "BACKWARD"))) {
    got <- .hzr_selection_spec(ops)
    expect_equal(got$direction, "backward", info = paste(ops, collapse = " "))
    expect_equal(got$slstay, 0.05, info = paste(ops, collapse = " "))
  }
})

test_that("out-of-range SELECTION thresholds fall back to SAS's defaults (#160)", {
  # stpwprc.c:27-35 clamps an SLE outside (0,1) to 0.3 and drops an
  # out-of-range SLS back to its default; przconc.c:25 clamps MOVE to 1
  # under NOSTEPWISE. Passing them through meant slstay = 0 (nothing is ever
  # removed) or slentry = 5 (everything enters), silently.
  got <- .hzr_selection_spec(c("SLENTRY=5", "SLSTAY=0"))
  expect_equal(got$slentry, 0.3)
  expect_equal(got$slstay, 0.2)
  expect_true(any(grepl("SLENTRY", got$untranslated$construct)))
  expect_true(any(grepl("SLSTAY", got$untranslated$construct)))
  got <- .hzr_selection_spec(c("NOSTEPWISE", "MOVE=3"))
  expect_equal(got$max_move, 1)
  expect_true(any(grepl("MOVE", got$untranslated$construct)))
})

test_that("each refusal names only what actually fired (#160)", {
  # The reason was one boilerplate string naming every refusable construct,
  # so a test asserting "FAST" passed for a MAXVARS job: the assertions
  # could not discriminate and a swap of the switch arms survived.
  reason <- function(stmts) {
    u <- .sel_job(stmts)$untranslated
    u$reason[u$construct == "SELECTION"][1L]
  }
  r <- reason("SELECTION FAST; EARLY A, B;")
  expect_match(r, "FAST")
  expect_no_match(r, "MAXVARS")
  expect_no_match(r, "ROBUST")
  r <- reason("SELECTION MAXVARS=2; EARLY A, B;")
  expect_match(r, "MAXVARS")
  expect_no_match(r, "FAST")
  expect_no_match(r, "ROBUST")
  r <- reason("SELECTION; EARLY A/ORDER=1, B;")
  expect_match(r, "ORDER")
  expect_no_match(r, "FAST")
})

test_that("a phase option PROC HAZARD rejects is refused, not screened (#160)", {
  # hazard_y.y is `phaseopt : MOVE '=' NUMBER`, so a bare /MOVE is a syntax
  # error. The refusal grepped for "/MOVE=" and missed it, and the job got a
  # running screen.
  job <- .sel_job("SELECTION; EARLY STRONG/MOVE, NOISE;")
  expect_identical(job$calls$fit[[3L]][[1L]], as.name("stop"))
})

test_that("force_in names only variables a built phase carries (#160)", {
  # sel_force_in was collected across all three phases while scope and
  # movable were keyed to BUILT phases, so a /I in a phase with no shape
  # operand reached force_in, where hzr_stepwise() ignores it silently.
  job <- .sel_job("SELECTION; EARLY A, B/I; LATE ZZZ/I;")
  fi <- eval(job$calls$fit[[3L]][["force_in"]])
  expect_equal(fi, "B")
})

test_that("the screen check names its own fit and avoids %||% (#160)", {
  # do.call(substitute) rewrites the symbol but not a string literal, so a
  # second block's check told the reader to look at `fit`, which in that
  # document is a DIFFERENT object. %||% is base R only since 4.4 and this
  # package supports 4.1, so an emitted chunk must not use it.
  f <- withr::local_tempfile(fileext = ".sas")
  one <- paste("PROC HAZARD DATA=D CONDITION=14; EVENT DEAD; TIME TT;",
               "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005; SELECTION SLE=0.2;",
               "EARLY STRONG, NOISE;")
  writeLines(paste0("%HAZARD( ", one, " );\n%HAZARD( ", one, " );"), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  chk <- grep("^screen_check", names(job$calls), value = TRUE)
  expect_length(chk, 2L)
  # Executed, because the label is pasted from pieces in the source: what
  # matters is the message the reader sees.
  env <- new.env(parent = baseenv())
  env$fit_2 <- list(criteria = list(n_uncomputable_scores = 2L))
  expect_warning(eval(job$calls[[chk[2L]]], env),
                 "fit_2\\$criteria\\$uncomputable_reasons")
  all_src <- paste(vapply(job$calls, function(c0) {
    paste(deparse(c0), collapse = " ")
  }, character(1)), collapse = " ")
  expect_no_match(all_src, "%||%", fixed = TRUE)
})

test_that("the listwise guard says what a candidate actually is (#160)", {
  # Candidates are outside the BASE model but the screen does see them, so
  # the guard's "hazard() never sees these variables" was false for them.
  job <- .sel_job("SELECTION; EARLY STRONG, NOISE;")
  msg <- paste(deparse(job$calls$status), collapse = " ")
  expect_no_match(msg, "never sees these variables")
  expect_match(msg, "candidate")
})
