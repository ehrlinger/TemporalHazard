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
  expect_equal(got$direction, "both")
  expect_equal(got$slentry, 0.05)
  expect_equal(got$slstay, 0.1)
})

test_that("a bare SELECTION with no options at all enables stepwise", {
  got <- .hzr_selection_spec(character(0))
  expect_equal(got$direction, "both")
})

test_that("NOSTEPWISE keeps its entry threshold, because it still screens", {
  got <- .hzr_selection_spec(c("NOSTEPWISE", "SLENTRY=0.05"))
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
  # max_move is deliberately NOT emitted: see the MOVE test below.
  expect_null(cl[["max_move"]])

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
  # The MOVE divergence is named too: a variable can be frozen here that
  # PROC HAZARD would still move.
  expect_match(txt, "still move")
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
  # PROC HAZARD's MAXSTEPS default is INT_MAX (stpwprc.c:73-74), not
  # hzr_stepwise()'s 50, so it is written out like the other defaults.
  expect_equal(cl(.sel_job("SELECTION; EARLY A, B;"))[["max_steps"]],
               .Machine$integer.max)
  d <- .sel_job("SELECTION; EARLY A, B;")
  expect_equal(cl(d)[["criterion"]], "score")
  expect_equal(cl(d)[["direction"]], "both")
  expect_equal(cl(d)[["slentry"]], 0.3)
  expect_equal(cl(d)[["slstay"]], 0.2)
  expect_null(cl(d)[["max_move"]])

  # /E under SELECTION: out of the base AND out of scope, but still in the
  # listwise guard, because PROC HAZARD deletes rows where it is missing.
  j <- .sel_job("SELECTION; EARLY A, B/E, C/S;")
  expect_equal(deparse(j$calls$fit_base[[3L]][["phases"]][[2L]][["formula"]]), "~C")
  expect_equal(deparse(cl(j)[["scope"]]), "list(phase_1 = ~A, phase_2 = NULL)")
  expect_true(any(grepl("\\bB\\b", deparse(j$calls$status))))
})

test_that("ROBUST translates, recorded as the optimizer choice it is (#160)", {
  skip_on_cran()
  # ROBUST and SEMIROBUST choose PROC HAZARD's OPTIMIZER for the stepwise
  # step (stpwprc.c:60-69 sets swnewt/swhess; hazrd2.c:31-34 swaps them in
  # around the step; cmpmeth.c: Newton vs quasi-Newton, Hessian vs
  # steepest-descent start). They do not change the variance or the Wald
  # tests that drive selection. Both earlier rulings -- refuse, then
  # translate "because it changes the variance" -- rested on the false
  # premise that they did. Executed, not just asserted on the row text.
  for (kw in c("ROBUST", "SEMIROBUST")) {
    job <- .sel_job(paste0("SELECTION ", kw, " SLE=0.2; EARLY STRONG, NOISE;"))
    expect_identical(job$calls$fit[[3L]][[1L]], as.name("hzr_stepwise"),
                     info = kw)
    u <- job$untranslated
    reason <- u$reason[u$construct == kw][1L]
    expect_match(reason, "optimizer", info = kw)
    # It must NOT claim a variance difference, which is the false mechanism
    # two rulings rested on. Pinned so it cannot be reinstated from a stale
    # memory of the decision thread.
    expect_no_match(reason, "variance", info = kw)
    # Nor may it promise the same optimum: the likelihood is multimodal.
    expect_no_match(reason, "not the converged estimates", info = kw)
    expect_match(reason, "optimum", info = kw)
    res <- suppressWarnings(render_sim(job, list(D = .sel_data())))
    expect_true(res$ok, info = paste(kw, paste(res$results, collapse = "; ")))
    expect_s3_class(res$env$fit, "hzr_stepwise")
  }
})

test_that("the callout does not claim ROBUST changes the variance (#160)", {
  # INVERTED from an earlier test that required the callout to name ROBUST
  # as a variance divergence. ROBUST is an optimizer choice: it does not
  # change the variance, though a different optimizer can reach a different
  # optimum on a multimodal likelihood and so change the selected model (the
  # $untranslated row says so). What this guards is the false VARIANCE
  # mechanism, which two rulings rested on.
  doc <- paste(TemporalHazard:::.hzr_render_qmd(
    .sel_job("SELECTION ROBUST SLE=0.2; EARLY STRONG, NOISE;")), collapse = " ")
  expect_no_match(doc, "ROBUST variance")
  expect_no_match(doc, "robust variance")
  expect_no_match(doc, "SEMIROBUST variance")
  # The callout is still there, with its real mechanisms.
  expect_match(doc, "may select a different model")
  expect_match(doc, "approximate variances")
  expect_no_match(doc, "reproduces PROC HAZARD")
  # The entry statistic reproduces SAS's approximate Q (score-test.R:15-18);
  # only the Wald removal tests use the full Hessian. Say drops, not both.
  expect_match(doc, "drop decisions", fixed = TRUE)
  expect_no_match(doc, "each enter and drop", fixed = TRUE)
})

test_that("MOVE= is recorded, not mapped onto max_move (#160)", {
  # The two count different things: PROC HAZARD counts DELETIONS only
  # (hazrd4.c:361-362) per (variable, PHASE) slot (setstat.c:15), while
  # hzr_stepwise()'s max_move counts entries AND exits keyed by variable
  # name across phases. A row is written whether or not the job sets MOVE=,
  # because PROC HAZARD's default of 1 applies either way.
  cases <- c("SELECTION MOVE=3; EARLY A, B;" = "MOVE=3",
             "SELECTION; EARLY A, B;" = "MOVE (PROC HAZARD default 1)")
  for (stmts in names(cases)) {
    job <- .sel_job(stmts)
    expect_null(job$calls$fit[[3L]][["max_move"]], info = stmts)
    u <- job$untranslated
    hit <- u$construct == cases[[stmts]]
    expect_equal(sum(hit), 1L, info = stmts)
    expect_match(u$reason[hit], "deletions, separately for each phase",
                 info = stmts)
    # At the default a removed variable can never re-enter (swvari.c:159-160);
    # hzr_stepwise() re-enters freely, which is the difference that bites.
    expect_match(u$reason[hit], "never re-enter", info = stmts)
  }
})

test_that("the callout says a removed variable may re-enter here (#160)", {
  doc <- paste(TemporalHazard:::.hzr_render_qmd(
    .sel_job("SELECTION; EARLY STRONG, NOISE;")), collapse = " ")
  expect_match(doc, "can never return to it")
  expect_match(doc, "can re-enter variables PROC HAZARD would have kept out")
})

test_that("a fractional MAXSTEPS truncates, as PROC HAZARD's int cast does (#160)", {
  # stpwprc.c:82 casts MAXSTEPS to int.
  cl <- function(j) j$calls$fit[[3L]]
  expect_identical(cl(.sel_job("SELECTION MAXSTEPS=2.9; EARLY A, B;"))[["max_steps"]], 2)
})

test_that("a negative MAXSTEPS is refused, as PROC HAZARD refuses the job (#160)", {
  # stpwprc.c:76-79 logs an ERROR and exits, so there is no run to
  # translate; passing it through gave hzr_stepwise() a budget that ends
  # the screen on its first test.
  job <- .sel_job("SELECTION MAXSTEPS=-3; EARLY A, B;")
  expect_identical(job$calls$fit[[3L]][[1L]], as.name("stop"))
  expect_error(eval(job$calls$fit), "MAXSTEPS")
})

test_that("a /I in a phase this job does not build is not a cross-phase pin (#160)", {
  # force_in was collected over all three phase statements while movable
  # covered only built phases, so a /I naming an unbuilt phase refused a
  # job that has no phase for it to conflict with.
  job <- .sel_job("SELECTION; EARLY STRONG; LATE STRONG/I;")
  expect_identical(job$calls$fit[[3L]][[1L]], as.name("hzr_stepwise"))
})

test_that("RESTRICT without a SELECTION stays a recorded gap (#160)", {
  # Recording it under a saw_restrict flag counted it as MAPPED and dropped
  # its $untranslated row, inflating coverage for a statement nothing
  # implements.
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste("%HAZARD( PROC HAZARD DATA=D CONDITION=14; EVENT DEAD;",
                   "TIME TT; PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005;",
                   "EARLY A, B; RESTRICT A; );"), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  expect_true(any(job$untranslated$construct == "RESTRICT"))
  expect_identical(job$calls$fit[[3L]][[1L]], as.name("hazard"))
})

test_that("printing options are recorded, not refused (#160)", {
  job <- .sel_job("SELECTION NOPRINTS NOPRINTQ; EARLY STRONG, NOISE;")
  expect_identical(job$calls$fit[[3L]][[1L]], as.name("hzr_stepwise"))
  expect_true(any(grepl("NOPRINTS", job$untranslated$construct)))
})

test_that("the emitted screen reports candidates it could not score (#160, #400)", {
  # This used to grep the rendered document for "uncomputable", a word the
  # callout body also contains, so it passed with the check chunk deleted.
  # Assert the chunk exists and EXECUTE it. The mocks carry the criteria a
  # real hzr_stepwise() fit carries since #399: an attempt count, the named
  # reasons, and whether the screen stopped on them.
  job <- .sel_job("SELECTION SLE=0.2; EARLY STRONG, NOISE;")
  expect_true("screen_check" %in% names(job$calls))
  env <- new.env(parent = baseenv())
  crit <- function(reasons, stopped = FALSE) {
    list(criteria = list(n_uncomputable_scores = sum(reasons),
                         uncomputable_reasons = reasons,
                         stopped_uncomputable = stopped))
  }
  # A completed screen that could not score two candidates says so, by reason.
  env$fit <- crit(c(nuisance_singular = 2L))
  msg <- tryCatch({
    eval(job$calls$screen_check, env)
    "none"
  }, warning = conditionMessage)
  expect_match(msg, "2 candidate score(s)", fixed = TRUE)
  expect_match(msg, "nuisance_singular = 2", fixed = TRUE)
  env$fit <- crit(stats::setNames(integer(0), character(0)))
  expect_no_warning(eval(job$calls$screen_check, env))
  # A Wald test without a variance is not an uncomputable SCORE, and
  # hzr_stepwise() reports it itself, naming the variables (#389, #400).
  env$fit <- crit(c(wald_no_variance = 3L))
  expect_no_warning(eval(job$calls$screen_check, env))
  # A screen that STOPPED on uncomputable candidates has already warned, with
  # every reason; the check must not say it a second time (#400).
  env$fit <- crit(c(nuisance_singular = 3L), stopped = TRUE)
  expect_no_warning(eval(job$calls$screen_check, env))
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
  env$fit_2 <- list(criteria = list(n_uncomputable_scores = 2L,
                                    uncomputable_reasons = c(nuisance_singular = 2L),
                                    stopped_uncomputable = FALSE))
  expect_warning(eval(job$calls[[chk[2L]]], env),
                 "fit_2\\$criteria\\$uncomputable_reasons")
  # The substitution renames the SYMBOL `fit`, so a component reached as
  # `$fit` would be renamed too (#160: fit_2$fit_2$vcov). A screen that
  # scored everything must not warn under the second block's name either.
  env$fit_2 <- list(criteria = list(n_uncomputable_scores = 0L,
                                    uncomputable_reasons = stats::setNames(integer(0), character(0)),
                                    stopped_uncomputable = FALSE))
  expect_no_warning(eval(job$calls[[chk[2L]]], env))
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

# --- SELECTION must not blind a refusal that reads the base model ----------
# Under a forward or two-way screen a bare phase variable is WITHHELD from the
# base model (setstat.c, H->sw == 1). A refusal whose condition reads the
# base model's covariates therefore cannot see a candidate. #311's no-DATA=
# refusal was one: with SELECTION it went quiet and the screen then failed
# on "`data` must be a data frame", naming neither DATA= nor the fix.

.nodata_job <- function(stmts, parms = "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005;",
                        data = "", env = parent.frame()) {
  f <- withr::local_tempfile(fileext = ".sas", .local_envir = env)
  writeLines(paste0("%HAZARD( PROC HAZARD", data, " CONDITION=14; EVENT DEAD; ",
                    "TIME TT; ", parms, " ", stmts, " );"), f)
  suppressWarnings(hzr_translate_sas(f))
}

test_that("no DATA= is refused whatever SELECTION withholds from the base (#160, #311)", {
  # Only the first and last cases were blind before the fix; the middle three
  # put a variable in the base model and are controls for the refusal itself.
  for (stmts in c("SELECTION; EARLY A, B;",           # candidates only
                  "SELECTION; EARLY A, B/S;",         # a candidate beside a movable
                  "SELECTION; EARLY A/I, B;",         # a candidate beside a /I
                  "SELECTION BACKWARD; EARLY A, B;",  # candidates start IN the base
                  "SELECTION; EARLY A/E;")) {         # only an excluded variable
    job <- .nodata_job(stmts)
    expect_identical(job$calls$fit[[3L]][[1L]], as.name("stop"), info = stmts)
    expect_error(eval(job$calls$fit), "names no DATA= dataset", info = stmts)
    expect_true("DATA=" %in% job$untranslated$construct, info = stmts)
    expect_false("fit_base" %in% names(job$calls), info = stmts)
  }
  # Control: the same jobs WITH DATA= are not refused for it.
  job <- .nodata_job("SELECTION; EARLY A, B;", data = " DATA=D")
  expect_identical(job$calls$fit[[3L]][[1L]], as.name("hzr_stepwise"))
  expect_false("DATA=" %in% job$untranslated$construct)
})

test_that("no DATA= with only an /E variable is refused, SELECTION or not (#311)", {
  # A DELIBERATE WIDENING of #311, decided by the release coordinator: this
  # job used to translate, and its listwise guard then read A from whatever
  # environment rendered the document. That is the lookup #311 closed for
  # base covariates. Narrowing the guard to SELECTION candidates removes it.
  job <- .nodata_job("EARLY A/E;")
  expect_identical(job$calls$fit[[3L]][[1L]], as.name("stop"))
  expect_error(eval(job$calls$fit), "names no DATA= dataset")
  expect_true("DATA=" %in% job$untranslated$construct)
  expect_null(job$calls$status)
})

test_that("no DATA= with only an unbuilt phase's covariates is refused, truthfully (#311)", {
  # The same widening as /E, found by review: these covariates are read only
  # by the missing-row guard, so the job used to translate and read them from
  # the rendering environment. The message must not claim "a phase has
  # covariates" -- no emitted phase does -- so it names the variables instead.
  cases <- c("PARMS MUE=0.2 THALF=0.15 NU=1;" = "LATE X;",       # no MUL
             "PARMS MUE=0.2 THALF=0.15 NU=1; " = "CONSTANT X;",  # no MUC
             "PARMS MUC=0.0005;" = "EARLY X;")                   # no MUE
  for (i in seq_along(cases)) {
    info <- paste(names(cases)[i], cases[[i]])
    job <- .nodata_job(cases[[i]], parms = names(cases)[i])
    expect_identical(job$calls$fit[[3L]][[1L]], as.name("stop"), info = info)
    msg <- tryCatch(eval(job$calls$fit), error = conditionMessage)
    expect_match(msg, "names no DATA= dataset, but its phase statements name X,",
                 fixed = TRUE, info = info)
    expect_no_match(msg, "a phase has covariates", info = info)
    expect_match(job$untranslated$reason[job$untranslated$construct == "DATA="],
                 "name X,", fixed = TRUE, info = info)
  }
  # And where an emitted phase DOES carry covariates, the message says so.
  job <- .nodata_job("EARLY A;")
  expect_error(eval(job$calls$fit), "a phase has covariates")
})

test_that("a no-DATA= job keeps its SELECTION refusal reason too (#160, #340 item 3)", {
  # The DATA= refusal returns first. The reader who adds DATA= as told must
  # not then meet a second refusal the first document never mentioned.
  job <- .nodata_job("SELECTION FAST; EARLY A, B;")
  u <- job$untranslated
  expect_true("DATA=" %in% u$construct)
  expect_equal(sum(u$construct == "SELECTION"), 1L)
  expect_match(u$reason[u$construct == "SELECTION"], "FAST is a different search")
  expect_error(eval(job$calls$fit), "names no DATA= dataset")
  # And a SELECTION this translator can run records no SELECTION refusal.
  ok <- .nodata_job("SELECTION; EARLY A, B;")
  expect_false("SELECTION" %in% ok$untranslated$construct)
})

test_that("candidates of a phase that is not built still leave a row (#160)", {
  # PROC HAZARD skips an unbuilt phase's covariates (setstat.c:9-12) and the
  # translator records them. Under SELECTION they were withheld from the list
  # the row is written from, so a candidate-only phase vanished silently.
  job <- .nodata_job("SELECTION; EARLY A; CONSTANT B, C/S;", data = " DATA=D",
                     parms = "PARMS MUE=0.2 THALF=0.15 NU=1;")
  u <- job$untranslated
  hit <- grepl("constant phase covariates with no active MUC", u$reason)
  expect_equal(sum(hit), 1L)
  expect_equal(u$construct[hit], "B C")
  job <- .nodata_job("SELECTION; EARLY A, B/S; CONSTANT C;", data = " DATA=D",
                     parms = "PARMS MUC=0.0005;")
  u <- job$untranslated
  hit <- grepl("early phase material with no active MUE", u$reason)
  expect_equal(sum(hit), 1L)
  expect_match(u$construct[hit], "\\bA B\\b")
  job <- .nodata_job("SELECTION; EARLY A; LATE X, Y/S;", data = " DATA=D",
                     parms = "PARMS MUE=0.2 THALF=0.15 NU=1;")
  u <- job$untranslated
  hit <- grepl("late phase material with no active MUL", u$reason)
  expect_equal(sum(hit), 1L)
  expect_match(u$construct[hit], "\\bX Y\\b")
})

test_that("every other refusal still fires when the job carries SELECTION (#160)", {
  # The census behind the no-DATA= fix above: each translator refusal whose
  # condition could depend on the model, run WITH a SELECTION statement. A
  # "no" needs a test as much as a "yes": the next refusal added should
  # show which side of the line it falls on. None of these read the base
  # model's covariates -- they read the censoring statements, the PARMS
  # operands or the SETG3 shape rules -- so SELECTION cannot blind them.
  # (The HAZPRED grid, %repeat, rewrite-step and INHAZ refusals are not in
  # the census: they read other blocks, never the PROC HAZARD model.)
  cases <- list(
    list(parms = "PARMS MUE=0.2 THALF=0.15 NU=1;", extra = "LCENSOR ST; ICENSOR C3 = CT;",
         msg = "combines LCENSOR"),
    # This one CRASHED the translation under SELECTION: the scope names were
    # paste0("phase_", seq_along(built)), which is "phase_" when no phase is
    # built, so setNames() failed and the whole file failed to translate.
    list(parms = "PARMS THALF=0.15 NU=1;", extra = "", msg = "selects no phase"),
    # Operands this parser cannot read build no phase and are not refused
    # (SAS expands the macro and runs the job). A MUE with no shape operand
    # used to stand here (it now builds on SAS's defaults, #345), then a
    # template's `?` (now a syntax stop, U1), then a spaced operand (joined
    # and read now, #421).
    list(parms = "PARMS &ALLPARMS;", extra = "",
         msg = "builds no phase this translator could use"))
  for (cs in cases) {
    for (sel in c("", "SELECTION;")) {
      info <- paste(cs$parms, cs$extra, sel)
      job <- .nodata_job(paste(cs$extra, sel, "EARLY A, B;"), parms = cs$parms,
                         data = " DATA=D")
      expect_identical(job$calls$fit[[3L]][[1L]] %||% job$calls$fit[[1L]],
                       as.name("stop"), info = info)
      expect_error(eval(job$calls$fit, new.env()), cs$msg, info = info)
    }
  }
  # The SETG3 entry refusals are $untranslated rows naming PROC HAZARD's
  # own refusal; SELECTION must not change which rows are written.
  setg3 <- "PARMS MUE=0.2 THALF=0.15 NU=1 MUL=0.01 TAU=0 FIXTAU GAMMA=2 ETA=2;"
  rows <- lapply(c("", "SELECTION;"), function(sel) {
    u <- .nodata_job(paste(sel, "EARLY A, B;"), parms = setg3, data = " DATA=D")$untranslated
    u$reason[grepl("SETG3", u$reason)]
  })
  expect_true(any(grepl("SETG3900", rows[[1L]])))
  expect_identical(rows[[1L]], rows[[2L]])
})

# --- Copilot review of c91c35b5, each finding reproduced before fixing -----

test_that("a SELECTION job with no DATA= is refused even with nothing to screen (#160)", {
  # No phase variable at all, so neither the covariate nor the listwise test
  # fired, and hzr_stepwise() then stopped on "`data` must be a data frame".
  # A screen always needs `data`, whatever it screens.
  job <- .nodata_job("SELECTION;")
  expect_identical(job$calls$fit[[3L]][[1L]], as.name("stop"))
  msg <- tryCatch(eval(job$calls$fit), error = conditionMessage)
  expect_match(msg, "names no DATA= dataset")
  expect_match(msg, "hzr_stepwise() needs `data`", fixed = TRUE)
  expect_false("fit_base" %in% names(job$calls))
})

test_that("a no-DATA= refusal mentions SELECTION only for a SELECTION job (#311)", {
  plain <- tryCatch(eval(.nodata_job("EARLY A/E;")$calls$fit), error = conditionMessage)
  expect_no_match(plain, "SELECTION")
  expect_match(plain, "PROC HAZARD deletes rows where any is missing")
  sel <- tryCatch(eval(.nodata_job("SELECTION; EARLY A, B;")$calls$fit),
                  error = conditionMessage)
  expect_match(sel, "hzr_stepwise() needs `data`", fixed = TRUE)
})

test_that("a negative MAXSTEPS is refused with its own reason, not force_in's (#160)", {
  job <- .sel_job("SELECTION MAXSTEPS=-3; EARLY A, B;")
  msg <- tryCatch(eval(job$calls$fit), error = conditionMessage)
  expect_match(msg, "a negative MAXSTEPS", fixed = TRUE)
  expect_match(msg, "stpwprc.c:76-79", fixed = TRUE)
  expect_no_match(msg, "force_in")
})

test_that("/I in a phase that is not built pins nothing (#160)", {
  # PROC HAZARD skips an unbuilt phase's variables (setstat.c:9-12), flags
  # included. force_in used to be filtered against every variable of the
  # built phases, so LATE A/I with no MUL pinned the early phase's movable A.
  job <- .nodata_job("SELECTION; EARLY A, B; LATE A/I;", data = " DATA=D",
                     parms = "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005;")
  cl <- job$calls$fit[[3L]]
  expect_identical(cl[[1L]], as.name("hzr_stepwise"))
  expect_null(cl[["force_in"]])
  expect_equal(deparse(cl[["scope"]]), "list(phase_1 = ~A + B, phase_2 = NULL)")
  # Control: the same /I in a BUILT phase is still a cross-phase pin, refused.
  built <- .nodata_job("SELECTION; EARLY A, B; LATE A/I;", data = " DATA=D",
                       parms = "PARMS MUE=0.2 THALF=0.15 NU=1 MUL=0.01 TAU=1 GAMMA=2 ETA=2;")
  expect_identical(built$calls$fit[[3L]][[1L]], as.name("stop"))
  expect_error(eval(built$calls$fit), "/I in one phase, movable in another", fixed = TRUE)
})

test_that("the screen refits under the job's own control, not hazard()'s defaults (#160)", {
  skip_on_cran()
  # hzr_stepwise() forwards only its `...` to each refit. Without the job's
  # control there, a NOCONSERVE job fitted its base without Conservation of
  # Events and then reported a selected model refitted WITH it (and at
  # hazard()'s default maxit): a different fitting contract, silently.
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=D NOCONSERVE MAXITER=77; EVENT DEAD; TIME TT;",
    "PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005;",
    "SELECTION SLE=0.2; EARLY STRONG, NOISE; );"), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  expect_identical(job$calls$fit[[3L]][["control"]],
                   job$calls$fit_base[[3L]][["control"]])
  res <- suppressWarnings(render_sim(job, list(D = .sel_data())))
  expect_true(res$ok, info = paste(res$results, collapse = "; "))
  final <- res$env$fit
  # A step was taken, so the reported model IS a refit.
  expect_true("phase_1.STRONG" %in% names(stats::coef(final)))
  expect_false(final$spec$control$conserve)
  expect_equal(final$spec$control$maxit, 77)
})

test_that("the screen check leaves the reader's own objects alone (#400)", {
  # The chunk runs in the reader's session, so each local it creates shares
  # a namespace with the job's data: a bare `rs` would have overwritten a
  # job's DATA=rs, and a bare `n_unscored` a variable of that name (Codex on
  # #419). Assert the property, not the name: a pre-existing object keeps
  # its value, and the chunk adds no non-dotted name to the frame, whether
  # it warns or not.
  job <- .sel_job("SELECTION SLE=0.2; EARLY STRONG, NOISE;")
  for (reasons in list(c(nuisance_singular = 2L),
                       stats::setNames(integer(0), character(0)))) {
    env <- new.env(parent = baseenv())
    env$fit <- list(criteria = list(n_uncomputable_scores = sum(reasons),
                                    uncomputable_reasons = reasons,
                                    stopped_uncomputable = FALSE))
    env$n_unscored <- "user data"
    env$rs <- "user data"
    before <- ls(env)
    suppressWarnings(eval(job$calls$screen_check, env))
    expect_identical(env$n_unscored, "user data")
    expect_identical(env$rs, "user data")
    expect_identical(ls(env), before)
  }
})

test_that("the screen check does not repeat hzr_stepwise()'s Wald report (#400)", {
  skip_on_cran() # two multiphase fits
  # Rebuilt from a REAL backward fit, not a vcov = NULL mock. A multiphase
  # ICENSOR fit has no variances without numDeriv, so no removal can be
  # tested. hzr_stepwise() says so itself since #389/#399; the emitted check
  # used to say it twice more, and called the Wald failure an uncomputable
  # score.
  orig <- base::requireNamespace
  local_mocked_bindings(
    requireNamespace = function(package, ...) {
      if (identical(package, "numDeriv")) FALSE else orig(package, ...)
    },
    .package = "base"
  )
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(paste(
    "%HAZARD( PROC HAZARD DATA=D CONDITION=14; EVENT DEAD; TIME TT;",
    "ICENSOR C3 = CT; PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005;",
    "SELECTION BACKWARD SLS=0.05; EARLY STRONG, NOISE, MAL; );"
  ), f)
  job <- suppressWarnings(hzr_translate_sas(f))
  d <- .sel_data()
  set.seed(7)
  ic <- sample(which(d$DEAD == 1), 60)
  d$C3 <- 0L
  d$C3[ic] <- 1L
  d$DEAD[ic] <- 0L
  d$CT <- d$TT * 0.6
  env <- new.env()
  env$D <- d
  chunk_warnings <- function(nm) {
    w <- character()
    withCallingHandlers(eval(job$calls[[nm]], env), warning = function(cw) {
      w <<- c(w, conditionMessage(cw))
      invokeRestart("muffleWarning")
    })
    w
  }
  for (nm in setdiff(names(job$calls), c("fit", "screen_check"))) {
    suppressWarnings(eval(job$calls[[nm]], env))
  }
  expect_false(any(is.finite(env$fit_base$fit$se)))
  w_fit <- chunk_warnings("fit")
  # The screen reported it, and it was the Wald test that failed.
  expect_identical(env$fit$criteria$uncomputable_reasons, c(wald_no_variance = 3L))
  expect_true(any(grepl("no remaining candidate could be tested for removal", w_fit, fixed = TRUE)))
  # The check adds nothing to that.
  expect_identical(chunk_warnings("screen_check"), character(0))
})

test_that("SELECTION options written with spaces around = keep their values (#160)", {
  # SAS accepts `SLE = 0.2`; splitting on whitespace turned it into three
  # tokens, recorded as untranslated, and the screen then ran at SAS's
  # DEFAULT thresholds -- a runnable screen with the wrong rule.
  job <- .sel_job("SELECTION SLE = 0.2 SLS = 0.1 MAXSTEPS = 5; EARLY A, B;")
  cl <- job$calls$fit[[3L]]
  expect_equal(cl[["slentry"]], 0.2)
  expect_equal(cl[["slstay"]], 0.1)
  expect_equal(cl[["max_steps"]], 5)
  expect_false(any(c("SLE", "SLS", "MAXSTEPS", "0.2", "") %in% job$untranslated$construct))
})

test_that("the callout describes only what this screen's direction does (#160)", {
  # The callout's own text: "re-enter" also appears, rightly, in the MOVE
  # row of the untranslated callout elsewhere in the document.
  note <- function(stmts) {
    .sel_job(stmts)$notes$fit$body
  }
  both <- note("SELECTION; EARLY STRONG, NOISE;")
  expect_match(both, "The entry statistic", fixed = TRUE)
  expect_match(both, "drop decisions", fixed = TRUE)
  expect_match(both, "can re-enter variables PROC HAZARD would have kept out", fixed = TRUE)
  # BACKWARD never enters, so neither the entry statistic nor re-entry applies.
  back <- note("SELECTION BACKWARD; EARLY STRONG, NOISE;")
  expect_no_match(back, "The entry statistic", fixed = TRUE)
  expect_no_match(back, "re-enter", fixed = TRUE)
  expect_match(back, "drop decisions", fixed = TRUE)
  # NOSTEPWISE is forward only: it never removes, so neither does re-entry.
  fwd <- note("SELECTION NOSW; EARLY STRONG, NOISE;")
  expect_match(fwd, "The entry statistic", fixed = TRUE)
  expect_no_match(fwd, "drop decisions", fixed = TRUE)
  expect_no_match(fwd, "re-enter", fixed = TRUE)
  for (d in list(both, back, fwd)) expect_match(d, "approximate variances", fixed = TRUE)
})

# --- N3 (1.2.12 review): a SELECTION value PROC HAZARD's lexer does not read -
# .hzr_selection_spec() reads a value with as.numeric(), which reads 1E-3,
# 2E-1, +0.1 and 5. -- none of them a NUMBER to hazard_l.l:33-38, so in the
# STEP state (:53) each is unexpected text or an unexpected char, and
# initprz.c:75-77 stops the job with SYNTAX. Measured on the HAZARD binary on
# avc (2026-09-25; PARMS MUE=0.2 THALF=1; SELECTION <value>; EARLY AGE;):
# SLE=0.2, .2, 0.20, 0.2E-1, MAXSTEPS=5 and 5.0 fit (4 of 4 markers, with a
# stepwise section); SLE=ABC, 1E-3, SLS=2E-1, SLE=+0.1 and MAXSTEPS=5. exit
# SYNTAX with no markers. The translation emitted a screen with no warning,
# and with no row for any but ABC. Under U1 it now warns and records the row,
# and still screens at the value as this translation reads it.
.n3_refusal <- function(job) {
  slot <- grep("^refusal", names(job$calls), value = TRUE)
  if (!length(slot)) return(NA_character_)
  job$calls[[slot[[1L]]]][[2L]]
}

test_that("a SELECTION value the lexer does not read warns and still screens (N3)", {
  bad <- list(
    list(op = "SLE=1E-3", arg = "slentry", val = 0.001),
    list(op = "SLS=2E-1", arg = "slstay", val = 0.2),
    list(op = "SLE=+0.1", arg = "slentry", val = 0.1),
    list(op = "MAXSTEPS=5.", arg = "max_steps", val = 5),
    list(op = "SLE=ABC", arg = "slentry", val = 0.3))
  for (b in bad) {
    job <- .sel_job(paste0("SELECTION ", b$op, "; EARLY A, B;"))
    cl <- job$calls$fit[[3L]]
    expect_identical(cl[[1L]], as.name("hzr_stepwise"), info = b$op)
    expect_equal(cl[[b$arg]], b$val, info = b$op)
    msg <- .n3_refusal(job)
    expect_match(msg, "PROC HAZARD does not run this job", fixed = TRUE,
                 info = b$op)
    expect_match(msg, paste0(b$op, ": not a number PROC HAZARD's lexer reads"),
                 fixed = TRUE, info = b$op)
    row <- job$untranslated$reason[job$untranslated$construct == b$op]
    expect_length(row, 1L)
    expect_match(row, "rejects this job with a syntax error", fixed = TRUE,
                 info = b$op)
  }
})

test_that("a SELECTION value the lexer reads raises nothing (N3 known positive)", {
  good <- list(c("SLE=0.2", "slentry", 0.2), c("SLE=.2", "slentry", 0.2),
               c("SLE=0.20", "slentry", 0.2), c("SLE=0.2E-1", "slentry", 0.02),
               c("MAXSTEPS=5", "max_steps", 5), c("MAXSTEPS=5.0", "max_steps", 5))
  for (g in good) {
    job <- .sel_job(paste0("SELECTION ", g[[1L]], "; EARLY A, B;"))
    expect_equal(job$calls$fit[[3L]][[g[[2L]]]], as.numeric(g[[3L]]),
                 info = g[[1L]])
    expect_identical(.n3_refusal(job), NA_character_, info = g[[1L]])
    expect_false(g[[1L]] %in% job$untranslated$construct, info = g[[1L]])
  }
  # A macro value is SAS's to expand, so it carries no verdict.
  expect_identical(.n3_refusal(.sel_job("SELECTION SLE=&SLE; EARLY A, B;")),
                   NA_character_)
})

test_that("a later `(` clears a SELECTION value refusal, as the binary does (N3, #461)", {
  # Measured: `SELECTION SLE=1E-3; EARLY LOG();` exits with no SYNTAX and
  # fits (4 of 4 markers) but prints no stepwise section, while the same job
  # at SLE=0.2 prints "Forward Stepwise Selection". The parser's recovery
  # drops the screen and the translation still emits one, so the job takes
  # the #461 cleared wording rather than "does not run".
  job <- .sel_job("SELECTION SLE=1E-3; EARLY A, B; EARLY LOG();")
  msg <- .n3_refusal(job)
  expect_match(msg, "does not stop this job for the syntax error in SLE=1E-3",
               fixed = TRUE)
  expect_match(msg, "cannot emit PROC HAZARD's model", fixed = TRUE)
  expect_no_match(msg, "does not run", fixed = TRUE)
  row <- job$untranslated$reason[job$untranslated$construct == "SLE=1E-3"]
  expect_length(row, 1L)
  expect_match(row, "not a number PROC HAZARD's lexer reads", fixed = TRUE)
  expect_no_match(row, "rejects this job", fixed = TRUE)
  # A `(` in the SELECTION statement itself does not clear it: measured,
  # `EARLY AGE; SELECTION SLE=1E-3 ();` exits SYNTAX.
  same <- .n3_refusal(.sel_job("EARLY A, B; SELECTION SLE=1E-3 ();"))
  expect_match(same, "PROC HAZARD does not run this job", fixed = TRUE)
})

test_that("a SELECTION value refusal emits a screen that runs (N3)", {
  skip_on_cran()
  job <- .sel_job("SELECTION SLE=1E-3; EARLY STRONG, NOISE;")
  expect_equal(job$calls$fit[[3L]][["slentry"]], 0.001)
  w <- character(0)
  res <- withCallingHandlers(
    render_sim(job, list(D = .sel_data())),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  expect_true(res$ok, info = paste(res$results, collapse = "; "))
  # The document raises the refusal when its chunk runs, not only at
  # translation.
  expect_true(any(grepl("SLE=1E-3: not a number", w, fixed = TRUE)))
  sw <- res$env$fit
  expect_s3_class(sw, "hzr_stepwise")
  steps <- as.data.frame(sw)
  # STRONG's effect clears p < 0.001 on these data, so the screen entered it.
  expect_true("STRONG" %in% steps$variable[steps$action == "enter"])
})

# --- #505: PROC HAZARD accumulates SELECTION statements across the job ------
# Measured on the HAZARD binary on avc with a noise column Z (2026-09-25,
# `PARMS MUE=0.2 THALF=1; ...; EARLY AGE, Z;`), reading the SLE from "No
# other variables met the <SLE> significance level for entry":
# `SLE=0.05` alone 0.05; `SLS=0.1` alone 0.3 (the default);
# `SLE=0.05; SELECTION SLS=0.1;` 0.05; `SLS=0.1; SELECTION SLE=0.05;` 0.05;
# `SLE=0.05; SELECTION SLE=0.07;` 0.07. A direction keyword in either
# statement holds too: `BACKWARD; SELECTION SLE=0.05;` and
# `SLE=0.05; SELECTION BACKWARD;` both print "Backward Stepwise Selection",
# and `NOSW; SELECTION;` prints "Forward Selection". The translation kept
# only the last statement, so the #505 job screened at slentry 0.3.
test_that("SELECTION statements accumulate, a repeated option last-wins (#505)", {
  cl <- function(sel) .sel_job(paste(sel, "EARLY A, B;"))$calls$fit[[3L]]
  rows <- list(
    list("SELECTION SLE=0.05;", 0.05, 0.2),
    list("SELECTION SLS=0.1;", 0.3, 0.1),
    list("SELECTION SLE=0.05; SELECTION SLS=0.1;", 0.05, 0.1),
    list("SELECTION SLS=0.1; SELECTION SLE=0.05;", 0.05, 0.1),
    list("SELECTION SLE=0.05; SELECTION SLE=0.07;", 0.07, 0.2))
  for (r in rows) {
    got <- cl(r[[1L]])
    expect_equal(got[["slentry"]], r[[2L]], info = r[[1L]])
    expect_equal(got[["slstay"]], r[[3L]], info = r[[1L]])
  }
  expect_identical(cl("SELECTION BACKWARD; SELECTION SLE=0.05;")[["direction"]],
                   "backward")
  expect_identical(cl("SELECTION SLE=0.05; SELECTION BACKWARD;")[["direction"]],
                   "backward")
  expect_identical(cl("SELECTION NOSW; SELECTION;")[["direction"]], "forward")
  # The N3 check still reads the first statement's values.
  job <- .sel_job("SELECTION SLE=1E-3; SELECTION SLS=0.1; EARLY A, B;")
  expect_match(.n3_refusal(job), "SLE=1E-3: not a number", fixed = TRUE)
  expect_equal(job$calls$fit[[3L]][["slentry"]], 0.001)
  expect_identical(sum(job$untranslated$construct == "SLE=1E-3"), 1L)
})

test_that("a SELECTION option PROC HAZARD rejects warns with its own reason (#504 review)", {
  # Measured on the binary (avc plus noise Z, `SELECTION <op>; EARLY AGE, Z;`,
  # 2026-09-25): BOGUS=1, BOGUS, BOGUS=ABC, NOPRINTS=1, NOPRINTS=ABC,
  # `SLE 0.2` and MOVE=ABC each exit SYNTAX with no markers; NOPRINTS and
  # SLE=0.2 fit. With `EARLY AGE, Z();` instead, every one of them fits
  # (4 of 4 markers). The lexer-number check read every OPTION=value, so an
  # unknown option was refused as "not a number", and NOPRINTS=1 and a bare
  # BOGUS got a row and no warning.
  cases <- list(
    list(sel = "BOGUS=1", what = "BOGUS=1", why = "unknown SELECTION option"),
    list(sel = "BOGUS", what = "BOGUS", why = "unknown SELECTION option"),
    list(sel = "BOGUS=ABC", what = "BOGUS=ABC", why = "unknown SELECTION option"),
    list(sel = "NOPRINTS=1", what = "NOPRINTS=1",
         why = "a value on an option that takes none"),
    list(sel = "NOPRINTS=ABC", what = "NOPRINTS=ABC",
         why = "a value on an option that takes none"),
    list(sel = "SLE 0.2", what = "SLE 0.2", why = "no value"),
    list(sel = "MOVE=ABC", what = "MOVE=ABC",
         why = "not a number PROC HAZARD's lexer reads"))
  for (cs in cases) {
    job <- .sel_job(paste0("SELECTION ", cs$sel, "; EARLY A, B;"))
    msg <- .n3_refusal(job)
    expect_match(msg, "PROC HAZARD does not run this job", fixed = TRUE,
                 info = cs$sel)
    expect_match(msg, paste0(cs$what, ": ", cs$why), fixed = TRUE, info = cs$sel)
    if (!identical(cs$why, "not a number PROC HAZARD's lexer reads")) {
      expect_no_match(msg, "not a number", fixed = TRUE, info = cs$sel)
    }
    u <- job$untranslated
    expect_identical(sum(u$construct == cs$what), 1L, info = cs$sel)
    expect_match(u$reason[u$construct == cs$what], cs$why, fixed = TRUE,
                 info = cs$sel)
    expect_identical(job$calls$fit[[3L]][[1L]], as.name("hzr_stepwise"),
                     info = cs$sel)
  }
  # One row each, not a second "unknown" or "non-numeric" row beside it.
  u <- .sel_job("SELECTION BOGUS=1; EARLY A, B;")$untranslated
  expect_false("BOGUS" %in% u$construct)
  u <- .sel_job("SELECTION SLE=ABC; EARLY A, B;")$untranslated
  expect_false("SLE" %in% u$construct)
  u <- .sel_job("SELECTION MOVE=ABC; EARLY A, B;")$untranslated
  expect_identical(setdiff(grep("^MOVE", u$construct, value = TRUE),
                           "MOVE (PROC HAZARD default 1)"), "MOVE=ABC")
  # `SLE 0.2` screens at the default, as the value was never assigned.
  expect_equal(.sel_job("SELECTION SLE 0.2; EARLY A, B;")$calls$fit[[3L]][["slentry"]],
               0.3)
  # A value on a direction keyword still keeps the direction it names.
  back <- .sel_job("SELECTION BACKWARD=1; EARLY A, B;")
  expect_match(.n3_refusal(back), "BACKWARD=1: a value on an option", fixed = TRUE)
  expect_identical(back$calls$fit[[3L]][["direction"]], "backward")
  # Known positives: the same keywords written as the grammar has them.
  for (ok in c("NOPRINTS", "SLE=0.2", "BACKWARD", "MOVE=2")) {
    expect_identical(.n3_refusal(.sel_job(paste0("SELECTION ", ok, "; EARLY A, B;"))),
                     NA_character_, info = ok)
  }
  # A later `(` clears each, as the binary fits them all (#461 wording).
  for (cs in c("BOGUS=1", "NOPRINTS=1", "SLE 0.2")) {
    msg <- .n3_refusal(.sel_job(paste0("SELECTION ", cs, "; EARLY A, B();")))
    expect_match(msg, paste0("does not stop this job for the syntax error in ", cs),
                 fixed = TRUE, info = cs)
  }
})

test_that("an accumulated SELECTION job screens at the first statement's SLE (#505)", {
  skip_on_cran()
  job <- .sel_job("SELECTION SLE=0.05; SELECTION SLS=0.1; EARLY STRONG, NOISE;")
  res <- suppressWarnings(render_sim(job, list(D = .sel_data())))
  expect_true(res$ok, info = paste(res$results, collapse = "; "))
  sw <- res$env$fit
  expect_s3_class(sw, "hzr_stepwise")
  expect_equal(sw$criteria$slentry, 0.05)
  expect_equal(sw$criteria$slstay, 0.1)
  expect_true("STRONG" %in% as.data.frame(sw)$variable)
})
