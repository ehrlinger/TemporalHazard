test_that("a canonical AVC-style block becomes a hazard() call", {
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=AVCS CONSERVE P OUTHAZ=EX.HZD",
    "STEEPEST QUASI CONDITION=14 MI=200;",
    "EVENT DEAD; TIME INT_DEAD;",
    "PARMS MUE=0.2361727 THALF=0.1512095 NU=1.438652 M=1 FIXM MUC=0.0005436977; );"
  ))
  b <- .hzr_sas_blocks(txt)[[1L]]
  got <- .hzr_parse_hazard(b)

  expect_equal(got$outhaz, "EX.HZD")
  expect_equal(got$call[["data"]], as.name("AVCS"))
  expect_equal(got$call[["time"]], as.name("INT_DEAD"))
  # Only what hazard() reads. CONDITION and QUASI used to be emitted as
  # `condition` and `method`, which nothing reads: the translator counted two
  # options as mapped while they did nothing (#384).
  expect_equal(got$call[["control"]],
               quote(list(maxit = 200, conserve = TRUE)))
  # STEEPEST has no R equivalent and must be surfaced, not dropped.
  expect_true("STEEPEST" %in% got$untranslated$construct)
})

test_that("NOCONSERVE is emitted explicitly rather than defaulted", {
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A NOCONSERVE CONDITION=14;",
    "EVENT D; TIME T; PARMS MUE=1 THALF=1 NU=1; );"
  ))
  got <- .hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]])
  expect_equal(got$call[["control"]][["conserve"]], FALSE)
})

test_that("a LATE statement with comma-separated VAR=VALUE operands parses", {
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A;",
    "EVENT D; TIME T;",
    "PARMS MUL=0.01 TAU=2 GAMMA=1.5;",
    "LATE NOPREVTE=2.653528, NOTEE=-0.3680751; );"
  ))
  got <- .hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]])
  expect_equal(
    got$call[["phases"]],
    quote(list(hzr_phase("g3", tau = 2, gamma = 1.5, alpha = 1, eta = 2,
                          formula = ~NOPREVTE + NOTEE)))
  )
  # log_mu, log_tau, gamma, alpha (defaulted to PROC HAZARD's 1), eta
  # (defaulted to PROC HAZARD's 2, NOT hzr_phase()'s 1), then the two covariate
  # starts in the order they appear on the LATE statement. Compare evaluated: a
  # parsed `-0.3680751` literal is a unary-minus call, not the plain double the
  # translator builds, so quote()-comparing the call would spuriously differ in
  # representation.
  expect_equal(
    eval(got$call[["theta"]]),
    c(log(0.01), log(2), 1.5, 1, 2, 2.653528, -0.3680751)
  )
})

test_that("a non-numeric CONDITION is recorded, not coerced to NA", {
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=ABC;",
    "EVENT D; TIME T; PARMS MUE=1 THALF=1 NU=1; );"
  ))
  got <- expect_silent(.hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]]))
  expect_true("CONDITION" %in% got$untranslated$construct)
  expect_null(got$call[["control"]][["condition"]])
})

test_that("a non-numeric MI is recorded, not coerced to NA", {
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A MI=XYZ CONDITION=14;",
    "EVENT D; TIME T; PARMS MUE=1 THALF=1 NU=1; );"
  ))
  got <- expect_silent(.hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]]))
  expect_true(any(grepl("MI|MAXITER", got$untranslated$construct)))
  expect_null(got$call[["control"]][["maxit"]])
  # The option alongside it is still handled: CONDITION has no R equivalent,
  # so it is recorded rather than dropped (#384).
  expect_true("CONDITION" %in% got$untranslated$construct)
})

test_that("CONDITION and QUASI are recorded, never emitted into control (#384)", {
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=14 QUASI MI=50 CONSERVE;",
    "EVENT D; TIME T; PARMS MUE=1 THALF=1 NU=1; );"
  ))
  got <- .hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]])
  ctl <- as.list(got$call[["control"]])[-1L]
  # Every emitted name is one the fitter reads (the census in #384).
  read_by_fitter <- c("maxit", "reltol", "abstol", "n_starts", "conserve",
                      "phase_share_tol", "start_seed")
  expect_true(all(names(ctl) %in% read_by_fitter), info = toString(names(ctl)))
  expect_null(ctl$condition)
  expect_null(ctl$method)
  u <- got$untranslated
  # CONDITION= stops SAS's optimizer as ill-conditioned once log10 of the
  # Hessian approximation's condition estimate exceeds it (setopt.c:452-456);
  # hazard() has no such stop, only a warning on the final Hessian.
  expect_equal(sum(u$construct == "CONDITION"), 1L)
  expect_match(u$reason[u$construct == "CONDITION"], "setopt.c:452-456",
               fixed = TRUE)
  # QUASI chooses SAS's optimizer; hazard() has no choice to make.
  expect_equal(sum(u$construct == "QUASINEWTON"), 1L)
  expect_match(u$reason[u$construct == "QUASINEWTON"], "BFGS", fixed = TRUE)
  # Not "L-BFGS-B when bounded": every .hzr_optim_generic() caller passes
  # use_bounds = FALSE, so no fit hazard() runs takes that branch.
  expect_no_match(u$reason[u$construct == "QUASINEWTON"], "L-BFGS-B", fixed = TRUE)
  expect_lte(got$tokens_mapped, got$tokens_seen)
})

test_that("tokens_mapped never exceeds tokens_seen when values are bad", {
  txt <- .hzr_sas_normalise(paste(
    "%HAZARD( PROC HAZARD DATA=A CONDITION=ABC MI=XYZ STEEPEST;",
    "EVENT D; TIME T; PARMS MUE=1 THALF=1 NU=1; );"
  ))
  got <- .hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]])
  expect_lte(got$tokens_mapped, got$tokens_seen)
  expect_gte(nrow(got$untranslated), 3L)
})

test_that("a CONDITION= PROC HAZARD ignores is not described as a stop (#384)", {
  # hazpprc.c:48-56 stores CONDITION only for 3 <= n <= 14; outside that the
  # limit stays at the 0 stmtprc.c:74 set, setopt.c:454 skips the test, and
  # the built-in thresholds apply. Saying "CONDITION=20 stops the optimizer"
  # would name a cause that never fires (the #387 message shape).
  reason_for <- function(val) {
    txt <- .hzr_sas_normalise(paste0(
      "%HAZARD( PROC HAZARD DATA=A CONDITION=", val, ";",
      "EVENT D; TIME T; PARMS MUE=1 THALF=1 NU=1; );"))
    u <- .hzr_parse_hazard(.hzr_sas_blocks(txt)[[1L]])$untranslated
    u$reason[u$construct == "CONDITION"]
  }
  inside <- reason_for(14)
  expect_match(inside, "setopt.c:452-456", fixed = TRUE)
  for (val in c(2, 20)) {
    r <- reason_for(val)
    expect_match(r, "outside the 3 to 14 PROC HAZARD accepts", fixed = TRUE, info = as.character(val))
    expect_match(r, "hazpprc.c:48-56", fixed = TRUE, info = as.character(val))
    expect_no_match(r, "stops PROC HAZARD's optimizer", info = as.character(val))
  }
})
