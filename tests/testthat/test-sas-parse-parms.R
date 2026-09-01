test_that("an early-plus-constant PARMS maps to two phases and a theta", {
  ops <- c("MUE=0.2361727", "THALF=0.1512095", "NU=1.438652", "M=1", "FIXM",
           "MUC=0.0005436977")
  got <- .hzr_parse_parms(ops)
  expect_equal(
    got$phases,
    quote(list(
      hzr_phase("cdf", t_half = 0.1512095, nu = 1.438652, m = 1,
                fixed = "m"),
      hzr_phase("constant")
    ))
  )
  expect_equal(
    got$theta,
    quote(c(log(0.2361727), log(0.1512095), 1.438652, 1, log(0.0005436977)))
  )
})

test_that("WEIBULL becomes a G3 constrained at alpha = 1, eta = 1", {
  # PARMS ... WEIBULL is setopt(6) -> SETG3_weibull, g3flag += 2. The R general
  # form (((t/tau)^gamma + 1)^(1/alpha) - 1)^eta collapses at alpha = eta = 1
  # to (t/tau)^gamma, a Weibull cumulative hazard. See spec 7.1.
  ops <- c("MUL=0.01", "TAU=2", "GAMMA=1.5", "WEIBULL", "FIXTAU", "FIXGAMMA")
  got <- .hzr_parse_parms(ops)
  expect_equal(
    got$phases,
    quote(list(
      hzr_phase("g3", tau = 2, gamma = 1.5, alpha = 1, eta = 1,
                fixed = c("tau", "gamma", "alpha", "eta"))
    ))
  )
})

test_that("phase covariates become a formula on the owning phase", {
  ops <- c("MUE=0.2", "THALF=1", "NU=1", "MUC=0.001")
  got <- .hzr_parse_parms(ops, covars = list(early = c("AGE", "SEX"),
                                             constant = "AGE"))
  expect_equal(
    got$phases,
    quote(list(
      hzr_phase("cdf", t_half = 1, nu = 1, formula = ~AGE + SEX),
      hzr_phase("constant", formula = ~AGE)
    ))
  )
})

test_that("an unrecognised PARMS token is recorded, never dropped", {
  got <- .hzr_parse_parms(c("MUE=0.2", "THALF=1", "NU=1", "FIXGAE2"))
  expect_equal(got$untranslated$construct, "FIXGAE2")
  expect_equal(nrow(got$untranslated), 1L)
})

test_that("a raw EARLY VAR=value operand list produces a formula from names only", {
  # Real syntax: phasevaropt : phasevar phaseval phaseoptspec -- comma-
  # separated VAR=startvalue pairs, not bare names.
  ops <- c("MUE=0.2", "THALF=1", "NU=1")
  got <- .hzr_parse_parms(
    ops, covars = list(early = "NYHA=1.121142, I_PATH=0.9513664")
  )
  expect_equal(
    got$phases,
    quote(list(hzr_phase("cdf", t_half = 1, nu = 1, formula = ~NYHA + I_PATH)))
  )
})

test_that("a non-numeric phase-statement value is untranslated, not guessed", {
  ops <- c("MUE=0.2", "THALF=1", "NU=1")
  got <- .hzr_parse_parms(ops, covars = list(early = "NOBS=NUM"))
  expect_equal(
    got$phases,
    quote(list(hzr_phase("cdf", t_half = 1, nu = 1)))
  )
  expect_true("NOBS=NUM" %in% got$untranslated$construct)
})

test_that("a phase '/ options' tail is untranslated, not parsed", {
  ops <- c("MUE=0.2", "THALF=1", "NU=1")
  got <- .hzr_parse_parms(ops, covars = list(early = "AGE=1.2 / EXCLUDE=(SEX)"))
  expect_equal(
    got$phases,
    quote(list(hzr_phase("cdf", t_half = 1, nu = 1, formula = ~AGE)))
  )
  expect_true(any(grepl("phase options", got$untranslated$reason, fixed = TRUE)))
})

test_that("phase covariate starting values map into theta, in covariate order", {
  ops <- c("MUE=0.2", "THALF=1", "NU=1")
  got <- .hzr_parse_parms(ops, covars = list(
    early = "NYHA=1.121142, I_PATH=0.9513664, INC_SURG=1.375285"
  ))
  # log_mu, log_t_half, nu, m (defaulted, PARMS gave none), then the three
  # covariate starts in the order they appear on the EARLY statement.
  expect_equal(
    got$theta,
    quote(c(log(0.2), log(1), 1, 0, 1.121142, 0.9513664, 1.375285))
  )
  expect_false(any(grepl(
    "not yet mapped to theta", got$untranslated$reason, fixed = TRUE
  )))
})

# DELTA is unimplemented, not "absorbed by decompos()" (#181). Two comments in
# the package said absorbed, which made the omission look deliberate and
# harmless. It is neither: DELTA enters rho, the time argument and the density
# Jacobian separately, and R computes the delta = 0 branch of all three, so a
# job with DELTA != 0 is reproduced against a DIFFERENT function.

test_that("PARMS DELTA = 0 is a faithful translation, not a gap", {
  # R implements exactly this branch, so flagging it would be a false alarm --
  # and a note that fires on the safe case teaches the reader to ignore it on
  # the unsafe one.
  r <- .hzr_parse_parms(c("MUE=0.2", "THALF=0.15", "DELTA=0"))
  expect_equal(nrow(r$untranslated), 0L)
  r2 <- .hzr_parse_parms(c("MUE=0.2", "THALF=0.15", "FIXDELTA"))
  expect_equal(nrow(r2$untranslated), 0L)
})

test_that("PARMS DELTA != 0 is flagged as a different model, not a missing keyword", {
  r <- .hzr_parse_parms(c("MUE=0.2", "THALF=0.15", "DELTA=0.5"))
  expect_equal(nrow(r$untranslated), 1L)
  expect_identical(r$untranslated$construct, "DELTA=0.5")
  # The reason must say the emitted call fits something else. Before this, both
  # DELTA=0 and DELTA=0.5 produced the identical generic string "PARMS keyword
  # has no phase target", so the note distinguished nothing.
  expect_match(r$untranslated$reason, "DIFFERENT model")
  expect_match(r$untranslated$reason, "0\\.5")
  expect_false(grepl("no phase target", r$untranslated$reason))
})

test_that("a non-zero DELTA is flagged whatever else the PARMS carries", {
  # Negative, and alongside a FIXDELTA that pins it there.
  r <- .hzr_parse_parms(c("MUE=0.2", "DELTA=-0.25", "FIXDELTA", "NU=1.4"))
  expect_equal(nrow(r$untranslated), 1L)
  expect_match(r$untranslated$reason, "DIFFERENT model")
  # The rest of the statement still translates -- this is a flag, not a refusal.
  expect_true(r$has_phases)
})
