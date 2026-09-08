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

test_that("WEIBULL leaves unspecified ALPHA and ETA at their defaults, free", {
  # PARMS ... WEIBULL is setopt(6) -> SETG3_weibull (setg3.c:427), the
  # generalized Weibull, which admits all positive parameter values. With
  # ALPHA and ETA absent from PARMS they take hzr_phase()'s defaults of 1 --
  # so they are omitted from the emitted call -- and stay estimated. Only the
  # explicit FIXTAU/FIXGAMMA pin anything.
  #
  # Note those are *R's* defaults, not SAS's: stmtprc.c:30-37 starts an
  # unspecified late phase at gamma = 1, alpha = 1, eta = 2 (and tau at
  # 2*Tmax/3, data-dependent). This test pins the translator's behaviour, not
  # start-value parity with PROC HAZARD -- a separate, pre-existing gap that
  # WEIBULL jobs now share with every other path.
  #
  # This test previously asserted alpha = eta = 1, fixed. That was the
  # translator's behaviour, not SAS's: it read WEIBULL as a constraint to the
  # alpha = eta = 1 special case. The G3-collapses-to-Weibull identity at
  # alpha = eta = 1 is real and is covered by
  # test-g3-weibull-correspondence.R, but it is not what the WEIBULL keyword
  # requests.
  ops <- c("MUL=0.01", "TAU=2", "GAMMA=1.5", "WEIBULL", "FIXTAU", "FIXGAMMA")
  got <- .hzr_parse_parms(ops)
  expect_equal(
    got$phases,
    quote(list(
      hzr_phase("g3", tau = 2, gamma = 1.5, fixed = c("tau", "gamma"))
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

test_that("the DELTA reason string does not move with OutDec", {
  # The reason lands in the untranslated frame and is grepped by callers and
  # by the tests above, so it is data rather than only a message. format()
  # honours getOption("OutDec"); a session with OutDec = "," would write
  # "DELTA = 0,5" and break every grep silently.
  old <- options(OutDec = ",")
  on.exit(options(old), add = TRUE)
  r <- .hzr_parse_parms(c("MUE=0.2", "DELTA=0.5"))
  # MUE here has no early shape operand, so it is recorded too; target the
  # DELTA row rather than asserting over the whole frame.
  delta_row <- r$untranslated[grepl("DELTA", r$untranslated$reason), ]
  expect_equal(nrow(delta_row), 1L)
  expect_match(delta_row$reason, "DELTA = 0\\.5", fixed = FALSE)
  expect_false(any(grepl(",", r$untranslated$construct, fixed = TRUE)))
  expect_false(any(grepl("0,5", r$untranslated$reason, fixed = TRUE)))
  expect_false(any(grepl("0,2", r$untranslated$construct, fixed = TRUE)))
})

test_that("WEIBULL keeps the ALPHA and ETA that PARMS specified, both free", {
  # setg3.c:427 SETG3_weibull() is the GENERALIZED Weibull: "we admit all
  # positive values of the parameters". It validates gamma > 0, eta > 0 and
  # alpha >= 0 and bumps g3flag; it never assigns 1 to alpha or eta and never
  # fixes either. Operands are from the production job
  # hz.ce_cardioversion_repeated.ehb.sas, whose listing prints ALPHA and ETA
  # as "Estimated? Yes" with exactly these starting values.
  ops <- c("MUE=0.3012686", "THALF=0.01038402", "NU=0.1708571", "M=5.818869",
           "MUL=0.2450542", "TAU=0.5433813", "ALPHA=2.501719",
           "GAMMA=6.448979", "ETA=0.1365255", "WEIBULL")
  got <- .hzr_parse_parms(ops)
  expect_equal(
    got$phases,
    quote(list(
      hzr_phase("cdf", t_half = 0.01038402, nu = 0.1708571, m = 5.818869),
      hzr_phase("g3", tau = 0.5433813, gamma = 6.448979, alpha = 2.501719,
                eta = 0.1365255)
    ))
  )
  # $phases is what a reader sees; $theta is what reaches the optimizer, and
  # .hzr_parms_theta_block() reads the shape list independently -- so assert
  # it too. These are the starting values the old WEIBULL branch discarded.
  expect_equal(
    got$theta,
    quote(c(log(0.3012686), log(0.01038402), 0.1708571, 5.818869,
            log(0.2450542), log(0.5433813), 6.448979, 2.501719, 0.1365255))
  )
})

test_that("WEIBULL alongside FIXALPHA fixes alpha and only alpha", {
  # The fix must not overshoot: an explicit FIX<param> still pins that one
  # parameter. Distinguishes "WEIBULL fixes nothing" from "nothing is ever
  # fixed on a WEIBULL phase".
  ops <- c("MUL=0.01", "TAU=2", "GAMMA=1.5", "ALPHA=3", "ETA=4",
           "WEIBULL", "FIXALPHA")
  got <- .hzr_parse_parms(ops)
  expect_equal(
    got$phases,
    quote(list(
      hzr_phase("g3", tau = 2, gamma = 1.5, alpha = 3, eta = 4,
                fixed = "alpha")
    ))
  )
})

test_that("a MUL with no late shape operand is recorded, not dropped", {
  # PARMS names a late phase by giving it a scale. Building the phase only
  # when a *shape* operand appeared meant MUL could vanish with the phase --
  # a one-phase R model against SAS's two, reported as fully translated. The
  # starting values SAS would default to here (stmtprc.c: tau = 2*Tmax/3,
  # gamma = 1, alpha = 1, eta = 2) are not this parser's defaults, so the MU
  # is recorded as untranslated rather than guessed at.
  got <- .hzr_parse_parms(c("MUE=0.2", "THALF=1", "NU=1", "MUL=0.05", "WEIBULL"))
  expect_equal(nrow(got$untranslated), 1L)
  expect_equal(got$untranslated$construct, "MUL=0.05")
  expect_match(got$untranslated$reason, "no late phase")
})

test_that("a MUE with no early shape operand is recorded, not dropped", {
  got <- .hzr_parse_parms(c("MUE=0.2", "MUL=0.05", "TAU=2", "GAMMA=1.5"))
  expect_equal(nrow(got$untranslated), 1L)
  expect_equal(got$untranslated$construct, "MUE=0.2")
  expect_match(got$untranslated$reason, "no early phase")
})

test_that("alpha = 1 fixed with GAMMA and ETA both free is recorded", {
  # setg3.c:312-314 routes alpha == 1 && FIXALPHA into SETG3_ignore_tau(),
  # which at setg3.c:405-406 does `if(hzr_parms_ge_estim()) set_fixed(ETA)`:
  # with both shape parameters free, SAS fixes ETA and estimates the product
  # gamma*eta as gamma alone. At alpha = 1, tau = 1 the G3 form collapses to
  # t^(gamma*eta), so gamma and eta are not separately identifiable and SAS is
  # resolving that. hazard() would leave both free and fit the ridge -- one
  # more estimated parameter than PROC HAZARD, with different standard errors.
  # A model difference, not a reporting one, so it is recorded.
  #
  # No corpus job reaches this today (0 of 38 live PARMS blocks, measured
  # 2026-09-08); the 16 that fire SETG3_ignore_tau() all fix GAMMA.
  got <- .hzr_parse_parms(c("MUL=0.01", "TAU=1", "ALPHA=1", "GAMMA=2",
                            "ETA=3", "FIXALPHA"))
  expect_equal(nrow(got$untranslated), 1L)
  expect_match(got$untranslated$reason, "estimates GAMMA\\*ETA")
})

test_that("alpha = 1 fixed with GAMMA fixed is not recorded", {
  # The dominant corpus shape -- gamma fixed, eta free. SAS reparameterises to
  # eta <- gamma*eta, gamma <- 1, which with gamma = 1 is an identity. Nothing
  # to report. Distinguishes the guard above from "any fixed alpha = 1".
  got <- .hzr_parse_parms(c("MUL=0.01", "TAU=1", "ALPHA=1", "GAMMA=1",
                            "ETA=1.32", "FIXALPHA", "FIXGAMMA", "WEIBULL"))
  expect_equal(nrow(got$untranslated), 0L)
})

test_that("ALPHA = 0 under WEIBULL without FIXALPHA is recorded", {
  # setg3.c:437: SETG3_weibull() rejects alpha == 0 when g3flag == 3, i.e.
  # ALPHA=0 under WEIBULL without FIXALPHA, with error SETG3980. The job does
  # not run. hzr_phase() accepts alpha = 0 as the limiting exponential, so
  # without this the translator would emit a fit for a job SAS refuses.
  # ALPHA=0 FIXALPHA is legal (g3flag == 4) and stays unflagged.
  got <- .hzr_parse_parms(c("MUL=0.01", "TAU=2", "GAMMA=1.5", "ALPHA=0",
                            "ETA=1", "WEIBULL"))
  expect_equal(nrow(got$untranslated), 1L)
  expect_match(got$untranslated$reason, "SETG3980")
})

test_that("ALPHA = 0 with FIXALPHA under WEIBULL is not recorded", {
  got <- .hzr_parse_parms(c("MUL=0.01", "TAU=2", "GAMMA=1.5", "ALPHA=0",
                            "ETA=1", "WEIBULL", "FIXALPHA"))
  expect_equal(nrow(got$untranslated), 0L)
})

test_that("the alpha = 1 guard does not fire on a job with no late phase", {
  # alpha_val falls back to hzr_phase()'s default of 1 when PARMS named no
  # ALPHA, which is right for a late phase that exists and wrong for one that
  # does not: an early-only or constant-only job carrying a stray FIXALPHA was
  # flagged about GAMMA and ETA it has no phase for. A false positive on the
  # untranslated frame is not harmless -- that frame is how a caller decides
  # whether a translation is trustworthy.
  expect_equal(nrow(.hzr_parse_parms(c("MUE=0.2", "THALF=1", "NU=1",
                                       "FIXALPHA"))$untranslated), 0L)
  expect_equal(nrow(.hzr_parse_parms(c("MUE=0.2", "THALF=1", "NU=1",
                                       "FIXALPHA", "WEIBULL"))$untranslated), 0L)
  expect_equal(nrow(.hzr_parse_parms(c("MUC=0.01", "FIXALPHA"))$untranslated), 0L)
})

test_that("the alpha = 1 guard still fires when ALPHA was left to default", {
  # The other half of the same boundary: a late phase exists, PARMS never named
  # ALPHA, and FIXALPHA pins it at the default of 1 with GAMMA and ETA free.
  # SAS fixes ETA here, so this must still be recorded -- the scope fix above
  # must not buy its way out of the guard by requiring an explicit ALPHA.
  got <- .hzr_parse_parms(c("MUL=0.01", "TAU=2", "GAMMA=2", "ETA=3", "FIXALPHA"))
  expect_equal(nrow(got$untranslated), 1L)
  expect_match(got$untranslated$reason, "estimates GAMMA\\*ETA")
})
