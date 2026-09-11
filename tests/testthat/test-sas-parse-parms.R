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
  # ALPHA and ETA absent from PARMS they take PROC HAZARD's defaults
  # (stmtprc.c:30-37: alpha = 1, eta = 2) and stay estimated. Only the
  # explicit FIXTAU/FIXGAMMA pin anything.
  #
  # Those defaults are emitted explicitly rather than left to hzr_phase(),
  # whose own eta default is 1. Omitting them would put the printed call and
  # the theta vector into disagreement the moment the two default tables
  # differ, which is what this test guards: assert $theta as well as $phases.
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
      hzr_phase("g3", tau = 2, gamma = 1.5, alpha = 1, eta = 2,
                fixed = c("tau", "gamma"))
    ))
  )
  expect_equal(got$theta, quote(c(log(0.01), log(2), 1.5, 1, 2)))
})

test_that("phase covariates become a formula on the owning phase", {
  ops <- c("MUE=0.2", "THALF=1", "NU=1", "MUC=0.001")
  got <- .hzr_parse_parms(ops, covars = list(early = c("AGE", "SEX"),
                                             constant = "AGE"))
  expect_equal(
    got$phases,
    quote(list(
      hzr_phase("cdf", t_half = 1, nu = 1, m = 1, formula = ~AGE + SEX),
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
    quote(list(hzr_phase("cdf", t_half = 1, nu = 1, m = 1,
                         formula = ~NYHA + I_PATH)))
  )
})

test_that("a non-numeric phase-statement value is untranslated, not guessed", {
  ops <- c("MUE=0.2", "THALF=1", "NU=1")
  got <- .hzr_parse_parms(ops, covars = list(early = "NOBS=NUM"))
  expect_equal(
    got$phases,
    quote(list(hzr_phase("cdf", t_half = 1, nu = 1, m = 1)))
  )
  expect_true("NOBS=NUM" %in% got$untranslated$construct)
})

test_that("a phase '/ options' tail is untranslated, not parsed", {
  ops <- c("MUE=0.2", "THALF=1", "NU=1")
  got <- .hzr_parse_parms(ops, covars = list(early = "AGE=1.2 / EXCLUDE=(SEX)"))
  expect_equal(
    got$phases,
    quote(list(hzr_phase("cdf", t_half = 1, nu = 1, m = 1, formula = ~AGE)))
  )
  expect_true(any(grepl("phase options", got$untranslated$reason, fixed = TRUE)))
})

test_that("phase covariate starting values map into theta, in covariate order", {
  ops <- c("MUE=0.2", "THALF=1", "NU=1")
  got <- .hzr_parse_parms(ops, covars = list(
    early = "NYHA=1.121142, I_PATH=0.9513664, INC_SURG=1.375285"
  ))
  # log_mu, log_t_half, nu, m (defaulted to PROC HAZARD's 1, PARMS gave none),
  # then the three covariate starts in the order they appear on the EARLY
  # statement.
  expect_equal(
    got$theta,
    quote(c(log(0.2), log(1), 1, 1, 1.121142, 0.9513664, 1.375285))
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

test_that("alpha = 1 fixed with GAMMA and ETA both free fixes ETA, as SAS does", {
  # setg3.c:312-314 routes alpha == 1 && FIXALPHA into SETG3_ignore_tau(),
  # which at setg3.c:405 does `if(hzr_parms_ge_estim()) set_fixed(ETA)`: with
  # both shape parameters free, SAS fixes ETA. At alpha = 1, tau = 1 the G3
  # form collapses to t^(gamma*eta), so the pair is not separately
  # identifiable and SAS is resolving that.
  #
  # This was previously RECORDED rather than mirrored, which left hazard()
  # fitting the flat ridge with one more free parameter than PROC HAZARD. It
  # is now mirrored: unlike SETG3_verify_ge_2()'s GAMMA rewrite, this is exact
  # algebra rather than a numerical-branch restriction, so mirroring it costs
  # none of the general G3 shape hzr_phase() carries -- that shape is
  # genuinely degenerate at alpha = 1.
  got <- .hzr_parse_parms(c("MUL=0.01", "TAU=1", "ALPHA=1", "GAMMA=2",
                            "ETA=3", "FIXALPHA"))
  expect_equal(
    got$phases,
    quote(list(hzr_phase("g3", tau = 1, gamma = 2, alpha = 1, eta = 3,
                         fixed = c("tau", "alpha", "eta"))))
  )
  expect_equal(nrow(got$untranslated), 0L)
})

test_that("fixing ETA at alpha = 1 is what makes GAMMA estimable", {
  # Two-sided, because "the standard error is finite" alone could come from
  # anything: fit the SAME phase with ETA left free and show the ridge. At
  # alpha = 1, tau = 1 only gamma * eta is identified, so the free fit must
  # reach the same log-likelihood and the same product, and its gamma must
  # be undetermined: a standard error at least 100 times the fixed fit's
  # (about 0.02), or none at all. It was 7.5 when BFGS stopped partway along
  # the ridge; an optimizer that follows the ridge further, as the SAS/C
  # acceptance test makes it do, reaches a singular Hessian and reports none.
  skip_on_cran()
  got <- .hzr_parse_parms(c("MUL=0.01", "TAU=1", "ALPHA=1", "GAMMA=2",
                            "ETA=3", "FIXALPHA"))
  set.seed(1)
  n <- 400
  d <- data.frame(t = stats::rexp(n, 0.3), s = rep(c(1, 0), length.out = n))
  fit_one <- function(phases) {
    suppressWarnings(hazard(
      time = d$t, status = d$s, dist = "multiphase",
      phases = phases, theta = eval(got$theta), fit = TRUE
    ))
  }
  se_gamma <- function(f) {
    v <- stats::vcov(f)
    if (is.matrix(v)) sqrt(diag(v))[["phase_1.gamma"]] else NA_real_
  }
  product <- function(f) {
    f$fit$theta[["phase_1.gamma"]] * f$fit$theta[["phase_1.eta"]]
  }
  fixed <- fit_one(eval(got$phases))
  free <- fit_one(list(hzr_phase("g3", tau = 1, gamma = 2, alpha = 1,
                                 eta = 3, fixed = c("tau", "alpha"))))
  expect_true(is.finite(se_gamma(fixed)))
  expect_true(!is.finite(se_gamma(free)) ||
                se_gamma(free) > 100 * se_gamma(fixed))
  expect_lt(abs(free$fit$objective - fixed$fit$objective), 1e-4)
  expect_lt(abs(product(free) / product(fixed) - 1), 1e-3)
})

test_that("alpha = 1 fixed with GAMMA fixed is not recorded", {
  # The dominant corpus shape -- gamma fixed, eta free. SAS reparameterises to
  # eta <- gamma*eta, gamma <- 1, which with gamma = 1 is an identity. Nothing
  # to report. Distinguishes the guard above from "any fixed alpha = 1".
  got <- .hzr_parse_parms(c("MUL=0.01", "TAU=1", "ALPHA=1", "GAMMA=1",
                            "ETA=1.32", "FIXALPHA", "FIXGAMMA", "WEIBULL"))
  expect_false(any(grepl("estimates GAMMA", got$untranslated$reason)))
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
  #
  # Asserts that no row blames SETG3_ignore_tau(), not that there are no rows:
  # a stray FIXALPHA on a job with no late phase IS dropped material and earns
  # its own row (PROC HAZARD clears its status flag at stmtprc.c:118-121). An
  # nrow == 0 proxy would fail on that unrelated and correct row.
  no_setg3 <- function(ops) {
    u <- .hzr_parse_parms(ops)$untranslated
    expect_false(any(grepl("SETG3|GAMMA\\*ETA|estimates GAMMA", u$reason)))
  }
  no_setg3(c("MUE=0.2", "THALF=1", "NU=1", "FIXALPHA"))
  no_setg3(c("MUE=0.2", "THALF=1", "NU=1", "FIXALPHA", "WEIBULL"))
  no_setg3(c("MUC=0.01", "FIXALPHA"))
})

test_that("ETA is fixed even when ALPHA was left to its default", {
  # The other half of the same boundary: a late phase exists, PARMS never named
  # ALPHA, and FIXALPHA pins it at the default of 1 with GAMMA and ETA free.
  # SAS fixes ETA here too, so the mirror must not buy its way out by
  # requiring an explicit ALPHA operand.
  got <- .hzr_parse_parms(c("MUL=0.01", "TAU=2", "GAMMA=2", "ETA=3", "FIXALPHA"))
  expect_true("eta" %in% eval(got$phases)[[1L]]$fixed)
  # TAU=2 is discarded by the same branch, which is the one row here.
  expect_equal(got$untranslated$construct, "TAU=2")
})

test_that("the alpha = 1 guard does not fire without MUL", {
  # A phase is active in PROC HAZARD iff its MU is specified and positive.
  # parmprc.c:19 registers MUL through setparmno(22, 7, 3, ...), and
  # setparmno.c:14 sets C->phase[3] = 1 only when the value is > 0. The four
  # shape operands go through setprmf (parmprc.c:20-23), which never touches
  # C->phase[]. With phase 3 off, stmtprc.c:113-122 zeroes TAU/GAMMA/ALPHA/ETA
  # and shape.c:31 never calls SETG3(), so SETG3_ignore_tau() cannot run and
  # there is nothing to warn about.
  #
  # Gating on shape operands instead of MUL got this wrong: TAU/GAMMA/ETA with
  # no MUL is not a late phase at all.
  # Since the MU gate landed, no MUL builds no late phase at all, so there is
  # no `fixed=` to inspect -- the stronger form of the same claim.
  got <- .hzr_parse_parms(c("TAU=2", "GAMMA=2", "ETA=3", "FIXALPHA"))
  expect_length(eval(got$phases), 0L)
})

test_that("the ALPHA = 0 WEIBULL guard does not fire without MUL", {
  # Same gate: with no MUL there is no phase 3, shape.c:31 never reaches
  # SETG3(), and SETG3_weibull() cannot raise SETG3980. Claiming PROC HAZARD
  # refuses the job would be wrong.
  got <- .hzr_parse_parms(c("TAU=2", "GAMMA=1.5", "ALPHA=0", "ETA=1", "WEIBULL"))
  expect_false(any(grepl("SETG3980", got$untranslated$reason)))
})

test_that("both late-phase guards still fire when MUL is present", {
  # The other side of the gate, so it cannot be satisfied by never firing.
  g1 <- .hzr_parse_parms(c("MUL=0.01", "TAU=2", "GAMMA=2", "ETA=3", "FIXALPHA"))
  expect_true("eta" %in% eval(g1$phases)[[1L]]$fixed)
  g2 <- .hzr_parse_parms(c("MUL=0.01", "TAU=2", "GAMMA=1.5", "ALPHA=0",
                           "ETA=1", "WEIBULL"))
  expect_true(any(grepl("SETG3980", g2$untranslated$reason)))
})

# ---------------------------------------------------------------------------
# Shape defaults for a partially specified PARMS block (PROC HAZARD's, not
# hzr_phase()'s). stmtprc.c:30-37 starts an unspecified early phase at
# thalf 1, nu 2, m 1 and an unspecified late phase at gamma 1, alpha 1,
# eta 2; hzr_phase() defaults nu 1, m 0, eta 1. theta is what reaches the
# optimizer, so mirroring R's defaults started a partially specified job
# somewhere PROC HAZARD would not have, and the multiphase likelihood is
# multimodal -- a different start is a different answer.
#
# No public-corpus PARMS block is partially specified (0 of 38 measured
# 2026-09-08; every live block gives MUE + THALF + NU + M together). That is
# NOT the same as "nothing exercises this": test-sas-translate-fits.R's
# end-to-end fit uses `MUE THALF NU MUC` with no M, so its early phase moved
# from m = 0 to m = 1 under this change. It needed no edit only because its
# assertions are on class, dist and theta LENGTH -- none of which can see a
# start-value change -- which is the reason these tests carry the burden of
# the behaviour, not the absence of a caller.
# ---------------------------------------------------------------------------

test_that("an unspecified early NU and M start at PROC HAZARD's values", {
  got <- .hzr_parse_parms(c("MUE=0.2", "THALF=0.5"))
  expect_equal(
    got$phases,
    quote(list(hzr_phase("cdf", t_half = 0.5, nu = 2, m = 1)))
  )
  expect_equal(got$theta, quote(c(log(0.2), log(0.5), 2, 1)))
})

test_that("an unspecified late GAMMA, ALPHA and ETA start at SAS's values", {
  got <- .hzr_parse_parms(c("MUL=0.01", "TAU=3"))
  expect_equal(
    got$phases,
    quote(list(hzr_phase("g3", tau = 3, gamma = 1, alpha = 1, eta = 2)))
  )
  expect_equal(got$theta, quote(c(log(0.01), log(3), 1, 1, 2)))
})

test_that("the emitted defaults are SAS's and differ from hzr_phase()'s", {
  # Without this the two tests above would keep passing if the values were
  # silently re-derived from hzr_phase()'s formals -- which is the state this
  # change exists to leave. Reads the formals rather than restating them, so
  # it fails loudly if hzr_phase()'s own defaults are ever changed to match
  # (at which point emitting them explicitly stops being load bearing and
  # this whole block should be revisited, not quietly deleted).
  f <- formals(hzr_phase)
  expect_false(identical(eval(f$nu), unname(.hzr_parms_early_sas_default[["nu"]])))
  expect_false(identical(eval(f$m), unname(.hzr_parms_early_sas_default[["m"]])))
  expect_false(identical(eval(f$eta), unname(.hzr_parms_late_sas_default[["eta"]])))
  expect_equal(.hzr_parms_early_sas_default, c(t_half = 1, nu = 2, m = 1))
  expect_equal(.hzr_parms_late_sas_default,
               c(tau = 1, gamma = 1, alpha = 1, eta = 2))
})

test_that("the emitted phases call and the theta vector agree, defaults included", {
  # The failure this change had to avoid: filling theta from SAS's defaults
  # while the hzr_phase() call kept omitting them, so the printed call would
  # describe a fit that did not happen. Asserting both literals side by side
  # would share the generator's assumptions, so EXECUTE both instead --
  # evaluate the emitted call and rebuild theta from the resulting phase
  # objects through .hzr_phase_start(), the same route hazard() takes when no
  # theta is supplied. A disagreement of any kind fails here.
  cases <- list(
    c("MUE=0.2", "THALF=0.5"),
    c("MUE=0.2", "THALF=0.5", "NU=1.4"),
    c("MUL=0.01", "TAU=3"),
    c("MUL=0.01", "TAU=3", "GAMMA=2.5", "ALPHA=0.5", "ETA=1.25"),
    c("MUE=0.3", "THALF=0.1", "M=4", "MUC=0.002", "MUL=0.05", "TAU=2",
      "ETA=0.5")
  )
  mu_of <- function(ops, key) {
    hit <- grep(paste0("^", key, "="), ops, value = TRUE)
    if (!length(hit)) return(NULL)
    as.numeric(sub("^[^=]+=", "", hit))
  }
  for (ops in cases) {
    got <- .hzr_parse_parms(ops)
    phases <- eval(got$phases)
    expect_gt(length(phases), 0L)
    mus <- Filter(Negate(is.null),
                  lapply(c("MUE", "MUC", "MUL"), function(k) mu_of(ops, k)))
    expect_equal(length(mus), length(phases))
    from_phases <- unlist(Map(
      function(ph, mu) .hzr_phase_start(ph, n_covariates = 0L, mu_start = mu),
      phases, mus
    ))
    expect_equal(eval(got$theta), from_phases, info = paste(ops, collapse = " "))
  }
})

test_that("an incomplete shape list fails loudly in the theta block", {
  # .hzr_parms_theta_block()'s roxygen claims completeness is asserted rather
  # than defaulted. Reachable only by calling it directly, since the one
  # in-package caller always passes a .hzr_parms_fill_shape() result -- so
  # without this the claim is untested and the stopifnot() has never been seen
  # to fire.
  expect_error(
    .hzr_parms_theta_block("late", 0.1, list(tau = 1), numeric(0))
  )
  expect_error(
    .hzr_parms_theta_block("early", 0.1, list(t_half = 1, nu = 2), numeric(0))
  )
})

test_that("a late phase with no TAU records the data-dependent SAS default", {
  # An UNSPECIFIED TAU is 0.75*Tmax (readobs.c:153-154), not setg3.c:317's
  # 2*Tmax/3 -- readobs() has already replaced it before SETG3 runs, so the
  # tau <= 0 branch there is never reached for an absent TAU. This assertion
  # named the wrong rule and the wrong number until that was checked; the
  # explicit-TAU=0 case below is the one setg3.c:317 governs.
  got <- .hzr_parse_parms(c("MUL=0.01", "GAMMA=2"))
  expect_equal(nrow(got$untranslated), 1L)
  expect_equal(got$untranslated$construct, "TAU (unspecified)")
  expect_match(got$untranslated$reason, "0\\.75\\*Tmax")
  expect_false(grepl("2*Tmax/3", got$untranslated$reason, fixed = TRUE))
  expect_equal(got$phases,
               quote(list(hzr_phase("g3", tau = 1, gamma = 2, alpha = 1,
                                    eta = 2))))
})

test_that("an explicit TAU = 0 takes SETG3's own branch and is recorded", {
  # setg3.c:316-317's predicate is tau <= 0, and only a TAU the job WROTE can
  # still be non-positive there -- readobs.c:153-154 has already replaced an
  # absent one. hzr_phase() would reject tau = 0 outright, so without this the
  # translator emitted a call that errored where SAS supplies a default and
  # runs. Different rule and different constant from the absent case above.
  got <- .hzr_parse_parms(c("MUL=0.01", "TAU=0", "GAMMA=2"))
  expect_equal(got$untranslated$construct, "TAU=0")
  expect_match(got$untranslated$reason, "2\\*Tmax/3")
  expect_equal(got$phases,
               quote(list(hzr_phase("g3", tau = 1, gamma = 2, alpha = 1,
                                    eta = 2))))
})

test_that("a positive TAU is not recorded", {
  got <- .hzr_parse_parms(c("MUL=0.01", "TAU=2", "GAMMA=2"))
  expect_equal(nrow(got$untranslated), 0L)
})

test_that("SETG3_ignore_tau() pins TAU at 1 AND fixes it, and so does the call", {
  # setg3.c:378-379 is two statements: Late.tau = ONE, then
  # hzr_parm_set_fixed(HZ_TAU). An earlier version of this test read only the
  # first, concluded "exactly what is emitted, so there is nothing to report",
  # and asserted ZERO untranslated rows over a fit that returns no standard
  # errors. Both halves are now mirrored, so the translation is faithful and
  # there is genuinely nothing to report -- but assert the emitted `fixed=`,
  # not just the row count, or this reverts to the assertion that could not
  # fail.
  # WEIBULL matches every late-phase block in the corpus, and keeps this test
  # on the tau pin alone: without it GAMMA*ETA = 1.32 <= 2 also trips the
  # SETG3_verify_ge_2 row, which is a separate divergence tested below.
  got <- .hzr_parse_parms(c("MUL=0.01", "ALPHA=1", "GAMMA=1", "ETA=1.32",
                            "FIXALPHA", "FIXGAMMA", "WEIBULL"))
  expect_equal(
    got$phases,
    quote(list(hzr_phase("g3", tau = 1, gamma = 1, alpha = 1, eta = 1.32,
                         fixed = c("tau", "gamma", "alpha"))))
  )
  expect_equal(nrow(got$untranslated), 0L)
})

test_that("pinning TAU restores standard errors on the aliased log_mu", {
  # The point of the pin, executed rather than argued. With TAU free this fit
  # converges onto the log_mu/log_tau ridge and vcov() gives NA for both; with
  # TAU fixed as SAS fixes it, log_mu is identified again. An assertion on the
  # emitted `fixed=` alone would not show that.
  skip_on_cran()
  got <- .hzr_parse_parms(c("MUL=0.01", "ALPHA=1", "GAMMA=1", "ETA=1.32",
                            "FIXALPHA", "FIXGAMMA", "WEIBULL"))
  set.seed(1)
  n <- 400
  d <- data.frame(t = stats::rexp(n, 0.3), s = rep(c(1, 0), length.out = n))
  fit <- suppressWarnings(hazard(
    time = d$t, status = d$s, dist = "multiphase",
    phases = eval(got$phases), theta = eval(got$theta), fit = TRUE
  ))
  se <- sqrt(diag(stats::vcov(fit)))
  expect_true(is.finite(se[["phase_1.log_mu"]]))
})

test_that("a TAU the ignore_tau branch discards is recorded", {
  # PROC HAZARD overrides the value, and so now does the translator -- but a
  # starting value the job wrote and neither program uses is worth saying out
  # loud, which is what $untranslated is for.
  got <- .hzr_parse_parms(c("MUL=0.01", "TAU=5", "ALPHA=1", "GAMMA=1",
                            "ETA=1.32", "FIXALPHA", "FIXGAMMA", "WEIBULL"))
  expect_equal(got$untranslated$construct, "TAU=5")
  expect_match(got$untranslated$reason, "discards this TAU")
  expect_equal(
    got$phases,
    quote(list(hzr_phase("g3", tau = 1, gamma = 1, alpha = 1, eta = 1.32,
                         fixed = c("tau", "gamma", "alpha"))))
  )
})

test_that("the ignore_tau branch suppresses the 2*Tmax/3 row, not just reorders it", {
  # setg3.c:313-316 is an if/else: a job that fixes ALPHA at 1 never reaches
  # the 2*Tmax/3 assignment, so an absent TAU there is not a data-dependent
  # default. The first cut of this fell through to that row and reported a
  # divergence PROC HAZARD does not have.
  got <- .hzr_parse_parms(c("MUL=0.01", "ALPHA=1", "GAMMA=1", "ETA=1.32",
                            "FIXALPHA", "FIXGAMMA", "WEIBULL"))
  expect_false(any(grepl("Tmax", got$untranslated$reason, fixed = TRUE)))
})

test_that("SETG3_verify_ge_2's GAMMA rewrite is recorded, not mirrored", {
  # setg3.c:907-921 pushes a late phase with GAMMA*ETA <= 2 clear of the
  # boundary -- gamma = 3/eta with neither fixed. The SAS defaults gamma = 1,
  # eta = 2 land on GAMMA*ETA = 2 exactly, so this fires on every defaulted
  # non-WEIBULL late phase: PROC HAZARD would start at gamma = 1.5.
  #
  # It is deliberately NOT mirrored. The constraint keeps PROC HAZARD inside a
  # numerical branch it can evaluate; hzr_decompos_g3() carries the general
  # form and needs no such restriction, and reaching a late shape SAS cannot
  # is a goal here. Recorded so a parity run knows why the starts differ.
  got <- .hzr_parse_parms(c("MUL=0.01", "TAU=3"))
  expect_equal(got$untranslated$construct, "gamma=1 alpha=1 eta=2")
  expect_match(got$untranslated$reason, "optimizes from gamma = 1\\.5")
  expect_equal(
    got$phases,
    quote(list(hzr_phase("g3", tau = 3, gamma = 1, alpha = 1, eta = 2)))
  )
})

test_that("the verify_ge_2 row does not fire above the boundary or on WEIBULL", {
  # GAMMA*ETA > 2 is untouched by setg3.c:907, and the WEIBULL path returns at
  # setg3.c:347 before the sign dispatch reaches SETG3_verify_ge_2() at all --
  # which is why no corpus job trips this (every late block there is WEIBULL).
  above <- .hzr_parse_parms(c("MUL=0.01", "TAU=3", "GAMMA=2", "ETA=2"))
  expect_equal(nrow(above$untranslated), 0L)
  weib <- .hzr_parse_parms(c("MUL=0.01", "TAU=3", "GAMMA=1", "ETA=2", "WEIBULL"))
  expect_equal(nrow(weib$untranslated), 0L)
})

test_that("the product SETG3_verify_ge_2 reads survives the ignore_tau swap", {
  # setg3.c:403-421 moves the exponent between GAMMA and ETA but preserves
  # their PRODUCT, which is all `gte` reads -- so the boundary test is the same
  # whether or not that branch ran.
  #
  # Only while the product is positive, though: setg3.c:411-412 falls back to
  # the other operand when gamma*eta <= 0, and the product then changes. That
  # case is covered by the test below, not here.
  # gamma * eta = 3, above the boundary, in both orderings: 0.5 * 6 straight,
  # and 1 * 3 after the swap. Neither may record the rewrite.
  a <- .hzr_parse_parms(c("MUL=0.01", "TAU=2", "GAMMA=0.5", "ETA=6"))
  expect_equal(nrow(a$untranslated), 0L)
  b <- .hzr_parse_parms(c("MUL=0.01", "TAU=2", "GAMMA=0.5", "ETA=6",
                          "ALPHA=1", "FIXALPHA", "FIXGAMMA"))
  expect_false(any(grepl("optimizes from", b$untranslated$reason)))
  # And the same product BELOW the boundary does record, so the pair above is
  # not passing merely because the guard never fires on these shapes.
  lo <- .hzr_parse_parms(c("MUL=0.01", "TAU=2", "GAMMA=0.5", "ETA=3"))
  expect_true(any(grepl("optimizes from gamma", lo$untranslated$reason)))
})

test_that("the aliasing that motivates pinning TAU is real, not asserted", {
  # The comment justifying the pin claims G3 collapses to (t/tau)^(gamma*eta)
  # at alpha = 1, which is why log_mu and log_tau cannot both be identified.
  # Compute it rather than restate it.
  tt <- c(0.5, 1, 2, 5, 10)
  g <- hzr_decompos_g3(tt, tau = 2, gamma = 1.3, alpha = 1, eta = 1.7)$G
  expect_equal(g, (tt / 2)^(1.3 * 1.7))
  # And it does NOT collapse away from alpha = 1, so the guard's scope is right.
  g2 <- hzr_decompos_g3(tt, tau = 2, gamma = 1.3, alpha = 1.5, eta = 1.7)$G
  expect_false(isTRUE(all.equal(g2, (tt / 2)^(1.3 * 1.7))))
})

test_that("the TAU guard does not fire without a positive MUL", {
  # Same gate as the other two SETG3 guards: with phase 3 off, shape.c never
  # calls SETG3(), so no TAU default is ever applied.
  #
  # MUL=0 deliberately does NOT appear here. It is SAS's own "phase off"
  # marker, but this parser still builds the phase and emits log(0) = -Inf as
  # its starting log_mu with no untranslated row -- a pre-existing gap in the
  # phase-activation asymmetry noted in .hzr_parse_parms(), not something this
  # guard should be read as blessing. Asserting the TAU row is absent there
  # would pass over the -Inf and quietly certify it.
  got <- .hzr_parse_parms(c("GAMMA=2", "ETA=3"))
  expect_false(any(grepl("Tmax", got$untranslated$reason, fixed = TRUE)))
})

test_that("a MUL with no late shape operand records the MU, not TAU as well", {
  # The whole phase is missing and the MUL row already says so; a second row
  # about that phase's TAU would be noise on $untranslated, which is how a
  # caller decides whether a translation can be trusted.
  got <- .hzr_parse_parms(c("MUE=0.2", "THALF=1", "NU=1", "MUL=0.05"))
  expect_equal(nrow(got$untranslated), 1L)
  expect_match(got$untranslated$reason, "no late phase")
})

# ---------------------------------------------------------------------------
# The rest of the SETG3() chain (.hzr_setg3_notes()).
#
# Two kinds of divergence, both recorded and neither mirrored, per the split
# documented on that function: a REFUSAL is a job PROC HAZARD will not run, and
# a REWRITE is SETG3 constraining the late shape to stay inside a numerical
# branch it can evaluate. Only the two alpha = 1 identifiability fixes are
# mirrored, and those are tested above.
#
# None of this fires on the public corpus: every late-phase PARMS block there
# carries WEIBULL, which returns at setg3.c:347 before the sign dispatch, and
# all of them write TAU=1 FIXTAU ALPHA=1 FIXALPHA with a positive GAMMA and
# ETA. These tests therefore carry the whole behaviour.
# ---------------------------------------------------------------------------

test_that("each SETG3 entry refusal is recorded with its own code", {
  # setg3.c:269-284, checked before anything else -- including
  # SETG3_ignore_tau(). A FIX* on an operand SAS reads as unspecified is
  # fatal, and the value SAS reads is the job's own, so FIXTAU with no TAU
  # operand refuses on the stmtprc.c initializer of 0.
  refusal <- function(ops) {
    got <- .hzr_parse_parms(ops)
    expect_equal(nrow(got$untranslated), 1L, info = paste(ops, collapse = " "))
    got$untranslated$reason
  }
  # TAU=0 explicitly, NOT a bare FIXTAU: an unspecified TAU never reaches
  # setg3.c:269 as 0 (see the readobs.c test below), so a bare FIXTAU cannot
  # raise SETG3900. This assertion said otherwise until that was checked.
  expect_match(refusal(c("MUL=0.01", "TAU=0", "FIXTAU", "GAMMA=2", "ETA=2")),
               "SETG3900")
  expect_match(refusal(c("MUL=0.01", "TAU=1", "GAMMA=0", "FIXGAMMA", "ETA=2")),
               "SETG3910")
  expect_match(refusal(c("MUL=0.01", "TAU=1", "GAMMA=2", "ETA=2",
                         "ALPHA=-1", "FIXALPHA")),
               "SETG3920")
  expect_match(refusal(c("MUL=0.01", "TAU=1", "GAMMA=2", "ETA=0", "FIXETA")),
               "SETG3930")
})

test_that("a refusal names the operand, not just the SAS message code", {
  # The code alone is greppable but opaque; a caller reading $untranslated has
  # to know which operand to change.
  got <- .hzr_parse_parms(c("MUL=0.01", "TAU=0", "FIXTAU", "GAMMA=2",
                            "ETA=2"))
  expect_match(got$untranslated$reason, "TAU is fixed at a non-positive value")
  expect_match(got$untranslated$construct, "fixed:tau")
})

test_that("an ENTRY refusal suppresses the TAU row but a later one does not", {
  # setg3.c:269-284 returns before the TAU rules at :309-323, so a TAU row
  # there would describe code PROC HAZARD never reaches. The other twelve
  # codes fire from :429 onwards, by which point the TAU rule HAS run --
  # suppressing the row for those would hide a divergence that happened.
  # Keying the suppression on "refused at all" gets the second case wrong.
  entry <- .hzr_parse_parms(c("MUL=0.01", "TAU=0", "FIXTAU", "GAMMA=2",
                              "ETA=2"))
  expect_equal(nrow(entry$untranslated), 1L)
  expect_false(any(grepl("Tmax", entry$untranslated$reason, fixed = TRUE)))
  late <- .hzr_parse_parms(c("MUL=0.01", "GAMMA=1", "FIXGAMMA", "ETA=2",
                             "FIXETA"))
  expect_true(any(grepl("SETG31020", late$untranslated$reason)))
  expect_true(any(grepl("Tmax", late$untranslated$reason, fixed = TRUE)))
})

test_that("an unspecified TAU is 0.75*Tmax, not the stmtprc.c initializer", {
  # readobs.c:153-154 replaces an unspecified TAU with 0.75*Tmax whenever the
  # late phase is active, and readobs() runs at hazard.c:276, before hzrg()
  # reaches SETG3() at :292. So the stmtprc.c initializer of 0 never arrives
  # at setg3.c:269 -- an absent TAU is POSITIVE there, cannot raise SETG3900,
  # and never reaches the 2*Tmax/3 assignment at :317 either.
  #
  # Conflating "absent" with "non-positive" put a refusal in $untranslated for
  # a job PROC HAZARD runs, which is this package's signature defect inverted:
  # a finding that looks like one and is not.
  absent <- .hzr_parse_parms(c("MUL=0.01", "FIXTAU", "GAMMA=2", "ETA=2"))
  expect_equal(nrow(absent$untranslated), 1L)
  expect_match(absent$untranslated$reason, "0\\.75\\*Tmax")
  expect_false(any(grepl("SETG3900", absent$untranslated$reason)))
  # An explicit non-positive TAU is the other rule, and a different number.
  written <- .hzr_parse_parms(c("MUL=0.01", "TAU=0", "GAMMA=2", "ETA=2"))
  expect_match(written$untranslated$reason, "2\\*Tmax/3")
})

test_that("WEIBULL's own refusals are recorded, including a negative ALPHA", {
  # setg3.c:429-440. The pre-existing guard covered only alpha == 0; a
  # negative ALPHA raises the same SETG3980 and was emitted as a runnable fit.
  neg <- .hzr_parse_parms(c("MUL=0.01", "TAU=2", "GAMMA=1.5", "ALPHA=-2",
                            "ETA=1", "WEIBULL"))
  expect_match(neg$untranslated$reason, "SETG3980")
  g <- .hzr_parse_parms(c("MUL=0.01", "TAU=2", "GAMMA=0", "ETA=1", "WEIBULL"))
  expect_match(g$untranslated$reason, "SETG3960")
  e <- .hzr_parse_parms(c("MUL=0.01", "TAU=2", "GAMMA=1.5", "ETA=0",
                          "WEIBULL"))
  expect_match(e$untranslated$reason, "SETG3970")
})

test_that("ALPHA = 0 is faithful when fixed and a rewrite when free", {
  # setg3.c:332-335: alpha == 0 with FIXALPHA sets g3flag = 2, the limiting
  # exponential, which is exactly what hzr_phase() fits at alpha = 0 -- so
  # that job translates cleanly. Left free, SETG3_alpha_gener() (:854) instead
  # derives alpha = gamma*eta/3, and the emitted call would fit the
  # exponential SAS did not choose. The pair distinguishes the two; asserting
  # only the second would pass for a guard that fired on every ALPHA = 0.
  fixed <- .hzr_parse_parms(c("MUL=0.01", "TAU=1", "GAMMA=2", "ETA=2",
                              "ALPHA=0", "FIXALPHA"))
  expect_equal(nrow(fixed$untranslated), 0L)
  free <- .hzr_parse_parms(c("MUL=0.01", "TAU=1", "GAMMA=2", "ETA=2",
                             "ALPHA=0"))
  expect_match(free$untranslated$reason, "optimizes from alpha = 1\\.33333")
})

test_that("a non-positive GAMMA or ETA is a rewrite SAS makes and R cannot", {
  # setg3.c:533-585 and :587-636 treat a non-positive operand as "derive one
  # for me" -- gamma = 3*alpha/eta, then 3/eta if the product still sits at or
  # below 2. hzr_phase() requires gamma > 0 and would reject the emitted call
  # outright, so this is the one place SAS is the more permissive of the two;
  # recording it says which operand to supply.
  g <- .hzr_parse_parms(c("MUL=0.01", "TAU=1", "GAMMA=0", "ETA=2", "ALPHA=1"))
  expect_match(g$untranslated$reason, "optimizes from gamma = 1\\.5")
  # The claim above, made load-bearing: the emitted call really cannot be
  # built, so the row must not describe the difference as a starting value.
  expect_error(eval(g$phases), "gamma must be a positive scalar")
  expect_match(g$untranslated$reason, "cannot be built at all")
  e <- .hzr_parse_parms(c("MUL=0.01", "TAU=1", "GAMMA=2", "ETA=0", "ALPHA=1"))
  expect_match(e$untranslated$reason, "optimizes from .*eta = 1\\.5")
})

test_that("SETG3_alpha_fixup's refusal and its substitution are distinguished", {
  # setg3.c:838-851: with GAMMA*ETA/ALPHA <= 2 the C either derives
  # alpha = gamma*eta/3, or refuses with SETG31040 when ALPHA is fixed and it
  # cannot. Same operands, one FIXALPHA apart.
  free <- .hzr_parse_parms(c("MUL=0.01", "TAU=2", "GAMMA=2", "ETA=2",
                             "ALPHA=3"))
  expect_match(free$untranslated$reason, "optimizes from alpha = 1\\.33333")
  pinned <- .hzr_parse_parms(c("MUL=0.01", "TAU=2", "GAMMA=2", "ETA=2",
                               "ALPHA=3", "FIXALPHA"))
  expect_match(pinned$untranslated$reason, "SETG31040")
})

test_that("a late shape SETG3 leaves alone records nothing", {
  # The other side of every guard above, so none of them can be satisfied by
  # firing unconditionally: GAMMA*ETA = 4 > 2 and GAMMA*ETA/ALPHA = 4 > 2, so
  # verify_ge_2 and alpha_fixup are both no-ops and the translation is exact.
  got <- .hzr_parse_parms(c("MUL=0.01", "TAU=3", "GAMMA=2", "ETA=2"))
  expect_equal(nrow(got$untranslated), 0L)
})

test_that("the SETG3 trace is gated on an active late phase", {
  # Same gate as every other SETG3 guard: with no positive MUL, shape.c never
  # calls SETG3(), so neither a refusal nor a rewrite can happen.
  #
  # These jobs are NOT silent -- the MU gate records the orphaned operands, and
  # the no-phase job is refused outright -- so this asserts the absence of the
  # SETG3 rows specifically. An nrow == 0 here would be asserting the MU gate's
  # rows away, which is the opposite of what either guard wants.
  setg3_rows <- function(ops) {
    r <- .hzr_parse_parms(ops)$untranslated$reason
    sum(grepl("SETG3|optimizes from", r))
  }
  expect_equal(setg3_rows(c("GAMMA=0", "FIXGAMMA", "ETA=2")), 0L)
  expect_equal(setg3_rows(c("MUL=0", "FIXTAU", "GAMMA=2")), 0L)
  # And the paired active case does produce one, so the helper can fail.
  expect_gt(setg3_rows(c("MUL=0.01", "TAU=2", "GAMMA=1", "ETA=2")), 0L)
})

test_that("the ignore_tau swap is NOT product-preserving when the product is <= 0", {
  # setg3.c:411-412: `if(Late.gamma<=ZERO) Late.gamma = HZRstr.l.eta;`. With
  # gamma = -2 and eta = 3 the C does not carry -6 across -- it keeps eta,
  # leaving gamma = 3, eta = 1 and a product of 3. The trace has to reproduce
  # the fallback rather than the arithmetic, and the comparison has to notice.
  # Stated as an invariant without this exception, "the swap preserves the
  # product" is the kind of half-read claim that put a wrong assertion in this
  # file once already.
  tr <- .hzr_setg3_notes(tau_raw = 1, gamma = -2, alpha = 1, eta = 3,
                         fixed = "alpha", weibull = FALSE)
  expect_equal(unname(tr$shape[["gamma"]]), 3)
  expect_equal(unname(tr$shape[["eta"]]), 1)
  expect_false(isTRUE(all.equal(
    tr$shape[["gamma"]] * tr$shape[["eta"]], -2 * 3
  )))
  # And the divergence is reported rather than absorbed by the product test.
  got <- .hzr_parse_parms(c("MUL=0.01", "TAU=1", "ALPHA=1", "GAMMA=-2",
                            "ETA=3", "FIXALPHA"))
  expect_match(got$untranslated$reason, "optimizes from gamma\\*eta = 3")
  # The positive case still preserves it, so the two halves are distinguished.
  pos <- .hzr_setg3_notes(tau_raw = 1, gamma = 2, alpha = 1, eta = 3,
                          fixed = "alpha", weibull = FALSE)
  expect_equal(unname(pos$shape[["gamma"]] * pos$shape[["eta"]]), 6)
})

test_that("the four remaining sign branches use SETG3's own constants", {
  # setg3.c:638-812. These are the branches with the distinctive substitutions,
  # and none was covered: mutating `eta <- 3` to 2 in all_le_0, or
  # `gamma <- 1.5 * alpha` to `alpha` in alpha_gt_0, produced no failure.
  # Asserted through .hzr_setg3_notes() directly so the values are visible
  # rather than inferred from a reason string.
  fx <- character(0)
  # eta_gt_0 (:638): alpha <= 0, gamma <= 0, eta > 0 -> gamma = 3/eta, then
  # alpha_gener derives alpha = gamma*eta/3 = 1.
  a <- .hzr_setg3_notes(1, gamma = 0, alpha = 0, eta = 3, fixed = fx,
                        weibull = FALSE)
  expect_equal(unname(a$shape[["gamma"]]), 1)
  expect_equal(unname(a$shape[["alpha"]]), 1)
  # gamma_gt_0 (:678): alpha <= 0, gamma > 0, eta <= 0 -> eta = 3/gamma.
  b <- .hzr_setg3_notes(1, gamma = 2, alpha = 0, eta = 0, fixed = fx,
                        weibull = FALSE)
  expect_equal(unname(b$shape[["eta"]]), 1.5)
  # alpha_gt_0 (:717): alpha > 0, gamma <= 0, eta <= 0 -> eta = 2,
  # gamma = 1.5*alpha.
  cc <- .hzr_setg3_notes(1, gamma = 0, alpha = 2, eta = 0, fixed = fx,
                         weibull = FALSE)
  expect_equal(unname(cc$shape[["eta"]]), 2)
  expect_equal(unname(cc$shape[["gamma"]]), 3)
  # all_le_0 (:772): gamma = 1, eta = 3.
  d <- .hzr_setg3_notes(1, gamma = 0, alpha = 0, eta = 0, fixed = fx,
                        weibull = FALSE)
  expect_equal(unname(d$shape[["gamma"]]), 1)
  expect_equal(unname(d$shape[["eta"]]), 3)
})

test_that("the reachable refusal codes are exactly the nine, by exhaustive search", {
  # Sixteen codes are mirrored from the C, but only NINE can fire. The other
  # seven each guard a condition the ENTRY checks at setg3.c:269-284 already
  # refused: SETG31090/32050/33020 want a non-positive GAMMA that is fixed
  # (SETG3910), SETG32020/32080/33010 a non-positive ETA that is fixed
  # (SETG3930), and SETG31070 a fixed ALPHA on a branch where alpha <= 0 --
  # negative refuses at SETG3920, and zero sets g3flag = 2 so alpha_gener is
  # never called. That is a property of PROC HAZARD, not of this port, and it
  # is why they have no individual tests: an input that reaches them does not
  # exist.
  #
  # Searched rather than argued. Reasoning about reachability through three
  # files is how the SETG3900 claim went wrong; a grid says what is true.
  skip_on_cran()
  parms <- c("tau", "gamma", "alpha", "eta")
  subsets <- unlist(lapply(0:4, function(k) utils::combn(parms, k, simplify = FALSE)),
                    recursive = FALSE)
  seen <- character(0)
  for (tau in c(NA, -1, 0, 1, 5)) {
    for (g in c(-2, 0, 1, 2, 3)) {
      for (a in c(-2, 0, 1, 2, 3)) {
        for (e in c(-2, 0, 1, 2, 3)) {
          for (fx in subsets) {
            for (w in c(TRUE, FALSE)) {
              tr <- .hzr_setg3_notes(tau, g, a, e, fx, w)
              if (!is.null(tr$refusal)) seen <- union(seen, tr$refusal)
            }
          }
        }
      }
    }
  }
  expect_setequal(seen, c("(SETG3900)", "(SETG3910)", "(SETG3920)",
                          "(SETG3930)", "(SETG3960)", "(SETG3970)",
                          "(SETG3980)", "(SETG31020)", "(SETG31040)"))
})

test_that("the two non-entry refusals that ARE reachable fire on named inputs", {
  # SETG31020 (verify_ge_2, both operands fixed at or below the boundary) and
  # SETG31040 (alpha_fixup, ALPHA fixed with GAMMA*ETA/ALPHA <= 2) are the only
  # refusals raised after the entry checks. Both need positive operands, which
  # is exactly why they survive where the other seven do not.
  both <- .hzr_setg3_notes(1, gamma = 1, alpha = 1, eta = 2,
                           fixed = c("gamma", "eta"), weibull = FALSE)
  expect_equal(both$refusal, "(SETG31020)")
  expect_false(isTRUE(both$entry))
  pinned <- .hzr_setg3_notes(1, gamma = 2, alpha = 3, eta = 2,
                             fixed = "alpha", weibull = FALSE)
  expect_equal(pinned$refusal, "(SETG31040)")
})

test_that("every code in the refusal table has its own gloss", {
  # The table is kept complete against the C even where a code is unreachable,
  # so nothing may fall through to the generic text.
  for (code in names(.hzr_setg3_refusal)) {
    expect_false(identical(.hzr_setg3_refusal_reason(code),
                           "the operand combination is rejected"),
                 info = code)
  }
})

test_that("an unknown refusal code degrades instead of erroring", {
  # `[[` on an unmatched name in a character vector throws rather than
  # returning NULL, so the is.null() fallback this replaced was hollow: right
  # shape, unreachable, and the parser would have errored if a code were ever
  # added to the trace but not the table.
  expect_equal(.hzr_setg3_refusal_reason("(SETG39999)"),
               "the operand combination is rejected")
})

# ---------------------------------------------------------------------------
# A phase is active iff its MU was specified and positive (setparmno.c:11-14).
# ---------------------------------------------------------------------------

test_that("late shape operands with no MUL build no phase and are recorded", {
  # parmprc.c:19-23 registers MUL through setparmno() and TAU/GAMMA/ALPHA/ETA
  # through setprmf(), and setprmf() never touches C->phase[]. With phase 3
  # off, stmtprc.c:113-122 zeroes all four operands and shape.c:31 never calls
  # SETG3(). So this is a one-phase job in PROC HAZARD, and the operands are
  # discarded -- recorded here rather than dropped silently.
  # M is unspecified, so it starts at PROC HAZARD's default of 1, not
  # hzr_phase()'s 0 -- see the shape-defaults block above.
  got <- .hzr_parse_parms(c("MUE=0.2", "THALF=1", "NU=1", "TAU=2", "GAMMA=1.5"))
  expect_equal(got$phases,
               quote(list(hzr_phase("cdf", t_half = 1, nu = 1, m = 1))))
  expect_equal(got$theta, quote(c(log(0.2), log(1), 1, 1)))
  expect_equal(nrow(got$untranslated), 1L)
  expect_equal(got$untranslated$construct, "TAU=2 GAMMA=1.5")
  expect_match(got$untranslated$reason, "no active MUL")
})

test_that("early shape operands with no MUE build no phase and are recorded", {
  # The mirror of the above through stmtprc.c:101-112 and shape.c:19.
  # GAMMA/ALPHA/ETA are unspecified and start at PROC HAZARD's 1/1/2. That
  # product is exactly 2, so SETG3_verify_ge_2() would rewrite GAMMA -- a
  # second row, deliberately not mirrored (see the SETG3 block above).
  got <- .hzr_parse_parms(c("THALF=1", "NU=1", "MUL=0.05", "TAU=2"))
  expect_equal(
    got$phases,
    quote(list(hzr_phase("g3", tau = 2, gamma = 1, alpha = 1, eta = 2)))
  )
  expect_equal(got$theta, quote(c(log(0.05), log(2), 1, 1, 2)))
  expect_equal(nrow(got$untranslated), 2L)
  expect_true("THALF=1 NU=1" %in% got$untranslated$construct)
  expect_true(any(grepl("no active MUE", got$untranslated$reason)))
})

test_that("a MU of zero is inactive, not merely absent", {
  # setparmno.c:11 gates on stmtfld(parmno) > ZERO, so MUE=0 registers the
  # keyword but leaves C->phase[1] at 0. Absence and non-positivity are the
  # same outcome, and testing only for absence would let this through.
  got <- .hzr_parse_parms(c("MUE=0", "THALF=1", "NU=1", "MUC=0.001"))
  expect_equal(got$phases, quote(list(hzr_phase("constant"))))
  expect_equal(got$theta, quote(c(log(0.001))))
  expect_true(any(grepl("no active MUE", got$untranslated$reason)))
})

test_that("MUC gates the constant phase on positivity, not presence", {
  # shape.c:23 reads Common.phase[2], set by setparmno(21, 6, 2, ...) under the
  # same > ZERO gate as the other two. A presence test would build this phase.
  off <- .hzr_parse_parms(c("MUE=0.2", "THALF=1", "NU=1", "MUC=0"))
  expect_equal(off$phases,
               quote(list(hzr_phase("cdf", t_half = 1, nu = 1, m = 1))))
  expect_equal(off$theta, quote(c(log(0.2), log(1), 1, 1)))
  expect_equal(off$untranslated$construct, "MUC=0")

  # The paired positive case, so the assertion above cannot pass by the phase
  # never being built at all.
  on <- .hzr_parse_parms(c("MUE=0.2", "THALF=1", "NU=1", "MUC=0.001"))
  expect_equal(
    on$phases,
    quote(list(hzr_phase("cdf", t_half = 1, nu = 1, m = 1),
               hzr_phase("constant")))
  )
  expect_equal(on$theta, quote(c(log(0.2), log(1), 1, 1, log(0.001))))
  expect_equal(nrow(on$untranslated), 0L)
})

test_that("theta drops the block of a phase that was not built", {
  # .hzr_phase_theta_names() labels theta by *position*, so a phase that is
  # not built must not leave its block behind: the late block would otherwise
  # be read as the early phase's shape. Asserts the whole vector, not its
  # length -- a length check passes on a block of the wrong contents.
  got <- .hzr_parse_parms(c("THALF=9", "NU=7", "M=5", "MUC=0.001",
                            "MUL=0.05", "TAU=2", "GAMMA=1.5"))
  expect_equal(
    got$phases,
    quote(list(hzr_phase("constant"),
               hzr_phase("g3", tau = 2, gamma = 1.5, alpha = 1, eta = 2)))
  )
  expect_equal(got$theta, quote(c(log(0.001), log(0.05), log(2), 1.5, 1, 2)))
})

test_that("a PARMS with no active MU builds nothing and says so", {
  # The hollow-object case from AGENTS.md: has_phases must be FALSE rather
  # than the caller receiving an empty list() that looks like a translation.
  got <- .hzr_parse_parms(c("THALF=1", "NU=1", "TAU=2"))
  expect_false(got$has_phases)
  expect_equal(got$phases, quote(list()))
  expect_setequal(
    got$untranslated$construct,
    c("THALF=1 NU=1", "TAU=2", "PARMS with no positive MUE, MUC or MUL")
  )
})

test_that("no active MU at all is recorded as the refusal PROC HAZARD raises", {
  # modterm.c:18-22: with phase[1], phase[2] and phase[3] all 0 the reference
  # prints "No phase selected", sets C->errorno = 1001 and the job does not
  # run. Without this row .hzr_parse_job() omits `dist` (sas-parse-job.R:604),
  # the emitted call falls through to hazard()'s default distribution, and
  # hzr_translate_sas() reports full coverage for a job SAS refuses -- an
  # output that looks like a result. Same shape as the SETG3980 guard.
  got <- .hzr_parse_parms(c("MUE=0", "THALF=1"))
  expect_false(got$has_phases)
  expect_true(any(grepl("modterm.c", got$untranslated$reason, fixed = TRUE)))

  # The paired case that must NOT be refused: MUE is positive, so PROC HAZARD
  # does select phase 1 and does run the job -- this parser declines to
  # translate it for a different reason (no shape operand), and saying
  # "no phase selected" there would be a false positive.
  ok <- .hzr_parse_parms(c("MUE=0.2"))
  expect_false(any(grepl("modterm.c", ok$untranslated$reason, fixed = TRUE)))
  expect_equal(ok$untranslated$construct, "MUE=0.2")
})

test_that("a job with no PARMS operands at all is refused, not defaulted", {
  # The wider half of the same modterm.c rule. A PROC HAZARD job carrying no
  # PARMS statement reaches here with operands = character(0), so no
  # setparmno() call ever fires (parmprc.c:13,18,19), all three C->phase[]
  # stay at their stmtprc.c:87 zero, and modterm.c:18-22 raises ERROR 1001.
  # modterm() is universal rather than multiphase-only: its single call site
  # is outmods.c:91, outmods() is unconditional in main (hazard.c:296), and
  # hazard.c:299-302 routes 1001 to hzfxit("SEMANTIC") BEFORE results(). The
  # job never fits, so a translation that emits one is a wrong answer.
  got <- .hzr_parse_parms(character(0))
  expect_false(got$has_phases)
  expect_true(got$refused)
  expect_equal(got$untranslated$construct,
               "no PARMS operands (no MUE, MUC or MUL)")
  expect_true(any(grepl("modterm.c", got$untranslated$reason, fixed = TRUE)))
})

test_that("operands this parser cannot read are recorded but never refused", {
  # `refused` drives a stop() in place of the fit, so it is a claim about what
  # PROC HAZARD does and is sound only when the whole statement was understood.
  # .hzr_parse_hazard() splits on " ", so `PARMS MUE = 0.2` arrives as separate
  # "MUE", "=", "0.2" operands and nothing parses -- while HAZARD's lexer drops
  # whitespace unconditionally (hazard_l.l:32, rule at :50) and RUNS that job
  # with an active early phase. Refusing it would stop a job the reference
  # accepts, so the gate must key on comprehension, not on has_phases.
  spaced <- .hzr_parse_parms(c("MUE", "=", "0.2", "THALF", "=", "1"))
  expect_false(spaced$has_phases)
  expect_false(spaced$refused)
  expect_false(any(grepl("modterm.c", spaced$untranslated$reason, fixed = TRUE)))
  # Not refusing is not the same as declaring it fine -- the operands are still
  # reported, so this cannot pass by the parser having silently accepted them.
  expect_true("MUE" %in% spaced$untranslated$construct)

  # The paired readable case, so the assertions above cannot pass merely
  # because nothing ever refuses: MUE=0 is understood AND selects no phase.
  read <- .hzr_parse_parms(c("MUE=0", "THALF=1"))
  expect_false(read$has_phases)
  expect_true(read$refused)
})

test_that("a PARMS template with `?` placeholders is not refused as no-phase", {
  # The SECOND %HAZARD block of dist/examples/hm.dthar.TGA.sas (line 110,
  # PARMS at 122) carries literal `PARMS MUE=? THALF=? NU=? M=1 FIXM MUC=?;`
  # for the reader to fill in from the stepwise output above it. The file's
  # FIRST block (line 70, PARMS at 74) is fully valid and activates phases 1
  # and 2, so grepping the file for a valid PARMS will find one -- the shape
  # under test here belongs to the second block alone. Every value is non-numeric, so no MU is
  # active -- but that is this parser failing to read the statement, not PROC
  # HAZARD selecting no phase, and the reference would reject `?` as a SYNTAX
  # error long before modterm.c's ERROR 1001. Claiming "No phase selected"
  # here would attribute the wrong refusal to the reference, and now that
  # `refused` emits a stop() it would also replace the fit on that basis.
  got <- .hzr_parse_parms(c("MUE=?", "THALF=?", "M=1", "FIXM", "MUC=?"))
  expect_false(got$has_phases)
  expect_false(got$refused)
  expect_false(any(grepl("modterm.c", got$untranslated$reason, fixed = TRUE)))
  # The unreadable operands are each still reported by name, so the assertions
  # above cannot pass by the whole statement having been silently dropped.
  expect_true(all(c("MUE=?", "THALF=?", "MUC=?") %in% got$untranslated$construct))
})

test_that("covariates of a phase that is not built are recorded, not dropped", {
  # setstat.c:9-12 returns early for a phase whose C->phase[] is 0, so PROC
  # HAZARD drops these too -- the values agree and only silence would be wrong.
  # This is the hole the MU gate opened: with MUL present the covariate
  # survives into the emitted formula, so it must not simply vanish without it.
  got <- .hzr_parse_parms(c("MUE=0.2", "THALF=1", "NU=1", "TAU=2"),
                          covars = list(late = "AGE=1.2, SEX=0.5"))
  expect_equal(got$phases, quote(list(hzr_phase("cdf", t_half = 1, nu = 1, m = 1))))
  expect_equal(got$untranslated$construct, "TAU=2 AGE SEX")

  # The paired built case, so the assertion above cannot pass by the covariates
  # never having been parsed at all.
  kept <- .hzr_parse_parms(c("MUE=0.2", "THALF=1", "NU=1", "MUL=0.05", "TAU=2"),
                           covars = list(late = "AGE=1.2, SEX=0.5"))
  expect_equal(
    kept$phases,
    quote(list(hzr_phase("cdf", t_half = 1, nu = 1, m = 1),
               hzr_phase("g3", tau = 2, gamma = 1, alpha = 1, eta = 2,
                         formula = ~AGE + SEX)))
  )
  # The defaulted GAMMA*ETA is exactly 2, so SETG3_verify_ge_2() would rewrite
  # GAMMA -- one row, and not about the covariates this test is pinning.
  expect_false(any(grepl("AGE|SEX", kept$untranslated$construct)))
})

test_that("constant phase covariates with no active MUC are recorded", {
  got <- .hzr_parse_parms(c("MUE=0.2", "THALF=1", "NU=1", "MUC=0"),
                          covars = list(constant = "AGE=1.2"))
  expect_equal(got$phases, quote(list(hzr_phase("cdf", t_half = 1, nu = 1, m = 1))))
  expect_setequal(got$untranslated$construct, c("MUC=0", "AGE"))
})

test_that("FIX tokens of a phase that is not built are recorded", {
  # A reader grepping $untranslated for FIXTAU must find it: the token pinned a
  # parameter in the SAS job and has no effect on the emitted R call.
  got <- .hzr_parse_parms(c("THALF=1", "NU=1", "FIXTHALF", "MUL=0.05",
                            "TAU=2", "FIXTAU"))
  expect_equal(got$phases,
               quote(list(hzr_phase("g3", tau = 2, gamma = 1, alpha = 1,
                                    eta = 2, fixed = "tau"))))
  expect_true("THALF=1 NU=1 FIXTHALF" %in% got$untranslated$construct)
})
