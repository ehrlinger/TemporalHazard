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
  # TAU=1 is what the branch pins anyway, so it contributes no second row.
  expect_true("tau" %in% eval(got$phases)[[1L]]$fixed)
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
  expect_true(any(grepl("estimates GAMMA\\*ETA", got$untranslated$reason)))
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
  got <- .hzr_parse_parms(c("TAU=2", "GAMMA=2", "ETA=3", "FIXALPHA"))
  expect_false(any(grepl("estimates GAMMA", got$untranslated$reason)))
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
  expect_true(any(grepl("estimates GAMMA", g1$untranslated$reason)))
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
  # setg3.c:317 starts a non-positive TAU at 2*Tmax/3, which is unreproducible
  # at parse time. The emitted phase starts at tau = 1, so the difference is
  # recorded rather than passed off as a translation.
  got <- .hzr_parse_parms(c("MUL=0.01", "GAMMA=2"))
  expect_equal(nrow(got$untranslated), 1L)
  expect_equal(got$untranslated$construct, "TAU (unspecified)")
  expect_match(got$untranslated$reason, "2\\*Tmax/3")
  expect_equal(got$phases,
               quote(list(hzr_phase("g3", tau = 1, gamma = 2, alpha = 1,
                                    eta = 2))))
})

test_that("an explicit TAU = 0 takes the same SETG3 branch and is recorded", {
  # The C predicate is tau <= 0, not "absent", and hzr_phase() would reject
  # tau = 0 outright -- so without this the translator emitted a call that
  # errored where SAS supplies a default and runs.
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
  got <- .hzr_parse_parms(c("MUL=0.01", "ALPHA=1", "GAMMA=1", "ETA=1.32",
                            "FIXALPHA", "FIXGAMMA"))
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
                            "FIXALPHA", "FIXGAMMA"))
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
                            "ETA=1.32", "FIXALPHA", "FIXGAMMA"))
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
                            "FIXALPHA", "FIXGAMMA"))
  expect_false(any(grepl("Tmax", got$untranslated$reason, fixed = TRUE)))
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
