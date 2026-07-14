# SAS HAZARD parity tests
#
# Compare TemporalHazard::hazard() output to reference .lst captures from
# ~/Documents/GitHub/hazard/examples/.  Tests skip on CRAN and skip when
# the fixtures directory is unavailable.
#
# Reference parser lives in helper-sas-parity.R.
# Inventory + gap log: inst/dev/PRE-CRAN-PARITY-INVENTORY.md

skip_if_no_sas_fixtures <- function() {
  dir <- .hzr_sas_fixture_dir() # nolint: object_usage_linter.
  testthat::skip_if(is.na(dir), "SAS HAZARD fixture directory not available")
  dir
}

skip_if_no_omc_raw <- function() {
  p <- .hzr_omc_raw_path() # nolint: object_usage_linter.
  testthat::skip_if(is.na(p), "OMC raw data file not available")
  p
}

test_that("hz.deadp.KUL: 3-phase fixed-shape Weibull-tail matches SAS LL/MLEs", {
  testthat::skip_on_cran()
  dir <- skip_if_no_sas_fixtures()
  set.seed(101)  # deterministic multi-start; independent of test order

  ref <- .hzr_parse_sas_lst(file.path(dir, "hz.deadp.KUL.lst"))$fits[[1]]
  expect_equal(ref$loglik,           -3740.52,  tolerance = 1e-2)
  expect_equal(ref$n_obs,            5880L)
  expect_equal(ref$n_events,         545L)

  data(cabgkul, package = "TemporalHazard")
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data   = cabgkul,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.2, nu = 1, m = 1,
                            fixed = "shapes"),
      constant = hzr_phase("constant"),
      late     = hzr_phase("g3", tau = 1, gamma = 3, alpha = 1, eta = 1,
                            fixed = "shapes")
    ),
    fit     = TRUE,
    control = list(n_starts = 3, maxit = 500, conserve = FALSE)
  )

  # Final log-likelihood matches SAS to print precision.
  expect_equal(fit$fit$objective, ref$loglik, tolerance = 1e-2,
               label = "R log-likelihood vs SAS reference")

  # Free parameters on the optim scale (log_mu for each phase).  The full
  # theta vector includes fixed shape params; mask down to the 3 free.
  free <- !fit$fit$fixed_mask
  theta_free <- fit$fit$theta[free]
  V_free     <- fit$fit$vcov[free, free, drop = FALSE]
  expect_equal(length(theta_free), 3L)

  # Natural-scale MUE / MUC / MUL match.
  nat <- ref$natural
  mue_ref <- nat$estimate[nat$name == "MUE"]
  muc_ref <- nat$estimate[nat$name == "MUC"]
  mul_ref <- nat$estimate[nat$name == "MUL"]
  expect_equal(unname(exp(theta_free[1])), mue_ref, tolerance = 1e-3,
               label = "MUE (early) natural-scale")
  expect_equal(unname(exp(theta_free[2])), muc_ref, tolerance = 1e-3,
               label = "MUC (constant) natural-scale")
  expect_equal(unname(exp(theta_free[3])), mul_ref, tolerance = 1e-3,
               label = "MUL (late) natural-scale")

  # Std errors from the param-summary table (log-mu scale).
  params_ref <- ref$params
  se_e0 <- params_ref$se[params_ref$name == "E0"]
  se_c0 <- params_ref$se[params_ref$name == "C0"]
  se_l0 <- params_ref$se[params_ref$name == "L0"]
  expect_equal(sqrt(V_free[1, 1]), se_e0, tolerance = 5e-3,
               label = "SE(E0) vs SAS")
  expect_equal(sqrt(V_free[2, 2]), se_c0, tolerance = 5e-3,
               label = "SE(C0) vs SAS")
  expect_equal(sqrt(V_free[3, 3]), se_l0, tolerance = 5e-3,
               label = "SE(L0) vs SAS")

  # Off-diagonal vcov sanity (covariance vs SAS Asymptotic V-cov Matrix).
  expect_equal(V_free[1, 2], ref$vcov["E0", "C0"], tolerance = 5e-3,
               label = "cov(E0, C0)")
  expect_equal(V_free[2, 3], ref$vcov["C0", "L0"], tolerance = 5e-3,
               label = "cov(C0, L0)")
  expect_equal(V_free[1, 3], ref$vcov["E0", "L0"], tolerance = 5e-3,
               label = "cov(E0, L0)")
})

test_that("hz.death.AVC: 2-phase free-shape Early matches SAS LL/MLEs", {
  testthat::skip_on_cran()
  dir <- skip_if_no_sas_fixtures()
  set.seed(102)  # deterministic multi-start; independent of test order

  ref <- .hzr_parse_sas_lst(file.path(dir, "hz.death.AVC.lst"))$fits[[1]]
  expect_equal(ref$loglik,    -210.501, tolerance = 1e-2)
  expect_equal(ref$n_obs,     310L)
  expect_equal(ref$n_events,  70L)

  data(avc, package = "TemporalHazard")
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data   = avc,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.1512, nu = 1.44, m = 1,
                            fixed = "m"),
      constant = hzr_phase("constant")
    ),
    fit     = TRUE,
    control = list(n_starts = 3, maxit = 500, conserve = TRUE)
  )

  expect_equal(fit$fit$objective, ref$loglik, tolerance = 1e-2,
               label = "R log-likelihood vs SAS reference")

  # Natural-scale parameters from "Estimates for Model Parameters" table.
  nat <- ref$natural
  thalf_ref <- nat$estimate[nat$name == "THALF" & nat$phase == "Early"]
  nu_ref    <- nat$estimate[nat$name == "NU"    & nat$phase == "Early"]
  mue_ref   <- nat$estimate[nat$name == "MUE"]
  muc_ref   <- nat$estimate[nat$name == "MUC"]

  th <- fit$fit$theta
  expect_equal(unname(exp(th["early.log_t_half"])), thalf_ref, tolerance = 1e-3,
               label = "THALF (Early) natural-scale")
  expect_equal(unname(th["early.nu"]),              nu_ref,    tolerance = 1e-3,
               label = "NU (Early) natural-scale")
  expect_equal(unname(exp(th["early.log_mu"])),     mue_ref,   tolerance = 1e-3,
               label = "MUE natural-scale")
  expect_equal(unname(exp(th["constant.log_mu"])),  muc_ref,   tolerance = 1e-3,
               label = "MUC natural-scale (CoE-constrained)")

  # SE comparison (log-mu / log-t_half scale). With the full-information vcov
  # for CoE fits, all three log-scale SEs match SAS -- including E0, the
  # CoE-conserved early log_mu, whose variance was previously dropped (the old
  # ~0.059-vs-0.133 gap; now resolved).
  se <- sqrt(diag(fit$fit$vcov))
  names(se) <- names(th)
  se_ref <- function(nm) ref$params$se[ref$params$name == nm]
  expect_equal(unname(se[["early.log_t_half"]]), se_ref("E2"), tolerance = 5e-3,
               label = "SE(log_t_half / E2)")
  expect_equal(unname(se[["early.log_mu"]]),     se_ref("E0"), tolerance = 5e-3,
               label = "SE(early log_mu / E0) -- conserved phase, full-info vcov")
  expect_equal(unname(se[["constant.log_mu"]]),  se_ref("C0"), tolerance = 5e-3,
               label = "SE(constant log_mu / C0)")
})

# ---------------------------------------------------------------------------
# hm.deadp.VALVES: nu = 0, m < 0 parametric-form correctness check
# ---------------------------------------------------------------------------
# The inventory P1 gap #5 originally claimed a ~17-unit LL gap at nu=0, m<0.
# Investigation showed R's null model matches SAS exactly (LL = -1864.76),
# confirming that hzr_decompos() Case 2L (m<0, nu=0) is correctly
# implemented.  The previously observed "gap" was a spurious comparison
# between R's covariate-model LL and the SAS stepwise null LL.
#
# The final covariate model (SELECTION SLS=0.1 MAXSTEPS=5, 7 free params)
# belongs to the stepwise harness (gap #3); R finds a better optimum
# (-1786) than SAS (-1804.78) because the stepwise initialization puts
# SAS in a different basin.  That is not a bug.
#
# This test verifies the parametric-form correctness only (null model).
test_that("hm.deadp.VALVES: null model nu=0 m<0 LL matches SAS initial (Case 2L check)", {
  testthat::skip_on_cran()
  skip_if_no_sas_fixtures()
  set.seed(103)  # deterministic multi-start; independent of test order

  # SAS Initial Summary LL (no covariates, before stepwise selection):
  # shapes THALF=0.6781774 FIXTHALF, NU=0 FIXNU, M=-2.15596 FIXM; free: MUE, MUC.
  data(valves, package = "TemporalHazard")
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data   = valves,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.6781774, nu = 0, m = -2.15596,
                             fixed = c("t_half", "nu", "m")),
      constant = hzr_phase("constant")
    ),
    fit     = TRUE,
    control = list(n_starts = 3, maxit = 500, conserve = FALSE)
  )

  # SAS initial LL (from the first "Log likelihood = " before any step):
  # grep confirmed: -1864.76.
  expect_equal(fit$fit$objective, -1864.76, tolerance = 1e-2,
               label = "R null LL vs SAS initial (Case 2L: m<0, nu=0)")

  # Natural-scale mu values should be positive and finite.
  th <- fit$fit$theta
  expect_true(is.finite(th["early.log_mu"]),    label = "MUE finite")
  expect_true(is.finite(th["constant.log_mu"]), label = "MUC finite")

  # NOTE: the final stepwise covariate model (LL=-1804.78 with NYHA, I_PATH,
  # AGE_COP, BLACK) is tested via hzr_stepwise() in the stepwise harness
  # (gap #3).  R finds a better optimum (-1786) for that model because the
  # multi-start optimizer avoids the basin SAS's stepwise procedure falls
  # into.  This is not a bug; both optimize the same likelihood.
})

# ---------------------------------------------------------------------------
# hm.death.patient: 2-phase Early(CDF, free THALF/NU, M fixed) + Constant,
# multivariable in BOTH phases with a covariate (AGE_COP) shared across phases.
# Primary VALVES cohort, per-patient death. CONSERVE on, no truncation.
# ---------------------------------------------------------------------------
# This exercises the case the gap list flagged as P1 #1 (a covariate name
# appearing in both Early and Constant). The `.lst` parser routes the duplicate
# AGE_COP labels to distinct slots, and `vcov()` now returns a named matrix, so
# parity can be checked by name.
#
# The fit is warm-started from the SAS estimates. Blind multi-start lands in a
# different (shrunk-beta) basin here -- a 12-parameter CoE model with widely
# scaled covariates is hard for BFGS to optimize from beta = 0 -- so this test
# validates likelihood/model-spec parity (R agrees the SAS point is the MLE and
# reports the same log-likelihood there), not blind optimizer convergence.
test_that("hm.death.patient: 2-phase both-phase covariate model matches SAS", {
  testthat::skip_on_cran()
  dir <- skip_if_no_sas_fixtures()

  ref <- .hzr_parse_sas_lst(file.path(dir, "hm.death.patient.lst"))$fits[[1]]
  expect_equal(ref$loglik,   -1792.07, tolerance = 1e-2)
  expect_equal(ref$n_obs,    1533L)
  expect_equal(ref$n_events, 338L)

  # Reference estimates, pulled from the parsed fixture (no hard-coded numbers).
  nat <- ref$natural
  mue <- nat$estimate[nat$name == "MUE"]
  muc <- nat$estimate[nat$name == "MUC"]
  thalf <- nat$estimate[nat$name == "THALF" & nat$phase == "Early"]
  nu    <- nat$estimate[nat$name == "NU"    & nat$phase == "Early"]
  # name first so vapply(names(cov), par_beta, ..., phase = "Early") reads
  # unambiguously: the positional covariate name fills `name`, `phase` is named.
  par_beta <- function(name, phase) {
    ref$params$estimate[ref$params$phase == phase & ref$params$name == name]
  }
  # SAS covariate labels paired with the R `valves` column names (DOUBLE maps
  # to `double_` because `double` is a reserved word in R).
  early_cov <- c(AGE_COP = "age_cop", NYHA = "nyha", MV_NYHA = "mv_nyha",
                 DOUBLE = "double_", AO_PINC = "ao_pinc")
  const_cov <- c(AGE_COP = "age_cop", BLACK = "black", I_PATH = "i_path")
  beta_early <- vapply(names(early_cov), par_beta, numeric(1), phase = "Early")
  beta_const <- vapply(names(const_cov), par_beta, numeric(1), phase = "Constant")

  data(valves, package = "TemporalHazard")
  d <- valves
  # SAS PROC STANDARD ... REPLACE fills missing NYHA with the variable mean.
  d$nyha[is.na(d$nyha)] <- mean(d$nyha, na.rm = TRUE)
  d$mv_nyha <- d$mitral * d$nyha

  # Internal-scale warm start: log_mu, log_t_half, nu, m (fixed), early betas,
  # then constant log_mu, constant betas -- in phase/formula order.
  theta0 <- c(log(mue), log(thalf), nu, 1, unname(beta_early),
              log(muc), unname(beta_const))

  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data   = d,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = thalf, nu = nu, m = 1, fixed = "m",
                           formula = ~ age_cop + nyha + mv_nyha + double_ + ao_pinc),
      constant = hzr_phase("constant", formula = ~ age_cop + black + i_path)
    ),
    theta   = theta0,
    fit     = TRUE,
    control = list(n_starts = 1, maxit = 2000, conserve = TRUE)
  )

  expect_true(isTRUE(fit$fit$converged))
  # R agrees the SAS estimates are the MLE: log-likelihood matches to print
  # precision and the fit does not drift away from the seed.
  expect_equal(fit$fit$objective, ref$loglik, tolerance = 1e-2,
               label = "R log-likelihood vs SAS reference")

  th <- fit$fit$theta
  expect_equal(unname(exp(th[["early.log_mu"]])),      mue,   tolerance = 1e-4,
               label = "MUE (Early) natural-scale")
  expect_equal(unname(exp(th[["early.log_t_half"]])),  thalf, tolerance = 5e-3,
               label = "THALF (Early) natural-scale")
  expect_equal(unname(th[["early.nu"]]),               nu,    tolerance = 5e-3,
               label = "NU (Early) natural-scale")
  expect_equal(unname(exp(th[["constant.log_mu"]])),   muc,   tolerance = 1e-4,
               label = "MUC (Constant) natural-scale")

  # Covariate coefficients (same scale in R and SAS), including AGE_COP which
  # enters both phases and must resolve to distinct, name-addressable slots.
  for (nm in names(early_cov)) {
    expect_equal(unname(th[[paste0("early.", early_cov[[nm]])]]),
                 unname(beta_early[[nm]]), tolerance = 5e-3,
                 label = paste0("early.", nm, " coefficient"))
  }
  for (nm in names(const_cov)) {
    expect_equal(unname(th[[paste0("constant.", const_cov[[nm]])]]),
                 unname(beta_const[[nm]]), tolerance = 5e-3,
                 label = paste0("constant.", nm, " coefficient"))
  }

  # vcov() is name-addressable for the shared covariate (P1 #1 fix): the two
  # AGE_COP coefficients occupy distinct, finite-variance slots.
  V <- vcov(fit)
  expect_true(is.matrix(V))
  expect_true(all(c("early.age_cop", "constant.age_cop") %in% rownames(V)))
})

# ---------------------------------------------------------------------------
# hm.death.AVC.deciles: 2-phase Early(CDF, free THALF/NU, M fixed) + Constant,
# 6 Early covariates + 3 Constant covariates, with STATUS shared across both
# phases. Plus a %DECILES goodness-of-fit table that hzr_deciles() reproduces.
# ---------------------------------------------------------------------------
# The fit starts from the SAS job's documented PARMS starting values (NOT the
# converged .lst estimates), with n_starts = 1: a deterministic, non-circular
# convergence test (start where SAS started, check R reaches where SAS
# finished). Blind multi-start lands in a different basin for this 12-parameter
# CoE model (extreme MUC ~ 4.4e-7 constant phase), so we seed deterministically.
#
# The %DECILES goodness-of-fit table IS captured in the .lst (under the spaced
# title "D E C I L E   A N A L Y S I S"); .hzr_parse_sas_deciles() reads it and
# we assert full CASES/EXPECTED/ACTUAL parity against hzr_deciles().
test_that("hm.death.AVC.deciles: 2-phase multivariable model matches SAS", {
  testthat::skip_on_cran()
  dir <- skip_if_no_sas_fixtures()

  # Seed for determinism. Every multi-start parity test in this file now seeds
  # itself, so the previous save/restore RNG-neutral workaround is unnecessary:
  # tests no longer depend on the RNG state left by whatever ran before them.
  set.seed(107)

  ref <- .hzr_parse_sas_lst(file.path(dir, "hm.death.AVC.deciles.lst"))$fits[[1]]
  expect_equal(ref$loglik,   -160.408, tolerance = 1e-2)
  expect_equal(ref$n_obs,    310L)
  expect_equal(ref$n_events, 70L)

  nat <- ref$natural
  mue <- nat$estimate[nat$name == "MUE"]
  muc <- nat$estimate[nat$name == "MUC"]
  thalf <- nat$estimate[nat$name == "THALF" & nat$phase == "Early"]
  nu    <- nat$estimate[nat$name == "NU"    & nat$phase == "Early"]
  # name first so vapply(names(cov), par_beta, ..., phase = "Early") reads
  # unambiguously: the positional covariate name fills `name`, `phase` is named.
  par_beta <- function(name, phase) {
    ref$params$estimate[ref$params$phase == phase & ref$params$name == name]
  }
  # SAS covariate labels -> R `avc` column names (lowercase, identical here).
  early_cov <- c(AGE = "age", COM_IV = "com_iv", MAL = "mal",
                 OPMOS = "opmos", OP_AGE = "op_age", STATUS = "status")
  const_cov <- c(INC_SURG = "inc_surg", ORIFICE = "orifice", STATUS = "status")
  beta_early <- vapply(names(early_cov), par_beta, numeric(1), phase = "Early")
  beta_const <- vapply(names(const_cov), par_beta, numeric(1), phase = "Constant")

  data(avc, package = "TemporalHazard")
  d <- avc
  # SAS PROC STANDARD ... REPLACE fills missing INC_SURG with the variable mean.
  d$inc_surg[is.na(d$inc_surg)] <- mean(d$inc_surg, na.rm = TRUE)

  # Deterministic start from the SAS job's PARMS statement (the documented
  # starting values, not the converged answer). Order: early log_mu, log_t_half,
  # nu, m(fixed), early betas; then constant log_mu, constant betas.
  start_thalf <- 0.1905077
  start_nu    <- 1.437416
  theta0 <- c(
    log(0.3504743), log(start_thalf), start_nu, 1,
    -0.03205774, 1.336675, 0.6872028, -0.01963377, 0.0002086689, 0.5169533,
    log(4.391673e-07), 1.375285, 3.11765, 1.054988
  )

  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data   = d,
    dist   = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = start_thalf, nu = start_nu, m = 1,
                           fixed = "m",
                           formula = ~ age + com_iv + mal + opmos + op_age + status),
      constant = hzr_phase("constant", formula = ~ inc_surg + orifice + status)
    ),
    theta   = theta0,
    fit     = TRUE,
    control = list(n_starts = 1, maxit = 2000, conserve = TRUE)
  )

  expect_true(isTRUE(fit$fit$converged))
  expect_equal(fit$fit$objective, ref$loglik, tolerance = 1e-2,
               label = "R log-likelihood vs SAS reference")

  th <- fit$fit$theta
  # Tolerances are looser than the previous (circular) warm-start-from-the-
  # answer version: starting from SAS's PARMS, R converges to its own optimum,
  # which agrees with SAS's reported MLEs to ~1e-4 on a flat likelihood.
  expect_equal(unname(exp(th[["early.log_mu"]])),     mue,   tolerance = 1e-3,
               label = "MUE (Early) natural-scale")
  expect_equal(unname(exp(th[["early.log_t_half"]])), thalf, tolerance = 5e-3,
               label = "THALF (Early) natural-scale")
  expect_equal(unname(th[["early.nu"]]),              nu,    tolerance = 5e-3,
               label = "NU (Early) natural-scale")
  # MUC ~ 4.4e-7: the constant phase is barely identified here (near-singular
  # Hessian), so compare on the log scale with a relative tolerance rather than
  # an absolute one on a tiny number.
  expect_equal(unname(th[["constant.log_mu"]]),       log(muc), tolerance = 5e-3,
               label = "log MUC (Constant) internal-scale")

  for (nm in names(early_cov)) {
    expect_equal(unname(th[[paste0("early.", early_cov[[nm]])]]),
                 unname(beta_early[[nm]]), tolerance = 5e-3,
                 label = paste0("early.", nm, " coefficient"))
  }
  for (nm in names(const_cov)) {
    expect_equal(unname(th[[paste0("constant.", const_cov[[nm]])]]),
                 unname(beta_const[[nm]]), tolerance = 5e-3,
                 label = paste0("constant.", nm, " coefficient"))
  }

  # STATUS enters both phases and must resolve to distinct named vcov slots.
  V <- vcov(fit)
  expect_true(is.matrix(V))
  expect_true(all(c("early.status", "constant.status") %in% rownames(V)))

  # %DECILES goodness-of-fit parity. SAS GROUPS=5 at TIME=120 over the same
  # cohort. hzr_deciles() labels group 1 = lowest risk, the reverse of SAS
  # _DECILE_ 0 = highest risk, so align by reversing R's group order.
  sas_dec <- .hzr_parse_sas_deciles(file.path(dir, "hm.death.AVC.deciles.lst"))
  expect_false(is.null(sas_dec),
               label = "SAS decile table parsed from the .lst")
  expect_equal(nrow(sas_dec), 6L)                        # overall row + 5 groups
  sas_grp <- sas_dec[!is.na(sas_dec$decile), ]
  sas_grp <- sas_grp[order(sas_grp$decile), ]           # deciles 0..4

  dec <- hzr_deciles(fit, time = 120, groups = 5)
  expect_s3_class(dec, "hzr_deciles")
  r <- dec[order(-dec$group), ]                         # group 5..1 == decile 0..4

  expect_equal(r$n, sas_grp$cases,
               label = "decile group sizes (CASES) vs SAS")
  expect_equal(r$events, sas_grp$actual,
               label = "observed events per decile (ACTUAL) vs SAS")
  expect_equal(r$expected, sas_grp$expected, tolerance = 1e-2,
               label = "expected events per decile vs SAS")

  # Conservation of events: total expected == total observed == SAS total (70).
  ov <- attr(dec, "overall")
  expect_equal(ov$total_expected, ov$total_events, tolerance = 1e-3)
  expect_equal(ov$total_events, 70L)
})

# ---------------------------------------------------------------------------
# hm.death.AVC fit 1: phase-aware forward stepwise selection (DOCUMENTED GAP)
# ---------------------------------------------------------------------------
# This fixture runs SAS's `SELECTION SLE=0.2 SLS=0.1` with shaping parameters
# fixed during selection, offering the same candidate variables to BOTH the
# Early and Constant phases with per-variable flags (/I include, /S start,
# /E exclude).  SAS converges in 10 steps to Early{AGE,STATUS,COM_IV,MAL,
# OPMOS,OP_AGE} + Constant{STATUS,ORIFICE,INC_SURG}, log-likelihood -160.64
# (hm.death.AVC.lst step trace, lines 104-650).  That FINAL selected model is
# the saved "HMDEATH" fit, whose estimates and predictions are already covered
# by the `hm.death.AVC.deciles`, `hp.death.AVC.hm1`, and `hp.death.AVC.hm2`
# parity tests -- so nothing about the converged model is left unverified here.
#
# What is NOT a parity target is the SELECTION PATH.  Investigation (2026-06-06)
# confirmed R's hzr_stepwise() cannot reproduce SAS's path or final set on this
# fixture, for fundamental (not fixable-by-tuning) reasons:
#   * SAS uses APPROXIMATE variances during selection ("variances are
#     approximate because shaping parameter covariances are ignored", lst:101);
#     R uses the full Hessian, which here is near-singular with non-positive
#     variance estimates (rcond ~1e-12), so the Wald statistics that drive the
#     enter/drop decisions are not comparable between engines.
#   * SAS's /I and /S flags are PHASE-LEVEL (force MAL into Early only); R's
#     force_in is variable-level (phase-blind), so the constraint sets differ.
#   * R oscillates at p ~ slstay (= 0.1) and freezes; seeded with the SAS /S
#     start set it lands at LL -185.8, intercept-only base already -182.8,
#     both far worse than SAS's -160.64 -- a different optimizer basin, the
#     same phenomenon already documented for hm.deadp.VALVES above.
#
# This block is therefore a REGRESSION GUARD for the multiphase phase-aware
# stepwise code path on a real fixture (it must run end-to-end and return a
# well-formed result), plus a recorded SAS reference -- NOT a parity assertion.
# See inst/dev/FIXTURE-GAP-LIST.md (Group A, hm.death.AVC fit 1).
test_that("hm.death.AVC stepwise: multiphase selection runs (SAS path is a documented gap)", {
  testthat::skip_on_cran()
  skip_if_no_sas_fixtures()
  set.seed(110)

  # Recorded SAS reference (hm.death.AVC.lst): final stepwise model.
  sas_final_vars <- list(
    early    = c("age", "status", "com_iv", "mal", "opmos", "op_age"),
    constant = c("status", "orifice", "inc_surg")
  )
  expect_length(sas_final_vars$early,    6L)
  expect_length(sas_final_vars$constant, 3L)

  data(avc, package = "TemporalHazard")
  d <- avc
  d$inc_surg[is.na(d$inc_surg)] <- mean(d$inc_surg, na.rm = TRUE)

  # Base model seeds the SAS /S (start) + /I (include) variables in their
  # phases; shapes fixed, as SAS fixes them during selection.
  base <- hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.1511909, nu = 1.438631, m = 1,
                           fixed = "shapes", formula = ~ com_iv + mal),
      constant = hzr_phase("constant",
                           formula = ~ orifice + mal + op_age + inc_surg)),
    fit = TRUE, control = list(n_starts = 1, conserve = TRUE))

  sw <- hzr_stepwise(
    base,
    scope = list(early    = ~ age + status + orifice + inc_surg + opmos + op_age,
                 constant = ~ age + status + opmos),
    data = d, direction = "both", criterion = "wald",
    slentry = 0.2, slstay = 0.1,
    force_in = c("mal", "inc_surg"),
    control = list(n_starts = 1, conserve = TRUE), trace = FALSE)

  # Regression guard only: the phase-aware multiphase stepwise path completes
  # and returns a well-formed result on real data.  We deliberately do NOT
  # assert anything about WHICH variables R selects or the attained
  # log-likelihood relative to SAS -- the selection path is a documented gap
  # (see the block comment above), and asserting divergence would perversely
  # fail if hzr_stepwise() later improved toward SAS.
  expect_s3_class(sw, "hzr_stepwise")
  expect_true(is.data.frame(sw$steps))
  expect_true(all(c("action", "variable", "phase") %in% names(sw$steps)))
  expect_true(is.finite(sw$fit$objective))
})

# ---------------------------------------------------------------------------
# ac.death.AVC: actuarial life tables -- Kaplan-Meier (%KAPLAN) and
# Nelson-Aalen (%NELSONT), overall and KM stratified by COM_IV.
# ---------------------------------------------------------------------------
# No model fit: hzr_kaplan()/hzr_nelson() are non-parametric, deterministic.
# SAS prints times to ~3 decimals, so rows are aligned by index (both tables
# are sorted by time with the same event times). Tolerances reflect SAS print
# rounding (~5e-6 in practice).
test_that("ac.death.AVC: Kaplan-Meier and Nelson-Aalen life tables match SAS", {
  testthat::skip_on_cran()
  dir <- skip_if_no_sas_fixtures()
  lst <- file.path(dir, "ac.death.AVC.lst")
  # An older local fixture set may predate this capture; skip cleanly rather
  # than hard-error in readLines() if this specific .lst is absent.
  testthat::skip_if_not(file.exists(lst), "ac.death.AVC.lst fixture not present")

  data(avc, package = "TemporalHazard")

  # --- Kaplan-Meier, overall (%KAPLAN _CATG=ALL) ---------------------------
  km <- hzr_kaplan(avc$int_dead, avc$dead)
  sk <- .hzr_parse_sas_lifetable(lst, "kaplan")
  expect_false(is.null(sk), label = "KAPLAN life table parsed")
  expect_equal(nrow(km), nrow(sk))
  expect_equal(km$n_risk,  sk$NUMBER)
  expect_equal(km$n_event, sk$DEAD)
  expect_equal(km$survival, sk$CUM_SURV, tolerance = 1e-4,
               label = "KM survival vs SAS CUM_SURV")
  expect_equal(km$std_err, sk$SE_EXACT, tolerance = 1e-4,
               label = "KM standard error vs SAS SE_EXACT")
  # CUM_HAZ = -log(CUM_SURV); the final all-events row is S=0 (NA in SAS).
  ok <- !is.na(sk$CUM_HAZ)
  expect_equal(km$cumhaz[ok], sk$CUM_HAZ[ok], tolerance = 1e-4,
               label = "KM cumulative hazard (-log S) vs SAS CUM_HAZ")

  # --- Nelson-Aalen, overall (%NELSONT _CATG=ALL) --------------------------
  ne <- hzr_nelson(avc$int_dead, avc$dead)
  sn <- .hzr_parse_sas_lifetable(lst, "nelson")
  expect_false(is.null(sn), label = "NELSONT life table parsed")
  expect_equal(nrow(ne), nrow(sn))
  ok_n <- !is.na(sn$CUM_HAZ)
  expect_equal(ne$cumhaz[ok_n], sn$CUM_HAZ[ok_n], tolerance = 1e-4,
               label = "Nelson-Aalen cumulative hazard vs SAS CUM_HAZ")

  # --- Kaplan-Meier stratified by COM_IV (%KAPLAN STRATIFY=COM_IV) ----------
  for (g in 0:1) {
    sub <- avc[avc$com_iv == g, ]
    kg  <- hzr_kaplan(sub$int_dead, sub$dead)
    sg  <- .hzr_parse_sas_lifetable(lst, paste0("catg", g))
    expect_false(is.null(sg), label = paste0("COM_IV=", g, " stratum parsed"))
    expect_equal(nrow(kg), nrow(sg))
    expect_equal(kg$n_risk, sg$NUMBER,
                 label = paste0("COM_IV=", g, " n at risk"))
    expect_equal(kg$survival, sg$CUM_SURV, tolerance = 1e-4,
                 label = paste0("COM_IV=", g, " KM survival"))
  }
})

# ---------------------------------------------------------------------------
# hp.death.AVC: HAZPRED predictions from the saved hz.death.AVC fit.
# ---------------------------------------------------------------------------
# %HAZPRED predicts survival and hazard at a "digital nomogram" of time points
# from a saved model (EXAMPLES.HZDEATH = the hz.death.AVC 2-phase Early+Constant
# fit). We refit that model deterministically from its .sas PARMS and predict at
# the same exact times (reconstructed from the SAS DO loop; the printed MONTHS
# are rounded to 3 decimals, so predicting at the rounded values would inflate
# the error). Survival and the instantaneous multiphase hazard both go through
# the public predict() path (type = "survival" / type = "hazard"), with their
# HAZPRED confidence limits also asserted: survival via conf.type = "logit"
# (HAZPRED's logit transform; R's default is complementary-log-log) and hazard
# via the log scale (which already matches HAZPRED). The full-information vcov
# for CoE fits supplies the se(H) basis that makes the CLs match to ~1e-5.
test_that("hp.death.AVC: HAZPRED survival/hazard nomogram matches SAS", {
  testthat::skip_on_cran()
  dir <- skip_if_no_sas_fixtures()
  lst <- file.path(dir, "hp.death.AVC.lst")
  testthat::skip_if_not(file.exists(lst), "hp.death.AVC.lst fixture not present")

  nom <- .hzr_parse_sas_nomogram(lst)
  expect_false(is.null(nom), label = "HAZPRED nomogram parsed")

  # Exact prediction times from hp.death.AVC.sas:
  #   DO MONTHS = 1*DTY..7*DTY, 14*DTY, 30*DTY, 1,2,3,6,12,18, 24 TO 180 BY 12
  # with DTY = 12 / 365.2425 (one day, in months).
  dty <- 12 / 365.2425
  months <- c((1:7) * dty, 14 * dty, 30 * dty, 1, 2, 3, 6, 12, 18,
              seq(24, 180, by = 12))
  expect_equal(length(months), nrow(nom))

  # Saved model = hz.death.AVC fit; refit deterministically from its .sas PARMS.
  data(avc, package = "TemporalHazard")
  theta0 <- c(log(0.2361727), log(0.1512095), 1.438652, 1, log(0.0005437099))
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1,
    data    = avc,
    dist    = "multiphase",
    phases  = list(
      early    = hzr_phase("cdf", t_half = 0.1512095, nu = 1.438652, m = 1,
                           fixed = "m"),
      constant = hzr_phase("constant")
    ),
    theta   = theta0,
    fit     = TRUE,
    control = list(n_starts = 1, maxit = 2000, conserve = TRUE)
  )
  expect_equal(fit$fit$objective, -210.501, tolerance = 1e-2,
               label = "saved-model LL vs SAS")

  # Survival predictions (public API).
  surv <- predict(fit, newdata = data.frame(time = months), type = "survival")
  expect_equal(unname(surv), nom$SURVIV, tolerance = 1e-4,
               label = "HAZPRED _SURVIV nomogram")

  # SAS HAZPRED uses a 1-SD level (CLEVEL = 0.68268948 -> z = 1).
  lvl <- 2 * stats::pnorm(1) - 1

  # Survival confidence limits: conf.type = "logit" reproduces HAZPRED's
  # logit-scale CLs (the default "log-log" intentionally does not).
  scl <- predict(fit, newdata = data.frame(time = months), type = "survival",
                 se.fit = TRUE, level = lvl, conf.type = "logit")
  expect_equal(scl$lower, nom$CLLSURV, tolerance = 1e-4,
               label = "HAZPRED _CLLSURV (logit CL)")
  expect_equal(scl$upper, nom$CLUSURV, tolerance = 1e-4,
               label = "HAZPRED _CLUSURV (logit CL)")

  # Instantaneous multiphase hazard via the public predict() path.
  haz <- predict(fit, newdata = data.frame(time = months), type = "hazard")
  expect_equal(unname(haz), nom$HAZARD, tolerance = 1e-3,
               label = "HAZPRED _HAZARD nomogram")

  # Hazard confidence limits (log scale -- HAZPRED and R agree here).
  hcl <- predict(fit, newdata = data.frame(time = months), type = "hazard",
                 se.fit = TRUE, level = lvl)
  expect_equal(hcl$lower, nom$CLLHAZ, tolerance = 1e-4,
               label = "HAZPRED _CLLHAZ")
  expect_equal(hcl$upper, nom$CLUHAZ, tolerance = 1e-4,
               label = "HAZPRED _CLUHAZ")
})

# ---------------------------------------------------------------------------
# bs.death.AVC: %HAZBOOT bootstrap + per-resample stepwise (DOCUMENTED GAP)
# ---------------------------------------------------------------------------
# SAS %HAZBOOT(RESAMPL=5, SLE=0.12, SLS=0.1) draws bootstrap resamples and runs
# a *fresh stepwise variable selection* on each, reporting one selected model
# per resample (bs.death.AVC.lst: a "." marks a covariate not selected in that
# resample).  The headline result is a variable-selection FREQUENCY per phase
# (e.g. Early STATUS/COM_IV selected in 5/5 resamples, OP_AGE 3/5, MAL 2/5).
#
# R has no equivalent procedure.  hzr_bootstrap() resamples and refits a FIXED
# model (covariate-stability bootstrap); it has no embedded-selection mode, and
# building one would call hzr_stepwise() per resample -- which already diverges
# from SAS on this exact AVC model (see the hm.death.AVC stepwise documented
# gap).  Combined with %HAZBOOT's intrinsic non-determinism (SEED=-1) and the
# tiny RESAMPL=5, selection-frequency parity is not a meaningful target today.
#
# This block therefore (1) captures and asserts the SAS reference frequencies
# in parseable form -- so the parity test is half-written for if/when R gains a
# bootstrap-with-selection capability -- and (2) regression-guards that R's
# (fixed-model) hzr_bootstrap() runs end-to-end on this dataset.  It does NOT
# assert R-vs-SAS frequency parity.  See inst/dev/FIXTURE-GAP-LIST.md.
test_that("bs.death.AVC: SAS bootstrap selection frequencies parse (parity is a documented gap)", {
  testthat::skip_on_cran()
  dir <- skip_if_no_sas_fixtures()
  lst <- file.path(dir, "bs.death.AVC.lst")
  testthat::skip_if_not(file.exists(lst), "bs.death.AVC.lst fixture not present")

  boot <- .hzr_parse_sas_bootstrap(lst)
  expect_false(is.null(boot), label = "SAS %HAZBOOT tables parsed")
  expect_true(all(c("early", "constant") %in% names(boot)))
  expect_equal(nrow(boot$early), 5L)          # RESAMPL = 5

  # SAS reference selection frequencies (out of 5 resamples).  Names are the
  # SAS column labels (upper case), preserved as parsed.
  freq <- attr(boot, "selection_freq")
  expect_equal(freq$early[["STATUS"]], 5L)
  expect_equal(freq$early[["COM_IV"]], 5L)
  expect_equal(freq$early[["AGE"]],    4L)
  expect_equal(freq$early[["OPMOS"]],  4L)
  expect_equal(freq$early[["OP_AGE"]], 3L)
  expect_equal(freq$early[["MAL"]],    2L)
  expect_equal(freq$constant[["OPMOS"]],   5L)
  expect_equal(freq$constant[["ORIFICE"]], 3L)
})

test_that("bs.death.AVC: R fixed-model bootstrap runs on the AVC cohort", {
  testthat::skip_on_cran()
  # No SAS .lst needed here -- this exercises R's bootstrap path on the shipped
  # `avc` dataset, so it must run in public CI (do NOT gate on fixtures).
  set.seed(111)

  data(avc, package = "TemporalHazard")
  # Null 2-phase AVC model (matches hz.death.AVC), bootstrapped.  This exercises
  # R's bootstrap path on the fixture cohort; it is NOT the %HAZBOOT procedure.
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.1512, nu = 1.44, m = 1, fixed = "m"),
      constant = hzr_phase("constant")),
    fit = TRUE, control = list(n_starts = 1, conserve = TRUE))

  bs <- hzr_bootstrap(fit, n_boot = 20L, seed = 111L)
  expect_s3_class(bs, "hzr_bootstrap")
  expect_gt(bs$n_success, 0L)
})

test_that("bs.death.AVC: R scope= embedded stepwise selection runs on the AVC cohort", {
  testthat::skip_on_cran()
  # No SAS .lst needed here -- exercises R's new scope= bootstrap-selection
  # path (Task 2 of inst/dev/PLAN-bootstrap-stepwise-selection.md) on the
  # shipped avc dataset, so it must run in public CI (do NOT gate on
  # fixtures). This intentionally does NOT assert R-vs-SAS selection
  # frequency parity -- see the "bs.death.AVC" SAS-frequency test above and
  # inst/dev/BOOTSTRAP-SELECTION-DESIGN.md for why that remains a
  # documented, accepted gap.
  set.seed(222)

  data(avc, package = "TemporalHazard")
  # Same null 2-phase AVC model as the fixed-refit test above, with shapes
  # fixed -- the R equivalent of hz.death.AVC's output feeding bh.death.AVC.
  base <- hazard(
    survival::Surv(int_dead, dead) ~ 1, data = avc, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.1512, nu = 1.44, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")),
    fit = TRUE, control = list(n_starts = 1, conserve = TRUE))

  bs <- hzr_bootstrap(base, n_boot = 10L, seed = 222L,
                       scope = list(early    = ~ age + mal,
                                    constant = ~ status + com_iv),
                       slentry = 0.3, slstay = 0.2,
                       control = list(n_starts = 1, conserve = TRUE))

  expect_s3_class(bs, "hzr_bootstrap")
  expect_identical(bs$mode, "select")
  expect_gt(bs$n_success, 0L)
  # Phase-qualified shape parameters are never dropped by stepwise
  # selection, so both phases must be represented in the summary.
  params <- bs$summary$parameter
  expect_true(any(grepl("^early\\.", params)))
  expect_true(any(grepl("^constant\\.", params)))
})

# ---------------------------------------------------------------------------
# hp.death.AVC.hm1 / hm2: patient-specific HAZPRED predictions from the saved
# multivariable both-phase model (hm.death.AVC final fit, "HMDEATH").
# ---------------------------------------------------------------------------
# Both jobs predict survival + instantaneous hazard (with logit survival CLs
# and log hazard CLs at the 1-SD level) for two covariate profiles printed as
# a BY-group "digital nomogram".  Each printed row carries its own covariate
# vector, so newdata is built directly from the parsed table -- the test is
# agnostic to the specific profiles (hm1 varies MAL 0/1; hm2 varies COM_IV 0/1
# across dates of repair).
#
# Tolerances are looser than the null-model nomogram (hp.death.AVC): this is a
# 9-coefficient, near-singular fit (rcond ~5e-12; see hm.death.AVC.deciles),
# so R converges to its own optimum a few e-4 from SAS's reported MLEs, which
# propagates to ~5e-4 in survival and (amplified at the steep early-phase
# times) ~8e-3 in the instantaneous hazard.  SAS 1-SD level: CLEVEL = z = 1.
.hzr_check_hmdeath_nomogram <- function(df, fit) {
  lvl <- 2 * stats::pnorm(1) - 1
  nd <- data.frame(
    time    = df$MONTHS,
    age     = df$AGE,
    com_iv  = df$COM_IV,
    mal     = df$MAL,
    opmos   = df$OPMOS,
    op_age  = df$OPMOS * df$AGE,
    status  = df$STATUS,
    inc_surg = df$INC_SURG,
    orifice = df$ORIFICE
  )

  surv <- predict(fit, newdata = nd, type = "survival")
  testthat::expect_equal(unname(surv), df$SURVIV, tolerance = 5e-4,
               label = "HAZPRED _SURVIV (multivariable nomogram)")
  scl <- predict(fit, newdata = nd, type = "survival", se.fit = TRUE,
                 level = lvl, conf.type = "logit")
  testthat::expect_equal(scl$lower, df$CLLSURV, tolerance = 5e-4,
               label = "HAZPRED _CLLSURV (logit CL)")
  testthat::expect_equal(scl$upper, df$CLUSURV, tolerance = 5e-4,
               label = "HAZPRED _CLUSURV (logit CL)")

  haz <- predict(fit, newdata = nd, type = "hazard")
  testthat::expect_equal(unname(haz), df$HAZARD, tolerance = 8e-3,
               label = "HAZPRED _HAZARD (multivariable nomogram)")
  hcl <- predict(fit, newdata = nd, type = "hazard", se.fit = TRUE, level = lvl)
  testthat::expect_equal(hcl$lower, df$CLLHAZ, tolerance = 8e-3,
               label = "HAZPRED _CLLHAZ")
  testthat::expect_equal(hcl$upper, df$CLUHAZ, tolerance = 8e-3,
               label = "HAZPRED _CLUHAZ")
}

test_that("hp.death.AVC.hm1: patient-specific predictions (MAL 0/1) match SAS", {
  testthat::skip_on_cran()
  dir <- skip_if_no_sas_fixtures()
  lst <- file.path(dir, "hp.death.AVC.hm1.lst")
  testthat::skip_if_not(file.exists(lst), "hp.death.AVC.hm1.lst fixture not present")

  df <- .hzr_parse_sas_nomogram_mv(lst)
  expect_false(is.null(df), label = "hm1 multivariable nomogram parsed")
  expect_setequal(unique(df$MAL), c(0, 1))     # two covariate profiles

  set.seed(108)
  fit <- .hzr_fit_avc_hmdeath()
  expect_equal(fit$fit$objective, -160.408, tolerance = 1e-2,
               label = "HMDEATH saved-model LL vs SAS")
  .hzr_check_hmdeath_nomogram(df, fit)
})

test_that("hp.death.AVC.hm2: predictions by date of repair (COM_IV 0/1) match SAS", {
  testthat::skip_on_cran()
  dir <- skip_if_no_sas_fixtures()
  lst <- file.path(dir, "hp.death.AVC.hm2.lst")
  testthat::skip_if_not(file.exists(lst), "hp.death.AVC.hm2.lst fixture not present")

  df <- .hzr_parse_sas_nomogram_mv(lst)
  expect_false(is.null(df), label = "hm2 multivariable nomogram parsed")
  expect_setequal(unique(df$COM_IV), c(0, 1))  # complete vs partial canal

  set.seed(109)
  fit <- .hzr_fit_avc_hmdeath()
  expect_equal(fit$fit$objective, -160.408, tolerance = 1e-2,
               label = "HMDEATH saved-model LL vs SAS")
  .hzr_check_hmdeath_nomogram(df, fit)
})

# ---------------------------------------------------------------------------
# hs.death.AVC.hm1: stratified observed-vs-expected calibration + mean survival
# ---------------------------------------------------------------------------
# Patient-specific HAZPRED predictions from the saved HMDEATH model, summarized
# BY type of AV canal (COM_IV).  Two SAS tables are reproduced:
#   1. "Predict number of deaths" -- per stratum, EXPECTED = sum of predicted
#      cumulative hazard at each subject's own follow-up, PEXPECT = sum of
#      predicted death probability (1 - S), ACTUAL = observed deaths.  Totals
#      conserve events (14.76 + 55.24 = 70 = ACTUAL total).
#   2. Per-stratum mean survival at a digital time grid (MSURVIV = mean over the
#      stratum of each subject's predicted survival; NSURVIV = stratum size).
# This is the population-averaged / calibration counterpart to the per-patient
# hp.death.AVC.hm1/hm2 nomograms.  Same fit; tolerances reflect the same
# near-singular 9-coefficient model (sums over 154/156 subjects, so absolute
# tolerances are scaled up accordingly).
test_that("hs.death.AVC.hm1: stratified observed-vs-expected calibration matches SAS", {
  testthat::skip_on_cran()
  dir <- skip_if_no_sas_fixtures()
  lst <- file.path(dir, "hs.death.AVC.hm1.lst")
  testthat::skip_if_not(file.exists(lst), "hs.death.AVC.hm1.lst fixture not present")

  cal <- .hzr_parse_sas_calibration(lst)
  expect_false(is.null(cal), label = "SAS observed-vs-expected table parsed")
  expect_equal(nrow(cal), 2L)                       # two COM_IV strata
  cal <- cal[order(cal$com_iv), ]

  set.seed(112)
  fit <- .hzr_fit_avc_hmdeath()
  expect_equal(fit$fit$objective, -160.408, tolerance = 1e-2)

  data(avc, package = "TemporalHazard")
  d <- avc
  d$inc_surg[is.na(d$inc_surg)] <- mean(d$inc_surg, na.rm = TRUE)
  nd <- data.frame(
    time = d$int_dead, age = d$age, com_iv = d$com_iv, mal = d$mal,
    opmos = d$opmos, op_age = d$opmos * d$age, status = d$status,
    inc_surg = d$inc_surg, orifice = d$orifice)
  ch <- predict(fit, newdata = nd, type = "cumulative_hazard")
  sv <- predict(fit, newdata = nd, type = "survival")

  for (k in seq_len(nrow(cal))) {
    g <- cal$com_iv[k]
    i <- d$com_iv == g
    expect_equal(sum(ch[i]),       cal$expected[k], tolerance = 5e-3,
                 label = paste0("EXPECTED (sum cumhaz), COM_IV=", g))
    expect_equal(sum(1 - sv[i]),   cal$pexpect[k],  tolerance = 5e-3,
                 label = paste0("PEXPECT (sum 1-S), COM_IV=", g))
    expect_equal(sum(d$dead[i]),   cal$actual[k],
                 label = paste0("ACTUAL (observed deaths), COM_IV=", g))
  }
  # Conservation of events across strata: total expected == total observed.
  expect_equal(sum(ch), sum(cal$actual), tolerance = 1e-2,
               label = "total expected == total observed (CoE)")
})

test_that("hs.death.AVC.hm1: per-stratum mean survival curve matches SAS", {
  testthat::skip_on_cran()
  dir <- skip_if_no_sas_fixtures()
  lst <- file.path(dir, "hs.death.AVC.hm1.lst")
  testthat::skip_if_not(file.exists(lst), "hs.death.AVC.hm1.lst fixture not present")

  ms <- .hzr_parse_sas_strata_survival(lst)
  expect_false(is.null(ms), label = "SAS per-stratum mean-survival table parsed")
  expect_setequal(unique(ms$com_iv), c(0, 1))

  set.seed(113)
  fit <- .hzr_fit_avc_hmdeath()

  data(avc, package = "TemporalHazard")
  d <- avc
  d$inc_surg[is.na(d$inc_surg)] <- mean(d$inc_surg, na.rm = TRUE)

  # NSURVIV (stratum size) must match the cohort.
  for (g in c(0, 1)) {
    n_sas <- ms$NSURVIV[ms$com_iv == g][1]
    expect_equal(sum(d$com_iv == g), n_sas,
                 label = paste0("stratum size NSURVIV, COM_IV=", g))
  }

  # MSURVIV = mean over the stratum of each subject's predicted survival at the
  # grid time (MONTHS = YEARS * 12).  Compare row by row.
  r_msurviv <- numeric(nrow(ms))
  for (k in seq_len(nrow(ms))) {
    g <- ms$com_iv[k]
    i <- which(d$com_iv == g)
    nd <- data.frame(
      time = rep(ms$YEARS[k] * 12, length(i)),
      age = d$age[i], com_iv = d$com_iv[i], mal = d$mal[i], opmos = d$opmos[i],
      op_age = d$opmos[i] * d$age[i], status = d$status[i],
      inc_surg = d$inc_surg[i], orifice = d$orifice[i])
    r_msurviv[k] <- mean(predict(fit, newdata = nd, type = "survival"))
  }
  expect_equal(r_msurviv, ms$MSURVIV, tolerance = 5e-4,
               label = "per-stratum mean survival (MSURVIV) vs SAS")
})

# ---------------------------------------------------------------------------
# hz.te123.OMC: 2-phase (Early CDF + Late Weibull) repeated TE events
# ---------------------------------------------------------------------------
# Two fits:
#   Fit 1 -- left-truncated (LCENSOR STARTTME), no covariates, CONSERVE
#   Fit 2 -- modulated renewal (INT_TE = INT_TE - STARTTME), late-phase
#            covariates NOPREVTE + NOTEE, CONSERVE
#
# Data derivation: PRIMISOL DATA step reproduced in .hzr_derive_primisol();
# raw fixed-width file from ~/Documents/GitHub/hazard/examples/data/omc.
#
# Known gaps:
#   P1 #6 (fit 1): CoE for early+late 2-phase differs from SAS CONSERVE.
#   SAS keeps all 5 mus free under Lagrange constraint; R fixes late.log_mu
#   via closed-form CoE, yielding LL offset ~0.14 and MUE shifted.
#   Shape params (THALF, NU, ETA) still match well; asserted below.
#   MUE NOT asserted (affected by P1 #6 gap).
test_that("hz.te123.OMC fit 1: left-truncated 2-phase shape params match SAS", {
  testthat::skip_on_cran()
  dir  <- skip_if_no_sas_fixtures()
  skip_if_no_omc_raw()
  set.seed(104)  # deterministic multi-start; independent of test order

  ref  <- .hzr_parse_sas_lst(file.path(dir, "hz.te123.OMC.lst"))$fits[[1]]
  expect_equal(ref$loglik,   -322.226, tolerance = 1e-2)
  expect_equal(ref$n_obs,    382L)
  expect_equal(ref$n_events, 44L)

  raw      <- .hzr_omc_raw_path()
  primisol <- .hzr_derive_primisol(raw)

  # Verify the R derivation reproduces SAS observation counts.
  expect_equal(nrow(primisol$te), 382L,   label = "PRIMISOL n_obs")
  expect_equal(sum(primisol$te$te), 44L,  label = "PRIMISOL n_events")

  fit <- hazard(
    survival::Surv(starttme, int_te, te) ~ 1,
    data   = primisol$te,
    dist   = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.06, nu = -1.74, m = 0, fixed = "m"),
      late  = hzr_phase("g3",  tau = 1, gamma = 1, alpha = 1, eta = 1.32,
                         fixed = c("tau", "gamma", "alpha"))
    ),
    fit     = TRUE,
    control = list(n_starts = 3, maxit = 500, conserve = TRUE)
  )

  # LL: P1 #6 gap produces ~0.14 offset; assert within 0.2.
  expect_equal(fit$fit$objective, ref$loglik, tolerance = 0.2,
               label = "R LL vs SAS (P1 #6 gap: CoE implementation differs)")

  # Shape params from "Estimates for Model Parameters" natural-scale table.
  nat       <- ref$natural
  thalf_ref <- nat$estimate[nat$name == "THALF" & nat$phase == "Early"]
  nu_ref    <- nat$estimate[nat$name == "NU"    & nat$phase == "Early"]
  eta_ref   <- nat$estimate[nat$name == "ETA"   & nat$phase == "Late"]

  th <- fit$fit$theta
  expect_equal(unname(exp(th["early.log_t_half"])), thalf_ref, tolerance = 5e-3,
               label = "THALF (Early) natural-scale")
  expect_equal(unname(th["early.nu"]),              nu_ref,    tolerance = 5e-3,
               label = "NU (Early) natural-scale")
  expect_equal(unname(th["late.eta"]),              eta_ref,   tolerance = 5e-3,
               label = "ETA (Late) natural-scale")

  # NOTE: MUE NOT asserted — affected by P1 #6 CoE gap.
  # R: MUE ~ 0.01601; SAS: MUE = 0.01740.
})

test_that("hz.te123.OMC fit 2: modulated renewal + late covariates matches SAS LL/MLEs", {
  testthat::skip_on_cran()
  dir  <- skip_if_no_sas_fixtures()
  skip_if_no_omc_raw()
  set.seed(105)  # deterministic multi-start; independent of test order

  ref  <- .hzr_parse_sas_lst(file.path(dir, "hz.te123.OMC.lst"))$fits[[2]]
  expect_equal(ref$loglik,   -311.597, tolerance = 1e-2)
  expect_equal(ref$n_obs,    382L)
  expect_equal(ref$n_events, 44L)

  raw      <- .hzr_omc_raw_path()
  primisol <- .hzr_derive_primisol(raw)

  # Fit 2: modulated renewal.  INT_TE = INT_TE - STARTTME (relative time);
  # no LCENSOR; late phase has covariates NOPREVTE and NOTEE.
  fit <- hazard(
    survival::Surv(int_te, te) ~ 1,
    data   = primisol$te_mod,
    dist   = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.053, nu = -1.738, m = 0, fixed = "m"),
      late  = hzr_phase("g3",  tau = 1, gamma = 1, alpha = 1, eta = 1,
                         fixed = "shapes",
                         formula = ~ noprevte + notee)
    ),
    fit     = TRUE,
    control = list(n_starts = 3, maxit = 500, conserve = TRUE)
  )

  expect_equal(fit$fit$objective, ref$loglik, tolerance = 1e-2,
               label = "R LL vs SAS (fit 2)")

  nat       <- ref$natural
  thalf_ref <- nat$estimate[nat$name == "THALF" & nat$phase == "Early"]
  nu_ref    <- nat$estimate[nat$name == "NU"    & nat$phase == "Early"]
  mue_ref   <- nat$estimate[nat$name == "MUE"]
  mul_ref   <- nat$estimate[nat$name == "MUL"]

  th <- fit$fit$theta
  expect_equal(unname(exp(th["early.log_t_half"])), thalf_ref, tolerance = 1e-3,
               label = "THALF (Early) natural-scale")
  expect_equal(unname(th["early.nu"]),              nu_ref,    tolerance = 1e-3,
               label = "NU (Early) natural-scale")
  expect_equal(unname(exp(th["early.log_mu"])),     mue_ref,   tolerance = 1e-3,
               label = "MUE (Early) natural-scale")
  expect_equal(unname(exp(th["late.log_mu"])),      mul_ref,   tolerance = 1e-3,
               label = "MUL (Late) natural-scale")

  # Late-phase covariate coefficients — compare to SAS Parameter Estimate Summary.
  params_ref <- ref$params
  np_ref  <- params_ref$estimate[params_ref$name == "NOPREVTE"]
  ne_ref  <- params_ref$estimate[params_ref$name == "NOTEE"]
  expect_equal(unname(th["late.noprevte"]), np_ref, tolerance = 1e-3,
               label = "NOPREVTE coefficient")
  expect_equal(unname(th["late.notee"]),    ne_ref, tolerance = 1e-3,
               label = "NOTEE coefficient")

  # SEs from the parameter-estimate summary table (log-scale where applicable).
  V <- fit$fit$vcov
  param_idx <- function(nm) which(names(th) == nm)
  se_e2_ref  <- params_ref$se[params_ref$name == "E2"]
  se_np_ref  <- params_ref$se[params_ref$name == "NOPREVTE"]
  se_ne_ref  <- params_ref$se[params_ref$name == "NOTEE"]
  expect_equal(sqrt(V[param_idx("early.log_t_half"), param_idx("early.log_t_half")]),
               se_e2_ref, tolerance = 5e-3, label = "SE(log_t_half / E2)")
  expect_equal(sqrt(V[param_idx("late.noprevte"), param_idx("late.noprevte")]),
               se_np_ref, tolerance = 5e-2, label = "SE(NOPREVTE)")
  expect_equal(sqrt(V[param_idx("late.notee"),    param_idx("late.notee")]),
               se_ne_ref, tolerance = 5e-2, label = "SE(NOTEE)")
})

# ---------------------------------------------------------------------------
# hz.tm123.OMC: 2-phase (Early CDF + Late Weibull) permanent morbidity
# ---------------------------------------------------------------------------
# Uses WEIGHT MORBID: each TE event is weighted by severity grade.
# SAS WEIGHT applies the weight to the event log-hazard term only;
# censored intervals contribute to the survival function with weight = 1.
# R: pass weights = ifelse(te==1, morbid, 1) to reproduce this semantics.
#
# Known gaps:
#   Same P1 #6 gap as fit 1 above (CoE for early+late, MUE shifted ~10%).
#   LL offset ~0.43.  THALF and ETA assert within 5e-3.
test_that("hz.tm123.OMC: morbidity-weighted 2-phase shape params match SAS", {
  testthat::skip_on_cran()
  dir  <- skip_if_no_sas_fixtures()
  skip_if_no_omc_raw()
  set.seed(106)  # deterministic multi-start; independent of test order

  ref  <- .hzr_parse_sas_lst(file.path(dir, "hz.tm123.OMC.lst"))$fits[[1]]
  expect_equal(ref$loglik,   -581.528, tolerance = 1e-2)
  expect_equal(ref$n_obs,    382L)
  expect_equal(ref$n_events, 44L)

  raw      <- .hzr_omc_raw_path()
  primisol <- .hzr_derive_primisol(raw)

  # WEIGHT MORBID: apply severity weight to event rows; censored rows use
  # weight = 1 so they contribute normally to the survival function.
  wt <- ifelse(primisol$tm$te == 1L, primisol$tm$morbid, 1)

  fit <- hazard(
    survival::Surv(starttme, int_te, te) ~ 1,
    data    = primisol$tm,
    dist    = "multiphase",
    phases  = list(
      early = hzr_phase("cdf", t_half = 0.0109, nu = 1, m = 1e-6,
                         fixed = c("nu", "m")),
      late  = hzr_phase("g3",  tau = 1, gamma = 1, alpha = 1, eta = 1.4112,
                         fixed = c("tau", "gamma", "alpha"))
    ),
    weights = wt,
    fit     = TRUE,
    control = list(n_starts = 3, maxit = 500, conserve = TRUE)
  )

  # LL: P1 #6 gap produces ~0.43 offset; assert within 0.5.
  expect_equal(fit$fit$objective, ref$loglik, tolerance = 0.5,
               label = "R LL vs SAS (P1 #6 gap: CoE implementation differs)")

  nat       <- ref$natural
  thalf_ref <- nat$estimate[nat$name == "THALF" & nat$phase == "Early"]
  eta_ref   <- nat$estimate[nat$name == "ETA"   & nat$phase == "Late"]

  th <- fit$fit$theta
  expect_equal(unname(exp(th["early.log_t_half"])), thalf_ref, tolerance = 5e-3,
               label = "THALF (Early) natural-scale")
  expect_equal(unname(th["late.eta"]),              eta_ref,   tolerance = 5e-3,
               label = "ETA (Late) natural-scale")

  # NOTE: MUE NOT asserted — P1 #6 CoE gap shifts MUE ~10% (R 0.01643 vs SAS 0.01828).
  # NOTE: NU=1 FIXNU (fixed); M=1e-6 FIXM (fixed) — not compared.
})
