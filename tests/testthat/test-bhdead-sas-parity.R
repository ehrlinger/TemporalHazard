# SAS parity for the bh.dead analysis: multiphase shape fit (%HAZARD) and the
# hzr_bootstrap(scope=) embedded stepwise screen (%HAZBOOT).
#
# THIS REPOSITORY IS PUBLIC. No study value appears in this file: the data and
# the SAS reference both live on a secure volume, and every expectation is read
# from the fixture at run time. The test skips unless both are readable, so it
# never runs on CRAN, on CI, or on any machine without the volume mounted.
#
# Enable by building the fixture first:
#   source("inst/dev/bhdead-parity/parse-bhdead-lst.R"); .hzr_write_bhdead_fixture()
# See inst/dev/bhdead-parity/README.md.

.bhdead_helpers <- system.file("dev", "bhdead-parity", "bhdead-fixture.R",
                                package = "TemporalHazard")

# Shape-fit tolerance is calibrated from observed output: R reproduces SAS's
# converged shape fit to 3.4e-04 max relative error (2026-07-16), so 1e-3 is a
# real guard with headroom, not a guess.
#
# There is deliberately NO statistical tolerance for the bootstrap screen. SAS
# seeds from the time of day (seed=-1), so its selection frequencies are one
# random realisation; and one R replicate over the 92-variable pool costs ~25
# minutes, making a SAS-scale n_boot infeasible until criterion = "score"
# lands. At the feasible n_boot the frequencies are too coarse to support a
# meaningful floor. See BHDEAD-SAS-PARITY-DESIGN.md, "Deferred: score-test
# selection". The bootstrap test below is a smoke check by design.
.bhdead_tolerance <- list(
  shape_rel = 1e-3
)

.bhdead_setup <- function() {
  testthat::skip_if_not_installed("haven")
  if (!nzchar(.bhdead_helpers) || !file.exists(.bhdead_helpers)) {
    testthat::skip("bh.dead parity helpers not installed")
  }
  # Source into the calling test_that() frame (pf) so the sourced helpers
  # (e.g. .hzr_bhdead_candidates) stay visible to the rest of the test body.
  # Call them via pf$... here because unqualified calls from inside this
  # function resolve lexically against this function's own environment, not
  # against pf, even though the assignments landed in pf.
  pf <- parent.frame()
  source(.bhdead_helpers, local = pf)
  fix <- pf$.hzr_load_bhdead_fixture()
  testthat::skip_if(is.null(fix), "bh.dead SAS fixture not available")
  data_path <- pf$.hzr_bhblt_path()
  testthat::skip_if_not(file.exists(data_path), "bhblt data not available")
  caa <- haven::read_sas(data_path)
  names(caa) <- tolower(names(caa))
  list(fix = fix, caa = caa)
}

# Build the phase list for a given `fixed` specification, seeding every start
# value from the fixture.
.bhdead_phases <- function(sp, fixed) {
  list(
    early    = hzr_phase("cdf", t_half = sp[["thalf"]], nu = sp[["nu"]],
                         m = sp[["m"]], fixed = fixed),
    constant = hzr_phase("constant")
  )
}

.bhdead_theta <- function(sp) {
  c(early.log_mu     = log(sp[["mue"]]),
    early.log_t_half = log(sp[["thalf"]]),
    early.nu         = sp[["nu"]],
    early.m          = sp[["m"]],
    constant.log_mu  = log(sp[["muc"]]))
}

test_that("multiphase shape fit reproduces the SAS shape fit", {
  env <- .bhdead_setup()
  sp <- env$fix$shape$specified
  cv <- env$fix$shape$converged

  # hz.dead fixes only M; THALF and NU are estimated. noconserve -> conserve = FALSE.
  fit <- hazard(
    survival::Surv(iv_dead, dead) ~ 1, data = env$caa, dist = "multiphase",
    phases  = .bhdead_phases(sp, fixed = "m"),
    theta   = .bhdead_theta(sp),
    control = list(conserve = FALSE),
    fit     = TRUE
  )
  expect_true(isTRUE(fit$fit$converged))

  th <- coef(fit)
  expect_equal(exp(th[["early.log_t_half"]]), cv[["thalf"]],
               tolerance = .bhdead_tolerance$shape_rel)
  expect_equal(th[["early.nu"]], cv[["nu"]],
               tolerance = .bhdead_tolerance$shape_rel)
  expect_equal(exp(th[["early.log_mu"]]), cv[["mue"]],
               tolerance = .bhdead_tolerance$shape_rel)
  expect_equal(exp(th[["constant.log_mu"]]), cv[["muc"]],
               tolerance = .bhdead_tolerance$shape_rel)
})

test_that("every SAS candidate exists in the data", {
  env <- .bhdead_setup()
  expect_identical(
    setdiff(.hzr_bhdead_candidates(env$fix), names(env$caa)),
    character(0)
  )
})

test_that("bootstrap screen wiring: fixture scope -> hzr_bootstrap -> summary", {
  env <- .bhdead_setup()
  sp <- env$fix$shape$specified

  # The full 92-variable SAS candidate pool costs ~25 minutes per bootstrap
  # replicate (hzr_stepwise refits once per candidate per step), so a
  # SAS-scale run belongs in the analysis notebook, not the test suite. See
  # inst/dev/BHDEAD-SAS-PARITY-DESIGN.md, "Deferred: score-test selection".
  # Here we slice to the first 3 candidates and n_boot = 1 -- enough to
  # exercise fixture scope -> hzr_bootstrap -> summary wiring without the
  # runtime cost.
  scope_vars <- head(.hzr_bhdead_candidates(env$fix), 3)

  # bh.dead fixes nu and m; THALF stays free.
  base_fit <- hazard(
    survival::Surv(iv_dead, dead) ~ 1, data = env$caa, dist = "multiphase",
    phases  = .bhdead_phases(sp, fixed = c("nu", "m")),
    theta   = .bhdead_theta(sp),
    control = list(conserve = FALSE),
    fit     = TRUE
  )

  bs <- hzr_bootstrap(
    base_fit, n_boot = 1, seed = 123,
    scope   = list(early = stats::reformulate(scope_vars),
                   constant = stats::reformulate(scope_vars)),
    slentry = env$fix$meta$sle, slstay = env$fix$meta$sls,
    control = list(n_starts = 1, conserve = FALSE)
  )

  expect_s3_class(bs, "hzr_bootstrap")
  expect_gt(bs$n_success, 0)
  expect_gt(nrow(bs$summary), 0)
  # The intercepts are never dropped by stepwise, so they appear in every
  # successful replicate.
  intercepts <- bs$summary[bs$summary$parameter %in%
                             c("early.log_mu", "constant.log_mu"), ]
  expect_true(all(intercepts$n == bs$n_success))
})
