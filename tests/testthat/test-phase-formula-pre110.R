# A multiphase fit saved before 1.1.0 stores neither its data frame nor a
# record of whether a phase formula was used (#324). On the vector interface
# with a call naming `data`, the call alone cannot say whether that `data`
# was a frame or a symbol bound to NULL, so such a fit whose phase looks
# inherited is refused rather than judged.

avc_324 <- function() {
  stats::na.omit(get("avc", envir = asNamespace("TemporalHazard"))[
    , c("int_dead", "dead", "age", "mal")])
}
ctl_324 <- list(n_starts = 1L, conserve = FALSE)
pre_110 <- function(f) {
  attr(f$fit$x_list, "from_formula") <- NULL
  f$data["frame"] <- NULL
  f
}
screen_324 <- function(fit, d) {
  hzr_stepwise(fit, scope = list(constant = ~ mal), data = d,
               direction = "forward", criterion = "wald", trace = FALSE)
}
undecidable <- "saved by a version before 1.1.0"

test_that("a pre-1.1.0 fit whose `data =` symbol was NULL is refused (#324)", {
  skip_on_cran() # multiphase fits
  d <- avc_324()
  dd <- NULL
  f <- pre_110(suppressWarnings(hzr_saved_before_299(
    list(early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1,
                           fixed = "m", formula = ~ mal),
         constant = hzr_phase("constant")),
    time = d$int_dead, status = d$dead, data = dd,
    dist = "multiphase", fit = TRUE, control = ctl_324)))
  # The helper forwards through `...`, so restore the call as it was written.
  f$call$data <- quote(dd)
  # The fixture is the silent route: no frame, no record, a call naming `dd`.
  expect_false("frame" %in% names(f$data))
  expect_null(attr(f$fit$x_list, "from_formula"))
  expect_identical(f$call$data, quote(dd))

  expect_error(screen_324(f, d), undecidable, fixed = TRUE)
  expect_error(screen_324(f, d), "phase 'early' has a formula, `~mal`",
               fixed = TRUE)
  expect_error(hzr_bootstrap(f, n_boot = 2L), undecidable, fixed = TRUE)
})

test_that("a pre-1.1.0 real-data fit whose phase columns are the inherited ones is refused too (#324)", {
  skip_on_cran() # a multiphase fit
  # The accepted false positive. This fit did use its phase formula, but its
  # columns are exactly the global `x` it would otherwise have inherited, so
  # without a frame or a record it is indistinguishable from the NULL case.
  # Refusing it costs a refit; passing it would let the NULL case through.
  d <- avc_324()
  f <- pre_110(hazard(
    time = d$int_dead, status = d$dead, x = cbind(mal = d$mal), data = d,
    dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1,
                                    fixed = "m", formula = ~ mal),
                  constant = hzr_phase("constant", formula = ~ age)),
    fit = TRUE, control = ctl_324))
  expect_true(.hzr_phase_inherits_global(f, "early"))
  # Every phase has its own formula, so no other blocker refuses it: before
  # #324 this screen ran.
  expect_error(
    hzr_stepwise(f, scope = list(early = ~ age), data = d,
                 direction = "forward", criterion = "wald", trace = FALSE),
    undecidable, fixed = TRUE
  )
})

test_that("the same fits with a stored frame, or on the formula interface, are not refused (#324)", {
  skip_on_cran() # multiphase fits
  d <- avc_324()
  phases <- list(early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1,
                                   fixed = "m", formula = ~ mal),
                 constant = hzr_phase("constant"))
  real <- hazard(time = d$int_dead, status = d$dead, x = cbind(mal = d$mal),
                 data = d, dist = "multiphase", phases = phases, fit = TRUE,
                 control = ctl_324)
  # Saved by 1.1.0 or later without the record: the frame decides.
  with_frame <- real
  attr(with_frame$fit$x_list, "from_formula") <- NULL
  expect_true("frame" %in% names(with_frame$data))
  expect_null(.hzr_ignored_phase_formula(with_frame))

  # A pre-1.1.0 vector fit whose phase columns are not the inherited ones is
  # judged by them: it used its formula, so it is not refused.
  no_x <- pre_110(hazard(
    time = d$int_dead, status = d$dead, data = d, dist = "multiphase",
    phases = phases, fit = TRUE, control = ctl_324))
  expect_false(.hzr_phase_inherits_global(no_x, "early"))
  expect_null(.hzr_ignored_phase_formula(no_x))

  # A pre-1.1.0 formula-interface fit always had `data`.
  formula_fit <- pre_110(hazard(
    survival::Surv(int_dead, dead) ~ mal, data = d, dist = "multiphase",
    phases = phases, fit = TRUE, control = ctl_324))
  expect_true(.hzr_phase_inherits_global(formula_fit, "early"))
  expect_null(.hzr_ignored_phase_formula(formula_fit))
})
