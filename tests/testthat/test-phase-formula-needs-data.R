# A phase formula with covariates needs `data` (#299).
#
# hazard() builds a phase's design from its own formula only when `data` is
# supplied. On the vector interface without `data` the formula was ignored:
# the phase took the global `x`, or no columns at all, and the fit was a
# different model with no warning and no message. It is refused now.

avc_299 <- function() {
  stats::na.omit(get("avc", envir = asNamespace("TemporalHazard"))[
    , c("int_dead", "dead", "age", "mal")])
}
phases_299 <- function(early_formula = ~ mal) {
  list(early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m",
                         formula = early_formula),
       constant = hzr_phase("constant"))
}
ctl_299 <- list(n_starts = 1L, conserve = FALSE)

test_that("a phase formula with covariates and no `data` is refused", {
  d <- avc_299()
  int_dead <- d$int_dead
  dead <- d$dead
  mal <- d$mal # visible from the calling frame, which is the trap
  for (fit in c(TRUE, FALSE)) {
    expect_error(
      hazard(time = int_dead, status = dead, dist = "multiphase",
             phases = phases_299(), fit = fit, control = ctl_299),
      "Phase 'early' has the formula `~mal`.*`data =`",
      label = paste("fit =", fit)
    )
  }
})

test_that("the refusal covers a phase formula beside a global `x`", {
  # The phase inherited the global x instead: `~ log(age)` fitted as `age`.
  d <- avc_299()
  expect_error(
    hazard(time = d$int_dead, status = d$dead, x = cbind(age = d$age),
           dist = "multiphase", phases = phases_299(~ log(age)),
           fit = FALSE),
    "Phase 'early' has the formula `~log\\(age\\)`"
  )
})

test_that("`~ 1` beside a global `x`, or a constant term, is refused too", {
  # Without `data` a `~ 1` phase took the global x (early.age, logLik
  # -196.44); with `data` it has no columns (-211.72). A constant term such
  # as log(2) builds a column only in `data`.
  d <- avc_299()
  expect_error(
    hazard(time = d$int_dead, status = d$dead, x = cbind(age = d$age),
           dist = "multiphase", phases = phases_299(~ 1), fit = FALSE),
    "Phase 'early' has the formula `~1`.*take the global `x`.*drop `x`"
  )
  expect_error(
    hazard(time = d$int_dead, status = d$dead, dist = "multiphase",
           phases = phases_299(~ log(2)), fit = FALSE),
    "Phase 'early' has the formula `~log\\(2\\)`"
  )
})

test_that("an intercept-only or absent phase formula needs no `data`", {
  d <- avc_299()
  for (form in list(~ 1, NULL)) {
    f <- hazard(time = d$int_dead, status = d$dead, dist = "multiphase",
                phases = phases_299(form), fit = FALSE)
    expect_s3_class(f, "hazard")
    expect_identical(f$spec$phases$early$formula, form)
  }
})

test_that("the same fit with `data` is unchanged", {
  skip_on_cran() # multiphase fits
  d <- avc_299()
  # Pinned against main 9ec83e7, where the issue's control gave these.
  fv <- hazard(time = int_dead, status = dead, data = d, dist = "multiphase",
               phases = phases_299(), fit = TRUE, control = ctl_299)
  ff <- hazard(survival::Surv(int_dead, dead) ~ 1, data = d,
               dist = "multiphase", phases = phases_299(), fit = TRUE,
               control = ctl_299)
  for (f in list(fv, ff)) {
    expect_identical(NCOL(f$fit$x_list$early), 1L)
    expect_equal(f$fit$objective, -202.172323, tolerance = 1e-8)
    expect_equal(coef(f)[["early.mal"]], 1.10524535, tolerance = 1e-6)
  }
})

test_that("a fit that uses its phase formula records it and predicts from it", {
  skip_on_cran() # multiphase fits
  d <- avc_299()
  fits <- list(
    vector = hazard(time = int_dead, status = dead, data = d,
                    dist = "multiphase", phases = phases_299(), fit = TRUE,
                    control = ctl_299),
    formula = hazard(survival::Surv(int_dead, dead) ~ 1, data = d,
                     dist = "multiphase", phases = phases_299(), fit = TRUE,
                     control = ctl_299)
  )
  nd <- data.frame(time = c(0.1, 0.1), mal = c(0, 1))
  for (nm in names(fits)) {
    f <- fits[[nm]]
    expect_identical(unname(attr(f$fit$x_list, "from_formula")["early"]),
                     TRUE, label = nm)
    p <- predict(f, newdata = nd, type = "cumulative_hazard",
                 decompose = TRUE)
    # mal moves the early phase by exactly exp(beta) and nothing else.
    expect_gt(p$total[2], p$total[1])
    expect_equal(p$early[2] / p$early[1], exp(coef(f)[["early.mal"]]),
                 tolerance = 1e-8, label = nm)
    expect_equal(p$constant[2], p$constant[1], tolerance = 1e-12, label = nm)
  }
})

test_that("the saved-fit helper rebuilds what hazard() returned before #299", {
  skip_on_cran() # multiphase fits
  # With the refusal mocked away, hazard() takes the old path: the design
  # builder gets the phase formula and no `data`, and ignores it. The helper
  # has to rebuild that object, or every test built on it stops simulating
  # the saved fits it claims to.
  d <- avc_299()
  for (x in list(NULL, cbind(age = d$age))) {
    ph <- phases_299(~ log(age))
    set.seed(1)
    rebuilt <- suppressWarnings(hzr_saved_before_299(
      ph, time = d$int_dead, status = d$dead, x = x, dist = "multiphase",
      fit = TRUE, control = ctl_299))
    set.seed(1)
    old <- local({
      local_mocked_bindings(
        .hzr_check_phase_formula_data = function(...) invisible(NULL))
      suppressWarnings(hazard(
        time = d$int_dead, status = d$dead, x = x, dist = "multiphase",
        phases = ph, fit = TRUE, control = ctl_299))
    })
    # The old object really did ignore the formula...
    expect_false("early.log(age)" %in% names(coef(old)))
    expect_identical(old$spec$phases$early$formula, ~ log(age))
    # ...and the rebuilt one is that object.
    expect_setequal(names(rebuilt), names(old))
    for (el in setdiff(names(old), c("call", "call_env"))) {
      expect_equal(rebuilt[[el]], old[[el]], label = el)
    }
  }
})

test_that("the blocker leaves alone a record-less fit that used its formula", {
  skip_on_cran() # a multiphase fit
  # A phase formula whose columns carry the inherited names looks inherited
  # to the record-less fallback (its documented limit). The call's `data`
  # is what clears it: that fit did evaluate its phase formulas.
  d <- avc_299()
  f <- hazard(survival::Surv(int_dead, dead) ~ age, data = d,
              dist = "multiphase",
              phases = list(
                early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1,
                                  fixed = "m", formula = ~ age),
                constant = hzr_phase("constant", formula = ~ age)),
              fit = TRUE, control = ctl_299)
  attr(f$fit$x_list, "from_formula") <- NULL
  expect_true(.hzr_phase_inherits_global(f, "early"))
  expect_null(.hzr_inherit_blocker(f))
})

test_that("hzr_stepwise() refuses a saved fit whose phase formula was ignored", {
  skip_on_cran() # multiphase fits
  # Given `data`, every refit would build the ignored formula's columns into
  # a base model that never had them: on main the screen reported "ENTER
  # age" over a final model that also carried mal.
  # Simulates fits saved before #299, which hazard() can no longer make.
  d <- avc_299()
  hollow <- suppressWarnings(hzr_saved_before_299(
    phases_299(~ mal), time = d$int_dead, status = d$dead,
    dist = "multiphase", fit = TRUE, control = ctl_299))
  beside_x <- suppressWarnings(hzr_saved_before_299(
    phases_299(~ log(age)), time = d$int_dead, status = d$dead,
    x = cbind(age = d$age), dist = "multiphase", fit = TRUE,
    control = ctl_299))
  # A fit saved before the from_formula record is judged by its call and
  # columns instead.
  no_record <- hollow
  attr(no_record$fit$x_list, "from_formula") <- NULL
  # `~ 1` beside the global x took x; a refit given `data` would drop it.
  one_beside_x <- suppressWarnings(hzr_saved_before_299(
    phases_299(~ 1), time = d$int_dead, status = d$dead,
    x = cbind(age = d$age), dist = "multiphase", fit = TRUE,
    control = ctl_299))
  for (f in list(hollow, beside_x, no_record, one_beside_x)) {
    expect_error(
      hzr_stepwise(f, scope = list(early = ~ age), data = d,
                   direction = "forward", criterion = "wald",
                   trace = FALSE),
      "phase 'early' has a formula, `.+`, that the fit ignored"
    )
  }
  # hzr_bootstrap() already refuses both: neither stores a `data` frame.
  expect_error(hzr_bootstrap(hollow, n_boot = 2L), "vector interface")
  expect_error(hzr_bootstrap(beside_x, n_boot = 2L), "`x`")
})
