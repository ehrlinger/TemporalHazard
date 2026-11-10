# tests/testthat/test-limit-trigger.R
# #683: which fits the restart path acts on. Hand-built objects: fast.

phases2 <- list(early = hzr_phase("cdf", t_half = 1, nu = 1, m = 1, fixed = "m"),
                constant = hzr_phase("constant"))
mk <- function(boundary = NULL, weak = NULL, fixed = NULL) {
  th <- c(early.log_mu = 0, early.log_t_half = 0, early.nu = 1, early.m = 1,
          constant.log_mu = 0, early.age = 0.1)
  list(spec = list(dist = "multiphase", phases = phases2),
       fit = list(theta = th, boundary = boundary, weak = weak,
                  fixed_mask = if (is.null(fixed)) NULL else names(th) %in% fixed))
}

test_that("a clean fit does not trigger", {
  expect_null(.hzr_limit_trigger(mk()))
  expect_null(.hzr_limit_trigger(mk(boundary = NA, weak = NA)))
})

test_that("each boundary mechanism maps to its parameter and limit", {
  cases <- list(
    phase_beyond_data   = c("t_half", "Inf (the rise lies past the data)"),
    unbounded_phase     = c("t_half", "0 (below the first observed time)"),
    phase_discontinuity = c("nu", "0 (a step at t_half)"))
  for (mech in names(cases)) {
    tr <- .hzr_limit_trigger(mk(boundary = list(list(mechanism = mech, phase = "early"))))
    expect_length(tr, 1L)
    expect_identical(tr[[1]]$phase, "early")
    expect_identical(tr[[1]]$parameter, cases[[mech]][[1]])
    expect_identical(tr[[1]]$limit, cases[[mech]][[2]])
    expect_identical(tr[[1]]$source, mech)
  }
})

test_that("records that are not limits do not trigger", {
  for (mech in c("g3_alpha_one", "g3_fixge2_alpha_start", "g3_corner_supremum",
                 "coe_no_events_left", "coe_phase_vanished")) {
    expect_null(.hzr_limit_trigger(mk(boundary = list(list(mechanism = mech,
                                                           phase = "early")))))
  }
})

test_that("a weak direction triggers only when it names a shape", {
  tr <- .hzr_limit_trigger(mk(weak = list(params = c("early.log_mu", "early.nu"))))
  expect_identical(vapply(tr, `[[`, "", "parameter"), "nu")
  expect_identical(tr[[1]]$source, "weak")
  # A covariate ridge is collinearity, not a limit.
  expect_null(.hzr_limit_trigger(mk(weak = list(params = c("early.age", "constant.log_mu")))))
})

test_that("a covariate NAMED like a shape is not read as one", {
  ph <- list(early = hzr_phase("cdf", t_half = 1, nu = 1, m = 1))
  obj <- list(spec = list(dist = "multiphase", phases = ph),
              fit = list(theta = c(early.log_mu = 0, early.log_t_half = 0,
                                   early.nu = 1, early.m = 1, early.nu_1 = 0),
                         weak = list(params = "early.nu_1"), boundary = NULL))
  expect_null(.hzr_limit_trigger(obj))
})

test_that("a FIXED parameter never triggers", {
  expect_null(.hzr_limit_trigger(mk(weak = list(params = "early.m"), fixed = "early.m")))
  expect_null(.hzr_limit_trigger(mk(
    boundary = list(list(mechanism = "phase_discontinuity", phase = "early")),
    fixed = "early.nu")))
})

test_that("a weak reading in a VANISHED phase does not trigger", {
  van <- list(list(mechanism = "coe_phase_vanished", phase = "early"))
  expect_null(.hzr_limit_trigger(mk(boundary = van,
                                    weak = list(params = c("early.log_t_half", "early.nu")))))
  # but a limit record in that phase still does
  tr <- .hzr_limit_trigger(mk(boundary = c(van, list(list(mechanism = "unbounded_phase",
                                                            phase = "early")))))
  expect_identical(tr[[1]]$source, "unbounded_phase")
})

test_that("one parameter reached two ways is reported once", {
  tr <- .hzr_limit_trigger(mk(
    boundary = list(list(mechanism = "phase_discontinuity", phase = "early")),
    weak = list(params = c("early.nu", "early.log_mu"))))
  expect_length(tr, 1L)
})
