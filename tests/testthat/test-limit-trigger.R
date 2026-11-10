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

test_that("the grid is 3 t_half x 3 nu, covariates zeroed, nothing else moved", {
  start <- c(early.log_mu = -1, early.log_t_half = log(10), early.nu = 0.2,
             early.m = 1, constant.log_mu = -3, early.age = 0.4, constant.age = 0.2)
  tr <- list(list(phase = "early", parameter = "nu", limit = "x", source = "weak"))
  et <- c(1, 2, 3, 4, 5, 6, 7, 8)
  g <- .hzr_limit_grid(start, tr, phases2, event_times = et, fixed = "early.m")
  expect_length(g, 9L)
  th <- vapply(g, `[[`, 0, "early.log_t_half")
  nu <- vapply(g, `[[`, 0, "early.nu")
  expect_setequal(unique(round(exp(th), 6)),
                  round(unname(stats::quantile(et, c(.25, .5, .75))), 6))
  expect_setequal(unique(nu), c(0.5, 1, 2))
  for (s in g) {
    expect_identical(names(s), names(start))
    expect_equal(unname(s[c("early.age", "constant.age")]), c(0, 0))
    expect_equal(unname(s[c("early.log_mu", "early.m", "constant.log_mu")]),
                 unname(start[c("early.log_mu", "early.m", "constant.log_mu")]))
  }
  expect_false(anyDuplicated(lapply(g, unname)) > 0L)
})

test_that("a fixed t_half collapses the grid to nu alone", {
  start <- c(early.log_mu = -1, early.log_t_half = 0, early.nu = 1, early.m = 1,
             constant.log_mu = -3)
  tr <- list(list(phase = "early", parameter = "nu", limit = "x", source = "weak"))
  g <- .hzr_limit_grid(start, tr, phases2, event_times = 1:8,
                       fixed = c("early.m", "early.log_t_half"))
  expect_length(g, 3L)
  expect_true(all(vapply(g, `[[`, 0, "early.log_t_half") == 0))
})

test_that("no events gives no grid", {
  start <- c(early.log_mu = -1, early.log_t_half = 0, early.nu = 1, early.m = 1,
             constant.log_mu = -3)
  tr <- list(list(phase = "early", parameter = "t_half", limit = "x", source = "weak"))
  expect_length(.hzr_limit_grid(start, tr, phases2, numeric(0), "early.m"), 0L)
})

test_that("a g3 trigger varies tau at the quartiles and nothing else", {
  ph <- list(late = hzr_phase("g3", tau = 1, gamma = 1, alpha = 1, eta = 1))
  start <- c(late.log_mu = 0, late.log_tau = 0, late.gamma = 2, late.alpha = 1,
             late.eta = 3)
  tr <- list(list(phase = "late", parameter = "tau", limit = "x", source = "weak"))
  g <- .hzr_limit_grid(start, tr, ph, event_times = 1:8, fixed = character(0))
  expect_length(g, 3L)
  expect_true(all(vapply(g, `[[`, 0, "late.gamma") == 2))
})

test_that("the limit line has a sentence for each state", {
  expect_identical(.hzr_limit_line(NULL), "  Limit check: none")
  expect_identical(.hzr_limit_line(NA), "  Limit check: not run")
  trig <- list(list(phase = "early", parameter = "nu",
                    limit = "0 (a step at t_half)", source = "weak"))
  lim <- list(trigger = trig, outcome = "restarted",
              grid = data.frame(start = 1:2, objective = c(-5, -4)), winner = 2L)
  expect_identical(.hzr_limit_line(lim),
    "  Limit check: restarted (early.nu ran to 0 (a step at t_half)); start 2 of 2 kept")
  kept <- list(trigger = trig, outcome = "kept",
               grid = data.frame(start = 1:2, objective = c(-5, -4)), winner = NA_integer_)
  expect_identical(.hzr_limit_line(kept),
    "  Limit check: restart tried (early.nu ran to 0 (a step at t_half)); none better, original kept")
})

test_that("print() shows the limit line for a multiphase object", {
  obj <- hazard(time = c(1, 2, 3), status = c(1, 0, 1), dist = "multiphase",
    phases = list(constant = hzr_phase("constant")), fit = FALSE)
  expect_match(paste(utils::capture.output(print(obj)), collapse = "\n"),
               "Limit check: not run", fixed = TRUE)
  expect_match(paste(utils::capture.output(print(summary(obj))), collapse = "\n"),
               "Limit check: not run", fixed = TRUE)
})
