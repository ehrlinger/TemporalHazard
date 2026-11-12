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

test_that("the grid's event times include left- and interval-censored events", {
  # Exact events at `time`; a left- or interval-censored row at its upper
  # bound (an interval row's `time` is its lower bound); right-censored rows
  # are not events. Exact events only left an all-interval fit no grid.
  d <- list(time = c(1, 2, 3, 4, 5), time_upper = c(NA, 2.5, 3.5, NA, NA),
            status = c(1, 2, -1, 0, 1))
  expect_identical(.hzr_limit_event_times(d), c(1, 2.5, 3.5, 5))
  d_int <- list(time = c(1, 2), time_upper = c(1.5, 2.5), status = c(2, 2))
  expect_identical(.hzr_limit_event_times(d_int), c(1.5, 2.5))
  expect_gt(length(.hzr_limit_grid(
    c(early.log_mu = 0, early.log_t_half = 0, early.nu = 1, early.m = 0),
    list(list(phase = "early", parameter = "t_half", limit = "x",
              source = "phase_beyond_data")),
    list(early = hzr_phase("cdf", t_half = 1, nu = 1, m = 0)),
    event_times = .hzr_limit_event_times(d_int), fixed = character(0))), 0L)
  # Without time_upper (the exact-only vector path), `time` is used.
  expect_identical(.hzr_limit_event_times(list(time = c(1, 2), status = c(1, 0))), 1)
})

test_that("a row of weight 0 does not place the grid (#737 review)", {
  # The likelihood ignores a zero-weight row, so its event must not move the
  # grid's quantiles. NA weights are left out with it.
  d <- list(time = c(1, 2, 3, 40), status = c(1, 1, 1, 1),
            weights = c(1, 1, 1, 0))
  expect_identical(.hzr_limit_event_times(d), c(1, 2, 3))
  d$weights[2] <- NA
  expect_identical(.hzr_limit_event_times(d), c(1, 3))
  d$weights <- NULL
  expect_identical(.hzr_limit_event_times(d), c(1, 2, 3, 40))
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
  lim$outcome <- "refused"
  expect_identical(.hzr_limit_line(lim), paste0(
    "  Limit check: REFUSED (early.nu ran to 0 (a step at t_half)); ",
    "this is the drifted fit, kept for inspection, not a result"))
  kept <- list(trigger = trig, outcome = "kept",
               grid = data.frame(start = 1:2, objective = c(-5, -4)), winner = NA_integer_)
  expect_identical(.hzr_limit_line(kept),
    "  Limit check: restart tried (early.nu ran to 0 (a step at t_half)); none better, original kept")
  skipped <- list(trigger = trig, outcome = "skipped", grid = NULL,
                  winner = NA_integer_)
  expect_identical(.hzr_limit_line(skipped),
    paste0("  Limit check: restart not tried (early.nu ran to 0 (a step at ",
           "t_half)); control$limit_restart = FALSE"))
})

test_that("a fit saved before the limit check prints \"not run\", not \"none\"", {
  # A 1.3.0 multiphase object has no $fit$limit element, which reads NULL:
  # "none" would claim it was examined and found clean.
  obj <- hazard(time = c(1, 2, 3), status = c(1, 0, 1), dist = "multiphase",
    phases = list(constant = hzr_phase("constant")), fit = FALSE)
  obj$fit$limit <- NULL
  expect_false("limit" %in% names(obj$fit))
  expect_match(paste(utils::capture.output(print(obj)), collapse = "\n"),
               "Limit check: not run", fixed = TRUE)
  expect_match(paste(utils::capture.output(print(summary(obj))), collapse = "\n"),
               "Limit check: not run", fixed = TRUE)
  # A fit that was examined and found clean still reads "none".
  obj$fit["limit"] <- list(NULL)
  expect_true("limit" %in% names(obj$fit))
  expect_match(paste(utils::capture.output(print(obj)), collapse = "\n"),
               "Limit check: none", fixed = TRUE)
})

test_that("print() shows the limit line for a multiphase object", {
  obj <- hazard(time = c(1, 2, 3), status = c(1, 0, 1), dist = "multiphase",
    phases = list(constant = hzr_phase("constant")), fit = FALSE)
  expect_match(paste(utils::capture.output(print(obj)), collapse = "\n"),
               "Limit check: not run", fixed = TRUE)
  expect_match(paste(utils::capture.output(print(summary(obj))), collapse = "\n"),
               "Limit check: not run", fixed = TRUE)
})

# #734: a fitted multiphase fit whose Hessian was unusable has $weak NA, and
# that NA names no parameter. It triggers weak-only, once per free G1 t_half.
mk_na <- function(phases = phases2, theta = NULL, fixed = NULL, boundary = NULL,
                  weak_cause = "covariance is not positive definite",
                  se_cause = NULL, converged = TRUE, objective = -10,
                  degraded = c("standard_errors", "weak_direction_check")) {
  th <- theta %||% c(early.log_mu = 0, early.log_t_half = 0, early.nu = 1,
                     early.m = 1, constant.log_mu = 0, early.age = 0.1)
  causes <- c(standard_errors = se_cause %||% "no finite positive variance for: early.m",
              weak_direction_check = weak_cause)[degraded]
  list(spec = list(dist = "multiphase", phases = phases),
       fit = list(theta = th, boundary = boundary, weak = NA,
                  converged = converged, objective = objective,
                  fixed_mask = if (is.null(fixed)) NULL else names(th) %in% fixed),
       degraded = degraded, degraded_causes = causes)
}
na_limit <- "no usable Hessian (identification not checked)"

test_that("an unusable Hessian on a fitted fit triggers weak-only (#734)", {
  for (cause in c("covariance is not positive definite",
                  "Hessian not positive-definite at the optimum",
                  "covariance has non-finite entries",
                  "standard errors unavailable")) {
    tr <- .hzr_limit_trigger(mk_na(weak_cause = cause))
    expect_length(tr, 1L)
    expect_identical(tr[[1]]$phase, "early")
    expect_identical(tr[[1]]$parameter, "t_half")
    expect_identical(tr[[1]]$limit, na_limit)
    expect_identical(tr[[1]]$source, "weak_na")
    expect_false(.hzr_limit_is_hard(tr))
  }
})

test_that("the NA trigger names every free G1 t_half and nothing else", {
  ph <- list(early = hzr_phase("cdf", t_half = 1, nu = 1, m = 1),
             mid = hzr_phase("hazard", t_half = 1, nu = 1, m = 1),
             constant = hzr_phase("constant"),
             late = hzr_phase("g3", tau = 1, gamma = 1, alpha = 1, eta = 1))
  th <- c(early.log_mu = 0, early.log_t_half = 0, early.nu = 1, early.m = 1,
          mid.log_mu = 0, mid.log_t_half = 0, mid.nu = 1, mid.m = 1,
          constant.log_mu = 0, late.log_mu = 0, late.log_tau = 0,
          late.gamma = 1, late.alpha = 1, late.eta = 1)
  tr <- .hzr_limit_trigger(mk_na(phases = ph, theta = th))
  expect_identical(vapply(tr, `[[`, "", "phase"), c("early", "mid"))
  expect_true(all(vapply(tr, `[[`, "", "parameter") == "t_half"))
  # A fixed t_half is the user's choice, and a vanished phase is a different
  # limit (#565): neither is named.
  tr <- .hzr_limit_trigger(mk_na(phases = ph, theta = th, fixed = "mid.log_t_half"))
  expect_identical(vapply(tr, `[[`, "", "phase"), "early")
  van <- list(list(mechanism = "coe_phase_vanished", phase = "early"))
  tr <- .hzr_limit_trigger(mk_na(phases = ph, theta = th, boundary = van))
  expect_identical(vapply(tr, `[[`, "", "phase"), "mid")
})

test_that("a fixed t_half with nu and m free gets no NA record", {
  # The NA reading names t_half, and a fixed t_half is the user's choice, so
  # a phase whose only fixed shape is t_half is not named even though nu and
  # m are free (#734).
  ph <- list(early = hzr_phase("cdf", t_half = 1, nu = 1, m = 1),
             constant = hzr_phase("constant"))
  expect_length(.hzr_limit_trigger(mk_na(phases = ph)), 1L)
  expect_null(.hzr_limit_trigger(mk_na(phases = ph, fixed = "early.log_t_half")))
})

test_that("a hard record leaves the NA reading out, so the hard path is unchanged", {
  # The NA reading adds records only to a fit nothing else triggers on: a
  # hard trigger's grid, refusal text and rescue are what they were (#734).
  tr <- .hzr_limit_trigger(mk_na(
    boundary = list(list(mechanism = "phase_discontinuity", phase = "early"))))
  expect_identical(vapply(tr, `[[`, "", "source"), "phase_discontinuity")
  expect_true(.hzr_limit_is_hard(tr))
})

test_that("unusable = FALSE reads an NA weak as before #734", {
  # How a hard or flat-direction trigger judges its grid starts: a start with
  # an unusable Hessian could rescue them before #734, and still can.
  expect_null(.hzr_limit_trigger(mk_na(), unusable = FALSE))
  expect_length(.hzr_limit_trigger(mk_na()), 1L)
})

test_that("an NA weak that is not an unusable Hessian does not trigger", {
  # fit = FALSE: no objective, and "fitting" is on the record.
  expect_null(.hzr_limit_trigger(mk_na(objective = NA_real_, converged = NA,
    degraded = c("fitting", "standard_errors", "weak_direction_check"))))
  # Imported from SAS: fitting on the record.
  expect_null(.hzr_limit_trigger(mk_na(
    degraded = c("fitting", "standard_errors", "weak_direction_check"))))
  # The fit did not converge: there is no optimum to take a Hessian at.
  expect_null(.hzr_limit_trigger(mk_na(converged = FALSE, objective = NA_real_)))
  # numDeriv not installed: the check could not be made, which says nothing
  # about the fit, so the estimates must not depend on an installed Suggests.
  expect_null(.hzr_limit_trigger(mk_na(weak_cause = "standard errors unavailable",
    se_cause = "numDeriv not installed and no analytic Hessian")))
  expect_null(.hzr_limit_trigger(mk_na(weak_cause = "standard errors unavailable",
    se_cause = paste0("Conservation-of-Events full-information variance could ",
                      "not be computed (numDeriv not installed)"))))
  # Every shape fixed: no free t_half to name.
  expect_null(.hzr_limit_trigger(mk_na(fixed = c("early.log_t_half", "early.nu", "early.m"))))
  # Not multiphase.
  expect_null(.hzr_limit_trigger(list(spec = list(dist = "weibull", phases = NULL),
    fit = list(theta = c(mu = 0, nu = 1), weak = NA, converged = TRUE, objective = -1),
    degraded = "weak_direction_check")))
})

test_that("the NA trigger reads sensibly in the limit line", {
  trig <- list(list(phase = "early", parameter = "t_half", limit = na_limit,
                    source = "weak_na"))
  kept <- list(trigger = trig, outcome = "kept",
               grid = data.frame(start = 1:2, objective = c(-5, -4)), winner = NA_integer_)
  expect_identical(.hzr_limit_line(kept), paste0(
    "  Limit check: restart tried (early.t_half: no usable Hessian ",
    "(identification not checked)); none better, original kept"))
})
