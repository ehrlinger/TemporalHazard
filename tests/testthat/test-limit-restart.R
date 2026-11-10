# tests/testthat/test-limit-restart.R
skip_on_cran()
# #683/#729: a shape parameter running to its limit is restarted, or refused;
# the drifted iterate is never returned. hm.death.patient blind start.

valves_729 <- function() {
  data(valves, package = "TemporalHazard", envir = environment())
  d <- valves
  d$nyha[is.na(d$nyha)] <- mean(d$nyha, na.rm = TRUE)
  d$mv_nyha <- d$mitral * d$nyha
  d
}
fit_729 <- function(n_starts = 1L, limit_restart = TRUE, d = valves_729()) {
  hazard(survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 10, nu = 0.2, m = 1, fixed = "m",
                    formula = ~ age_cop + nyha + mv_nyha + double_ + ao_pinc),
                  constant = hzr_phase("constant", formula = ~ age_cop + black + i_path)),
    fit = TRUE,
    control = list(n_starts = n_starts, maxit = 2000, conserve = TRUE,
                   limit_restart = limit_restart))
}
sas_ll_729 <- -1792.07   # hm.death.patient.lst, asserted in test-sas-parity.R

for (ns in c(1L, 5L)) {
  test_that(paste0("#729 blind start, n_starts = ", ns, ", reaches SAS's optimum"), {
    expect_warning(fit <- fit_729(ns), class = "hzr_limit_restarted")
    expect_equal(fit$fit$objective / sas_ll_729, 1, tolerance = 1e-5)
    lim <- fit$fit$limit
    expect_type(lim, "list")
    expect_identical(lim$outcome, "restarted")
    # The drift it recovered from is on the record, and is the bad one.
    expect_lt(lim$original$objective, sas_ll_729 - 10)
    expect_true(any(vapply(lim$trigger, `[[`, "", "phase") == "early"))
    # The grid ran, varied, and the winner is a grid row.
    expect_gte(nrow(lim$grid), 2L)
    expect_gt(length(unique(round(lim$grid$objective[is.finite(lim$grid$objective)], 4))), 1L)
    expect_false(lim$grid$triggered[[lim$winner]])
    expect_equal(lim$grid$objective[[lim$winner]], fit$fit$objective)
    # A rescued fit is an ordinary fit, and is not itself at a limit.
    expect_null(.hzr_limit_trigger(fit))
    expect_true(isTRUE(fit$fit$converged))
    # The fixed shape (m) has no variance; every free parameter does.
    free <- !as.logical(fit$fit$fixed_mask)
    expect_true(any(free))
    expect_true(all(is.finite(sqrt(diag(vcov(fit))[free]))))
    # Its call is the user's, so hzr_bootstrap() refits from the user's start.
    expect_null(fit$call$theta)
    # The grid forces n_starts = 1; the stored control is the user's, which
    # hzr_bootstrap() and hzr_stepwise() refit from.
    expect_identical(fit$spec$control$n_starts, ns)
    expect_identical(.hzr_inherit_control(list(), fit)$control$n_starts, ns)
  })
}

test_that("limit_restart = FALSE refuses at once, naming the parameter", {
  # n_starts = 5 lands on phase_beyond_data, a hard trigger.
  err <- tryCatch(suppressWarnings(fit_729(5L, limit_restart = FALSE)),
                  hzr_limit_refusal = function(e) e)
  expect_s3_class(err, "hzr_limit_refusal")
  expect_identical(err$phase, "early")
  expect_identical(err$parameter, "t_half")
  expect_true(all(vapply(err$trigger, `[[`, "", "source") != "weak"))
  expect_true(any(vapply(err$trigger, `[[`, "", "source") == "phase_beyond_data"))
  expect_match(conditionMessage(err), "hzr_phase\\(")
  expect_match(conditionMessage(err), "FIXTHALF", fixed = TRUE)
  expect_lt(err$best_objective, sas_ll_729 - 10)
})

test_that("limit_restart = FALSE never refits a weak-only trigger", {
  fit <- suppressWarnings(fit_729(1L, limit_restart = FALSE))
  expect_lt(fit$fit$objective, sas_ll_729 - 100)   # the drifted fit, about -2036.66
  expect_null(fit$fit$limit)
})

weak_only_fit <- function(limit_restart = TRUE) {
  set.seed(42)
  n <- 100
  time   <- rexp(n, rate = 0.5) + 0.01
  status <- sample(0:1, n, replace = TRUE, prob = c(0.2, 0.8))
  phases <- list(early = hzr_phase("cdf", t_half = 1, nu = 1.5, m = 0),
                 const = hzr_phase("constant"))
  hazard(time = time, status = status, dist = "multiphase",
         phases = phases, fit = TRUE,
         control = list(n_starts = 5, maxit = 500, limit_restart = limit_restart))
}

test_that("a weak-only trigger tries the grid and keeps the original", {
  # A ridge (weak flag only): the grid finds nothing better, so the original
  # fit is returned, with the tried grid on the record. Never a refusal.
  fit <- NULL
  w <- testthat::capture_warnings(fit <- weak_only_fit())
  expect_s3_class(fit, "hazard")
  expect_false(any(grepl("Restarted from", w)))
  lim <- fit$fit$limit
  expect_type(lim, "list")
  expect_identical(lim$outcome, "kept")
  expect_true(is.na(lim$winner))
  expect_gte(nrow(lim$grid), 2L)
  expect_true(all(vapply(lim$trigger, `[[`, "", "source") == "weak"))
  # The same fit as with the layer off.
  off <- suppressWarnings(weak_only_fit(limit_restart = FALSE))
  expect_null(off$fit$limit)
  expect_identical(fit$fit$objective, off$fit$objective)
  expect_identical(fit$fit$theta, off$fit$theta)
})

test_that("#448: when the best objective is at the step, the fit refuses", {
  # Task 1 grid: best -21.66 at the step, best clean start -208.31. Returning
  # the clean start would be 130 units worse than today's -78.15.
  data(avc, package = "TemporalHazard", envir = environment())
  tt <- avc$int_dead
  st <- avc$dead
  k <- is.finite(tt) & is.finite(st) & tt > 0
  err <- tryCatch(suppressWarnings(hazard(
    time = tt[k], status = st[k], dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
                  late = hzr_phase("g3", tau = 5, gamma = 1, alpha = 1, eta = 1,
                                   constraint = "alpha_gamma_eta")),
    fit = TRUE, control = list(n_starts = 1L))),
    hzr_limit_refusal = function(e) e)
  expect_s3_class(err, "hzr_limit_refusal")
  expect_identical(err$parameter, "nu")
  # The best objective it reports is the step's, above today's returned fit.
  expect_gt(err$best_objective, -78.15)
})

test_that("a clean fit is unchanged and records $limit NULL", {
  data(avc, package = "TemporalHazard", envir = environment())
  k <- is.finite(avc$int_dead) & avc$int_dead > 0
  fit <- hazard(time = avc$int_dead[k], status = avc$dead[k], dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "m"),
                  constant = hzr_phase("constant")),
    fit = TRUE, control = list(n_starts = 1L))
  expect_true("limit" %in% names(fit$fit))
  expect_null(fit$fit$limit)
})

test_that("an unfitted object records $limit NA", {
  fit <- hazard(time = c(1, 2, 3), status = c(1, 0, 1), dist = "multiphase",
    phases = list(constant = hzr_phase("constant")), fit = FALSE)
  expect_true(is.na(fit$fit$limit))
})

test_that("limit_restart must be TRUE or FALSE", {
  expect_error(fit_729(1L, limit_restart = NA), "limit_restart")
  expect_error(fit_729(1L, limit_restart = "yes"), "limit_restart")
})

test_that("the refusal message counts restarts that ended in an error", {
  trig <- list(list(phase = "early", parameter = "t_half",
                    limit = "0 (below the first observed time)",
                    source = "unbounded_phase"))
  with_err <- .hzr_limit_refusal(trig, -100, 9L, n_errored = 3L)
  expect_match(conditionMessage(with_err), "3 ended in an error", fixed = TRUE)
  without <- .hzr_limit_refusal(trig, -100, 9L, n_errored = 0L)
  expect_no_match(conditionMessage(without), "ended in an error")
})

# The layer re-evaluates the user's call. It must accept every call hazard()
# accepted without it, and give the same answer: through `...` wrappers,
# do.call(), local variables, and with weights named as a column of `data`.
wrap_data <- function() {
  set.seed(1)
  data.frame(stop = stats::rexp(60), event = stats::rbinom(60, 1, 0.7),
             w = stats::runif(60, 0.5, 1.5))
}
wrap_const <- function(...) {
  hazard(..., dist = "multiphase",
         phases = list(constant = hzr_phase("constant")), fit = TRUE)
}

test_that("a `...` wrapper with weights named as a column fits as the direct call", {
  d <- wrap_data()
  for (lr in c(TRUE, FALSE)) {
    direct <- hazard(survival::Surv(stop, event) ~ 1, data = d, weights = w,
                     dist = "multiphase",
                     phases = list(constant = hzr_phase("constant")),
                     fit = TRUE, control = list(limit_restart = lr))
    wrapped <- wrap_const(survival::Surv(stop, event) ~ 1, data = d,
                          weights = w, control = list(limit_restart = lr))
    expect_equal(wrapped$fit$objective / direct$fit$objective, 1,
                 tolerance = 1e-10)
    # The weights were used: an unweighted fit differs.
    unweighted <- wrap_const(survival::Surv(stop, event) ~ 1, data = d)
    expect_gt(abs(unweighted$fit$objective - direct$fit$objective), 1e-3)
  }
})

test_that("do.call() and a wrapper's local variables fit as the direct call", {
  d <- wrap_data()
  direct <- hazard(survival::Surv(stop, event) ~ 1, data = d, weights = w,
                   dist = "multiphase",
                   phases = list(constant = hzr_phase("constant")), fit = TRUE)
  same <- function(fit) {
    expect_equal(fit$fit$objective / direct$fit$objective, 1, tolerance = 1e-10)
  }
  args <- list(formula = survival::Surv(stop, event) ~ 1, data = d,
               weights = quote(w))
  same(do.call(hazard, c(args, list(dist = "multiphase",
    phases = list(constant = hzr_phase("constant")), fit = TRUE))))
  same(do.call(wrap_const, args))
  # Data local to the calling function, weights named as its column. (A
  # weights VECTOR local to that function is not found through a wrapper,
  # with or without the layer: hazard() evaluates it in its own caller.)
  local_wrap <- function() {
    dd <- wrap_data()
    wrap_const(survival::Surv(stop, event) ~ 1, data = dd, weights = w)
  }
  same(local_wrap())
  # Data and weights local to the wrapper itself.
  inner_wrap <- function(...) {
    dd <- wrap_data()
    hazard(..., data = dd, weights = w, dist = "multiphase",
           phases = list(constant = hzr_phase("constant")), fit = TRUE)
  }
  same(inner_wrap(survival::Surv(stop, event) ~ 1))
})

test_that("the vector interface fits through a `...` wrapper", {
  d <- wrap_data()
  direct <- hazard(time = d$stop, status = d$event, weights = d$w,
                   dist = "multiphase",
                   phases = list(constant = hzr_phase("constant")), fit = TRUE)
  vec_wrap <- function(dd) {
    wrap_const(time = stop, status = event, weights = w, data = dd)
  }
  expect_equal(vec_wrap(d)$fit$objective / direct$fit$objective, 1,
               tolerance = 1e-10)
})

test_that("a hard trigger through a `...` wrapper refuses as the direct call", {
  data(avc, package = "TemporalHazard", envir = environment())
  k <- is.finite(avc$int_dead) & is.finite(avc$dead) & avc$int_dead > 0
  dd <- data.frame(tt = avc$int_dead[k], st = avc$dead[k])
  ph <- list(early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
             late = hzr_phase("g3", tau = 5, gamma = 1, alpha = 1, eta = 1,
                              constraint = "alpha_gamma_eta"))
  wrap_448 <- function(...) {
    hazard(..., dist = "multiphase", phases = ph, fit = TRUE,
           control = list(n_starts = 1L))
  }
  refuse <- function(expr) {
    tryCatch(suppressWarnings(expr), hzr_limit_refusal = function(e) e)
  }
  direct <- refuse(hazard(time = tt, status = st, data = dd,
                          dist = "multiphase", phases = ph, fit = TRUE,
                          control = list(n_starts = 1L)))
  wrapped <- refuse(wrap_448(time = tt, status = st, data = dd))
  expect_s3_class(direct, "hzr_limit_refusal")
  expect_s3_class(wrapped, "hzr_limit_refusal")
  expect_identical(wrapped$parameter, direct$parameter)
  expect_equal(wrapped$best_objective / direct$best_objective, 1,
               tolerance = 1e-10)
})

test_that("a #729 restart through a `...` wrapper returns the direct call's fit", {
  d <- valves_729()
  d$one <- 1
  ph <- list(early = hzr_phase("cdf", t_half = 10, nu = 0.2, m = 1, fixed = "m",
                               formula = ~ age_cop + nyha + mv_nyha + double_ +
                                 ao_pinc),
             constant = hzr_phase("constant", formula = ~ age_cop + black + i_path))
  ctl <- list(n_starts = 1L, maxit = 2000, conserve = TRUE)
  wrap_729 <- function(...) {
    hazard(..., dist = "multiphase", phases = ph, fit = TRUE, control = ctl)
  }
  # Weights named as a column, so every grid start re-evaluates it too.
  expect_warning(direct <- hazard(survival::Surv(int_dead, dead) ~ 1, data = d,
                                  weights = one, dist = "multiphase",
                                  phases = ph, fit = TRUE, control = ctl),
                 class = "hzr_limit_restarted")
  expect_warning(wrapped <- wrap_729(survival::Surv(int_dead, dead) ~ 1,
                                     data = d, weights = one),
                 class = "hzr_limit_restarted")
  expect_identical(direct$fit$limit$outcome, "restarted")
  expect_identical(wrapped$fit$limit$outcome, "restarted")
  expect_identical(wrapped$fit$limit$winner, direct$fit$limit$winner)
  expect_equal(wrapped$fit$objective / direct$fit$objective, 1,
               tolerance = 1e-10)
  expect_equal(direct$fit$objective / sas_ll_729, 1, tolerance = 1e-5)
})
