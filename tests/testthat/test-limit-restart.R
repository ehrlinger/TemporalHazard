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

# Where n_starts = 5 drifted on macOS: early.t_half at e^275.8, a rise past
# the data. Starting there puts the fit on that limit on any platform; where
# a blind start drifts is not (Linux R-devel did not reach it).
drift_729 <- c(113.563, 275.830, -2.301, 1, 0.021, 0.557, 0.205, 0.447, 0.624,
               7.752, -1.481, 1.003, 0.409)

test_that("limit_restart = FALSE refuses at once, naming the parameter", {
  d <- valves_729()
  err <- tryCatch(suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(early = hzr_phase("cdf", t_half = 10, nu = 0.2, m = 1, fixed = "m",
                    formula = ~ age_cop + nyha + mv_nyha + double_ + ao_pinc),
                  constant = hzr_phase("constant", formula = ~ age_cop + black + i_path)),
    theta = drift_729, fit = TRUE,
    control = list(n_starts = 1L, maxit = 2000, conserve = TRUE,
                   limit_restart = FALSE))),
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
  # Not refitted, but not clean either: the trigger is on the record and the
  # line says the restart was not tried.
  lim <- fit$fit$limit
  expect_identical(lim$outcome, "skipped")
  expect_true(all(vapply(lim$trigger, `[[`, "", "source") == "weak"))
  expect_null(lim$grid)
  expect_match(paste(utils::capture.output(print(fit)), collapse = "\n"),
               "Limit check: restart not tried (", fixed = TRUE)
})

weak_only_fit <- function(limit_restart = TRUE, fit = TRUE) {
  set.seed(42)
  n <- 100
  time   <- rexp(n, rate = 0.5) + 0.01
  status <- sample(0:1, n, replace = TRUE, prob = c(0.2, 0.8))
  phases <- list(early = hzr_phase("cdf", t_half = 1, nu = 1.5, m = 0),
                 const = hzr_phase("constant"))
  hazard(time = time, status = status, dist = "multiphase",
         phases = phases, fit = fit,
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
  expect_identical(off$fit$limit$outcome, "skipped")
  expect_identical(off$fit$limit$trigger, lim$trigger)
  expect_identical(fit$fit$objective, off$fit$objective)
  expect_identical(fit$fit$theta, off$fit$theta)
})

test_that("fit = 1 is checked as fit = TRUE is", {
  # The body fits on `if (fit && ...)`, so fit = 1 fits; it must not skip the
  # check and print "Limit check: none" over an unexamined fit.
  f1 <- suppressWarnings(weak_only_fit(fit = 1))
  ft <- suppressWarnings(weak_only_fit())
  expect_identical(f1$fit$limit$outcome, "kept")
  expect_identical(f1$fit$limit$trigger, ft$fit$limit$trigger)
  expect_identical(f1$fit$theta, ft$fit$theta)
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

# #683 review: the layer used to re-evaluate the user's call for every inner
# fit, so an argument written in the call ran once per grid start (ten times
# here) and a resampling `data =` gave each start a different dataset. Every
# argument is now evaluated once, and every inner fit sees those values. The
# #448 fixture is a hard trigger that runs the whole grid and refuses.
avc_448 <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  k <- is.finite(avc$int_dead) & is.finite(avc$dead) & avc$int_dead > 0
  data.frame(tt = avc$int_dead[k], st = avc$dead[k], one = 1)
}
# alpha = 1 is replaced by the value the constraint derives, with a warning.
ph_448 <- suppressWarnings(list(
  early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m"),
  late = hzr_phase("g3", tau = 5, gamma = 1, alpha = 1, eta = 1,
                   constraint = "alpha_gamma_eta")))
refuse_448 <- function(expr) {
  tryCatch(suppressWarnings(expr), hzr_limit_refusal = function(e) e)
}
# The grid ran: a refusal that tried no restart would make "once" trivial.
expect_grid_ran <- function(err) {
  expect_s3_class(err, "hzr_limit_refusal")
  expect_match(conditionMessage(err), "None of [1-9][0-9]* restarts")
}

test_that("an argument written in the call is evaluated once (#683 review)", {
  dd <- avc_448()
  n <- 0L
  err <- refuse_448(hazard(time = tt, status = st,
                           data = {
                             n <- n + 1L
                             dd
                           },
                           dist = "multiphase", phases = ph_448, fit = TRUE,
                           control = list(n_starts = 1L)))
  expect_grid_ran(err)
  expect_identical(n, 1L)
  # A column-name argument holding user code runs once too, in the data mask.
  m <- 0L
  err <- refuse_448(hazard(time = tt, status = st, data = dd,
                           weights = {
                             m <<- m + 1L
                             one
                           },
                           dist = "multiphase", phases = ph_448, fit = TRUE,
                           control = list(n_starts = 1L)))
  expect_grid_ran(err)
  expect_identical(m, 1L)
})

test_that("a resampling data argument draws one resample (#683 review)", {
  dd <- avc_448()
  set.seed(1)
  one_draw <- dd[sample(nrow(dd), replace = TRUE), ]
  seed_after_one <- get(".Random.seed", envir = globalenv())
  ref <- refuse_448(hazard(time = tt, status = st, data = one_draw,
                           dist = "multiphase", phases = ph_448, fit = TRUE,
                           control = list(n_starts = 1L)))
  set.seed(1)
  got <- refuse_448(hazard(time = tt, status = st,
                           data = dd[sample(nrow(dd), replace = TRUE), ],
                           dist = "multiphase", phases = ph_448, fit = TRUE,
                           control = list(n_starts = 1L)))
  # The layer drew nothing: the stream is where one sample() leaves it.
  # (expect_true: testthat's diff of two seeds overflows an integer.)
  expect_true(identical(get(".Random.seed", envir = globalenv()),
                        seed_after_one))
  # Every inner fit saw that one resample: the outcome is the fixed data's.
  expect_identical(class(got), class(ref))
  if (inherits(ref, "hzr_limit_refusal")) {
    expect_grid_ran(got)
    expect_identical(got$best_objective, ref$best_objective)
  } else {
    expect_identical(got$fit$objective, ref$fit$objective)
    expect_identical(got$data$time, ref$data$time)
  }
})

test_that("the vector interface is evaluated once on a triggering fit", {
  dd <- avc_448()
  n <- 0L
  err <- refuse_448(hazard(time = {
                             n <- n + 1L
                             dd$tt
                           },
                           status = dd$st, weights = dd$one,
                           dist = "multiphase", phases = ph_448, fit = TRUE,
                           control = list(n_starts = 1L)))
  expect_grid_ran(err)
  expect_identical(n, 1L)
  direct <- refuse_448(hazard(time = dd$tt, status = dd$st,
                              dist = "multiphase", phases = ph_448, fit = TRUE,
                              control = list(n_starts = 1L)))
  expect_identical(err$best_objective, direct$best_objective)
})

# #683 review: user code inside the formula ran once per inner fit, because
# each fit parsed the formula into its own response and design: cnt() below
# ran ten times on the #448 grid against once without the layer. The parse is
# now done once per hazard() call, and every inner fit reads it.
fit_cnt_448 <- function(dd, limit_restart = TRUE) {
  n <- 0L
  cnt <- function(v) {
    n <<- n + 1L
    v
  }
  err <- refuse_448(hazard(survival::Surv(cnt(tt), st) ~ 1, data = dd,
                           dist = "multiphase", phases = ph_448, fit = TRUE,
                           control = list(n_starts = 1L,
                                          limit_restart = limit_restart)))
  list(err = err, n = n)
}

# #737 review: the layer's own inner fit is not checked again, but a
# multiphase hazard() that fit runs (here inside a formula term) is the
# user's, and is checked. With limit_restart = FALSE the nested #448 fit
# refuses at once; returning a fit means it skipped the check.
test_that("a multiphase hazard() run inside a fit is checked (#737 review)", {
  dd <- avc_448()
  seen <- character(0)
  nest <- function(v) {
    got <- refuse_448(hazard(time = dd$tt, status = dd$st,
                             dist = "multiphase", phases = ph_448, fit = TRUE,
                             control = list(n_starts = 1L,
                                            limit_restart = FALSE)))
    seen <<- c(seen, if (inherits(got, "hzr_limit_refusal")) "refused" else "returned")
    v
  }
  outer <- hazard(survival::Surv(nest(tt), st) ~ 1, data = dd,
                  dist = "multiphase",
                  phases = list(constant = hzr_phase("constant")), fit = TRUE)
  expect_null(outer$fit$limit)
  expect_gt(length(seen), 0L)
  expect_identical(unique(seen), "refused")
})

test_that("a refusal carries the restart grid it tried (#737 review)", {
  dd <- avc_448()
  err <- refuse_448(hazard(time = dd$tt, status = dd$st, dist = "multiphase",
                           phases = ph_448, fit = TRUE,
                           control = list(n_starts = 1L)))
  expect_grid_ran(err)
  expect_s3_class(err$grid, "data.frame")
  expect_identical(names(err$grid), c("start", "objective", "triggered", "error"))
  n <- as.integer(sub(".*None of ([0-9]+) restarts.*", "\\1", conditionMessage(err)))
  expect_identical(nrow(err$grid), n)
  off <- refuse_448(hazard(time = dd$tt, status = dd$st, dist = "multiphase",
                           phases = ph_448, fit = TRUE,
                           control = list(n_starts = 1L, limit_restart = FALSE)))
  expect_null(off$grid)
})

test_that("user code in the formula runs as often as without the layer", {
  dd <- avc_448()
  got <- fit_cnt_448(dd)
  expect_grid_ran(got$err)
  off <- fit_cnt_448(dd, limit_restart = FALSE)
  expect_s3_class(off$err, "hzr_limit_refusal")
  expect_identical(off$n, 1L)
  expect_identical(got$n, off$n)
  # Rows at time 0 are dropped and the formula parsed again on the kept rows:
  # that second parse is done once per call too.
  d0 <- rbind(dd[1:3, ], dd)
  d0$tt[1:3] <- 0
  got <- fit_cnt_448(d0)
  expect_grid_ran(got$err)
  off <- fit_cnt_448(d0, limit_restart = FALSE)
  expect_identical(off$n, 2L)
  expect_identical(got$n, off$n)
})

test_that("a random term in the formula draws once (#683 review)", {
  dd <- avc_448()
  # The draw is multiplied by 0: a jitter of even 1e-12 moves this fixture
  # off its trigger, and without one the grid does not run. So the response
  # is `tt`, and only the stream records how many draws were made.
  set.seed(3)
  invisible(stats::runif(nrow(dd)))
  seed_after_one <- get(".Random.seed", envir = globalenv())
  ref <- refuse_448(hazard(survival::Surv(tt, st) ~ 1, data = dd,
                           dist = "multiphase", phases = ph_448, fit = TRUE,
                           control = list(n_starts = 1L)))
  set.seed(3)
  got <- refuse_448(hazard(
    survival::Surv(tt + 0 * stats::runif(length(tt)), st) ~ 1,
    data = dd, dist = "multiphase", phases = ph_448, fit = TRUE,
    control = list(n_starts = 1L)))
  expect_grid_ran(got)
  # The stream is where one runif() leaves it, as without the layer.
  expect_true(identical(get(".Random.seed", envir = globalenv()),
                        seed_after_one))
  expect_grid_ran(ref)
  expect_identical(got$best_objective, ref$best_objective)
})

test_that("#734: an unusable Hessian in the m < 0 joint limit is checked, and kept", {
  # The #551 generator, seed 15: the early phase runs to t_half and nu near
  # infinity with m < 0 (log t_half 512.9, nu 317.8, m -0.312), a power law.
  # The Hessian is not positive definite there, so $weak is NA and, before
  # #734, nothing triggered. Every grid start drifts to the same limit, so
  # no clean start can win and the original is kept (measured 2026-11-11:
  # 9 starts, all triggered, objectives -399.586 to -399.638).
  set.seed(15)
  n <- 300
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  early <- stats::runif(n) < 0.3
  tt <- ifelse(early, stats::rexp(n, 8), stats::rexp(n, 0.2 * exp(1.5 * x1)))
  cens <- stats::runif(n, 2, 15)
  d <- data.frame(tt = pmin(tt, cens) + 0.001, ev = as.numeric(tt <= cens))
  fit_551 <- function(limit_restart = TRUE) {
    hazard(time = d$tt, status = d$ev, dist = "multiphase",
           phases = list(early = hzr_phase("cdf", t_half = 1, nu = 1.5, m = 0),
                         const = hzr_phase("constant")),
           fit = TRUE, control = list(limit_restart = limit_restart))
  }
  fit <- suppressWarnings(fit_551())
  expect_identical(fit$fit$weak, NA)
  expect_gt(fit$fit$theta[["early.log_t_half"]], 100)
  expect_lt(fit$fit$theta[["early.m"]], 0)
  lim <- fit$fit$limit
  expect_type(lim, "list")
  expect_identical(lim$outcome, "kept")
  expect_identical(lim$trigger, list(list(
    phase = "early", parameter = "t_half",
    limit = "no usable Hessian (identification not checked)", source = "weak_na")))
  expect_identical(nrow(lim$grid), 9L)
  expect_true(all(lim$grid$triggered))
  expect_true(all(is.finite(lim$grid$objective)))
  expect_match(paste(utils::capture.output(print(fit)), collapse = "\n"), paste0(
    "Limit check: restart tried (early.t_half: no usable Hessian ",
    "(identification not checked)); none better, original kept"), fixed = TRUE)
  # Weak-only: with the layer off, the same fit, its trigger recorded and
  # the restart marked not tried.
  off <- suppressWarnings(fit_551(limit_restart = FALSE))
  expect_identical(off$fit$limit$outcome, "skipped")
  expect_identical(off$fit$limit$trigger, fit$fit$limit$trigger)
  expect_identical(off$fit$theta, fit$fit$theta)
})

# #683 review: user code in a PHASE formula ran once per design built: in
# each inner fit, in the layer's own setup, and again when a grid start
# counted the phases' columns, 20 times on this fixture against once with
# limit_restart = FALSE. The #448 data with the early phase on `mal` is a
# weak-only trigger whose grid runs and keeps the original.
avc_mal <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  k <- is.finite(avc$int_dead) & is.finite(avc$dead) & avc$int_dead > 0
  data.frame(tt = avc$int_dead[k], st = avc$dead[k], mal = avc$mal[k])
}
fit_mal <- function(dd, early_formula, limit_restart = TRUE) {
  ph <- suppressWarnings(list(
    early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m",
                      formula = early_formula),
    late = hzr_phase("g3", tau = 5, gamma = 1, alpha = 1, eta = 1,
                     constraint = "alpha_gamma_eta")))
  suppressWarnings(hazard(survival::Surv(tt, st) ~ 1, data = dd,
                          dist = "multiphase", phases = ph, fit = TRUE,
                          control = list(n_starts = 1L,
                                         limit_restart = limit_restart)))
}
expect_kept_grid <- function(fit) {
  expect_identical(fit$fit$limit$outcome, "kept")
  expect_identical(nrow(fit$fit$limit$grid), 9L)
}

test_that("user code in a phase formula runs as often as without the layer", {
  dd <- avc_mal()
  n <- 0L
  cnt <- function(v) {
    n <<- n + 1L
    v
  }
  got <- fit_mal(dd, ~ cnt(mal))
  expect_kept_grid(got)
  n_on <- n
  n <- 0L
  off <- fit_mal(dd, ~ cnt(mal), limit_restart = FALSE)
  expect_identical(off$fit$limit$outcome, "skipped")
  expect_identical(n, 1L)
  expect_identical(n_on, n)
  # The design every fit used is the one the user's code returned.
  expect_identical(got$fit$theta, off$fit$theta)
})

test_that("a phase formula's warning reaches the user as often as without the layer", {
  # The layer's setup rebuilt the first fit's phases and replayed the cached
  # warning beside the first fit's own: twice, against once (#683 review).
  dd <- avc_mal()
  wf <- function(v) {
    warning("planted phase-formula warning")
    v
  }
  count <- function(limit_restart) {
    k <- 0L
    ph <- suppressWarnings(list(
      early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m",
                        formula = ~ wf(mal)),
      late = hzr_phase("g3", tau = 5, gamma = 1, alpha = 1, eta = 1,
                       constraint = "alpha_gamma_eta")))
    fit <- withCallingHandlers(
      hazard(survival::Surv(tt, st) ~ 1, data = dd, dist = "multiphase",
             phases = ph, fit = TRUE,
             control = list(n_starts = 1L, limit_restart = limit_restart)),
      warning = function(w) {
        if (grepl("planted phase-formula warning", conditionMessage(w),
                  fixed = TRUE)) k <<- k + 1L
        invokeRestart("muffleWarning")
      })
    list(fit = fit, k = k)
  }
  on <- count(TRUE)
  expect_kept_grid(on$fit)
  off <- count(FALSE)
  expect_identical(off$k, 1L)
  expect_identical(on$k, off$k)
})

test_that("the phase-design cache rebuilds for another formula or other rows", {
  # The cache is keyed on the phase name, and an entry is reused only for an
  # identical formula and data; anything else is built fresh, as without it.
  dd <- avc_mal()
  n <- 0L
  cnt <- function(v) {
    n <<- n + 1L
    v
  }
  old <- .hzr_limit_state$cache
  .hzr_limit_state$cache <- new.env(parent = emptyenv())
  on.exit(.hzr_limit_state$cache <- old, add = TRUE)
  f1 <- ~ cnt(mal)
  a <- .hzr_phase_formula_design("early", f1, dd)
  b <- .hzr_phase_formula_design("early", f1, dd)
  expect_identical(n, 1L)
  expect_identical(b, a)
  f2 <- ~ I(2 * cnt(mal))
  d2 <- .hzr_phase_formula_design("early", f2, dd)
  expect_identical(n, 2L)
  expect_equal(unname(d2$x[, 1L]), 2 * unname(a$x[, 1L]))
  short <- .hzr_phase_formula_design("early", f1, dd[1:50, ])
  expect_identical(n, 3L)
  expect_identical(nrow(short$x), 50L)
})

test_that("a random term in a phase formula draws once (#683 review)", {
  dd <- avc_mal()
  # Multiplied by 0, so the design is `mal` and the fit stays on its trigger;
  # only the stream records how many draws were made.
  set.seed(5)
  invisible(stats::runif(nrow(dd)))
  seed_after_one <- get(".Random.seed", envir = globalenv())
  set.seed(5)
  got <- fit_mal(dd, ~ I(mal + 0 * stats::runif(length(mal))))
  expect_kept_grid(got)
  expect_true(identical(get(".Random.seed", envir = globalenv()),
                        seed_after_one))
  set.seed(5)
  off <- fit_mal(dd, ~ I(mal + 0 * stats::runif(length(mal))),
                 limit_restart = FALSE)
  expect_true(identical(get(".Random.seed", envir = globalenv()),
                        seed_after_one))
  expect_identical(got$fit$theta, off$fit$theta)
})
