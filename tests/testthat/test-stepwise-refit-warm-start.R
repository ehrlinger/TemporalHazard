# A multiphase stepwise refit started from the phase specs' default values,
# not from the model it extends (#551). The candidate model contains the base
# (the base's estimates with the new coefficient at 0 reproduce the base
# log-likelihood exactly), so a refit started there cannot end below the
# base. From the cold start it did, by up to 25 log-likelihood units, and
# reported converged = TRUE: every multiphase entry was tested against a
# model that was not at its optimum.

.ws_avc <- local({
  data(avc, package = "TemporalHazard")
  avc
})

# The base's estimates with the named coefficient(s) of `cand` absent from
# the base set to 0: the candidate model's copy of the base.
.ws_base_plus_zero <- function(base, cand) {
  cb <- coef(base)
  th <- coef(cand)
  th[] <- 0
  common <- intersect(names(cb), names(th))
  th[common] <- cb[common]
  th
}

.ws_tol <- function(ll) 1e-6 * max(1, abs(ll))

# Each base in the grid, with the entries to try on it. Only the bases a
# block names are built, so a block builds only what it uses.
.ws_grid <- function(which = NULL) {
  a <- .ws_avc
  a284 <- stats::na.omit(a[, c("int_dead", "dead", "age", "mal")])
  a284$grp <- cut(a284$age, 3)
  d_ref <- a
  d_ref$inc_surg[is.na(d_ref$inc_surg)] <- mean(d_ref$inc_surg, na.rm = TRUE)
  set.seed(1)
  w <- stats::runif(nrow(a), 0.5, 2)
  phase_fixed <- function() {
    list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes"),
         constant = hzr_phase("constant"))
  }
  builders <- list(
    inherits = function() {
      list(
      data = a284, ctl = list(n_starts = 1L, conserve = FALSE),
      base = suppressWarnings(hazard(
                 survival::Surv(int_dead, dead) ~ grp, data = a284,
                 dist = "multiphase",
                 phases = list(
                   early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1,
                                     fixed = "m"),
                   constant = hzr_phase("constant", formula = ~ age)),
                 control = list(n_starts = 1L, conserve = FALSE),
                 fit = TRUE)),
      adds = list(c("mal", "constant"))
      )
    },
    reference_coe = function() {
      list(
      data = d_ref, ctl = list(n_starts = 1, conserve = TRUE),
      base = suppressWarnings(hazard(
                 survival::Surv(int_dead, dead) ~ 1, data = d_ref,
                 dist = "multiphase",
                 phases = list(
                   early = hzr_phase("cdf", t_half = 0.1511909,
                                     nu = 1.438631, m = 1, fixed = "shapes",
                                     formula = ~ com_iv + mal),
                   constant = hzr_phase("constant", formula = ~ orifice +
                                          mal + op_age + inc_surg)),
                 control = list(n_starts = 1, conserve = TRUE),
                 fit = TRUE)),
      adds = list(c("age", "early"), c("opmos", "early"),
                  c("age", "constant"), c("status", "constant"))
      )
    },
    intercepts = function() {
      list(
      data = a, ctl = list(n_starts = 1L),
      base = suppressWarnings(hazard(
                 survival::Surv(int_dead, dead) ~ 1, data = a,
                 dist = "multiphase", phases = phase_fixed(),
                 control = list(n_starts = 1L),
                 fit = TRUE)),
      adds = list(c("opmos", "constant"), c("age", "early"),
                  c("mal", "constant"))
      )
    },
    weighted = function() {
      list(
      data = a, ctl = list(n_starts = 1L),
      base = suppressWarnings(hazard(
                 survival::Surv(int_dead, dead) ~ 1, data = a,
                 weights = w, dist = "multiphase", phases = phase_fixed(),
                 control = list(n_starts = 1L),
                 fit = TRUE)),
      adds = list(c("age", "early"), c("opmos", "constant"))
      )
    },
    windows = function() {
      list(
      data = a284, ctl = list(n_starts = 1L, conserve = FALSE),
      base = suppressWarnings(hazard(
                 survival::Surv(int_dead, dead) ~ age, data = a284,
                 dist = "multiphase", time_windows = 1,
                 phases = list(
                   early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1,
                                     fixed = "m"),
                   constant = hzr_phase("constant", formula = ~ age)),
                 control = list(n_starts = 1L, conserve = FALSE),
                 fit = TRUE)),
      adds = list(c("mal", "constant"))
      )
    }
  )
  if (is.null(which)) which <- names(builders)
  lapply(builders[which], function(build) build())
}

# Every start the multiphase optimizer is handed during `expr`, in order.
.ws_starts_handed <- function(expr) {
  orig <- .hzr_optim_multiphase
  seen <- list()
  local_mocked_bindings(.hzr_optim_multiphase = function(...) {
    # Wrapped, so a NULL start (the default) is recorded, not dropped.
    seen[[length(seen) + 1L]] <<- list(start = list(...)$theta_start)
    orig(...)
  })
  value <- expr
  list(value = value, seen = lapply(seen, `[[`, "start"))
}

test_that("a multiphase refit is fitted from two starts, warm and default", {
  g <- .ws_grid("inherits")$inherits
  out <- .ws_starts_handed(suppressWarnings(.hzr_refit_with_scope(
    g$base, action = "add", var = "mal", phase = "constant",
    data = g$data, control = g$ctl
  )))
  cand <- out$value
  expected <- .ws_base_plus_zero(g$base, cand)
  # Known positive: the candidate model at that start reproduces the base.
  ev <- suppressWarnings(hzr_evaluate(cand, theta = expected))
  expect_equal(ev$logLik, g$base$fit$objective, tolerance = 1e-10)
  # Two fits: the first from exactly that start, the second from the phase
  # specs' default values -- log(0.1) for each rate, 0 for each coefficient,
  # the spec's free shapes -- with the fixed shape at the base's value.
  expect_length(out$seen, 2L)
  expect_equal(unname(out$seen[[1]]), unname(expected), tolerance = 0)
  expect_equal(names(out$seen[[1]]), names(expected))
  dflt <- out$seen[[2]]
  expect_equal(names(dflt), names(expected))
  expect_equal(unname(dflt[c("early.log_mu", "constant.log_mu")]),
               rep(log(0.1), 2), tolerance = 1e-12)
  covs <- c("early.grp(264,528]", "early.grp(528,792]", "constant.age",
            "constant.mal")
  expect_equal(unname(dflt[covs]), rep(0, 4), tolerance = 0)
  # The free shapes keep the spec's values, not the base's estimates, or the
  # second start is nearly the first and no longer a second optimum search.
  expect_equal(unname(dflt[c("early.log_t_half", "early.nu")]),
               c(log(0.15), 1.4), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(
    unname(dflt[c("early.log_t_half", "early.nu")]),
    unname(coef(g$base)[c("early.log_t_half", "early.nu")])
  )))
  expect_equal(dflt[["early.m"]], coef(g$base)[["early.m"]], tolerance = 0)
  # The kept fit is the better of the two, and says which it was.
  obj <- cand$fit$refit_objectives
  expect_named(obj, c("warm", "default"))
  expect_true(all(is.finite(obj)))
  expect_identical(cand$fit$refit_start, names(which.max(obj))[1])
  expect_equal(cand$fit$objective, max(obj), tolerance = 0)

  # A drop starts from the base's estimates with the dropped slot removed.
  g2 <- .ws_grid("reference_coe")$reference_coe
  out <- .ws_starts_handed(suppressWarnings(.hzr_refit_with_scope(
    g2$base, action = "drop", var = "op_age", phase = "constant",
    data = g2$data, control = g2$ctl
  )))
  expected <- coef(g2$base)[names(coef(out$value))]
  expect_false(anyNA(expected))
  expect_false("constant.op_age" %in% names(expected))
  expect_equal(unname(out$seen[[1]]), unname(expected), tolerance = 0)
})

test_that("the better start is kept, and a failed start never wins", {
  fit_like <- function(obj, converged = TRUE) {
    structure(list(fit = list(objective = obj, converged = converged)),
              class = "hazard")
  }
  pick <- function(w, d) {
    .hzr_refit_best_start(list(value = w, warnings = list()),
                          list(value = d, warnings = list()))
  }
  expect_identical(pick(fit_like(-10), fit_like(-12))$fit$refit_start,
                   "warm")
  expect_identical(pick(fit_like(-12), fit_like(-10))$fit$refit_start,
                   "default")
  # A tie keeps the warm start, which cannot end below the base.
  expect_identical(pick(fit_like(-10), fit_like(-10))$fit$refit_start,
                   "warm")
  # Not converged, non-finite, or an error: the other start wins.
  expect_identical(pick(fit_like(-1, FALSE), fit_like(-10))$fit$refit_start,
                   "default")
  expect_identical(pick(fit_like(NaN), fit_like(-10))$fit$refit_start,
                   "default")
  expect_identical(
    pick(simpleError("warm failed"), fit_like(-10))$fit$refit_start,
    "default")
  expect_identical(
    pick(fit_like(-10), simpleError("default failed"))$fit$refit_start,
    "warm")
  # Both unusable: the default start's outcome, as before #551.
  expect_error(pick(simpleError("w"), simpleError("default failed")),
               "default failed")
  both_bad <- pick(fit_like(-1, FALSE), fit_like(-2, FALSE))
  expect_false(both_bad$fit$converged)
  expect_equal(both_bad$fit$objective, -2)
  # Only the kept fit's warnings reach the caller.
  w <- testthat::capture_warnings(.hzr_refit_best_start(
    list(value = fit_like(-10), warnings = list(simpleWarning("from warm"))),
    list(value = fit_like(-12),
         warnings = list(simpleWarning("from default")))
  ))
  expect_identical(w, "from warm")
})

test_that("no multiphase entry refit ends below its base (#551)", {
  skip_on_cran()
  grid <- .ws_grid()
  n_checked <- 0L
  for (nm in names(grid)) {
    g <- grid[[nm]]
    ll0 <- g$base$fit$objective
    expect_true(is.finite(ll0), label = paste(nm, "base logLik"))
    for (a in g$adds) {
      cand <- suppressWarnings(.hzr_refit_with_scope(
        g$base, action = "add", var = a[[1]], phase = a[[2]],
        data = g$data, control = g$ctl
      ))
      lab <- paste0(nm, ": ", a[[1]], "@", a[[2]])
      expect_true(isTRUE(cand$fit$converged), label = lab)
      expect_gte(cand$fit$objective, ll0 - .ws_tol(ll0), label = lab)
      # The warm start alone already meets the property; the kept fit is
      # the better of the two.
      obj <- cand$fit$refit_objectives
      expect_gte(obj[["warm"]], ll0 - .ws_tol(ll0), label = lab)
      expect_equal(cand$fit$objective, max(obj, na.rm = TRUE),
                   tolerance = 0, label = lab)
      n_checked <- n_checked + 1L
    }
  }
  # The grid ran every entry it lists.
  expect_equal(n_checked, sum(lengths(lapply(grid, `[[`, "adds"))))
})

test_that("a drop refit cannot end below its own start (#551)", {
  skip_on_cran()
  # Without Conservation of Events, so the objective at the start is the
  # plain log-likelihood hzr_evaluate() returns.
  grid <- .ws_grid(c("inherits", "windows"))
  drops <- list(list(g = grid$inherits, v = "age", phase = "constant"),
                list(g = grid$windows, v = "age", phase = "constant"))
  for (d in drops) {
    cand <- suppressWarnings(.hzr_refit_with_scope(
      d$g$base, action = "drop", var = d$v, phase = d$phase,
      data = d$g$data, control = d$g$ctl
    ))
    # The start is the base's estimates with the dropped slot removed.
    start <- coef(d$g$base)[names(coef(cand))]
    expect_false(anyNA(start))
    ll_start <- suppressWarnings(hzr_evaluate(cand, theta = start))$logLik
    expect_gte(cand$fit$objective, ll_start - .ws_tol(ll_start),
               label = paste("drop", d$v))
  }
})

test_that("fixed shapes stay fixed and start 1 is the warm start (#551)", {
  skip_on_cran()
  g <- .ws_grid("reference_coe")$reference_coe
  cand <- suppressWarnings(.hzr_refit_with_scope(
    g$base, action = "add", var = "age", phase = "early",
    data = g$data, control = list(n_starts = 3, conserve = TRUE)
  ))
  fixed <- c("early.log_t_half", "early.nu", "early.m")
  expect_equal(coef(cand)[fixed], coef(g$base)[fixed], tolerance = 0)
  ll0 <- g$base$fit$objective
  expect_gte(cand$fit$objective, ll0 - .ws_tol(ll0))

  # Start 1 of several IS the warm start, unperturbed: the same model fitted
  # from base+0 with one start reaches exactly what start 1 of three reaches.
  warm <- .ws_base_plus_zero(g$base, cand)
  fit_from <- function(n) {
    suppressWarnings(hazard(
      survival::Surv(int_dead, dead) ~ 1, data = g$data,
      dist = "multiphase", phases = cand$spec$phases, theta = warm,
      control = list(n_starts = n, conserve = TRUE), fit = TRUE
    ))
  }
  one <- fit_from(1)
  three <- fit_from(3)
  expect_equal(nrow(three$fit$starts), 3L)
  expect_equal(three$fit$starts$objective[1], one$fit$objective,
               tolerance = 1e-12)
  expect_gte(three$fit$starts$objective[1], ll0 - .ws_tol(ll0))
})

test_that("both starts hold a fixed shape at the base's value (#551)", {
  skip_on_cran()
  # A fixed shape is held at its start value. This base was fitted with a
  # user theta whose t_half (0.5) differs from the spec's (0.15), so a start
  # taking the spec's value fits another model -- one that won by about 15
  # log-likelihood units on the moved shape, credited to the candidate.
  a <- stats::na.omit(.ws_avc)
  base <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = a, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.15, nu = 1, m = 1,
                        fixed = "shapes"),
      constant = hzr_phase("constant")),
    theta = c(log(0.1), log(0.5), 1, 1, log(0.1)),
    control = list(n_starts = 1L), fit = TRUE
  ))
  fixed <- c("early.log_t_half", "early.nu", "early.m")
  # The premise: the base holds t_half at 0.5, not the spec's 0.15.
  expect_equal(coef(base)[["early.log_t_half"]], log(0.5), tolerance = 1e-12)
  out <- .ws_starts_handed(suppressWarnings(.hzr_refit_with_scope(
    base, action = "add", var = "age", phase = "constant", data = a,
    control = list(n_starts = 1L)
  )))
  expect_length(out$seen, 2L)
  for (s in out$seen) {
    expect_equal(unname(s[fixed]), unname(coef(base)[fixed]), tolerance = 0)
  }
  cand <- out$value
  expect_equal(unname(coef(cand)[fixed]), unname(coef(base)[fixed]),
               tolerance = 0)
  # age's own gain, whichever start won, is small here (0.016 on this data);
  # the moved shape was worth about 15.
  expect_lt(cand$fit$objective - base$fit$objective, 1)
})

test_that("a global covariate pool with phase formulas refits (#551)", {
  skip_on_cran()
  # The global formula lists the candidates while every phase has its own
  # formula, so theta has fewer entries than the global design has columns.
  a <- stats::na.omit(.ws_avc)
  base <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ age + mal + com_iv + opmos + orifice +
      op_age + status + inc_surg,
    data = a, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                        fixed = "shapes", formula = ~ 1),
      constant = hzr_phase("constant", formula = ~ 1)),
    control = list(n_starts = 1L), fit = TRUE
  ))
  expect_lt(length(base$fit$theta), ncol(base$data$x))
  cand <- suppressWarnings(.hzr_refit_with_scope(
    base, action = "add", var = "age", phase = "constant", data = a,
    control = list(n_starts = 1L)
  ))
  expect_true("constant.age" %in% names(coef(cand)))
  expect_gte(cand$fit$objective,
             base$fit$objective - .ws_tol(base$fit$objective))
  # The warm fit itself ran. A failed warm start is masked by the default
  # one, so the result alone cannot show that it was refused.
  expect_true(is.finite(cand$fit$refit_objectives[["warm"]]))
})

test_that("a theta forwarded to a refit is refused by name", {
  g <- .ws_grid("inherits")$inherits
  expect_error(
    .hzr_refit_with_scope(g$base, action = "add", var = "mal",
                          phase = "constant", data = g$data,
                          theta = c(1, 2)),
    "`theta` cannot be passed to a stepwise refit"
  )
})
