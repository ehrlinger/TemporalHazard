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

# Each base in the grid, with the entries to try on it. Built lazily so a
# skipped block builds nothing.
.ws_grid <- function() {
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
  list(
    inherits = list(
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
    ),
    reference_coe = list(
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
    ),
    intercepts = list(
      data = a, ctl = list(n_starts = 1L),
      base = suppressWarnings(hazard(
                 survival::Surv(int_dead, dead) ~ 1, data = a,
                 dist = "multiphase", phases = phase_fixed(),
                 control = list(n_starts = 1L),
                 fit = TRUE)),
      adds = list(c("opmos", "constant"), c("age", "early"),
                  c("mal", "constant"))
    ),
    weighted = list(
      data = a, ctl = list(n_starts = 1L),
      base = suppressWarnings(hazard(
                 survival::Surv(int_dead, dead) ~ 1, data = a,
                 weights = w, dist = "multiphase", phases = phase_fixed(),
                 control = list(n_starts = 1L),
                 fit = TRUE)),
      adds = list(c("age", "early"), c("opmos", "constant"))
    ),
    windows = list(
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
  )
}

test_that("a multiphase refit starts from the base's estimates (#551)", {
  g <- .ws_grid()$inherits
  ns <- asNamespace("TemporalHazard")
  orig <- get(".hzr_optim_multiphase", envir = ns)
  seen <- NULL
  local_mocked_bindings(.hzr_optim_multiphase = function(...) {
    seen <<- list(...)$theta_start
    orig(...)
  })
  cand <- suppressWarnings(.hzr_refit_with_scope(
    g$base, action = "add", var = "mal", phase = "constant",
    data = g$data, control = g$ctl
  ))
  expected <- .ws_base_plus_zero(g$base, cand)
  # Known positive: the candidate model at that start reproduces the base.
  ev <- suppressWarnings(hzr_evaluate(cand, theta = expected))
  expect_equal(ev$logLik, g$base$fit$objective, tolerance = 1e-10)
  # The optimizer received exactly that start, not NULL.
  expect_false(is.null(seen))
  expect_equal(unname(seen), unname(expected), tolerance = 0)
  expect_equal(names(seen), names(expected))

  # A drop starts from the base's estimates with the dropped slot removed.
  seen <- NULL
  g2 <- .ws_grid()$reference_coe
  cand <- suppressWarnings(.hzr_refit_with_scope(
    g2$base, action = "drop", var = "op_age", phase = "constant",
    data = g2$data, control = g2$ctl
  ))
  expected <- coef(g2$base)[names(coef(cand))]
  expect_false(anyNA(expected))
  expect_false("constant.op_age" %in% names(expected))
  expect_false(is.null(seen))
  expect_equal(unname(seen), unname(expected), tolerance = 0)
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
  grid <- .ws_grid()
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
  g <- .ws_grid()$reference_coe
  cand <- suppressWarnings(.hzr_refit_with_scope(
    g$base, action = "add", var = "age", phase = "early",
    data = g$data, control = list(n_starts = 3, conserve = TRUE)
  ))
  fixed <- c("early.log_t_half", "early.nu", "early.m")
  expect_equal(coef(cand)[fixed], coef(g$base)[fixed], tolerance = 0)
  # Perturbed starts cannot make it worse than the unperturbed first start.
  ll0 <- g$base$fit$objective
  expect_gte(cand$fit$objective, ll0 - .ws_tol(ll0))
  expect_equal(nrow(cand$fit$starts), 3L)
})
