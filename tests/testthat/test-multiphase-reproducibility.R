# test-multiphase-reproducibility.R
# The multiphase multi-start must be reproducible and must not disturb the
# caller's RNG stream.  Starts after the first perturb the initial values; when
# those perturbations were drawn from the ambient stream, the identical call run
# twice returned different coefficients and a different objective, and every fit
# silently advanced the caller's stream.

mp_repro_data <- function(n = 150L) {
  set.seed(24)
  data.frame(
    t  = stats::rexp(n, 0.2),
    ev = rep(c(1, 0, 0), length.out = n)
  )
}

mp_repro_fit <- function(d, ...) {
  suppressWarnings(hazard(
    data = d, time = t, status = ev, fit = TRUE,
    dist = "multiphase",
    phases = list(hzr_phase("cdf", t_half = 0.15, nu = 1),
                  hzr_phase("constant")),
    theta = c(log(0.2), log(0.15), 1, 0, log(5e-04)),
    ...
  ))
}

# ---------------------------------------------------------------------------
# The defect: two identical calls, no set.seed() in between
# ---------------------------------------------------------------------------

test_that("two identical unseeded multiphase fits return identical theta", {
  d <- mp_repro_data()

  a <- mp_repro_fit(d, control = list(condition = 14))
  b <- mp_repro_fit(d, control = list(condition = 14))

  expect_identical(a$fit$theta, b$fit$theta)
  expect_identical(a$fit$objective, b$fit$objective)
})

test_that("a multiphase fit leaves the caller's RNG stream untouched", {
  d <- mp_repro_data()

  set.seed(99)
  before <- .Random.seed
  invisible(mp_repro_fit(d, control = list(condition = 14)))

  expect_identical(.Random.seed, before)
})

test_that("interleaved RNG use does not change the multiphase fit", {
  d <- mp_repro_data()

  a <- mp_repro_fit(d, control = list(condition = 14))
  invisible(stats::rnorm(17))          # move the ambient stream
  b <- mp_repro_fit(d, control = list(condition = 14))

  expect_identical(a$fit$theta, b$fit$theta)
})

# ---------------------------------------------------------------------------
# control$start_seed selects the start ensemble
# ---------------------------------------------------------------------------

test_that("the same start_seed reproduces a fit and a different one may not", {
  d <- mp_repro_data()

  a  <- mp_repro_fit(d, control = list(condition = 14, start_seed = 7L))
  a2 <- mp_repro_fit(d, control = list(condition = 14, start_seed = 7L))
  b  <- mp_repro_fit(d, control = list(condition = 14, start_seed = 8L))

  expect_identical(a$fit$theta, a2$fit$theta)
  # A different ensemble must actually be a different ensemble.  The optimum it
  # reaches may coincide, so assert on the draws rather than on the fit.
  expect_false(identical(
    .hzr_start_perturbations(5L, 5L, seed = 7L),
    .hzr_start_perturbations(5L, 5L, seed = 8L)
  ))
})

test_that("start_seed is rejected when it is not a single whole number", {
  d <- mp_repro_data()

  expect_error(mp_repro_fit(d, control = list(start_seed = c(1L, 2L))),
               "start_seed")
  expect_error(mp_repro_fit(d, control = list(start_seed = NA_integer_)),
               "start_seed")
  expect_error(mp_repro_fit(d, control = list(start_seed = Inf)),
               "start_seed")
  expect_error(mp_repro_fit(d, control = list(start_seed = "3")),
               "start_seed")

  # Out of integer range. Left to set.seed() this raises "supplied seed is not
  # a valid integer", after a coercion warning, and names neither the argument
  # nor the fit it came from.
  expect_error(mp_repro_fit(d, control = list(start_seed = 3e9)),
               "start_seed")

  # Fractional seeds are rejected, not truncated: set.seed() truncates, so 3.9
  # and 3 would silently select the same ensemble and a seed sweep over
  # 3.1/3.5/3.9 would report three fits having tried one set of starts.
  expect_error(mp_repro_fit(d, control = list(start_seed = 3.9)),
               "start_seed")
  set.seed(3.9)
  from_fractional <- stats::rnorm(3)
  set.seed(3L)
  from_whole <- stats::rnorm(3)
  expect_identical(from_fractional, from_whole)
})

test_that("a negative start_seed is accepted and reproducible", {
  # set.seed() takes negative integers; rejecting them would remove half the
  # seed space for no reason.
  d <- mp_repro_data()

  a <- mp_repro_fit(d, control = list(condition = 14, start_seed = -7L))
  b <- mp_repro_fit(d, control = list(condition = 14, start_seed = -7L))

  expect_identical(a$fit$theta, b$fit$theta)
  expect_false(identical(
    .hzr_start_perturbations(5L, 5L, seed = -7L),
    .hzr_start_perturbations(5L, 5L, seed = 7L)
  ))
})

# ---------------------------------------------------------------------------
# .hzr_start_perturbations(): the one place the fitting path touches the RNG
# ---------------------------------------------------------------------------

test_that(".hzr_start_perturbations is deterministic and correctly shaped", {
  p1 <- .hzr_start_perturbations(5L, 3L, seed = 1L)
  p2 <- .hzr_start_perturbations(5L, 3L, seed = 1L)

  expect_identical(p1, p2)
  expect_length(p1, 4L)                       # one per start after the first
  expect_true(all(vapply(p1, length, 1L) == 3L))
  expect_true(all(vapply(p1, function(z) all(is.finite(z)), TRUE)))
  # The draws must differ from one another, or the "multi" in multi-start is a lie.
  expect_false(identical(p1[[1L]], p1[[2L]]))
})

test_that(".hzr_start_perturbations returns nothing for a single start", {
  expect_identical(.hzr_start_perturbations(1L, 3L, seed = 1L), list())
  expect_identical(.hzr_start_perturbations(0L, 3L, seed = 1L), list())
})

test_that(".hzr_start_perturbations restores an existing .Random.seed", {
  set.seed(4321)
  before <- .Random.seed

  invisible(.hzr_start_perturbations(5L, 3L, seed = 1L))

  expect_identical(.Random.seed, before)
})

test_that(".hzr_start_perturbations leaves no .Random.seed where there was none", {
  # A fresh session has no .Random.seed; the helper must not create one, or the
  # first post-fit draw in that session would follow from our seed, not R's.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    saved <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", saved, envir = globalenv()), add = TRUE)
    rm(".Random.seed", envir = globalenv())
  }

  invisible(.hzr_start_perturbations(5L, 3L, seed = 1L))

  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

# ---------------------------------------------------------------------------
# The default start_seed is a behavioral choice, not an arbitrary constant
# ---------------------------------------------------------------------------

test_that("the default start_seed converges on the reference fixture", {
  skip_on_cran()

  # The fixture from test-multiphase-gradient.R.  That test guards itself with
  # skip_if(), so a bad default would turn it into a permanent skip rather than
  # a failure.  Assert the default here, where it fails loudly instead.  Its
  # assembled start used to fail on its own -- see the n_starts = 1 test below,
  # which pins the fix that made it converge.
  set.seed(42)
  n <- 100
  time   <- stats::rexp(n, rate = 0.5) + 0.01
  status <- sample(0:1, n, replace = TRUE, prob = c(0.2, 0.8))
  phases <- list(early = hzr_phase("cdf", t_half = 1, nu = 1.5, m = 0),
                 const = hzr_phase("constant"))

  fit <- suppressWarnings(hazard(
    time = time, status = status, dist = "multiphase",
    phases = phases, fit = TRUE,
    control = list(n_starts = 3, maxit = 500)   # default start_seed
  ))

  expect_true(fit$fit$converged)
})


# ---------------------------------------------------------------------------
# Per-start outcomes and the diagnostics that name a discarded error
# ---------------------------------------------------------------------------

test_that("the assembled start converges on its own", {
  skip_on_cran()

  # The regression that motivated `starts`: this fixture's start 1 raised
  # "$ operator is invalid for atomic vectors" from the CoE adjustment, the
  # loop discarded it, and the fit was reported as a convergence failure.  At
  # n_starts = 1 there is nothing to hide behind.
  set.seed(42)
  n <- 100
  time   <- stats::rexp(n, rate = 0.5) + 0.01
  status <- sample(0:1, n, replace = TRUE, prob = c(0.2, 0.8))
  phases <- list(early = hzr_phase("cdf", t_half = 1, nu = 1.5, m = 0),
                 const = hzr_phase("constant"))

  fit <- suppressWarnings(hazard(
    time = time, status = status, dist = "multiphase",
    phases = phases, fit = TRUE,
    control = list(n_starts = 1, maxit = 500)
  ))

  expect_true(fit$fit$converged)
  expect_identical(fit$fit$starts$status, "ok")
})


test_that("every start_seed converges on the reference fixture", {
  skip_on_cran()

  # Before the cumhaz guard was shape-corrected, 12 of these 50 seeds failed
  # outright, so the default seed was load-bearing.  Sweep rather than trust
  # the default: a single passing seed is not evidence.
  set.seed(42)
  n <- 100
  time   <- stats::rexp(n, rate = 0.5) + 0.01
  status <- sample(0:1, n, replace = TRUE, prob = c(0.2, 0.8))
  phases <- list(early = hzr_phase("cdf", t_half = 1, nu = 1.5, m = 0),
                 const = hzr_phase("constant"))

  failures <- vapply(1:50, function(s) {
    fit <- tryCatch(
      suppressWarnings(hazard(
        time = time, status = status, dist = "multiphase",
        phases = phases, fit = TRUE,
        control = list(n_starts = 3, maxit = 500, start_seed = s)
      )),
      error = function(e) NULL
    )
    is.null(fit) || !isTRUE(fit$fit$converged)
  }, logical(1))

  expect_equal(sum(failures), 0L)
})


test_that("starts records every start and marks the one that won", {
  skip_on_cran()

  set.seed(42)
  n <- 100
  time   <- stats::rexp(n, rate = 0.5) + 0.01
  status <- sample(0:1, n, replace = TRUE, prob = c(0.2, 0.8))
  phases <- list(early = hzr_phase("cdf", t_half = 1, nu = 1.5, m = 0),
                 const = hzr_phase("constant"))

  fit <- suppressWarnings(hazard(
    time = time, status = status, dist = "multiphase",
    phases = phases, fit = TRUE,
    control = list(n_starts = 5, maxit = 500)
  ))

  starts <- fit$fit$starts
  expect_s3_class(starts, "data.frame")
  expect_identical(nrow(starts), 5L)
  expect_named(starts, c("start", "status", "objective", "convergence",
                         "best", "message"))
  expect_identical(starts$start, 1:5)
  expect_true(all(starts$convergence == 0L))

  # Exactly one winner, and it is the best objective among the usable starts.
  expect_identical(sum(starts$best), 1L)
  expect_equal(starts$objective[starts$best], max(starts$objective, na.rm = TRUE))
  expect_equal(starts$objective[starts$best], fit$fit$objective)

  # A start that errored carries its message; one that did not is NA.
  expect_true(all(is.na(starts$message[starts$status == "ok"])))
})


test_that("a start that errors is named, not silently discarded", {
  skip_on_cran()

  # Force one start to throw by making the objective error on a marked theta,
  # so the discard path is exercised without waiting for a real defect.
  set.seed(42)
  n <- 100
  time   <- stats::rexp(n, rate = 0.5) + 0.01
  status <- sample(0:1, n, replace = TRUE, prob = c(0.2, 0.8))
  phases <- list(early = hzr_phase("cdf", t_half = 1, nu = 1.5, m = 0),
                 const = hzr_phase("constant"))

  calls <- 0L
  orig <- TemporalHazard:::.hzr_optim_generic
  stub <- function(...) {
    calls <<- calls + 1L
    if (calls == 1L) stop("synthetic start failure", call. = FALSE)
    orig(...)
  }

  fit <- testthat::with_mocked_bindings(
    suppressWarnings(hazard(
      time = time, status = status, dist = "multiphase",
      phases = phases, fit = TRUE,
      control = list(n_starts = 3, maxit = 500)
    )),
    .hzr_optim_generic = stub,
    .package = "TemporalHazard"
  )

  starts <- fit$fit$starts
  expect_identical(starts$status[1], "error")
  expect_match(starts$message[1], "synthetic start failure")
  expect_false(starts$best[1])
  expect_true(any(starts$best))

  # And it is warned about rather than absorbed.
  expect_warning(
    testthat::with_mocked_bindings(
      {
        calls <- 0L
        hazard(time = time, status = status, dist = "multiphase",
               phases = phases, fit = TRUE,
               control = list(n_starts = 3, maxit = 500))
      },
      .hzr_optim_generic = stub,
      .package = "TemporalHazard"
    ),
    "synthetic start failure"
  )
})


test_that("failing every start names the error instead of calling it convergence", {
  skip_on_cran()

  set.seed(42)
  n <- 100
  time   <- stats::rexp(n, rate = 0.5) + 0.01
  status <- sample(0:1, n, replace = TRUE, prob = c(0.2, 0.8))
  phases <- list(early = hzr_phase("cdf", t_half = 1, nu = 1.5, m = 0),
                 const = hzr_phase("constant"))

  stub <- function(...) stop("synthetic total failure", call. = FALSE)

  expect_error(
    testthat::with_mocked_bindings(
      suppressWarnings(hazard(
        time = time, status = status, dist = "multiphase",
        phases = phases, fit = TRUE,
        control = list(n_starts = 3, maxit = 500)
      )),
      .hzr_optim_generic = stub,
      .package = "TemporalHazard"
    ),
    "synthetic total failure"
  )
})


test_that("a start that stops at maxit is not reported as ok", {
  skip_on_cran()

  # A finite objective is not a converged one. optim() returns code 1 when it
  # stops at maxit and attaches a perfectly ordinary value, so a status keyed
  # only on is.finite() would call every one of these "ok" -- and one of them
  # can carry the best objective and become the reported fit.
  #
  # conserve = FALSE is load-bearing: with CoE on, a parameter is fixed, which
  # turns on the Nelder-Mead warm-up. That warm-up runs at its own hardcoded
  # maxit, reaches the optimum, and leaves BFGS nothing to do, so `maxit` here
  # would not bite at all.
  set.seed(42)
  n <- 100
  time   <- stats::rexp(n, rate = 0.5) + 0.01
  status <- sample(0:1, n, replace = TRUE, prob = c(0.2, 0.8))
  phases <- list(early = hzr_phase("cdf", t_half = 1, nu = 1.5, m = 0),
                 const = hzr_phase("constant"))

  fit <- suppressWarnings(hazard(
    time = time, status = status, dist = "multiphase",
    phases = phases, fit = TRUE,
    control = list(n_starts = 3, maxit = 2, conserve = FALSE)
  ))

  starts <- fit$fit$starts
  expect_identical(starts$status, rep("nonconverged", 3L))
  expect_true(all(starts$convergence == 1L))
  expect_true(all(is.finite(starts$objective)))   # finite, and still not "ok"
  expect_false(fit$fit$converged)

  # The winner is still chosen on the objective, unchanged -- the point is that
  # you can now see it did not converge.
  expect_identical(sum(starts$best), 1L)
  expect_equal(starts$objective[starts$best], max(starts$objective))
})
