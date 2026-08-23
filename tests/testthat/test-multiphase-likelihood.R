# test-multiphase-likelihood.R — Tests for multiphase likelihood engine
#
# Tests theta splitting, cumhaz/hazard computation, log-likelihood
# evaluation at known parameters, and optimizer convergence on synthetic data.

# ============================================================================
# Helpers
# ============================================================================

# Simple synthetic data: exponential-ish with known rate
make_test_data <- function(n = 200, seed = 42) {
  set.seed(seed)
  time   <- rexp(n, rate = 0.3) + 0.01
  status <- rbinom(n, 1, prob = 0.7)
  age    <- rnorm(n, mean = 60, sd = 10)
  data.frame(time = time, status = status, age = age)
}


# ============================================================================
# .hzr_split_theta()
# ============================================================================

test_that(".hzr_split_theta partitions correctly for 2-phase model", {
  phases <- .hzr_validate_phases(list(
    early = hzr_phase("cdf"),
    late  = hzr_phase("hazard")
  ))
  cov_counts <- c(early = 0L, late = 0L)

  # early: log_mu + log_t_half + nu + m = 4

  # late:  log_mu + log_t_half + nu + m = 4
  theta <- 1:8
  split <- .hzr_split_theta(theta, phases, cov_counts)

  expect_equal(split$early, 1:4)
  expect_equal(split$late, 5:8)
})

test_that(".hzr_split_theta handles constant phase (1 param)", {
  phases <- .hzr_validate_phases(list(
    early = hzr_phase("cdf"),
    bg    = hzr_phase("constant"),
    late  = hzr_phase("hazard")
  ))
  cov_counts <- c(early = 0L, bg = 0L, late = 0L)

  # early: 4, bg: 1 (log_mu only), late: 4 = total 9
  theta <- 1:9
  split <- .hzr_split_theta(theta, phases, cov_counts)

  expect_equal(split$early, 1:4)
  expect_equal(split$bg, 5)
  expect_equal(split$late, 6:9)
})

test_that(".hzr_split_theta includes covariates", {
  phases <- .hzr_validate_phases(list(
    early = hzr_phase("cdf"),
    late  = hzr_phase("hazard")
  ))
  cov_counts <- c(early = 2L, late = 1L)

  # early: 4 + 2 = 6, late: 4 + 1 = 5, total = 11
  theta <- seq(0.1, 1.1, by = 0.1)
  split <- .hzr_split_theta(theta, phases, cov_counts)

  expect_equal(length(split$early), 6)
  expect_equal(length(split$late), 5)
})


# ============================================================================
# .hzr_unpack_phase_theta()
# ============================================================================

test_that("unpack_phase_theta works for cdf phase", {
  phase <- hzr_phase("cdf")
  theta_j <- c(-2.3, log(2), 1.5, 0.3)  # log_mu, log_t_half, nu, m
  pars <- .hzr_unpack_phase_theta(theta_j, phase)

  expect_equal(pars$log_mu, -2.3)
  expect_equal(pars$log_t_half, log(2))
  expect_equal(pars$nu, 1.5)
  expect_equal(pars$m, 0.3)
  expect_equal(length(pars$beta), 0)
})

test_that("unpack_phase_theta works for constant phase", {
  phase <- hzr_phase("constant")
  theta_j <- c(-1.5, 0.02)  # log_mu, beta_1
  pars <- .hzr_unpack_phase_theta(theta_j, phase)

  expect_equal(pars$log_mu, -1.5)
  expect_true(is.na(pars$log_t_half))
  expect_true(is.na(pars$nu))
  expect_true(is.na(pars$m))
  expect_equal(pars$beta, 0.02)
})


# ============================================================================
# .hzr_multiphase_cumhaz() — known-parameter evaluation
# ============================================================================

test_that("single constant phase reduces to mu * t", {
  phases <- .hzr_validate_phases(list(bg = hzr_phase("constant")))
  cov_counts <- c(bg = 0L)
  x_list <- list(bg = NULL)

  # log_mu = log(0.5) => mu = 0.5
  theta <- log(0.5)
  time <- c(1, 2, 5)

  cumhaz <- .hzr_multiphase_cumhaz(time, theta, phases, cov_counts, x_list)
  expect_equal(cumhaz, 0.5 * time, tolerance = 1e-12)
})

test_that("single cdf phase produces mu * G(t)", {
  phases <- .hzr_validate_phases(list(early = hzr_phase("cdf")))
  cov_counts <- c(early = 0L)
  x_list <- list(early = NULL)

  # log_mu = 0 => mu = 1; log_t_half = log(3); nu = 2; m = 0
  theta <- c(0, log(3), 2, 0)
  time <- seq(0.5, 10, by = 0.5)

  cumhaz <- .hzr_multiphase_cumhaz(time, theta, phases, cov_counts, x_list)
  expected_G <- hzr_decompos(time, t_half = 3, nu = 2, m = 0)$G
  expect_equal(cumhaz, expected_G, tolerance = 1e-10)
})

test_that("two-phase additive: constant + cdf", {
  phases <- .hzr_validate_phases(list(
    bg    = hzr_phase("constant"),
    early = hzr_phase("cdf")
  ))
  cov_counts <- c(bg = 0L, early = 0L)
  x_list <- list(bg = NULL, early = NULL)

  # bg: log_mu = log(0.1) => mu_bg = 0.1
  # early: log_mu = log(0.5), log_t_half = log(2), nu = 1, m = 0
  theta <- c(log(0.1), log(0.5), log(2), 1, 0)
  time <- c(1, 5, 10)

  cumhaz <- .hzr_multiphase_cumhaz(time, theta, phases, cov_counts, x_list)
  G <- hzr_decompos(time, t_half = 2, nu = 1, m = 0)$G
  expected <- 0.1 * time + 0.5 * G
  expect_equal(cumhaz, expected, tolerance = 1e-10)
})

test_that("per_phase = TRUE returns decomposed contributions", {
  phases <- .hzr_validate_phases(list(
    early = hzr_phase("cdf"),
    late  = hzr_phase("hazard")
  ))
  cov_counts <- c(early = 0L, late = 0L)
  x_list <- list(early = NULL, late = NULL)

  theta <- c(log(0.5), log(2), 2, 0,   # early
             log(0.3), log(5), 1, 0)    # late

  time <- c(1, 3, 7)
  result <- .hzr_multiphase_cumhaz(time, theta, phases, cov_counts, x_list,
                                     per_phase = TRUE)

  expect_true(is.list(result))
  expect_true("total" %in% names(result))
  expect_true("early" %in% names(result))
  expect_true("late" %in% names(result))
  expect_equal(result$total, result$early + result$late, tolerance = 1e-12)
})


test_that("per_phase = TRUE keeps its shape at an infeasible decomposition", {
  # m < 0 together with nu < 0 is not a usable decomposition, so the cumhaz
  # short-circuits to Inf and lets the caller penalise the step.  It has to
  # short-circuit in the shape the caller asked for: the CoE wrapper reads
  # decomp$total on every objective evaluation, and a bare vector there throws
  # out of optim() and takes the whole optimization start with it.
  phases <- .hzr_validate_phases(list(
    early = hzr_phase("cdf"),
    const = hzr_phase("constant")
  ))
  cov_counts <- c(early = 0L, const = 0L)
  x_list <- list(early = NULL, const = NULL)

  theta <- c(log(0.5), log(2), -0.5, -0.5,   # early: nu < 0 and m < 0
             log(0.3))                        # const
  time <- c(1, 3, 7)

  result <- .hzr_multiphase_cumhaz(time, theta, phases, cov_counts, x_list,
                                     per_phase = TRUE)
  expect_true(is.list(result))
  expect_named(result, c("early", "const", "total"))
  expect_equal(result$total, rep(Inf, length(time)))
  expect_equal(result$early, rep(Inf, length(time)))

  # per_phase = FALSE is unchanged: the bare vector it always returned.
  expect_equal(
    .hzr_multiphase_cumhaz(time, theta, phases, cov_counts, x_list),
    rep(Inf, length(time))
  )
})


test_that("the CoE adjustment declines at an infeasible decomposition", {
  # .hzr_conserve_events() runs inside the objective on every evaluation, so it
  # meets whatever theta BFGS proposes -- including shapes the cumhaz refuses.
  # It must decline (return theta untouched) rather than error.
  phases <- .hzr_validate_phases(list(
    early = hzr_phase("cdf"),
    const = hzr_phase("constant")
  ))
  cov_counts <- c(early = 0L, const = 0L)
  x_list <- list(early = NULL, const = NULL)

  theta <- c(log(0.5), log(2), -0.5, -0.5, log(0.3))
  time   <- c(1, 3, 7)
  status <- c(1, 0, 1)

  expect_identical(
    .hzr_conserve_events(theta, "const", 5L, time, status,
                         phases, cov_counts, x_list, total_events = 2),
    theta
  )
})


# ============================================================================
# .hzr_multiphase_hazard()
# ============================================================================

test_that("constant phase hazard is mu * 1 = mu", {
  phases <- .hzr_validate_phases(list(bg = hzr_phase("constant")))
  cov_counts <- c(bg = 0L)
  x_list <- list(bg = NULL)

  theta <- log(0.5)
  time <- c(1, 2, 5)

  haz <- .hzr_multiphase_hazard(time, theta, phases, cov_counts, x_list)
  expect_equal(haz, rep(0.5, 3), tolerance = 1e-12)
})


# ============================================================================
# .hzr_logl_multiphase() — basic properties
# ============================================================================

test_that("logl is finite for valid parameters", {
  dat <- make_test_data(n = 50)
  phases <- .hzr_validate_phases(list(
    early = hzr_phase("cdf", t_half = 2, nu = 1, m = 0)
  ))
  cov_counts <- c(early = 0L)
  x_list <- list(early = NULL)

  theta <- c(log(0.1), log(2), 1, 0)

  ll <- .hzr_logl_multiphase(
    theta = theta, time = dat$time, status = dat$status,
    phases = phases, covariate_counts = cov_counts, x_list = x_list
  )

  expect_true(is.finite(ll))
  expect_true(ll < 0)  # log-likelihood should be negative
})

test_that("logl returns -Inf for m < 0 and nu < 0", {
  dat <- make_test_data(n = 50)
  phases <- .hzr_validate_phases(list(
    bad = hzr_phase("cdf", t_half = 2, nu = 1, m = 0)
  ))
  cov_counts <- c(bad = 0L)
  x_list <- list(bad = NULL)

  # Force invalid: nu = -1 and m = -1 in the theta vector
  theta <- c(log(0.1), log(2), -1, -1)

  ll <- .hzr_logl_multiphase(
    theta = theta, time = dat$time, status = dat$status,
    phases = phases, covariate_counts = cov_counts, x_list = x_list
  )

  expect_equal(ll, -Inf)
})

test_that("logl is higher at better parameters (constant-rate sanity check)", {
  # Generate data from a known constant-rate model: H(t) = lambda * t
  set.seed(99)
  true_lambda <- 0.3
  n <- 200
  time <- rexp(n, rate = true_lambda) + 0.01
  status <- rep(1L, n)  # all events

  phases <- .hzr_validate_phases(list(bg = hzr_phase("constant")))
  cov_counts <- c(bg = 0L)
  x_list <- list(bg = NULL)

  # True parameter
  ll_true <- .hzr_logl_multiphase(
    theta = log(true_lambda), time = time, status = status,
    phases = phases, covariate_counts = cov_counts, x_list = x_list
  )

  # Wrong parameter (much higher rate)
  ll_wrong <- .hzr_logl_multiphase(
    theta = log(5.0), time = time, status = status,
    phases = phases, covariate_counts = cov_counts, x_list = x_list
  )

  expect_true(ll_true > ll_wrong)
})


# ============================================================================
# .hzr_multiphase_theta_user() — back-transformation
# ============================================================================

test_that("theta_user back-transforms correctly", {
  phases <- .hzr_validate_phases(list(
    early = hzr_phase("cdf"),
    bg    = hzr_phase("constant")
  ))
  cov_counts <- c(early = 0L, bg = 0L)

  # internal: [log(0.5), log(3), 2, 0.1, log(0.2)]
  theta_int <- c(log(0.5), log(3), 2, 0.1, log(0.2))
  theta_usr <- .hzr_multiphase_theta_user(theta_int, phases, cov_counts)

  expect_equal(unname(theta_usr[1]), 0.5, tolerance = 1e-10)  # mu early
  expect_equal(unname(theta_usr[2]), 3.0, tolerance = 1e-10)  # t_half early
  expect_equal(unname(theta_usr[3]), 2.0)                     # nu early
  expect_equal(unname(theta_usr[4]), 0.1)                     # m early
  expect_equal(unname(theta_usr[5]), 0.2, tolerance = 1e-10)  # mu bg
})


# ============================================================================
# .hzr_optim_multiphase() — convergence on simple data
# ============================================================================

test_that("optimizer recovers constant rate from exponential data", {
  skip_on_cran()  # optimization test — may be slow

  set.seed(123)
  true_lambda <- 0.4
  n <- 300
  dat <- data.frame(
    time   = rexp(n, rate = true_lambda) + 0.01,
    status = rep(1L, n)
  )

  phases <- list(bg = hzr_phase("constant"))

  result <- .hzr_optim_multiphase(
    time = dat$time, status = dat$status,
    phases = phases,
    control = list(n_starts = 3, maxit = 500)
  )

  expect_equal(result$convergence, 0)

  # Recover mu ~ true_lambda (on internal scale, par is log_mu)
  mu_hat <- unname(exp(result$par[1]))
  expect_equal(mu_hat, true_lambda, tolerance = 0.15)
})

test_that("optimizer handles two phases (cdf + constant)", {
  skip_on_cran()

  # Simulate: H(t) = 0.3 * G(t; t_half=2, nu=1, m=0) + 0.1 * t
  set.seed(456)
  n <- 400
  t_grid <- seq(0.01, 15, length.out = n)

  # Compute survival under the model and simulate events
  G <- hzr_decompos(t_grid, t_half = 2, nu = 1, m = 0)$G
  cumhaz_true <- 0.3 * G + 0.1 * t_grid
  surv_true <- exp(-cumhaz_true)

  # Simple simulation: event if uniform < 1 - S(t)
  set.seed(789)
  u <- runif(n)
  status <- as.integer(u > surv_true)

  dat <- data.frame(time = t_grid, status = status)

  phases <- list(
    early = hzr_phase("cdf", t_half = 2, nu = 1, m = 0),
    bg    = hzr_phase("constant")
  )

  result <- suppressWarnings(
    .hzr_optim_multiphase(
      time = dat$time, status = dat$status,
      phases = phases,
      control = list(n_starts = 3, maxit = 1000)
    )
  )

  expect_true(is.finite(result$value))
  expect_equal(result$convergence, 0)

  # At minimum, the log-likelihood should be better than starting values
  theta_start <- unlist(lapply(names(.hzr_validate_phases(phases)), function(nm) {
    .hzr_phase_start(.hzr_validate_phases(phases)[[nm]])
  }))
  ll_start <- .hzr_logl_multiphase(
    theta = theta_start, time = dat$time, status = dat$status,
    phases = .hzr_validate_phases(phases),
    covariate_counts = c(early = 0L, bg = 0L),
    x_list = list(early = NULL, bg = NULL)
  )
  expect_true(result$value >= ll_start)
})


# ============================================================================
# Regression: 3+ covariate phase must not emit vector-recycling warnings
# ============================================================================

test_that("3-covariate phase does not emit mu_j * phi_j recycling warning", {
  # Reproduces the issue where model.matrix() silently drops NA rows in a
  # phase-specific formula, leaving x_list[[nm]] shorter than time and
  # triggering "longer object length is not a multiple of shorter object
  # length" on every mu_j * phi_j multiplication inside
  # .hzr_multiphase_cumhaz() / .hzr_multiphase_hazard().  The avc dataset
  # has 5 NAs in inc_surg, so a 3-covariate phase referencing it would
  # produce a 305-row design matrix against a 310-row time vector.
  skip_if_not(
    tryCatch({
      data("avc", package = "TemporalHazard", envir = environment())
      exists("avc", inherits = FALSE)
    }, error = function(e) FALSE),
    "avc dataset not available"
  )

  w <- testthat::capture_warnings(
    fit <- hazard(
      survival::Surv(int_dead, dead) ~ 1, data = avc,
      dist = "multiphase",
      phases = list(
        early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                              fixed = "shapes",
                              formula = ~ age + mal + inc_surg),
        constant = hzr_phase("constant")
      ),
      fit = TRUE,
      control = list(n_starts = 1L, maxit = 10L)
    )
  )
  # The bug under test produced a vector-recycling warning; assert it is absent.
  # (An incidental non-PD-Hessian warning at maxit = 10 is unrelated and allowed.)
  expect_false(any(grepl("not a multiple|longer object length", w)))
})
