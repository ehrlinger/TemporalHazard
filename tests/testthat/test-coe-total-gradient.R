# ---------------------------------------------------------------------------
# The gradient the optimizer is given under Conservation of Events (#565).
#
# Under CoE the objective is L(phi, c(phi)): the conserved phase's log_mu, c,
# is re-solved from the other parameters at every evaluation. The optimizer
# was handed the partial score at the conserved theta, which leaves out how c
# moves with phi, and on the avc six-variable model it stopped 1.80
# log-likelihood units short of SAS while reporting convergence.
#
# The oracle is a finite-difference gradient of the very objective the
# optimizer receives, so the comparison is between the two functions that
# have to agree, not between two analytic derivations.
# ---------------------------------------------------------------------------

# What .hzr_optim_multiphase() hands the optimizer, captured without fitting.
.coe_optimizer_inputs <- function(...) {
  seen <- NULL
  local_mocked_bindings(.hzr_optim_generic = function(...) {
    seen <<- list(...)
    stop("captured")
  })
  res <- suppressWarnings(
    try(hazard(..., dist = "multiphase", fit = TRUE), silent = TRUE)
  )
  if (is.null(seen)) stop("the optimizer was never reached: ", res)
  seen
}

.coe_grad_pair <- function(inp, theta) {
  obj <- function(p) {
    inp$logl_fn(p, inp$time, inp$status, inp$time_lower, inp$time_upper,
                inp$x, weights = inp$weights)
  }
  list(
    analytic = unname(inp$gradient_fn(theta, inp$time, inp$status,
                                      inp$time_lower, inp$time_upper, inp$x,
                                      weights = inp$weights)),
    oracle = numDeriv::grad(obj, theta, method.args = list(d = 1e-3))
  )
}

.coe_sim <- function(n = 300, seed = 11) {
  set.seed(seed)
  z <- rbinom(n, 1, 0.4)
  age <- round(rnorm(n, 0, 1), 2)
  u <- runif(n)
  t_event <- ifelse(u < 0.3, rexp(n, 4 * exp(0.5 * z)),
                    rexp(n, 0.15 * exp(0.3 * age)))
  cens <- runif(n, 0.5, 8)
  data.frame(
    stop = pmin(t_event, cens) + 0.01,
    event = as.integer(t_event <= cens),
    start = 0, z = z, age = age,
    w = sample(c(0.5, 1, 2), n, replace = TRUE)
  )
}

# One cell of the grid: the optimizer's gradient against the oracle at the
# start and at two displaced points, where the score is far from zero.
.expect_coe_gradient <- function(inp, label, conserved = NULL) {
  env <- environment(inp$gradient_fn)
  expect_true(isTRUE(env$use_conserve), label = paste(label, "runs under CoE"))
  if (!is.null(conserved)) {
    expect_identical(env$fixmu_phase, conserved, label = label)
  }
  expect_true(isTRUE(inp$gradient_exact),
              label = paste(label, "declares its gradient exact"))
  p <- length(inp$theta_start)
  # Displaced from the start, which the warm-up may have left near an optimum.
  offsets <- list(0.15 * cos(seq_len(p)), -0.1 * sin(2 * seq_len(p)),
                  0.2 * sin(seq_len(p) + 1))
  for (off in offsets) {
    theta <- inp$theta_start + off
    pair <- .coe_grad_pair(inp, theta)
    # The comparison is only evidence where there is a score to compare.
    expect_gt(max(abs(pair$oracle)), 1e-2, label = paste(label, "oracle size"))
    expect_equal(pair$analytic, pair$oracle, tolerance = 1e-4, label = label)
  }
}

test_that("the CoE gradient is the gradient of the CoE objective", {
  skip_if_not_installed("numDeriv")
  d <- .coe_sim()

  # Two phases, free shapes, covariates in both.
  two <- list(
    early = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1, formula = ~ z),
    constant = hzr_phase("constant", formula = ~ age)
  )
  inp <- .coe_optimizer_inputs(survival::Surv(stop, event) ~ 1, data = d, phases = two,
                               control = list(n_starts = 1))
  .expect_coe_gradient(inp, "two phases")

  # Fixed shapes.
  fixed <- list(
    early = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1, fixed = "shapes",
                      formula = ~ z + age),
    constant = hzr_phase("constant", formula = ~ age)
  )
  inp <- .coe_optimizer_inputs(survival::Surv(stop, event) ~ 1, data = d, phases = fixed,
                               control = list(n_starts = 1))
  .expect_coe_gradient(inp, "fixed shapes")

  # Weights.
  inp <- .coe_optimizer_inputs(survival::Surv(stop, event) ~ 1, data = d, phases = fixed,
                               weights = w, control = list(n_starts = 1))
  expect_false(all(inp$weights == 1))
  .expect_coe_gradient(inp, "weights")

  # Left truncation: the conserved quantity is H(stop) - H(start).
  set.seed(5)
  dl <- d
  dl$start <- pmin(runif(nrow(d), 0, 0.5), 0.6 * d$stop)
  inp <- .coe_optimizer_inputs(survival::Surv(start, stop, event) ~ 1, data = dl,
                               phases = fixed, control = list(n_starts = 1))
  expect_true(any(inp$time_lower > 0))
  .expect_coe_gradient(inp, "left truncation")

  # Weights and left truncation together.
  inp <- .coe_optimizer_inputs(survival::Surv(start, stop, event) ~ 1, data = dl,
                               phases = fixed, weights = w,
                               control = list(n_starts = 1))
  expect_true(any(inp$time_lower > 0) && !all(inp$weights == 1))
  .expect_coe_gradient(inp, "weights and truncation")

  # A derived shape, whose score is folded into its sources' afterwards.
  tied <- list(
    early = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1, fixed = "shapes",
                      formula = ~ z),
    late = hzr_phase("g3", tau = 4, gamma = 2, alpha = 1, eta = 1,
                     constraint = "eta_gamma", formula = ~ age)
  )
  inp <- .coe_optimizer_inputs(survival::Surv(stop, event) ~ 1, data = d,
                               phases = tied, control = list(n_starts = 1))
  .expect_coe_gradient(inp, "derived shape")
})

test_that("the CoE gradient is exact whichever phase is conserved", {
  skip_if_not_installed("numDeriv")
  d <- .coe_sim()
  three <- list(
    early = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1, fixed = "shapes",
                      formula = ~ z),
    constant = hzr_phase("constant", formula = ~ age),
    late = hzr_phase("g3", tau = 4, gamma = 2, alpha = 1, eta = 1,
                     fixed = "shapes", formula = ~ age)
  )
  seen <- character()
  for (big in names(three)) {
    # Which phase is conserved is chosen from the data; name it here so that
    # every phase takes the role.
    local_mocked_bindings(.hzr_select_fixmu_phase = function(...) big)
    inp <- .coe_optimizer_inputs(survival::Surv(stop, event) ~ 1, data = d,
                                 phases = three, control = list(n_starts = 1))
    .expect_coe_gradient(inp, paste("conserving", big), conserved = big)
    seen <- c(seen, environment(inp$gradient_fn)$fixmu_phase)
  }
  expect_setequal(seen, names(three))
})

test_that("where CoE has nothing to solve, the partial score is the gradient", {
  skip_if_not_installed("numDeriv")
  d <- .coe_sim()
  ph <- list(
    early = hzr_phase("cdf", t_half = 0.3, nu = 1, m = 1, fixed = "shapes",
                      formula = ~ z),
    constant = hzr_phase("constant", formula = ~ age)
  )
  inp <- .coe_optimizer_inputs(survival::Surv(stop, event) ~ 1, data = d,
                               phases = ph, control = list(n_starts = 1))
  env <- environment(inp$gradient_fn)
  # Raise the free phase's scale until it alone predicts more events than
  # were observed: the conserved phase would need a negative share, so the
  # solve declines and its log_mu stays where it was.
  free_mu <- setdiff(unlist(env$log_mu_positions), env$fixmu_pos)
  theta <- inp$theta_start
  theta[match(free_mu, env$free_idx)] <- theta[match(free_mu, env$free_idx)] + 4
  coe <- .hzr_conserve_events(
    env$expand_theta(theta), env$fixmu_phase, env$fixmu_pos, inp$time,
    inp$status, env$phases, env$covariate_counts, env$x_list,
    env$total_events, weights = inp$weights, time_lower = inp$time_lower,
    details = TRUE
  )
  expect_false(coe$solved)
  pair <- .coe_grad_pair(inp, theta)
  expect_gt(max(abs(pair$oracle)), 1)
  expect_equal(pair$analytic, pair$oracle, tolerance = 1e-4)
})

test_that("the default start reaches SAS's optimum on the avc model (#565)", {
  skip_on_cran()
  data(avc, package = "TemporalHazard", envir = environment())
  a <- stats::na.omit(avc)
  fit <- hazard(
    survival::Surv(int_dead, dead) ~ 1, data = a, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.1512095, nu = 1.438652, m = 1,
                        fixed = c("t_half", "nu", "m"),
                        formula = ~ age + mal + orifice + com_iv),
      constant = hzr_phase("constant", formula = ~ age + inc_surg)
    ),
    fit = TRUE, control = list(n_starts = 1)
  )
  expect_true(fit$spec$control$conserve_applied)
  # SAS prints -182.659; the partial-score fit stopped at -184.462.
  expect_equal(fit$fit$objective, -182.659, tolerance = 1e-5)
  expect_lte(fit$fit$rel_gradient, .Machine$double.eps^(1 / 3))
})
