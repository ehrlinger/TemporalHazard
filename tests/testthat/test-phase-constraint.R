# test-phase-constraint.R
# hzr_phase(constraint = ): a g3 shape derived from the others, as SAS/C's
# FIXGAE2 (alpha = gamma*eta/2) and FIXGE2 (eta = 2/gamma) do (#325).
#
# The derived slot keeps its place in theta and is masked like a fixed one,
# so what these tests guard is the calculus a fixed mask alone gets wrong:
# the score, the Hessian and the covariance must all carry the derived
# parameter's dependence on its sources.

# One g3 phase's survival times, drawn by inverting its cumulative hazard.
constraint_data <- function(theta, n = 1500) {
  withr::local_seed(20260916)
  ph <- hzr_phase("g3", tau = exp(theta[2]), gamma = theta[3],
                  alpha = theta[4], eta = theta[5])
  grid <- seq(0.001, 60, length.out = 20000)
  cumhaz <- .hzr_multiphase_cumhaz(grid, theta, list(late = ph), c(late = 0L),
                                   list(late = NULL))
  event_time <- stats::approx(cumhaz, grid, xout = -log(stats::runif(n)),
                              rule = 2, ties = "ordered")$y
  censor <- stats::runif(n, 2, 25)
  list(time = pmin(event_time, censor),
       status = as.numeric(event_time <= censor))
}

# Positions of the searched parameters in the one-phase theta
# [log_mu, log_tau, gamma, alpha, eta].
constraint_free <- function(constraint) {
  if (constraint == "alpha_gamma_eta") c(1L, 2L, 3L, 5L) else 1:4
}

# The log-likelihood as a function of the searched parameters only, written
# from the public-facing pieces rather than from the optimizer's wrappers.
constraint_reduced_ll <- function(d, phases, constraint,
                                  free = constraint_free(constraint),
                                  base = numeric(5)) {
  counts <- c(late = 0L)
  function(par) {
    theta <- base
    theta[free] <- par
    theta <- .hzr_apply_constraints(theta, phases, counts)
    haz <- .hzr_multiphase_hazard(d$time[d$status == 1], theta, phases, counts,
                                  list(late = NULL))
    sum(log(haz)) -
      sum(.hzr_multiphase_cumhaz(d$time, theta, phases, counts,
                                 list(late = NULL)))
  }
}

# ---------------------------------------------------------------------------
# The constructor
# ---------------------------------------------------------------------------

test_that("the derived start is computed, and a differing supplied value warns", {
  ph <- hzr_phase("g3", tau = 14, gamma = 22, eta = 0.18,
                  constraint = "alpha_gamma_eta")
  expect_equal(ph$alpha, 22 * 0.18 / 2)
  expect_identical(ph$constraint, "alpha_gamma_eta")

  expect_warning(
    ph2 <- hzr_phase("g3", gamma = 4, alpha = 3, eta = 0.5,
                     constraint = "alpha_gamma_eta"),
    "alpha = 3 was replaced by 1"
  )
  expect_equal(ph2$alpha, 1)
  # Supplying the value the rule gives is not a replacement.
  expect_no_warning(hzr_phase("g3", gamma = 4, alpha = 1, eta = 0.5,
                              constraint = "alpha_gamma_eta"))

  expect_equal(hzr_phase("g3", gamma = 4, constraint = "eta_gamma")$eta, 0.5)
})

test_that("a derived shape cannot be fixed, and 'shapes' leaves it out", {
  expect_error(
    hzr_phase("g3", constraint = "alpha_gamma_eta", fixed = "alpha"),
    "alpha cannot be fixed"
  )
  expect_error(
    hzr_phase("g3", constraint = "eta_gamma", fixed = c("tau", "eta")),
    "eta cannot be fixed"
  )
  ph <- hzr_phase("g3", constraint = "eta_gamma", fixed = "shapes")
  expect_setequal(ph$fixed, c("tau", "gamma", "alpha"))
})

test_that("a constraint on a non-g3 phase is refused", {
  expect_error(hzr_phase("cdf", constraint = "alpha_gamma_eta"),
               "applies only to g3 phases")
  expect_identical(hzr_phase("cdf")$constraint, NULL)
})

test_that("the free mask drops exactly the derived slot", {
  counts <- c(early = 0L, late = 0L)
  early <- hzr_phase("cdf")
  expect_equal(
    .hzr_phase_free_mask(
      list(early = early,
           late = hzr_phase("g3", constraint = "alpha_gamma_eta")), counts),
    c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE)
  )
  expect_equal(
    .hzr_phase_free_mask(
      list(early = early, late = hzr_phase("g3", constraint = "eta_gamma")),
      counts),
    c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
  )
})

test_that("a phase saved before the argument existed reads as unconstrained", {
  ph <- hzr_phase("g3", gamma = 4, alpha = 3, eta = 0.5)
  ph$constraint <- NULL
  expect_identical(.hzr_phase_constraint(ph), "none")
  theta <- c(0, 0, 4, 3, 0.5)
  expect_identical(.hzr_apply_constraints(theta, list(late = ph), c(late = 0L)),
                   theta)
})

test_that("print shows the derivation", {
  expect_output(print(hzr_phase("g3", constraint = "eta_gamma")),
                "derived: eta = 2 / gamma")
})

# ---------------------------------------------------------------------------
# The calculus, against numDeriv on the reduced likelihood
# ---------------------------------------------------------------------------

# Points across the shape range. Beyond gamma of about 30 the analytic g3
# Hessian this builds on is itself inaccurate against numDeriv (a separate
# defect in the unconstrained Hessian), so the range stops there; the score,
# which does not depend on it, is exact well beyond.
calculus_points <- list(
  interior  = c(log(0.25), log(9), 5, 1.5, 0.9),
  low_gamma = c(log(0.25), log(9), 0.5, 0.7, 3),
  moderate  = c(log(0.25), log(9), 1.2, 1.5, 2),
  high_gamma = c(log(0.25), log(9), 15, 1.5, 0.1)
)

for (constraint in c("alpha_gamma_eta", "eta_gamma")) {
  test_that(paste("score and Hessian through", constraint, "match numDeriv"), {
    skip_if_not_installed("numDeriv")
    d <- constraint_data(c(log(0.3), log(10), 6, 2.4, 0.8))
    phases <- list(late = hzr_phase("g3", constraint = constraint))
    counts <- c(late = 0L)
    free <- constraint_free(constraint)
    ll <- constraint_reduced_ll(d, phases, constraint)
    curvature_seen <- 0

    for (point in names(calculus_points)) {
      # Arbitrary points, not optima: the curvature term is weighted by the
      # derived slot's own score, which none of these zero.
      theta <- .hzr_apply_constraints(calculus_points[[point]], phases, counts)

      score <- .hzr_gradient_multiphase(theta, d$time, d$status,
                                        phases = phases,
                                        covariate_counts = counts,
                                        x_list = list(late = NULL))
      folded <- .hzr_constraint_score(theta, score, phases, counts)[free]
      expect_equal(folded, numDeriv::grad(ll, theta[free]), tolerance = 1e-5,
                   label = paste("score at", point))

      hess <- .hzr_constrained_hessian(theta, d$time, d$status,
                                       phases = phases,
                                       covariate_counts = counts,
                                       x_list = list(late = NULL))[free, free]
      expect_equal(unname(hess), -numDeriv::hessian(ll, theta[free]),
                   tolerance = 1e-3, label = paste("Hessian at", point))

      plain <- .hzr_hessian_multiphase(theta, d$time, d$status,
                                       phases = phases,
                                       covariate_counts = counts,
                                       x_list = list(late = NULL))
      jac <- .hzr_constraint_jacobian(theta, phases, counts)
      first_order <- crossprod(jac, plain %*% jac)[free, free]
      curvature_seen <- max(curvature_seen,
                            max(abs(first_order - hess)) / max(abs(hess)))
    }
    # The second-order term is what a fixed mask would miss, and it is not
    # negligible here, so the Hessian comparisons above can fail without it.
    expect_gt(curvature_seen, 1e-3)
  })
}

# ---------------------------------------------------------------------------
# The covariance expansion
# ---------------------------------------------------------------------------

test_that("the derived variance is the delta method, and an NA does not spread", {
  phases <- list(late = hzr_phase("g3", constraint = "alpha_gamma_eta"))
  counts <- c(late = 0L)
  theta <- .hzr_apply_constraints(c(-1, 2, 6, 0, 0.8), phases, counts)
  free <- c(1L, 2L, 3L, 5L)
  v <- matrix(c(4, 1, 0.5, 0.2,
                1, 3, 0.3, 0.1,
                0.5, 0.3, 2, -0.4,
                0.2, 0.1, -0.4, 1), 4L, 4L)

  out <- .hzr_expand_vcov(v, free, theta, phases, counts)
  d1 <- c(0.8 / 2, 6 / 2)                       # d(alpha)/d(gamma, eta)
  expect_equal(out[free, free], v)
  expect_equal(out[4, 4], as.numeric(d1 %*% v[3:4, 3:4] %*% d1))
  expect_equal(out[4, free], as.vector(d1 %*% v[3:4, ]))
  expect_equal(out[free, 4], out[4, free])

  # One unresolved entry stays where it is.
  v_na <- v
  v_na[2, 2] <- NA_real_
  out_na <- .hzr_expand_vcov(v_na, free, theta, phases, counts)
  expect_equal(sum(is.na(out_na)), 1L)

  # Derived only from fixed parameters, there is nothing to carry.
  out_fixed <- .hzr_expand_vcov(v[1:2, 1:2], 1:2, theta, phases, counts)
  expect_true(all(is.na(out_fixed[4, ])))
})

# ---------------------------------------------------------------------------
# End to end
# ---------------------------------------------------------------------------

# eta_gamma is fitted with alpha held: with alpha free as well, simulated data
# of this size let gamma run to a ridge (gamma -> Inf, eta -> 0) where the
# fit says nothing about the constraint.
e2e_cases <- list(
  alpha_gamma_eta = list(truth = c(log(0.3), log(10), 6, 2.4, 0.8),
                         start = list(tau = 8, gamma = 5, eta = 1),
                         free = c(1L, 2L, 3L, 5L)),
  eta_gamma = list(truth = c(log(0.5), log(8), 2, 2, 1),
                   start = list(tau = 6, gamma = 3, alpha = 2,
                                fixed = "alpha"),
                   free = 1:3)
)

for (constraint in names(e2e_cases)) {
  test_that(paste("a fit under", constraint, "is the constrained optimum"), {
    skip_on_cran()
    skip_if_not_installed("numDeriv")
    case <- e2e_cases[[constraint]]
    d <- constraint_data(case$truth)
    phases <- list(late = do.call(hzr_phase, c(list("g3", constraint = constraint),
                                               case$start)))
    fit <- suppressWarnings(hazard(
      time = d$time, status = d$status, dist = "multiphase", fit = TRUE,
      phases = phases, control = list(conserve = FALSE, n_starts = 1L)
    ))
    theta <- unname(coef(fit))
    free <- case$free
    counts <- c(late = 0L)

    expect_true(fit$fit$converged)
    # SAS/C's acceptance test reads the score over the searched parameters,
    # which must include the derived slot's share; without it the constrained
    # optimum reads as unconverged.
    expect_lt(fit$fit$rel_gradient, .Machine$double.eps^(1 / 3))
    # The derived slot obeys its rule exactly, at the returned estimate.
    expect_identical(theta, unname(.hzr_apply_constraints(theta, phases,
                                                          counts)))

    # Independent optimum: nlm over the reduced likelihood, from the fit.
    ll <- constraint_reduced_ll(d, phases, constraint, free = free,
                                base = theta)
    polish <- stats::nlm(function(p) -ll(p), theta[free], gradtol = 1e-8)
    expect_equal(fit$fit$objective, -polish$minimum, tolerance = 1e-6)

    # The covariance is the inverse reduced information, carried onto the
    # derived slot by the Jacobian. Compared over every row that is not held
    # fixed, the derived one included.
    derived <- if (constraint == "alpha_gamma_eta") 4L else 5L
    rows <- sort(c(free, derived))
    jac <- .hzr_constraint_jacobian(theta, phases, counts)[, free]
    oracle <- jac %*% solve(-numDeriv::hessian(ll, theta[free])) %*% t(jac)
    expect_equal(unname(vcov(fit))[rows, rows], oracle[rows, rows],
                 tolerance = 1e-3)
    expect_gt(vcov(fit)[derived, derived], 0)

    # The derived standard error, with the derivative taken numerically too:
    # nothing here comes from the constraint helpers.
    derive <- function(p) {
      full <- theta
      full[free] <- p
      if (constraint == "alpha_gamma_eta") full[3] * full[5] / 2 else 2 / full[3]
    }
    grad_derived <- numDeriv::grad(derive, theta[free])
    v_numeric <- solve(-numDeriv::hessian(ll, theta[free]))
    expect_equal(sqrt(vcov(fit)[derived, derived]),
                 sqrt(as.numeric(grad_derived %*% v_numeric %*% grad_derived)),
                 tolerance = 1e-3)
  })
}

test_that("the constraint binds: freeing the derived shape fits a different model", {
  skip_on_cran()
  d <- constraint_data(c(log(0.3), log(10), 6, 2.4, 0.8))
  fit_one <- function(ph) {
    suppressWarnings(hazard(time = d$time, status = d$status,
                            dist = "multiphase", fit = TRUE,
                            phases = list(late = ph),
                            control = list(conserve = FALSE, n_starts = 1L)))
  }
  constrained <- fit_one(hzr_phase("g3", tau = 8, gamma = 5, eta = 1,
                                   constraint = "alpha_gamma_eta"))
  free <- fit_one(hzr_phase("g3", tau = 8, gamma = 5, alpha = 2.5, eta = 1))
  theta_c <- coef(constrained)
  theta_f <- coef(free)
  expect_equal(theta_c[["late.alpha"]],
               theta_c[["late.gamma"]] * theta_c[["late.eta"]] / 2)
  expect_gt(abs(theta_f[["late.alpha"]] -
                  theta_f[["late.gamma"]] * theta_f[["late.eta"]] / 2), 1e-3)
  expect_gte(free$fit$objective, constrained$fit$objective)
})

test_that("Conservation of Events runs through the constraint", {
  skip_on_cran()
  withr::local_seed(2)
  truth <- c(log(0.01), log(0.3), log(10), 6, 2.4, 0.8)
  phases0 <- list(const = hzr_phase("constant"),
                  late = hzr_phase("g3", tau = 10, gamma = 6, alpha = 2.4,
                                   eta = 0.8))
  counts <- c(const = 0L, late = 0L)
  nulls <- list(const = NULL, late = NULL)
  grid <- seq(0.001, 60, length.out = 20000)
  cumhaz <- .hzr_multiphase_cumhaz(grid, truth, phases0, counts, nulls)
  n <- 4000
  event_time <- stats::approx(cumhaz, grid, xout = -log(stats::runif(n)),
                              rule = 2, ties = "ordered")$y
  censor <- stats::runif(n, 2, 25)
  time <- pmin(event_time, censor)
  status <- as.numeric(event_time <= censor)

  phases <- list(const = hzr_phase("constant"),
                 late = hzr_phase("g3", tau = 8, gamma = 5, eta = 1,
                                  constraint = "alpha_gamma_eta"))
  fits <- lapply(c(FALSE, TRUE), function(conserve) {
    suppressWarnings(hazard(time = time, status = status, dist = "multiphase",
                            fit = TRUE, phases = phases,
                            theta = c(log(0.02), log(0.2), log(8), 1, 1, 1),
                            control = list(conserve = conserve,
                                           n_starts = 1L)))
  })
  expect_false(fits[[1]]$spec$control$conserve_applied)
  expect_true(fits[[2]]$spec$control$conserve_applied)
  theta <- coef(fits[[2]])
  expect_equal(theta[["late.alpha"]],
               theta[["late.gamma"]] * theta[["late.eta"]] / 2)
  # CoE's solution is the unconstrained-mu MLE, so both land together, and
  # the full-information recompute carries the derived variance too.
  expect_equal(fits[[2]]$fit$objective, fits[[1]]$fit$objective,
               tolerance = 1e-6)
  expect_equal(unname(vcov(fits[[2]])), unname(vcov(fits[[1]])),
               tolerance = 1e-2)
  expect_true(is.finite(vcov(fits[[2]])[5, 5]))
})

test_that("the derived shape is kept out of the weak-direction check and z tests", {
  skip_on_cran()
  d <- constraint_data(c(log(0.3), log(10), 6, 2.4, 0.8))
  seen <- NULL
  # The check fires only on an ill-conditioned fit, so capture what it is
  # given rather than wait for one: the derived row, an exact function of
  # gamma and eta, must not reach it as though it were estimated.
  local_mocked_bindings(.hzr_weak_direction_impl = function(vcov, ...) {
    seen <<- vcov
    list(weak = NULL, reason = NA_character_)
  })
  fit <- suppressWarnings(hazard(
    time = d$time, status = d$status, dist = "multiphase", fit = TRUE,
    phases = list(late = hzr_phase("g3", tau = 8, gamma = 5, eta = 1,
                                   constraint = "alpha_gamma_eta")),
    control = list(conserve = FALSE, n_starts = 1L)
  ))
  expect_true(is.finite(vcov(fit)[4, 4]))
  expect_true(all(is.na(seen[4, ])))
  expect_true(all(is.finite(diag(seen)[-4])))

  tab <- summary(fit)$coefficients
  expect_true(is.finite(tab[4, "std_error"]))
  expect_true(is.na(tab[4, "z_stat"]) && is.na(tab[4, "p_value"]))
  expect_true(all(is.finite(tab[-4, "z_stat"])))
})

test_that("a translated FIXGAE2 job runs as a constrained fit", {
  skip_on_cran()
  f <- withr::local_tempfile(fileext = ".sas")
  writeLines(c(
    "%HAZARD( PROC HAZARD DATA=BUILT NOCONSERVE; EVENT EV; TIME IV;",
    "  PARMS MUL=0.1 TAU=8 ALPHA=2 GAMMA=5 ETA=1 FIXGAE2 WEIBULL; );"
  ), f)
  out <- withr::local_tempdir()
  job <- hzr_translate_sas(f, out_dir = out)
  expect_false("FIXGAE2" %in% job$untranslated$construct)
  qmd <- readLines(file.path(out, sub("[.]sas$", ".qmd", basename(f))))
  expect_true(any(grepl("constraint = \"alpha_gamma_eta\"", qmd, fixed = TRUE)))

  # Execute the emitted chunks, not their text.
  chunk <- function(label) {
    start <- which(qmd == paste("#| label:", label))
    end <- which(qmd == "```")
    end <- end[end > start][1L]
    paste(qmd[(start + 1L):(end - 1L)], collapse = "\n")
  }
  d <- constraint_data(c(log(0.3), log(10), 6, 2.4, 0.8))
  env <- new.env(parent = asNamespace("TemporalHazard"))
  env$BUILT <- data.frame(IV = d$time, EV = d$status)
  eval(parse(text = chunk("status")), env)
  # The emitted theta sits on the constraint, so nothing may be replaced;
  # other warnings (conditioning, starts) are not this test's business.
  withCallingHandlers(
    eval(parse(text = chunk("fit")), env),
    warning = function(w) {
      if (grepl("was replaced by", conditionMessage(w), fixed = TRUE)) {
        stop("translated theta was off its constraint: ", conditionMessage(w))
      }
      invokeRestart("muffleWarning")
    }
  )
  theta <- coef(env$fit)
  expect_true(env$fit$fit$converged)
  expect_equal(unname(theta[4]), unname(theta[3] * theta[5] / 2))
})

test_that("a theta off the constraint is replaced, with a warning, fitted or not", {
  phases <- list(late = hzr_phase("g3", tau = 8, gamma = 5, eta = 1,
                                  constraint = "alpha_gamma_eta"))
  d <- constraint_data(c(log(0.3), log(10), 6, 2.4, 0.8), n = 200)
  off <- c(log(0.2), log(8), 5, 9, 1)            # alpha should be 2.5

  expect_warning(
    unfitted <- hazard(time = d$time, status = d$status, dist = "multiphase",
                       phases = phases, theta = off, fit = FALSE),
    "was replaced by 2.5"
  )
  expect_equal(unname(coef(unfitted))[4], 2.5)

  # With a phase formula the slot cannot be located unfitted: say so.
  dat <- data.frame(time = d$time, status = d$status,
                    z = seq_along(d$time) %% 2)
  covaried <- list(late = hzr_phase("g3", tau = 8, gamma = 5, eta = 1,
                                    formula = ~ z,
                                    constraint = "alpha_gamma_eta"))
  expect_warning(
    hazard(time = time, status = status, data = dat, dist = "multiphase",
           phases = covaried, theta = c(off, 0), fit = FALSE),
    "applied only under fit = TRUE"
  )

  skip_on_cran()
  expect_warning(
    fitted <- hazard(time = d$time, status = d$status, dist = "multiphase",
                     phases = phases, theta = off, fit = TRUE,
                     control = list(conserve = FALSE, n_starts = 1L)),
    "was replaced by 2.5"
  )
  theta <- unname(coef(fitted))
  expect_equal(theta[4], theta[3] * theta[5] / 2)
})
