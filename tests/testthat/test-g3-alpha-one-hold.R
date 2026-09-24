# test-g3-alpha-one-hold.R
# A g3 phase with alpha FIXED at 1 is re-expressed as PROC HAZARD's
# SETG3_ignore_tau() does (#415), and FIXGE2's alpha rules follow SETG3 (#418).
# hazard/src/model/setg3.c:313-315, 380-425 and 843-850.

# Weibull data: at alpha = 1, G3 = (t / tau)^(gamma * eta), so the g3 phase IS
# a Weibull cumulative hazard mu * t^p once tau is held at 1.
hold_data <- function(n = 300, shape = 2, scale_mu = 0.5, seed = 415) {
  withr::local_seed(seed)
  ev <- (stats::rexp(n) / scale_mu)^(1 / shape)
  cen <- stats::runif(n, 0, 3)
  list(time = pmin(ev, cen), status = as.numeric(ev <= cen))
}

fit_g3 <- function(d, phase, theta = NULL) {
  hazard(time = d$time, status = d$status, dist = "multiphase", fit = TRUE,
         phases = list(late = phase), theta = theta,
         control = list(n_starts = 1L))
}

boundary_records <- function(fit, mechanism) {
  Filter(function(r) identical(r$mechanism, mechanism), fit$fit$boundary)
}

test_that("alpha fixed at 1 holds tau and eta, and fits the Weibull", {
  skip_if_not_installed("numDeriv")
  d <- hold_data()
  w <- NULL
  fit <- withCallingHandlers(
    fit_g3(d, hzr_phase("g3", tau = 2, gamma = 3, alpha = 1, eta = 1,
                        fixed = "alpha")),
    hzr_g3_alpha_one = function(e) {
      w <<- e
      invokeRestart("muffleWarning")
    }
  )
  expect_s3_class(w, "hzr_boundary")
  expect_match(conditionMessage(w), "held at setup", fixed = TRUE)

  th <- coef(fit)
  # The held values are the ones SAS uses, and they stay there.
  expect_equal(unname(th[c("late.log_tau", "late.alpha", "late.eta")]),
               c(0, 1, 1))
  expect_equal(unname(fit$fit$fixed_mask), c(FALSE, TRUE, FALSE, TRUE, TRUE))

  # A known positive the package did not compute through g3: the held fit is
  # the Weibull MLE, so gamma-hat is its shape and the likelihoods agree.
  wb <- survival::survreg(survival::Surv(d$time, d$status) ~ 1,
                          dist = "weibull")
  expect_equal(unname(th[["late.gamma"]]), 1 / wb$scale, tolerance = 1e-4)
  expect_equal(fit$fit$objective, as.numeric(stats::logLik(wb)),
               tolerance = 1e-6)

  rec <- boundary_records(fit, "g3_alpha_one")
  expect_length(rec, 1L)
  expect_identical(rec[[1]]$phase, "late")
  expect_setequal(rec[[1]]$parameter, c("tau", "eta"))
  # The unheld fit leaves this ridge to the detector; the held one has none.
  expect_null(fit$fit$weak)
})

test_that("the hold keeps the likelihood: held and unheld maxima agree", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  d <- hold_data(seed = 7)
  held <- suppressWarnings(
    fit_g3(d, hzr_phase("g3", tau = 2, gamma = 3, alpha = 1, eta = 1,
                        fixed = "alpha")))
  # tau and eta pinned by the caller at the values the hold would choose: no
  # hold fires, so this is the same model reached without it.
  pinned <- fit_g3(d, hzr_phase("g3", tau = 1, gamma = 3, alpha = 1, eta = 1,
                                fixed = c("tau", "alpha", "eta")))
  expect_null(boundary_records(pinned, "g3_alpha_one")[[1]])
  expect_equal(held$fit$objective, pinned$fit$objective, tolerance = 1e-8)
  expect_equal(coef(held), coef(pinned), tolerance = 1e-5)
})

test_that("the hold follows SETG3's branches for which of gamma/eta carries", {
  hold <- function(fixed, gamma = 3, eta = 0.5, constraint = "none") {
    args <- list("g3", tau = 2, gamma = gamma, alpha = 1, fixed = fixed,
                 constraint = constraint)
    if (constraint == "none") args$eta <- eta
    ph <- do.call(hzr_phase, args)
    th <- c(late.log_mu = 0, late.log_tau = log(2), late.gamma = ph$gamma,
            late.alpha = 1, late.eta = ph$eta)
    .hzr_g3_alpha_one_hold(th, list(late = ph))
  }
  # Both free: eta is fixed, gamma carries the product (setg3.c:406-413).
  h <- hold("alpha")
  expect_equal(unname(h$theta[c("late.gamma", "late.eta")]), c(1.5, 1))
  expect_setequal(h$phases$late$fixed, c("alpha", "tau", "eta"))
  # Eta fixed by the caller: the same branch.
  h <- hold(c("alpha", "eta"))
  expect_equal(unname(h$theta[c("late.gamma", "late.eta")]), c(1.5, 1))
  # Gamma fixed, eta free: eta carries it, gamma = 1 (setg3.c:414-418).
  h <- hold(c("alpha", "gamma"))
  expect_equal(unname(h$theta[c("late.gamma", "late.eta")]), c(1, 1.5))
  expect_false("eta" %in% h$phases$late$fixed)
  # FIXGE2: gamma = 2, eta = 1, both held (setg3.c:393-403) ...
  h <- hold("alpha", constraint = "eta_gamma")
  expect_equal(unname(h$theta[c("late.gamma", "late.eta")]), c(2, 1))
  expect_true("gamma" %in% h$phases$late$fixed)
  # ... or eta = 2, gamma = 1 when ETA = 2 was given, i.e. gamma = 1 here.
  h <- hold("alpha", gamma = 1, constraint = "eta_gamma")
  expect_equal(unname(h$theta[c("late.gamma", "late.eta")]), c(1, 2))
  expect_equal(unname(h$theta[["late.log_tau"]]), 0)
})

test_that("the hold keys on alpha FIXED at exactly 1, as SAS does", {
  ph_free <- hzr_phase("g3", tau = 2, gamma = 3, alpha = 1, eta = 0.5)
  ph_off <- hzr_phase("g3", tau = 2, gamma = 3, alpha = 0.9, eta = 0.5,
                      fixed = "alpha")
  th <- function(a) {
    c(late.log_mu = 0, late.log_tau = log(2), late.gamma = 3,
      late.alpha = a, late.eta = 0.5)
  }
  expect_length(.hzr_g3_alpha_one_hold(th(1), list(late = ph_free))$records,
                0L)
  expect_length(.hzr_g3_alpha_one_hold(th(0.9), list(late = ph_off))$records,
                0L)
  # Supplied theta decides, not the phase's starting value: a phase built at
  # alpha = 0.9 and fitted from alpha = 1 is held.
  expect_length(.hzr_g3_alpha_one_hold(th(1), list(late = ph_off))$records,
                1L)
})

test_that("a hold that changes nothing is not announced", {
  ph <- hzr_phase("g3", tau = 1, gamma = 3, alpha = 1, eta = 1,
                  fixed = "shapes")
  th <- c(late.log_mu = 0, late.log_tau = 0, late.gamma = 3, late.alpha = 1,
          late.eta = 1)
  expect_length(.hzr_g3_alpha_one_hold(th, list(late = ph))$records, 0L)
})

test_that("FIXGE2 refuses a fixed alpha above 1 (SETG31040)", {
  d <- hold_data()
  for (a in c(1.0001, 2)) {
    expect_error(
      fit_g3(d, hzr_phase("g3", tau = 2, gamma = 3, alpha = a,
                          fixed = "alpha", constraint = "eta_gamma")),
      "SETG31040"
    )
  }
  # Below 1 is SAS's accepted region and fits.
  fit <- suppressWarnings(
    fit_g3(d, hzr_phase("g3", tau = 2, gamma = 3, alpha = 0.5,
                        fixed = "alpha", constraint = "eta_gamma")))
  expect_null(boundary_records(fit, "g3_fixge2_alpha_start")[[1]])
})

test_that("FIXGE2 moves a free alpha start at 1 or above to 2/3, and says so", {
  ph <- function(a) {
    hzr_phase("g3", tau = 2, gamma = 3, alpha = a, constraint = "eta_gamma")
  }
  th <- function(a) {
    c(late.log_mu = 0, late.log_tau = log(2), late.gamma = 3,
      late.alpha = a, late.eta = 2 / 3)
  }
  for (a in c(1, 2.5)) {
    out <- .hzr_g3_fixge2_alpha(th(a), list(late = ph(a)))
    expect_equal(unname(out$theta[["late.alpha"]]), 2 / 3)
    expect_identical(out$records[[1]]$mechanism, "g3_fixge2_alpha_start")
  }
  out <- .hzr_g3_fixge2_alpha(th(0.99), list(late = ph(0.99)))
  expect_equal(unname(out$theta[["late.alpha"]]), 0.99)
  expect_length(out$records, 0L)
  # Without FIXGE2 nothing moves.
  none <- hzr_phase("g3", tau = 2, gamma = 3, alpha = 2, eta = 1)
  expect_length(.hzr_g3_fixge2_alpha(th(2), list(late = none))$records, 0L)

  # End to end: the warning is classed and the record reaches the fit.
  skip_if_not_installed("numDeriv")
  d <- hold_data()
  classes <- character(0)
  fit <- withCallingHandlers(
    fit_g3(d, ph(1)),
    warning = function(w) {
      classes <<- c(classes, class(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true("hzr_g3_fixge2_alpha_start" %in% classes)
  expect_length(boundary_records(fit, "g3_fixge2_alpha_start"), 1L)
})

test_that("each boundary record keeps its own mechanism's lead-in", {
  recs <- list(list(mechanism = "g3_alpha_one", detail = "a."),
               list(mechanism = "unbounded_phase", detail = "b."))
  msg <- .hzr_boundary_message(recs)
  expect_match(msg, "^shape parameters held at setup: a\\.")
  expect_match(msg, "fitted outside the observed support: b.", fixed = TRUE)
})
