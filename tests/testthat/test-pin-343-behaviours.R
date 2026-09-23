# Three behaviours that are correct and were untested (#417).
#
# Found while verifying #343 item by item: none of them is a silent wrong
# answer, and none was pinned. An untested correct behaviour regresses
# without anyone noticing, which is most of this package's history, so each
# one below is pinned with the input that exercises it. Measured unchanged on
# 71277ff8, after #425 and #427 moved the same paths.

pin_data_417 <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  stats::na.omit(avc[, c("int_dead", "dead", "age", "mal", "com_iv")])
}

# A refit that returns `value` instead of a fit, as #343 items 1 and 2 do.
pin_refit_417 <- function(fit, value) {
  env <- new.env(parent = fit$call_env %||% globalenv())
  assign("odd_refit", function(...) value, envir = env)
  fit$call[[1L]] <- as.name("odd_refit")
  fit$call_env <- env
  fit
}

test_that("select mode counts a refit with no estimates as failed (#417)", {
  skip_on_cran() # bootstrap replicates with a screen each
  # #343 item 2. Refit mode is pinned in test-bootstrap-vector-interface.R
  # with the specific reason "refit returned no parameter estimates". Select
  # mode was not pinned at all. It is correct -- the run completes and
  # nothing counts as a success -- but its reason is the generic one from
  # hzr_stepwise()'s input check, so this records what it is rather than
  # what it might become (issue #417 item 1(b), for the string pass).
  d <- pin_data_417()
  base <- hazard(survival::Surv(int_dead, dead) ~ 1, data = d,
                 dist = "weibull", theta = c(0.1, 1), fit = TRUE)
  odd <- pin_refit_417(base, list(fit = list(objective = 1, theta = NULL)))
  bs <- suppressWarnings(hzr_bootstrap(odd, n_boot = 3L, seed = 1L,
                                       scope = ~ age, criterion = "wald"))
  expect_identical(bs$n_success, 0L)
  expect_identical(bs$n_failed, 3L)
  expect_identical(sum(bs$failure_reasons), 3L)
  expect_identical(names(bs$failure_reasons), "`fit` must be a `hazard` object.")
  # Known positive for the fixture: the same base bootstraps normally.
  ok <- suppressWarnings(hzr_bootstrap(base, n_boot = 3L, seed = 1L,
                                       scope = ~ age, criterion = "wald"))
  expect_identical(ok$n_success, 3L)
})

test_that("a current fit with `data$frame` removed is not refused (#417)", {
  skip_on_cran() # a multiphase fit
  # #343 item 6. Removing the stored frame, for example to drop data before
  # saving, used to make a fit look pre-1.1.0. Since the fit records whether
  # each phase used its formula (attr(x_list, "from_formula")), that record
  # decides, and a fit made by this version is NOT refused. Only the
  # pre-record case is pinned elsewhere, by a helper that strips both.
  d <- pin_data_417()[, c("int_dead", "dead", "age", "mal")]
  f <- suppressWarnings(hazard(
    time = d$int_dead, status = d$dead, x = cbind(mal = d$mal), data = d,
    dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1, fixed = "m",
                        formula = ~ mal),
      constant = hzr_phase("constant", formula = ~ age)
    ),
    fit = TRUE, control = list(n_starts = 1L, conserve = FALSE)
  ))
  # The mechanism: the fit carries the record that decides.
  expect_false(is.null(attr(f$fit$x_list, "from_formula")))
  f$data["frame"] <- NULL
  expect_null(.hzr_ignored_phase_formula(f))
  # Known positive: strip the record too, as a fit saved before it existed,
  # and the refusal fires.
  g <- f
  attr(g$fit$x_list, "from_formula") <- NULL
  expect_match(.hzr_ignored_phase_formula(g), "phase 'early' has a formula",
               fixed = TRUE)
})

test_that("the drop refusal names the reduced design on both paths (#417)", {
  skip_on_cran() # two fits and two backward steps
  # #343 item 3. A drop that removes no column is refused with "the reduced
  # design (...)", not "the refit's design", because on the single-
  # distribution path that design is built BEFORE any refit. Both paths
  # share the string deliberately; only the single-distribution one was
  # pinned, so a change to the multiphase path could have diverged silently.
  set.seed(5)
  n <- 400
  dd <- data.frame(z = stats::rnorm(n),
                   f = factor(sample(c("a", "b"), n, TRUE)))
  dd$time <- stats::rexp(n) * exp(-0.5 * dd$z * (dd$f == "b"))
  dd$status <- 1L
  single <- suppressWarnings(hazard(
    survival::Surv(time, status) ~ z + z:f, data = dd, dist = "weibull",
    theta = c(0.5, 1, 0, 0), fit = TRUE
  ))
  multi <- suppressWarnings(hazard(
    survival::Surv(time, status) ~ 1, data = dd, dist = "multiphase",
    phases = list(constant = hzr_phase("constant", formula = ~ z + z:f)),
    fit = TRUE
  ))
  step_of <- function(fit) {
    suppressWarnings(.hzr_stepwise_backward_step(fit, data = dd,
                                                 criterion = "wald",
                                                 slstay = 0.2))
  }
  s_reason <- step_of(single)$refit_failure_reasons[["z"]]
  m_reason <- step_of(multi)$refit_failure_reasons[[1L]]
  for (reason in list(s_reason, m_reason)) {
    expect_match(reason, "removes no column: the reduced design (",
                 fixed = TRUE)
    expect_no_match(reason, "refit's design", fixed = TRUE)
  }
  # The two paths share the string: same wording, same named columns.
  expect_identical(m_reason, s_reason)
})
