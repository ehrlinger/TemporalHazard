# The score path absorbs NUMERICAL failures on purpose: a Hessian it cannot
# invert, a non-finite gradient. It used to absorb EVERY error, so a data
# defect raised inside the likelihood came back as "the information matrix
# could not be inverted" and sent the reader looking at conditioning (#407).
#
# Both directions are tested here, because narrowing a catch can turn a
# handled condition into a user-visible error just as easily as it can fix
# a mislabelled one. The second test pins what must NOT change.

sp_fit_sas <- function(d) {
  suppressWarnings(hazard(
    time = d$int_dead, status = d$st, time_lower = d$tl, time_upper = d$tu,
    data = d, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    objective = "sas", fit = TRUE,
    control = list(n_starts = 2L, maxit = 300L)
  ))
}

sp_data <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc[, c("int_dead", "dead", "age", "mal")])
  d$tl <- d$int_dead
  d$tu <- d$int_dead + 1
  d$st <- 2L # every row interval-censored, so the sas bound guards are live
  d
}

test_that("a data defect reaching the score path is not reported as numerical (#407)", {
  skip_on_cran() # a multiphase fit plus a screen
  d <- sp_data()
  fit <- sp_fit_sas(d)
  # Known positive: the guard really does fire on this data, as a CLASSED
  # error, when the likelihood is called directly. Without this the test
  # below could pass because nothing was ever wrong.
  bad <- fit
  bad$data$time_lower[which(bad$data$status == 2L)[1L]] <- NA_real_
  expect_error(
    .hzr_logl_multiphase(
      theta = bad$fit$theta, time = bad$data$time, status = bad$data$status,
      time_lower = bad$data$time_lower, time_upper = bad$data$time_upper,
      x = bad$data$x, weights = bad$data$weights, phases = bad$spec$phases,
      covariate_counts = bad$fit$covariate_counts, x_list = bad$fit$x_list,
      objective = "sas"
    ),
    class = "hzr_data_error"
  )
  # The screen must surface that, not relabel it. Before #407 it returned
  # normally with "no remaining candidate could be tested for entry".
  expect_error(
    suppressWarnings(hzr_stepwise(
      bad, data = d, scope = list(early = ~ mal, constant = NULL),
      direction = "forward", criterion = "score", slentry = 0.99,
      trace = FALSE
    )),
    class = "hzr_data_error"
  )
})

test_that("a numerical failure is still absorbed by the score path (#407)", {
  skip_on_cran() # a multiphase fit
  # The half that must NOT change, and it needs an input that genuinely
  # ERRORS inside the wrapped expression. A first version used a constant
  # candidate column and passed without ever reaching the catch, so it could
  # not fail: a mutant that re-raised EVERY error survived it.
  d <- sp_data()
  fit <- sp_fit_sas(d)
  zero <- rep(0, length(fit$fit$theta))
  # Known positive: this theta really does raise, and it is NOT a data error.
  err <- tryCatch(
    numDeriv::hessian(
      function(par) {
        -.hzr_logl_multiphase(
          par, time = fit$data$time, status = fit$data$status,
          time_lower = fit$data$time_lower, time_upper = fit$data$time_upper,
          x = fit$data$x, weights = fit$data$weights,
          phases = fit$spec$phases,
          covariate_counts = fit$fit$covariate_counts,
          x_list = fit$fit$x_list, objective = "sas"
        )
      }, zero
    ),
    error = function(e) e
  )
  expect_s3_class(err, "error")
  expect_false(inherits(err, "hzr_data_error"))
  # The score path must swallow it and report NULL, as it always has.
  expect_null(.hzr_score_multiphase_hessian(fit, theta = zero))
  expect_null(.hzr_score_information(fit, theta = zero))
})
