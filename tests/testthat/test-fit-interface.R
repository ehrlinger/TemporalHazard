# hazard() records the interface it used in spec$interface (#432). The stored
# call cannot answer the question afterwards: a wrapper writes
# `formula = fml` whether `fml` held a formula or NULL, and a formula given
# beside `time =` is the one used.

fi_data_432 <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  stats::na.omit(avc[, c("int_dead", "dead", "age", "mal", "com_iv")])
}

fi_msg_432 <- function(expr) tryCatch(expr, error = conditionMessage)

# A wrapper forwarding a NULL formula by variable, with a design matrix.
fi_wrap_x_432 <- function(dat) {
  fml <- NULL
  hazard(formula = fml, time = dat$int_dead, status = dat$dead,
         x = cbind(age = dat$age), data = dat, dist = "weibull",
         theta = c(0.1, 1, 0), fit = TRUE)
}

test_that("hazard() records the interface each route took (#432)", {
  d <- fi_data_432()
  f <- survival::Surv(int_dead, dead) ~ age
  by_formula <- hazard(f, data = d, dist = "weibull", theta = c(0.1, 1, 0))
  by_vector <- hazard(time = d$int_dead, status = d$dead, dist = "weibull",
                      theta = c(0.1, 1))
  # A formula beside `time =` is the one used.
  beside <- hazard(formula = f, time = d$int_dead, data = d,
                   dist = "weibull", theta = c(0.1, 1, 0))
  expect_identical(by_formula$spec$interface, "formula")
  expect_identical(by_vector$spec$interface, "vector")
  expect_identical(beside$spec$interface, "formula")
  expect_identical(fi_wrap_x_432(d)$spec$interface, "vector")
})

test_that("the recorded interface decides; legacy only without it (#432)", {
  d <- fi_data_432()
  fit <- hazard(time = d$int_dead, status = d$dead, dist = "weibull",
                theta = c(0.1, 1))
  # The legacy answer is not evaluated when the field is present.
  expect_false(.hzr_is_formula_fit(fit, stop("legacy was evaluated")))
  fit$spec$interface <- NULL
  expect_true(.hzr_is_formula_fit(fit, TRUE))
  expect_false(.hzr_is_formula_fit(fit, FALSE))
})

test_that("a vector fit's formula name rebound later is not adopted (#432)", {
  d <- fi_data_432()
  plain <- hazard(time = d$int_dead, status = d$dead,
                  x = cbind(age = d$age), data = d, dist = "weibull",
                  theta = c(0.1, 1, 0), fit = TRUE)
  want <- fi_msg_432(hzr_bootstrap(plain, n_boot = 3L, seed = 1L))
  # Known positive: the plain fit is refused, so matching it is a refusal.
  expect_match(want, "passed directly as `x`", fixed = TRUE)

  wf <- fi_wrap_x_432(d)
  expect_identical(wf$call$formula, quote(fml))
  # At top level the name stays live; rebinding it where the call is
  # evaluated is that later assignment. Same response, other covariate.
  expect_true(is.environment(wf$call_env))
  assign("fml", survival::Surv(int_dead, dead) ~ mal, envir = wf$call_env)

  expect_identical(fi_msg_432(hzr_bootstrap(wf, n_boot = 3L, seed = 1L)),
                   want)
  expect_null(.hzr_stored_formula(wf))
  expect_match(
    fi_msg_432(hzr_stepwise(wf, scope = ~ com_iv, data = d,
                            direction = "forward", criterion = "wald",
                            trace = FALSE)),
    "built via the vector interface", fixed = TRUE
  )
})

test_that("a formula passed by variable beside `time =` screens (#432)", {
  skip_on_cran() # two forward screens
  d <- fi_data_432()
  f <- survival::Surv(int_dead, dead) ~ age
  screen <- function(fit) {
    sw <- suppressWarnings(hzr_stepwise(
      fit, scope = ~ age + mal, data = d, direction = "forward",
      criterion = "wald", slentry = 0.99, trace = FALSE
    ))
    list(steps = sw$steps[, c("action", "variable")], theta = sw$fit$theta)
  }
  plain <- hazard(f, data = d, dist = "weibull", theta = c(0.1, 1, 0),
                  fit = TRUE)
  beside <- hazard(formula = f, time = d$int_dead, data = d,
                   dist = "weibull", theta = c(0.1, 1, 0), fit = TRUE)
  want <- screen(plain)
  # Known positive: the screen takes a step, so a refit ran.
  expect_identical(nrow(want$steps), 1L)
  expect_identical(screen(beside), want)
})

test_that("predict() reads a wrapper fit's newdata as a vector fit's (#432)", {
  d <- fi_data_432()
  plain <- hazard(time = d$int_dead, status = d$dead,
                  x = cbind(age = d$age), data = d, dist = "weibull",
                  theta = c(0.1, 1, 0), fit = TRUE)
  wf <- fi_wrap_x_432(d)
  assign("fml", survival::Surv(int_dead, dead) ~ mal, envir = wf$call_env)
  # newdata carrying more than the design column, as a data frame row does.
  nd <- transform(d[1:3, ], time = 1)
  want <- predict(plain, newdata = nd, type = "survival")
  expect_length(want, 3L)
  expect_identical(predict(wf, newdata = nd, type = "survival"), want)
})

test_that("a multiphase wrapper fit screens as the plain fit does (#432)", {
  skip_on_cran() # multiphase screens
  d <- fi_data_432()
  ph <- function() {
    list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes", formula = ~ age),
         constant = hzr_phase("constant"))
  }
  ctl <- list(n_starts = 1L, maxit = 300L)
  plain <- suppressWarnings(hazard(
    time = d$int_dead, status = d$dead, data = d, dist = "multiphase",
    phases = ph(), fit = TRUE, control = ctl
  ))
  wrap_m <- function(dat) {
    fml <- NULL
    suppressWarnings(hazard(
      formula = fml, time = dat$int_dead, status = dat$dead, data = dat,
      dist = "multiphase", phases = ph(), fit = TRUE, control = ctl
    ))
  }
  wf <- wrap_m(d)
  assign("fml", survival::Surv(int_dead, dead) ~ com_iv, envir = wf$call_env)
  screen <- function(f, direction) {
    sw <- suppressWarnings(hzr_stepwise(
      f, scope = list(constant = ~ mal), data = d, direction = direction,
      criterion = "wald", slentry = 0.99, slstay = 1e-12, trace = FALSE
    ))
    list(steps = sw$steps[, c("action", "variable", "phase")],
         theta = sw$fit$theta)
  }
  want <- screen(plain, "both")
  # Known positive: the screen enters `mal` and then drops a term, so both
  # an entry refit and a drop refit ran.
  expect_identical(want$steps$action[1], "enter")
  expect_true("drop" %in% want$steps$action)
  expect_identical(screen(wf, "both"), want)
})

test_that("a wrapper fit without `x` resamples as the plain fit does (#432)", {
  d <- fi_data_432()
  plain <- hazard(time = d$int_dead, status = d$dead, data = d,
                  dist = "weibull", theta = c(0.1, 1), fit = TRUE)
  want <- hzr_bootstrap(plain, n_boot = 3L, seed = 1L)
  # Known positive: real replicates, so identical ones are evidence.
  expect_true(all(want$summary$sd > 0))
  wrap <- function(dat) {
    fml <- NULL
    hazard(formula = fml, time = dat$int_dead, status = dat$dead,
           data = dat, dist = "weibull", theta = c(0.1, 1), fit = TRUE)
  }
  wf <- wrap(d)
  # Each replicate re-evaluates the call; `fml` must not be read again.
  assign("fml", survival::Surv(int_dead, dead) ~ mal, envir = wf$call_env)
  got <- hzr_bootstrap(wf, n_boot = 3L, seed = 1L)
  expect_identical(got$replicates, want$replicates)
})

test_that("a multiphase wrapper fit is checked as the plain fit is (#432)", {
  skip_on_cran() # multiphase fits
  d <- fi_data_432()
  ph <- function(fm) {
    list(early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes", formula = fm),
         constant = hzr_phase("constant"))
  }
  ctl <- list(n_starts = 1L, maxit = 300L)
  fit_m <- function(dat, ...) {
    suppressWarnings(hazard(time = dat$int_dead, status = dat$dead,
                            data = dat, dist = "multiphase", fit = TRUE,
                            control = ctl, ...))
  }
  wrap_m <- function(dat, ...) {
    fml <- NULL
    suppressWarnings(hazard(formula = fml, time = dat$int_dead,
                            status = dat$dead, data = dat,
                            dist = "multiphase", fit = TRUE, control = ctl,
                            ...))
  }
  screen <- function(f, dd) {
    fi_msg_432(suppressWarnings(hzr_stepwise(
      f, scope = list(constant = ~ mal), data = dd, direction = "forward",
      criterion = "wald", slentry = 0.99, trace = FALSE
    )))
  }
  # Rows reordered: the vector fit's stored vectors are read by position.
  set.seed(3)
  shuffled <- d[sample(nrow(d)), ]
  want <- screen(fit_m(d, phases = ph(~ age)), shuffled)
  expect_match(want, "does not hold the rows the model was fitted on",
               fixed = TRUE)
  expect_identical(screen(wrap_m(d, phases = ph(~ age)), shuffled), want)

  # Phases inheriting a design matrix passed as `x`.
  x <- cbind(age = d$age)
  want <- screen(fit_m(d, x = x, phases = ph(NULL)), d)
  expect_match(want, "inherit a design matrix passed directly as `x`",
               fixed = TRUE)
  expect_identical(screen(wrap_m(d, x = x, phases = ph(NULL)), d), want)
})

test_that("a formula fit's status column is its Surv() event (#432)", {
  d <- fi_data_432()
  # `status =` beside a formula is not what the fit used.
  beside <- hazard(formula = survival::Surv(int_dead, dead) ~ age,
                   status = d$com_iv, data = d, dist = "weibull",
                   theta = c(0.1, 1, 0))
  expect_identical(.hzr_fit_row_inputs(beside)$status$column, "dead")
  vec <- hazard(time = d$int_dead, status = d$dead, dist = "weibull",
                theta = c(0.1, 1))
  expect_identical(.hzr_fit_row_inputs(vec)$status$column, "dead")
})
