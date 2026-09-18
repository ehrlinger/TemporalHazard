# On the formula interface, `weights = <name>` is looked up in `data` first,
# then the calling frame, as the vector interface has done since #151 and as
# stats::lm() does (#392). It used to skip `data`: a column-only name was not
# found, and a name bound both as a column and in the calling frame silently
# read the calling frame's vector.

fw_data <- function() {
  set.seed(1)
  data.frame(t = rexp(40), s = rbinom(40, 1, 0.7), x = rnorm(40),
             wc = runif(40, 0.5, 1.5))
}

test_that("a weights name that is a column of data is read from data", {
  d <- fw_data()
  want <- hazard(survival::Surv(t, s) ~ x, data = d, weights = d$wc,
                 dist = "weibull", theta = c(1, 1, 0), fit = TRUE)
  expect_no_warning(
    got <- hazard(survival::Surv(t, s) ~ x, data = d, weights = wc,
                  dist = "weibull", theta = c(1, 1, 0), fit = TRUE)
  )
  expect_identical(got$fit$theta, want$fit$theta)
  expect_identical(got$fit$objective, want$fit$objective)
})

test_that("a name that is both a column and a calling-frame vector reads the column, and warns", {
  d <- fw_data()
  wc <- rep(1, 40) # disagrees with d$wc, so the reading is observable
  column <- hazard(survival::Surv(t, s) ~ x, data = d, weights = d$wc,
                   dist = "weibull", theta = c(1, 1, 0), fit = TRUE)
  frame <- hazard(survival::Surv(t, s) ~ x, data = d, weights = rep(1, 40),
                  dist = "weibull", theta = c(1, 1, 0), fit = TRUE)
  # The two candidates give different fits, so the check below can fail.
  expect_false(isTRUE(all.equal(column$fit$objective, frame$fit$objective)))
  expect_warning(
    got <- hazard(survival::Surv(t, s) ~ x, data = d, weights = wc,
                  dist = "weibull", theta = c(1, 1, 0), fit = TRUE),
    "'wc' \\(weights\\).*both a column of 'data'.*The column was used"
  )
  expect_identical(got$fit$objective, column$fit$objective)
  expect_identical(got$fit$theta, column$fit$theta)
})

test_that("a calling-frame name that is not a column still resolves there, without a warning", {
  d <- fw_data()
  w_local <- d$wc * 2
  want <- hazard(survival::Surv(t, s) ~ x, data = d, weights = d$wc * 2,
                 dist = "weibull", theta = c(1, 1, 0), fit = TRUE)
  expect_no_warning(
    got <- hazard(survival::Surv(t, s) ~ x, data = d, weights = w_local,
                  dist = "weibull", theta = c(1, 1, 0), fit = TRUE)
  )
  expect_identical(got$fit$objective, want$fit$objective)
})

test_that("both interfaces give the same warning on the same ambiguous input", {
  # Interface parity is why the formula path now shares the vector path's
  # rule, so pin that the message is one message, not two that drift.
  d <- fw_data()
  wc <- rep(1, 40)
  catch <- function(expr) {
    msg <- NULL
    withCallingHandlers(expr, warning = function(w) {
      if (grepl("both a column of 'data'", conditionMessage(w))) {
        msg <<- conditionMessage(w)
      }
      invokeRestart("muffleWarning")
    })
    msg
  }
  via_formula <- catch(hazard(survival::Surv(t, s) ~ x, data = d,
                              weights = wc, dist = "weibull",
                              theta = c(1, 1, 0), fit = TRUE))
  via_vector <- catch(hazard(data = d, time = t, status = s, x = as.matrix(d["x"]),
                             weights = wc, dist = "weibull",
                             theta = c(1, 1, 0), fit = TRUE))
  expect_false(is.null(via_formula))
  expect_identical(via_formula, via_vector)
})
