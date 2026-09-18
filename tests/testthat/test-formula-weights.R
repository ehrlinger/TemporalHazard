# On the formula interface, `weights = <name>` is looked up in `data` first,
# then the calling frame, as the vector interface has done since #151 and as
# stats::lm() does (#392). It used to skip `data`: a column-only name was not
# found, and a name bound both as a column and in the calling frame silently
# read the calling frame's vector.

# The weights column is named fw_wcol, a name no caller plausibly has: the
# column-only test asserts no ambiguity warning, which is correct only while
# no variable of that name is visible from the calling frame.
fw_data <- function() {
  set.seed(1)
  data.frame(t = rexp(40), s = rbinom(40, 1, 0.7), x = rnorm(40),
             fw_wcol = runif(40, 0.5, 1.5))
}

test_that("a weights name that is a column of data is read from data", {
  d <- fw_data()
  want <- hazard(survival::Surv(t, s) ~ x, data = d, weights = d$fw_wcol,
                 dist = "weibull", theta = c(1, 1, 0), fit = TRUE)
  expect_no_warning(
    got <- hazard(survival::Surv(t, s) ~ x, data = d, weights = fw_wcol,
                  dist = "weibull", theta = c(1, 1, 0), fit = TRUE)
  )
  expect_identical(got$fit$theta, want$fit$theta)
  expect_identical(got$fit$objective, want$fit$objective)
})

test_that("a name that is both a column and a calling-frame vector reads the column, and warns", {
  d <- fw_data()
  fw_wcol <- rep(1, 40) # disagrees with d$fw_wcol, so the reading is observable
  column <- hazard(survival::Surv(t, s) ~ x, data = d, weights = d$fw_wcol,
                   dist = "weibull", theta = c(1, 1, 0), fit = TRUE)
  frame <- hazard(survival::Surv(t, s) ~ x, data = d, weights = rep(1, 40),
                  dist = "weibull", theta = c(1, 1, 0), fit = TRUE)
  # The two candidates give different fits, so the check below can fail.
  expect_false(isTRUE(all.equal(column$fit$objective, frame$fit$objective)))
  expect_warning(
    got <- hazard(survival::Surv(t, s) ~ x, data = d, weights = fw_wcol,
                  dist = "weibull", theta = c(1, 1, 0), fit = TRUE),
    "'fw_wcol' \\(weights\\).*both a column of 'data'.*The column was used"
  )
  expect_identical(got$fit$objective, column$fit$objective)
  expect_identical(got$fit$theta, column$fit$theta)
})

test_that("a calling-frame name that is not a column still resolves there, without a warning", {
  d <- fw_data()
  w_local <- d$fw_wcol * 2
  want <- hazard(survival::Surv(t, s) ~ x, data = d, weights = d$fw_wcol * 2,
                 dist = "weibull", theta = c(1, 1, 0), fit = TRUE)
  expect_no_warning(
    got <- hazard(survival::Surv(t, s) ~ x, data = d, weights = w_local,
                  dist = "weibull", theta = c(1, 1, 0), fit = TRUE)
  )
  expect_identical(got$fit$objective, want$fit$objective)
})

test_that("both interfaces read the same vector and warn on the same ambiguous input", {
  # Interface parity is about behaviour, which is what the design decided:
  # the same column read and a warning raised, with one diagnosis. The
  # remedy clause differs on purpose: only the vector interface can drop
  # `data` to reach the calling frame.
  d <- fw_data()
  fw_wcol <- rep(1, 40)
  catch <- function(expr) {
    msg <- NULL
    val <- withCallingHandlers(expr, warning = function(w) {
      if (grepl("both a column of 'data'", conditionMessage(w))) {
        msg <<- conditionMessage(w)
      }
      invokeRestart("muffleWarning")
    })
    list(fit = val, msg = msg)
  }
  via_formula <- catch(hazard(survival::Surv(t, s) ~ x, data = d,
                              weights = fw_wcol, dist = "weibull",
                              theta = c(1, 1, 0), fit = TRUE))
  via_vector <- catch(hazard(data = d, time = t, status = s,
                             x = as.matrix(d["x"]), weights = fw_wcol,
                             dist = "weibull", theta = c(1, 1, 0), fit = TRUE))
  expect_false(is.null(via_formula$msg))
  expect_false(is.null(via_vector$msg))
  diagnosis <- function(m) sub("The column was used\\..*$", "", m)
  expect_identical(diagnosis(via_formula$msg), diagnosis(via_vector$msg))
  expect_identical(via_formula$fit$fit$objective, via_vector$fit$fit$objective)
  expect_match(via_vector$msg, "or omit 'data' to use the calling frame's value\\.$")
  expect_match(via_formula$msg, "not a column of 'data'\\.$")
})

test_that("a namespace-qualified call does not raise a false ambiguity warning", {
  # `::` names a namespace and an export; neither is looked up in data, so a
  # column that happens to share either name is not ambiguous (#401 review).
  # Found on the vector interface too, where it had shipped since #151.
  d <- fw_data()
  names(d)[names(d) == "fw_wcol"] <- "base"
  base <- "a caller variable named base"
  w_ok <- d$base
  want <- hazard(survival::Surv(t, s) ~ x, data = d, weights = w_ok,
                 dist = "weibull", theta = c(1, 1, 0), fit = TRUE)
  expect_no_warning(
    got <- hazard(survival::Surv(t, s) ~ x, data = d,
                  weights = base::abs(w_ok), dist = "weibull",
                  theta = c(1, 1, 0), fit = TRUE)
  )
  expect_identical(got$fit$objective, want$fit$objective)
  expect_no_warning(
    hazard(data = d, time = t, status = s, x = as.matrix(d["x"]),
           weights = base::abs(w_ok), dist = "weibull", theta = c(1, 1, 0),
           fit = TRUE)
  )
})

test_that("an ambiguous argument inside a namespace-qualified call still warns", {
  # The control: skipping the namespace and export names must not skip the
  # call's own arguments, which ARE looked up in data.
  d <- fw_data()
  fw_wcol <- rep(1, 40)
  base <- "a caller variable named base"
  d$base <- 1
  expect_warning(
    hazard(survival::Surv(t, s) ~ x, data = d, weights = base::abs(fw_wcol),
           dist = "weibull", theta = c(1, 1, 0), fit = TRUE),
    "'fw_wcol' \\(weights\\): the name is both a column"
  )
})
