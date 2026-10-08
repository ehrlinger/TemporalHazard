# The likelihood's data-defect guards raise a classed condition,
# `hzr_data_error`, so a caller that absorbs numerical failures can let a
# data defect through instead of relabelling it. The score path caught every
# error at six sites and reported it as "information matrix could not be
# inverted" (#407); narrowing those sites needs a class to narrow on.
# Messages are unchanged: only the class is new.

test_that("a left-censored row under objective = 'sas' is a data error (#407)", {
  expect_error(.hzr_check_sas_status(c(1, -1), "sas"),
               "does not support left-censored rows",
               class = "hzr_data_error")
})

test_that("an incomplete status is a data error (#407)", {
  expect_error(.hzr_check_sas_data(c(1, NA), c(1, 2), NULL, NULL,
                                   "likelihood"),
               "'status' must be complete", class = "hzr_data_error")
})

test_that("a bound of the wrong length is a data error (#407)", {
  expect_error(.hzr_check_sas_data(c(1, 2, 2), c(1, 2, 3), c(0.5, 1), NULL,
                                   "sas"),
               "time_lower has length 2, but status has length 3",
               class = "hzr_data_error")
})

test_that("an NA interval bound is a data error (#407)", {
  expect_error(.hzr_check_sas_data(c(1, 2), c(1, 2), c(0, NA), c(1, 2),
                                   "sas"),
               "requires both bounds on every interval-censored row",
               class = "hzr_data_error")
})

test_that("an interval with upper <= lower is a data error at entry (#407)", {
  expect_error(.hzr_check_sas_data(c(1, 2), c(1, 2), c(0, 3), c(1, 2),
                                   "sas"),
               "requires upper > lower on every", class = "hzr_data_error")
})

test_that("an interval with upper <= lower is a data error in the likelihood (#407)", {
  expect_error(.hzr_logl_interval(0.1, 0.2, lower = 2, upper = 1,
                                  weights = 1, objective = "sas"),
               "requires upper > lower on every", class = "hzr_data_error")
})

test_that("the class survives a full multiphase likelihood evaluation (#407)", {
  # The route #407 names: the score path evaluates .hzr_logl_multiphase(),
  # so the class must reach its caller, not only the guard's.
  ph <- .hzr_validate_phases(list(e = hzr_phase("cdf", t_half = 3, nu = 1,
                                                m = 0),
                                  c = hzr_phase("constant")))
  time <- c(1, 2, 3, 4)
  status <- c(1, 0, 2, 1)
  lo <- c(0, 0, 3, 0)
  expect_error(
    .hzr_logl_multiphase(c(log(0.1), log(3), 1, 0, log(0.05)), time, status,
                         time_lower = lo, time_upper = time,
                         phases = ph,
                         covariate_counts = c(e = 0L, c = 0L),
                         x_list = list(e = NULL, c = NULL),
                         objective = "sas"),
    "requires upper > lower on every", class = "hzr_data_error"
  )
})
