# A single-distribution refit warm-starts from the base fit's theta. It used
# to append a zero for an added term, but terms() orders an interaction after
# every main effect, and hazard() keeps a named theta's names, so the new
# coefficient was reported under the interaction's name (#489). Each refit
# here is compared, names and values, with an independent fit of the model
# the refit ends at.

refit_489_data <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc[, c("int_dead", "dead", "age", "mal", "com_iv")])
  d$g <- factor(cut(d$age, c(-Inf, 6, 24, Inf), labels = c("lo", "mid", "hi")))
  d
}

refit_489_fit <- function(rhs, theta, data, ...) {
  f <- stats::as.formula(paste("survival::Surv(int_dead, dead) ~", rhs))
  hazard(f, data = data, dist = "weibull", theta = theta, fit = TRUE, ...)
}

# Names identical, values within the optimizer's tolerance. The mislabelled
# coefficient differed from the one bearing its name by a factor of ~600.
expect_same_fit_489 <- function(refit, independent) {
  expect_true(isTRUE(refit$fit$converged))
  expect_true(isTRUE(independent$fit$converged))
  expect_identical(names(coef(refit)), names(coef(independent)))
  expect_equal(unname(coef(refit)), unname(coef(independent)),
               tolerance = 1e-3)
}

test_that("an added term under an interaction keeps its own name (#489)", {
  d <- refit_489_data()
  base <- refit_489_fit(
    "age * mal",
    c(mu = 0.1, nu = 1, beta_age = 0, beta_mal = 0, beta_agexmal = 0), d
  )
  refit <- .hzr_refit_with_scope(base, action = "add", var = "com_iv",
                                 data = d)
  independent <- refit_489_fit(
    "age * mal + com_iv",
    c(mu = 0.1, nu = 1, beta_age = 0, beta_mal = 0, com_iv = 0,
      beta_agexmal = 0), d
  )
  expect_same_fit_489(refit, independent)
  # The two coefficients the defect swapped, read by name.
  expect_equal(coef(refit)[["com_iv"]], coef(independent)[["com_iv"]],
               tolerance = 1e-3)
  expect_equal(coef(refit)[["beta_agexmal"]],
               coef(independent)[["beta_agexmal"]], tolerance = 1e-3)
})

test_that("a stepwise entry under an interaction reports it by name (#489)", {
  d <- refit_489_data()
  base <- refit_489_fit(
    "age * mal",
    c(mu = 0.1, nu = 1, beta_age = 0, beta_mal = 0, beta_agexmal = 0), d
  )
  sw <- hzr_stepwise(base, scope = "com_iv", data = d,
                     direction = "forward", slentry = 0.99, trace = FALSE)
  expect_identical(sw$steps$variable, "com_iv")
  independent <- refit_489_fit(
    "age * mal + com_iv",
    c(mu = 0.1, nu = 1, beta_age = 0, beta_mal = 0, com_iv = 0,
      beta_agexmal = 0), d
  )
  # The stepwise result is itself the final fit.
  expect_false(is.null(sw$fit$theta))
  expect_same_fit_489(sw, independent)
})

test_that("dropping a term after a factor removes its own coefficient (#489)", {
  # The drop used the term's position among the terms, but a factor before it
  # spans several columns: dropping `mal` removed `ghi`'s slot and carried
  # `beta_mal` onto the `ghi` column.
  d <- refit_489_data()
  base <- refit_489_fit(
    "age + g + mal",
    c(mu = 0.1, nu = 1, beta_age = 0, gmid = 0, ghi = 0, beta_mal = 0), d
  )
  refit <- .hzr_refit_with_scope(base, action = "drop", var = "mal", data = d)
  independent <- refit_489_fit(
    "age + g",
    c(mu = 0.1, nu = 1, beta_age = 0, gmid = 0, ghi = 0), d
  )
  expect_same_fit_489(refit, independent)
})

test_that("a factor enters and leaves with one coefficient per column (#489)", {
  d <- refit_489_data()
  base <- refit_489_fit("age + mal",
                        c(mu = 0.1, nu = 1, beta_age = 0, beta_mal = 0), d)
  added <- .hzr_refit_with_scope(base, action = "add", var = "g", data = d)
  expect_same_fit_489(added, refit_489_fit(
    "age + mal + g",
    c(mu = 0.1, nu = 1, beta_age = 0, beta_mal = 0, gmid = 0, ghi = 0), d
  ))

  with_g <- refit_489_fit(
    "age + g + mal",
    c(mu = 0.1, nu = 1, beta_age = 0, gmid = 0, ghi = 0, beta_mal = 0), d
  )
  dropped <- .hzr_refit_with_scope(with_g, action = "drop", var = "g",
                                   data = d)
  expect_same_fit_489(dropped, refit_489_fit(
    "age + mal", c(mu = 0.1, nu = 1, beta_age = 0, beta_mal = 0), d
  ))
})

test_that("a windowed refit matches each window's block by name (#489)", {
  d <- refit_489_data()
  base <- refit_489_fit(
    "age + mal",
    c(mu = 0.1, nu = 1, age_w1 = 0, mal_w1 = 0, age_w2 = 0, mal_w2 = 0), d,
    time_windows = 1
  )
  refit <- .hzr_refit_with_scope(base, action = "add", var = "com_iv",
                                 data = d)
  independent <- refit_489_fit(
    "age + mal + com_iv",
    c(mu = 0.1, nu = 1, age_w1 = 0, mal_w1 = 0, com_iv_w1 = 0,
      age_w2 = 0, mal_w2 = 0, com_iv_w2 = 0), d,
    time_windows = 1
  )
  expect_same_fit_489(refit, independent)
})

test_that(".hzr_refit_warm_start() places coefficients by column", {
  theta <- c(mu = 1, nu = 2, a = 3, b = 4, ab = 5)
  expect_identical(
    .hzr_refit_warm_start(theta, 2L, c("x", "y", "x:y"),
                          c("x", "y", "z", "x:y"), NULL, NULL),
    c(mu = 1, nu = 2, a = 3, b = 4, z = 0, ab = 5)
  )
  # Unnamed stays unnamed: the rest of the package names it positionally.
  expect_identical(
    .hzr_refit_warm_start(unname(theta), 2L, c("x", "y", "x:y"),
                          c("x", "y", "z", "x:y"), NULL, NULL),
    c(1, 2, 3, 4, 0, 5)
  )
  # Within each window's block, not across blocks.
  expect_identical(
    .hzr_refit_warm_start(c(s = 9, x_w1 = 1, y_w1 = 2, x_w2 = 3, y_w2 = 4),
                          1L, c("x", "y"), c("y", "z"), 1, 1),
    c(s = 9, y_w1 = 2, z_w1 = 0, y_w2 = 4, z_w2 = 0)
  )
  # Different windows: no block corresponds to another, so all start at 0.
  expect_identical(
    .hzr_refit_warm_start(c(s = 9, x_w1 = 1, x_w2 = 3), 1L, "x", "x", 1, 2),
    c(s = 9, x_w1 = 0, x_w2 = 0)
  )
  # A theta that does not fit the base design cannot be matched.
  expect_error(
    .hzr_refit_warm_start(c(1, 2, 3, 4), 1L, c("x", "y"), "x", NULL, NULL),
    "cannot be matched"
  )
})
