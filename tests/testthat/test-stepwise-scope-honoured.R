# A `scope` the screen cannot honour is refused, not ignored (#343).
#
# A backward screen never read `scope`: from `~ age + mal` with
# `scope = ~ age` it dropped `mal`, the variable the caller had left out, and
# a scope variable absent from the base was never tested. A two-sided scope
# formula lost its left-hand side: `com_iv ~ age + mal` screened only `age`
# and `mal`. Each ran to a result with no message.

avc_343 <- function() {
  stats::na.omit(avc[, c("int_dead", "dead", "age", "mal", "com_iv")])
}
weibull_343 <- function(d, formula = survival::Surv(int_dead, dead) ~ age + mal,
                        theta = c(0.1, 1, 0, 0)) {
  hazard(formula, data = d, dist = "weibull", theta = theta, fit = TRUE)
}

test_that("a backward screen with a scope is refused, naming the remedy (#343)", {
  d <- avc_343()
  base <- weibull_343(d)
  for (sc in list(~ age, ~ age + mal + com_iv, "age")) {
    expect_error(
      hzr_stepwise(base, scope = sc, data = d, direction = "backward",
                   criterion = "wald", slstay = 1e-300, trace = FALSE),
      "`scope` has no effect when `direction = \"backward\"`", fixed = TRUE
    )
  }
  msg <- tryCatch(
    hzr_stepwise(base, scope = ~ age, data = d, direction = "backward",
                 trace = FALSE),
    error = conditionMessage
  )
  expect_match(msg, "force_in", fixed = TRUE)
})

test_that("a backward screen without a scope, and a two-way screen with one, still run (#343)", {
  d <- avc_343()
  base <- weibull_343(d)
  back <- suppressWarnings(
    hzr_stepwise(base, data = d, direction = "backward", criterion = "wald",
                 slstay = 1e-300, trace = FALSE)
  )
  expect_identical(back$steps$variable, c("age", "mal"))

  b0 <- weibull_343(d, survival::Surv(int_dead, dead) ~ 1, c(0.1, 1))
  both <- suppressWarnings(
    hzr_stepwise(b0, scope = ~ mal, data = d, direction = "both",
                 criterion = "wald", slentry = 0.9999, trace = FALSE)
  )
  expect_identical(both$steps$variable, "mal")
})

test_that("a multiphase backward screen with a scope list is refused (#343)", {
  skip_on_cran() # a multiphase fit
  d <- avc_343()
  fit <- hazard(survival::Surv(int_dead, dead) ~ 1, data = d,
                dist = "multiphase",
                phases = list(constant = hzr_phase("constant",
                                                   formula = ~ age + mal)),
                fit = TRUE, control = list(n_starts = 1L, conserve = FALSE))
  expect_error(
    hzr_stepwise(fit, scope = list(constant = ~ age), data = d,
                 direction = "backward", trace = FALSE),
    "`scope` has no effect when `direction = \"backward\"`", fixed = TRUE
  )
})

test_that("hzr_bootstrap() refuses a backward selection screen before seeding (#343)", {
  d <- avc_343()
  b0 <- weibull_343(d, survival::Surv(int_dead, dead) ~ 1, c(0.1, 1))
  set.seed(99)
  before <- .Random.seed
  expect_error(
    hzr_bootstrap(b0, n_boot = 2L, seed = 1L, scope = ~ age,
                  direction = "backward"),
    "`scope` has no effect when `direction = \"backward\"`", fixed = TRUE
  )
  expect_identical(.Random.seed, before)
})

test_that("a two-sided scope formula is refused, naming its left-hand side (#343)", {
  d <- avc_343()
  b0 <- weibull_343(d, survival::Surv(int_dead, dead) ~ 1, c(0.1, 1))
  expect_error(
    hzr_stepwise(b0, scope = com_iv ~ age + mal, data = d,
                 direction = "forward", criterion = "wald", slentry = 0.9999,
                 trace = FALSE),
    "`scope` must be one-sided: its left-hand side (`com_iv`) would be ignored",
    fixed = TRUE
  )
})

test_that("a two-sided multiphase scope element is refused, naming its phase (#343)", {
  skip_on_cran() # a multiphase fit
  d <- avc_343()
  fit <- hazard(survival::Surv(int_dead, dead) ~ 1, data = d,
                dist = "multiphase",
                phases = list(constant = hzr_phase("constant")),
                fit = TRUE, control = list(n_starts = 1L, conserve = FALSE))
  expect_error(
    hzr_stepwise(fit, scope = list(constant = com_iv ~ age), data = d,
                 direction = "forward", trace = FALSE),
    "`scope$constant` must be one-sided: its left-hand side (`com_iv`) would be ignored",
    fixed = TRUE
  )
})
