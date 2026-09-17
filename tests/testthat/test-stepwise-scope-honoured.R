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

test_that("hzr_bootstrap()'s backward refusal points at a screen it can run (#343)", {
  d <- avc_343()
  b0 <- weibull_343(d, survival::Surv(int_dead, dead) ~ 1, c(0.1, 1))
  msg <- tryCatch(hzr_bootstrap(b0, n_boot = 2L, scope = ~ age,
                                direction = "backward"),
                  error = conditionMessage)
  expect_match(msg, "pass an empty `scope` such as `~ 1`", fixed = TRUE)
  expect_no_match(msg, "leave `scope` unset", fixed = TRUE)
})

test_that("hzr_bootstrap() refuses selection arguments without a scope (#343)", {
  # Without `scope` there is no screen: these were ignored, and every term
  # came back at pct = 100.
  d <- avc_343()
  full <- weibull_343(d, survival::Surv(int_dead, dead) ~ age + mal + com_iv,
                      c(0.1, 1, 0, 0, 0))
  set.seed(7)
  before <- .Random.seed
  expect_error(
    hzr_bootstrap(full, n_boot = 3L, seed = 1L, direction = "backward",
                  force_in = "age", slstay = 1e-300),
    paste0("hzr_bootstrap(): `direction`, `slstay`, `force_in` only take ",
           "effect in a selection screen, which needs `scope`. Either pass ",
           "`scope` to screen on each replicate, or omit `direction`, ",
           "`slstay`, `force_in` to refit the fit's exact model"),
    fixed = TRUE
  )
  expect_identical(.Random.seed, before)
  for (a in list(list(criterion = "wald"), list(slentry = 0.5),
                 list(max_steps = 2L), list(max_move = 1L),
                 list(force_out = "mal"))) {
    expect_error(do.call(hzr_bootstrap, c(list(full, n_boot = 2L), a)),
                 paste0("`", names(a), "` only takes effect"), fixed = TRUE,
                 label = names(a))
  }
  # The fixed-model bootstrap itself is unchanged.
  b <- suppressWarnings(hzr_bootstrap(full, n_boot = 2L, seed = 1L))
  expect_equal(b$n_success, 2L)
})

test_that("a scope list naming a phase twice is refused (#343)", {
  expect_error(
    .hzr_refuse_unhonoured_scope(list(early = ~ age, early = ~ mal), "both"),
    "`scope` names `early` more than once", fixed = TRUE
  )
  expect_null(.hzr_refuse_unhonoured_scope(list(early = ~ age, late = ~ mal),
                                           "both"))
})

test_that("under direction = both, scope limits entry, not drops, as documented (#343)", {
  # SAS STEPWISE re-tests every term in the model. `mal` is outside the
  # scope and still leaves; only force_in keeps a term.
  d <- avc_343()
  base <- weibull_343(d)
  sw <- suppressWarnings(
    hzr_stepwise(base, scope = ~ age, data = d, direction = "both",
                 criterion = "wald", slstay = 1e-300, trace = FALSE)
  )
  expect_true("mal" %in% sw$steps$variable[sw$steps$action == "drop"])
  kept <- suppressWarnings(
    hzr_stepwise(base, scope = ~ age, data = d, direction = "both",
                 criterion = "wald", slstay = 1e-300, force_in = "mal",
                 trace = FALSE)
  )
  expect_false("mal" %in% kept$steps$variable[kept$steps$action == "drop"])
})

test_that("a wrapper forwarding the selection defaults without a scope still bootstraps (#343)", {
  # Every selection argument is unread without `scope`, so its default value
  # changes nothing and asks for nothing. Only a value that differs from the
  # default is a selection setting that would be ignored.
  d <- avc_343()
  full <- weibull_343(d)
  my_boot <- function(fit, direction = c("both", "forward", "backward"),
                      criterion = "score", slentry = 0.3, slstay = 0.2,
                      max_steps = 50, max_move = 4L,
                      force_in = character(), force_out = NULL) {
    hzr_bootstrap(fit, n_boot = 2L, seed = 1L, direction = direction,
                  criterion = criterion, slentry = slentry, slstay = slstay,
                  max_steps = max_steps, max_move = max_move,
                  force_in = force_in, force_out = force_out)
  }
  b <- suppressWarnings(my_boot(full))
  expect_equal(b$n_success, 2L)
  expect_identical(b$mode, "refit")
  expect_error(my_boot(full, direction = "forward"),
               "`direction` only takes effect", fixed = TRUE)
  expect_error(my_boot(full, slstay = 0.05),
               "`slstay` only takes effect", fixed = TRUE)
  expect_error(my_boot(full, force_out = "mal"),
               "`force_out` only takes effect", fixed = TRUE)
})

test_that("an empty scope under direction = backward is honoured, not refused (#343)", {
  # A backward screen enters nothing, and an empty scope offers nothing to
  # enter, so the two agree.
  d <- avc_343()
  base <- weibull_343(d)
  for (sc in list(~ 1, ~ 0, character(0), list())) {
    sw <- suppressWarnings(
      hzr_stepwise(base, scope = sc, data = d, direction = "backward",
                   criterion = "wald", slstay = 1e-300, trace = FALSE)
    )
    expect_identical(sw$steps$variable, c("age", "mal"),
                     label = deparse(sc))
  }
  expect_null(.hzr_refuse_unhonoured_scope(list(early = NULL, late = ~ 1),
                                           "backward"))
  expect_error(
    .hzr_refuse_unhonoured_scope(list(early = NULL, late = ~ age), "backward"),
    "`scope` has no effect", fixed = TRUE
  )
})
