# A `...` name hazard() does not declare is refused at the entry of
# hzr_stepwise() and hzr_bootstrap() (#386).
#
# Both forward `...` to every candidate refit, and hazard()'s own `...` is
# legacy pass-through that accepts any name. So a misspelled argument,
# `slentyr = 1e-6` for `slentry`, vanished into the refits and the screen ran
# at the default slentry = 0.30: a different model, with no warning.

dots_data_386 <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  stats::na.omit(avc[, c("int_dead", "dead", "age", "mal", "com_iv")])
}

dots_base_386 <- function(d) {
  hazard(survival::Surv(int_dead, dead) ~ 1, data = d, dist = "weibull",
         theta = c(0.1, 1), fit = TRUE)
}

screen_386 <- function(base, d, ...) {
  hzr_stepwise(base, scope = ~ age + mal + com_iv, data = d,
               direction = "forward", criterion = "wald", trace = FALSE, ...)
}

test_that("the typo is live: the two spellings select different models (#386)", {
  skip_on_cran() # three stepwise screens
  # Known positive for the fixture: slentry = 1e-6 admits one variable and
  # the default admits three, so a typo that falls back to the default is a
  # different answer, not a harmless no-op.
  d <- dots_data_386()
  base <- dots_base_386(d)
  strict <- screen_386(base, d, slentry = 1e-6)
  default <- screen_386(base, d)
  expect_identical(strict$steps$variable, "com_iv")
  expect_identical(sort(default$steps$variable), c("age", "com_iv", "mal"))
})

test_that("hzr_stepwise() refuses a misspelled argument and names it (#386)", {
  d <- dots_data_386()
  base <- dots_base_386(d)
  msg <- tryCatch(screen_386(base, d, slentyr = 1e-6),
                  error = conditionMessage)
  expect_match(msg, "^hzr_stepwise\\(\\): `slentyr`")
  expect_match(msg, "Did you mean `slentry`?", fixed = TRUE)
})

test_that("an unnamed argument in `...` is refused too (#386)", {
  # An unnamed value reaches `...` only once every formal before it is
  # filled; earlier, R matches it positionally to the next formal.
  d <- dots_data_386()
  base <- dots_base_386(d)
  expect_error(
    hzr_stepwise(base, ~ age, d, "forward", "wald", 0.30, 0.20, 50L, 4L,
                 character(), character(), FALSE, 1e-6),
    "an unnamed argument", fixed = TRUE
  )
})

test_that("a name hazard() declares is still forwarded (#386)", {
  skip_on_cran() # a stepwise screen
  # `control` is the documented use of `...`; maxit is read by every dist,
  # so it draws no control warning either.
  d <- dots_data_386()
  base <- dots_base_386(d)
  sw <- expect_no_warning(
    screen_386(base, d, slentry = 1e-6, control = list(maxit = 500L))
  )
  expect_identical(sw$steps$variable, "com_iv")
})

test_that("hzr_bootstrap() refuses it at entry, not per replicate (#386)", {
  # Inside a replicate, hzr_stepwise()'s refusal would be caught and tallied
  # as a replicate failure. The bootstrap must refuse before resampling, and
  # before seeding, so the caller's random number stream is left alone.
  d <- dots_data_386()
  base <- dots_base_386(d)
  set.seed(42)
  before <- .Random.seed
  msg <- tryCatch(
    hzr_bootstrap(base, n_boot = 3L, seed = 1L, scope = ~ age + mal,
                  criterion = "wald", slentyr = 1e-6),
    error = conditionMessage
  )
  expect_match(msg, "^hzr_bootstrap\\(\\): `slentyr`")
  expect_match(msg, "Did you mean `slentry`?", fixed = TRUE)
  expect_identical(.Random.seed, before)
})

test_that("hzr_bootstrap() still forwards `control` and `trace` (#386)", {
  skip_on_cran() # bootstrap replicates with a screen each
  d <- dots_data_386()
  base <- dots_base_386(d)
  bs <- hzr_bootstrap(base, n_boot = 2L, seed = 1L, scope = ~ com_iv,
                      criterion = "wald", control = list(maxit = 500L),
                      trace = TRUE)
  expect_identical(bs$n_success, 2L)
  expect_true("com_iv" %in% bs$summary$parameter)
})
