# hzr_stepwise() and hzr_bootstrap() check `...` against what a candidate
# refit reads (#386).
#
# Both forward `...` to every candidate refit, and hazard()'s own `...` is
# legacy pass-through that stores any name unread. So a misspelled argument,
# `slentyr = 1e-6` for `slentry`, vanished into the refits and the screen ran
# at the default slentry = 0.30: a different model, with no warning. The
# refit reads only control, weights, time_windows and objective; every other
# hazard() argument it sets itself, so a forwarded one was either ignored
# (time_lower on a formula refit) or failed every candidate (dist).

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

refused_386 <- function(expr) {
  tryCatch({
    expr
    "NO ERROR"
  }, error = conditionMessage)
}

test_that("the typo is live: the two spellings select different models (#386)", {
  skip_on_cran() # two stepwise screens
  # Known positive for the fixture: slentry = 1e-6 admits one variable and
  # the default admits three, so a typo that falls back to the default is a
  # different answer, not a harmless no-op.
  d <- dots_data_386()
  base <- dots_base_386(d)
  expect_identical(screen_386(base, d, slentry = 1e-6)$steps$variable,
                   "com_iv")
  expect_identical(sort(screen_386(base, d)$steps$variable),
                   c("age", "com_iv", "mal"))
})

test_that("hzr_stepwise() refuses a misspelled argument and names it (#386)", {
  d <- dots_data_386()
  base <- dots_base_386(d)
  msg <- refused_386(screen_386(base, d, slentyr = 1e-6))
  expect_match(msg, "^hzr_stepwise\\(\\): `slentyr` is not an argument")
  expect_match(msg, "Did you mean `slentry`?", fixed = TRUE)
  # An abbreviation is spelled out before the check, so a differing
  # objective is refused at entry rather than colliding with the refit's own
  # `objective` inside every candidate. (Not in hzr_bootstrap(), where R
  # binds `objec` to its own `object` formal first.)
  expect_match(refused_386(screen_386(base, d, objec = "sas")),
               "^hzr_stepwise\\(\\): `objective = \"sas\"` differs")
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

test_that("an argument the refit sets itself is refused, not ignored (#386)", {
  # `time_lower` on a formula refit was ignored with no message; `dist`
  # collided with the refit's own and failed every candidate.
  d <- dots_data_386()
  base <- dots_base_386(d)
  for (arg in list(list(time_lower = d$int_dead * 0.99),
                   list(dist = "exponential"))) {
    msg <- refused_386(do.call(screen_386, c(list(base, d), arg)))
    expect_match(msg, paste0("^hzr_stepwise\\(\\): `", names(arg),
                             "` cannot be passed through `...`"),
                 info = names(arg))
  }
})

test_that("`control` is forwarded, and so is an abbreviation of it (#386)", {
  skip_on_cran() # three stepwise screens
  # maxit = 1 stops every candidate refit short, so no variable enters;
  # the default admits three. The forwarding is therefore observable. An
  # abbreviation hazard() applied by partial matching still works.
  d <- dots_data_386()
  base <- dots_base_386(d)
  full <- suppressWarnings(screen_386(base, d, control = list(maxit = 1L)))
  abbr <- suppressWarnings(screen_386(base, d, contr = list(maxit = 1L)))
  expect_identical(nrow(full$steps), 0L)
  expect_identical(abbr$steps, full$steps)
  expect_identical(nrow(screen_386(base, d)$steps), 3L)
})

test_that("hzr_bootstrap() refuses at entry, not per replicate (#386)", {
  # Inside a replicate, hzr_stepwise()'s refusal would be caught and tallied
  # as a replicate failure. The bootstrap must refuse before resampling, and
  # before seeding, so the caller's random number stream is left alone.
  d <- dots_data_386()
  base <- dots_base_386(d)
  boot <- function(...) {
    hzr_bootstrap(base, n_boot = 3L, seed = 1L, scope = ~ age + mal,
                  criterion = "wald", ...)
  }
  set.seed(42)
  before <- .Random.seed
  msg <- refused_386(boot(slentyr = 1e-6))
  expect_match(msg, "^hzr_bootstrap\\(\\): `slentyr` is not an argument")
  expect_match(msg, "Did you mean `slentry`?", fixed = TRUE)
  expect_match(refused_386(boot(dist = "exponential")),
               "^hzr_bootstrap\\(\\): `dist` cannot be passed")
  expect_match(refused_386(boot(objective = "sas")),
               "^hzr_bootstrap\\(\\): `objective = \"sas\"` differs")
  expect_identical(.Random.seed, before)
})

test_that("hzr_bootstrap() still forwards `control` and `trace` (#386)", {
  skip_on_cran() # bootstrap replicates with a screen each
  d <- dots_data_386()
  base <- dots_base_386(d)
  boot <- function(...) {
    hzr_bootstrap(base, n_boot = 2L, seed = 1L, scope = ~ com_iv,
                  criterion = "wald", trace = TRUE, ...)
  }
  bs <- boot()
  expect_identical(bs$n_success, 2L)
  expect_true("com_iv" %in% bs$summary$parameter)
  # maxit = 1 reaches the screen's candidate refits, so com_iv never enters.
  short <- suppressWarnings(boot(control = list(maxit = 1L)))
  expect_false("com_iv" %in% short$summary$parameter)
  # The base fit's own objective, restated, is accepted.
  expect_identical(boot(objective = "likelihood")$n_success, 2L)
})
