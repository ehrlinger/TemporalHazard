# hzr_stepwise() and hzr_bootstrap() check `...` against what a candidate
# refit may take from it (#386).
#
# Both forward `...` to every candidate refit, and hazard()'s own `...` is
# legacy pass-through that stores any name unread. So a misspelled argument,
# `slentyr = 1e-6` for `slentry`, vanished into the refits and the screen ran
# at the default slentry = 0.30: a different model, with no warning. Only
# `control`, and an `objective` equal to the base fit's, may pass: every
# other hazard() argument describes the model a candidate is compared with,
# and a forwarded one was ignored (time_lower on a formula refit), failed
# every candidate (dist), or changed the candidates' estimand alone
# (weights, time_windows).

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

test_that("an argument the refit takes from the base fit is refused (#386)", {
  # time_lower was ignored by a formula refit; dist collided with the
  # refit's own and failed every candidate; weights and time_windows changed
  # the candidates' likelihood but not the base model's.
  d <- dots_data_386()
  base <- dots_base_386(d)
  for (arg in list(list(time_lower = d$int_dead * 0.99),
                   list(dist = "exponential"),
                   list(weights = rep(2, nrow(d))),
                   list(time_windows = 1))) {
    msg <- refused_386(do.call(screen_386, c(list(base, d), arg)))
    expect_match(msg, paste0("^hzr_stepwise\\(\\): `", names(arg),
                             "` cannot be passed through `...`"),
                 info = names(arg))
  }
})

test_that("ambiguous, repeated and malformed arguments are refused (#386)", {
  d <- dots_data_386()
  base <- dots_base_386(d)
  # pmatch() returns NA for an ambiguous prefix as for no match; the two
  # need different messages.
  expect_match(refused_386(screen_386(base, d, ti = 1)),
               "`ti` abbreviates more than one hazard() argument",
               fixed = TRUE)
  # After spelling out, `cont` is a second `control`.
  expect_match(refused_386(screen_386(base, d, control = list(),
                                      cont = list(maxit = 1L))),
               "`control` is given more than once", fixed = TRUE)
  # hazard() requires a list; NULL would fail every candidate.
  expect_match(refused_386(screen_386(base, d, control = NULL)),
               "`control` must be a list", fixed = TRUE)
})

test_that("`objective` passes only when it equals the base fit's (#386)", {
  skip_on_cran() # a stepwise screen
  d <- dots_data_386()
  base <- dots_base_386(d)
  expect_match(refused_386(screen_386(base, d, objective = "sas")),
               "^hzr_stepwise\\(\\): `objective = \"sas\"` differs")
  # An abbreviated name is spelled out first, so it is checked too.
  expect_match(refused_386(screen_386(base, d, objec = "sas")),
               "^hzr_stepwise\\(\\): `objective = \"sas\"` differs")
  # hazard() match.arg()s the value, so "lik" is the base fit's objective.
  sw <- screen_386(base, d, objective = "lik")
  expect_identical(sort(sw$steps$variable), c("age", "com_iv", "mal"))
})

test_that("`control` reaches the refits, and so does an abbreviation (#386)", {
  skip_on_cran() # two stepwise screens
  # The final model is the last accepted refit, and a fit records the
  # control it was given, so this reads the forwarded value itself rather
  # than inferring it from a screen that could be empty for other reasons.
  d <- dots_data_386()
  base <- dots_base_386(d)
  full <- screen_386(base, d, control = list(maxit = 77L))
  abbr <- screen_386(base, d, contr = list(maxit = 77L))
  expect_identical(nrow(full$steps), 3L)
  expect_identical(full$spec$control$maxit, 77L)
  expect_identical(abbr$spec$control$maxit, 77L)
  expect_identical(abbr$steps, full$steps)
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
  # Unresampled weights would be misaligned with every replicate's rows.
  expect_match(refused_386(boot(weights = rep(1, nrow(d)))),
               "^hzr_bootstrap\\(\\): `weights` cannot be passed")
  expect_match(refused_386(boot(objective = "sas")),
               "^hzr_bootstrap\\(\\): `objective = \"sas\"` differs")
  expect_identical(.Random.seed, before)
})

test_that("hzr_bootstrap() still forwards `control` and `trace` (#386)", {
  skip_on_cran() # bootstrap replicates with a screen each
  # reltol = 1e-2 stops each refit sooner but still converged, so com_iv
  # enters every replicate either way and only its estimate moves. A run
  # whose candidates all failed would not select com_iv at all.
  d <- dots_data_386()
  base <- dots_base_386(d)
  boot <- function(...) {
    hzr_bootstrap(base, n_boot = 2L, seed = 1L, scope = ~ com_iv,
                  criterion = "wald", trace = TRUE, ...)
  }
  est <- function(bs) bs$summary$mean[bs$summary$parameter == "com_iv"]
  tight <- boot()
  loose <- boot(control = list(reltol = 1e-2))
  expect_identical(c(tight$n_success, loose$n_success), c(2L, 2L))
  expect_identical(tight$summary$n[tight$summary$parameter == "com_iv"], 2L)
  expect_identical(loose$summary$n[loose$summary$parameter == "com_iv"], 2L)
  expect_false(identical(est(tight), est(loose)))
})
