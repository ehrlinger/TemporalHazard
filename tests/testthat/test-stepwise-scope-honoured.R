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

test_that("a two-sided scope names its left-hand side under every direction (#343)", {
  # NEWS promises the error names the left-hand side. The backward refusal
  # ran first, so a two-sided scope under backward got the generic "no
  # effect" message instead (Copilot on #367). The structural checks now
  # come first, for both callers, and an empty two-sided scope is caught too.
  lhs <- "`scope` must be one-sided: its left-hand side (`com_iv`) would be ignored"
  d <- avc_343()
  b0 <- weibull_343(d, survival::Surv(int_dead, dead) ~ 1, c(0.1, 1))
  for (dir in c("backward", "forward", "both")) {
    for (sc in list(com_iv ~ age + mal, com_iv ~ 1)) {
      expect_error(
        hzr_stepwise(b0, scope = sc, data = d, direction = dir,
                     criterion = "wald", trace = FALSE),
        lhs, fixed = TRUE, info = paste(dir, deparse(sc))
      )
    }
  }
  expect_error(
    hzr_bootstrap(b0, n_boot = 2L, seed = 1L, scope = com_iv ~ 1,
                  direction = "backward"),
    lhs, fixed = TRUE
  )
  # A duplicated phase is refused under backward too, empty entries or not.
  expect_error(
    .hzr_refuse_unhonoured_scope(list(early = ~ 1, early = ~ 1), "backward"),
    "`scope` names `early` more than once", fixed = TRUE
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
  for (sc in list(~ 1, ~ 0, character(0), list(), list(character(0)),
                  list(early = NULL, late = character(0)))) {
    sw <- suppressWarnings(
      hzr_stepwise(base, scope = sc, data = d, direction = "backward",
                   criterion = "wald", slstay = 1e-300, trace = FALSE)
    )
    expect_identical(sw$steps$variable, c("age", "mal"),
                     label = deparse(sc))
  }
  expect_null(.hzr_refuse_unhonoured_scope(list(early = NULL, late = ~ 1),
                                           "backward"))
  # An offset is not empty: "both" refuses it, so backward must not pass it.
  expect_error(.hzr_refuse_unhonoured_scope(~ offset(age), "backward"),
               "`scope` has no effect", fixed = TRUE)
  expect_error(
    .hzr_refuse_unhonoured_scope(list(early = NULL, late = ~ age), "backward"),
    "`scope` has no effect", fixed = TRUE
  )
})

test_that("hzr_bootstrap() with an empty scope under backward runs a real backward screen (#343)", {
  skip_on_cran() # eight replicate screens
  d <- avc_343()
  full <- weibull_343(d, survival::Surv(int_dead, dead) ~ age + mal + com_iv,
                      c(0.1, 1, 0, 0, 0))
  b <- suppressWarnings(
    hzr_bootstrap(full, n_boot = 8L, seed = 2L, scope = ~ 1,
                  direction = "backward", criterion = "wald", slstay = 0.05)
  )
  expect_identical(b$mode, "select")
  expect_equal(b$n_success, 8L)
  pct <- b$summary$pct[b$summary$parameter %in% c("age", "mal", "com_iv")]
  # A screen, not a refit: some term was dropped in some replicate.
  expect_length(pct, 3L)
  expect_true(any(pct < 100))
})

test_that("an empty-scope backward bootstrap reaches #389's untested-removal count", {
  # Interaction of #343 (an empty scope runs a backward screen per replicate,
  # and selection arguments without a scope are refused before seeding) with
  # #389 (hzr_bootstrap() counts replicates that decided a variable without a
  # Wald test). x3's removal test is NA while x2 is in the model: a replicate
  # that drops x2 then tests x3; one that keeps x2 leaves x3 untested.
  obj <- .fit_overfitted()
  orig <- .hzr_candidate_score
  local_mocked_bindings(
    .hzr_candidate_score = function(...) {
      a <- list(...)
      s <- orig(...)
      if (identical(a$mode, "drop")) {
        cols <- colnames(a$current$data$x)
        var <- cols[match(a$names, paste0("beta", seq_along(cols)))]
        if (identical(var, "x3") && "x2" %in% cols) {
          s$score <- NA_real_
          s$p_value <- NA_real_
          s$stat <- NA_real_
        }
      }
      s
    }
  )
  screens <- list()
  orig_sw <- hzr_stepwise
  local_mocked_bindings(
    hzr_stepwise = function(...) {
      r <- orig_sw(...)
      screens[[length(screens) + 1L]] <<- list(
        stopped = isTRUE(r$criteria$stopped_uncomputable),
        listed = length(c(r$criteria$wald_untested_removals,
                          r$criteria$wald_untested_entries)) > 0L
      )
      r
    }
  )
  w <- testthat::capture_warnings(
    boot <- hzr_bootstrap(obj$fit, n_boot = 4, seed = 1, scope = ~ 1,
                          direction = "backward", slstay = 0.20)
  )
  expect_identical(boot$n_failed, 0L)
  reps <- utils::tail(screens, boot$n_success)
  expect_length(reps, boot$n_success)
  stopped <- vapply(reps, `[[`, logical(1L), "stopped")
  listed <- vapply(reps, `[[`, logical(1L), "listed")
  expect_identical(boot$n_uncomputable_replicates, sum(stopped))
  expect_gt(boot$uncomputable_reasons[["wald_no_variance"]], 0L)
  hit <- grepl("successful replicates decided a variable without a Wald test",
               w)
  if (any(listed)) {
    expect_true(any(grepl(paste0("^", sum(listed), " of ", boot$n_success,
                                 " successful replicates decided"), w)))
  } else {
    expect_false(any(hit))
  }

  # The refusals come first: no replicate runs, so nothing is counted.
  n_before <- length(screens)
  expect_error(hzr_bootstrap(obj$fit, n_boot = 4, seed = 1,
                             direction = "backward", slstay = 0.10),
               "only take")
  expect_error(hzr_bootstrap(obj$fit, n_boot = 4, seed = 1, scope = ~ x3,
                             direction = "backward"),
               "backward")
  expect_identical(length(screens), n_before)
})
