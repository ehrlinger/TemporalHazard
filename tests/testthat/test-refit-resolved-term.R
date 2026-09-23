# A candidate enters as the variable it resolved to, not as the text of its
# name (#449).
#
# The refit used to paste the candidate's SPELLING into the formula text, so
# a column whose name reads as a different term entered as that term: a
# column `age ` (trailing space) beside `age` refit as `age`, and a literal
# column `age:mal` refit as the interaction. The screen reported the column
# it had chosen and fitted a different one, with no warning. The refit now
# writes the candidate's term label, which terms() itself produced, so it
# parses back to the column that was resolved.
#
# Every test runs a real screen to completion and reads the FINAL model. Two
# columns that hold different values give different fits, so "the right
# column was used" is checked against a direct fit on that column, and
# against a direct fit on the look-alike, which must differ.

rrt_avc <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  stats::na.omit(avc[, c("int_dead", "dead", "age", "mal")])
}

# `age` and a column `age ` whose values are NOT age's: a strong predictor
# built from `mal`, which is left out of the frame.
rrt_lookalike <- function(nm = "age ") {
  withr::local_seed(449L)
  d0 <- rrt_avc()
  d <- data.frame(d0[, c("int_dead", "dead", "age")],
                  x = d0$mal * 3 + stats::rnorm(nrow(d0), sd = 0.1),
                  check.names = FALSE)
  names(d)[names(d) == "x"] <- nm
  d
}

rrt_fit <- function(d, rhs, theta) {
  f <- stats::as.formula(paste("survival::Surv(int_dead, dead) ~", rhs))
  suppressWarnings(do.call(hazard, list(
    formula = f, data = d, dist = "weibull", theta = theta, fit = TRUE
  )))
}

rrt_final_terms <- function(sw) {
  attr(stats::terms(stats::formula(sw$call$formula)), "term.labels")
}

rrt_screen <- function(fit, d, ...) {
  suppressWarnings(hzr_stepwise(fit, data = d, direction = "forward",
                                slentry = 0.99, trace = FALSE, ...))
}

test_that("the fixture separates the two columns (#449)", {
  skip_on_cran() # fits
  # Known positive for every test below: the columns give different fits,
  # so a screen that fitted the wrong one would show a different logLik.
  d <- rrt_lookalike()
  both <- rrt_fit(d, "age + `age `", c(0.1, 1, 0, 0))
  age2 <- rrt_fit(d, "age + age", c(0.1, 1, 0))
  expect_gt(abs(both$fit$objective - age2$fit$objective), 1)
})

test_that("scope = NULL enters the column `age `, not `age` (#449)", {
  skip_on_cran() # full screens under three criteria
  # `age` is already in the model, so `age ` is the only candidate: pasted,
  # it read back as `age` and added nothing.
  d <- rrt_lookalike()
  base <- rrt_fit(d, "age", c(0.1, 1, 0))
  direct <- rrt_fit(d, "age + `age `", c(0.1, 1, 0, 0))
  for (crit in c("wald", "aic", "score")) {
    sw <- rrt_screen(base, d, criterion = crit)
    expect_setequal(rrt_final_terms(sw), c("age", "`age `"))
    expect_identical(sw$criteria$n_refit_failures, 0L, info = crit)
    # The DATA of `age ` was used: the final model is the direct fit.
    expect_equal(sw$fit$objective, direct$fit$objective, tolerance = 1e-6,
                 info = crit)
  }
})

test_that("a character scope enters the column it names (#449)", {
  skip_on_cran() # full screens
  d <- rrt_lookalike("age # x")
  base <- rrt_fit(d, "1", c(0.1, 1))
  direct <- rrt_fit(d, "`age # x`", c(0.1, 1, 0))
  for (crit in c("wald", "score")) {
    sw <- rrt_screen(base, d, criterion = crit, scope = "age # x")
    expect_identical(rrt_final_terms(sw), "`age # x`", info = crit)
    expect_equal(sw$fit$objective, direct$fit$objective, tolerance = 1e-6,
                 info = crit)
  }
})

test_that("score and refit agree for a literal `age:mal` column (#449, a)", {
  skip_on_cran() # full screens
  # The score read the COLUMN while the refit entered the INTERACTION, so the
  # p-value belonged to one variable and the model to another.
  withr::local_seed(1L)
  d0 <- rrt_avc()
  d <- data.frame(d0, `age:mal` = d0$mal * 2 + stats::rnorm(nrow(d0), sd = 0.5),
                  check.names = FALSE)
  base <- rrt_fit(d, "age + mal", c(0.1, 1, 0, 0))
  direct <- rrt_fit(d, "age + mal + `age:mal`", c(0.1, 1, 0, 0, 0))
  inter <- rrt_fit(d, "age + mal + age:mal", c(0.1, 1, 0, 0, 0))
  expect_gt(abs(direct$fit$objective - inter$fit$objective), 1e-3) # known +
  # The same column under a name nobody could misread: its score p-value.
  d2 <- d
  names(d2)[names(d2) == "age:mal"] <- "litcol"
  base2 <- rrt_fit(d2, "age + mal", c(0.1, 1, 0, 0))
  p_col <- rrt_screen(base2, d2, criterion = "score",
                      scope = "litcol")$steps$p_value[1L]
  for (sc in list(NULL, "age:mal")) {
    sw <- rrt_screen(base, d, criterion = "score", scope = sc)
    expect_identical(sw$steps$p_value[1L], p_col)
    expect_setequal(rrt_final_terms(sw), c("age", "mal", "`age:mal`"))
    expect_equal(sw$fit$objective, direct$fit$objective, tolerance = 1e-6)
  }
})

test_that("a multiphase default scope enters the column `age ` (#449)", {
  skip_on_cran() # multiphase screens
  d <- rrt_lookalike()
  mp <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(constant = hzr_phase("constant")),
    fit = TRUE, control = list(n_starts = 1L, maxit = 500L)
  ))
  sw <- suppressWarnings(hzr_stepwise(mp, data = d, direction = "forward",
                                      criterion = "wald", slentry = 0.99,
                                      trace = FALSE,
                                      control = list(n_starts = 1L)))
  terms_c <- attr(stats::terms(sw$spec$phases$constant$formula),
                  "term.labels")
  expect_setequal(terms_c, c("age", "`age `"))
  expect_identical(sw$criteria$n_refit_failures, 0L)
})

test_that("a multiphase score screen scores and refits the column `x2 ` (#449)", {
  skip_on_cran() # multiphase screens
  # The multiphase score built the candidate's phase formula from TEXT too,
  # so a strong column named `x2 ` was scored as the noise column `x2`; and
  # the Wald fallback that rescues a strong candidate refit `x2` as well.
  # Same data as test-score-wald-fallback.R's planted effect, with the
  # strong column renamed.
  withr::local_seed(11L)
  n <- 400L
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  x3 <- stats::rnorm(n)
  d <- data.frame(tt = stats::rexp(n, rate = 0.3 * exp(0.9 * x1)) + 0.01,
                  ev = stats::rbinom(n, 1, 0.85), x1 = x1, x2 = x2, x3 = x3)
  screen_on <- function(dd) {
    mp <- suppressWarnings(hazard(
      survival::Surv(tt, ev) ~ 1, data = dd, dist = "multiphase",
      phases = list(early = hzr_phase("cdf", t_half = 1, nu = 1.5, m = 0),
                    const = hzr_phase("constant")),
      fit = TRUE
    ))
    suppressWarnings(hzr_stepwise(mp, data = dd, direction = "forward",
                                  criterion = "score", max_steps = 1L,
                                  trace = FALSE))
  }
  ref <- screen_on(d) # the syntactic name: the known positive
  expect_identical(ref$steps$variable, "x1")
  expect_gt(ref$criteria$n_wald_fallbacks, 0L)
  d2 <- d
  names(d2)[names(d2) == "x1"] <- "x2 "
  sw <- screen_on(d2)
  expect_identical(sw$steps$variable, "`x2 `")
  expect_identical(sw$steps$phase, ref$steps$phase)
  expect_identical(sw$steps$p_value, ref$steps$p_value)
  expect_identical(sw$criteria$n_wald_fallbacks, ref$criteria$n_wald_fallbacks)
  expect_identical(attr(stats::terms(sw$spec$phases$const$formula),
                        "term.labels"), "`x2 `")
  expect_equal(sw$fit$objective, ref$fit$objective, tolerance = 1e-8)
})

test_that("score reads a formula-scope `_X1` from its column (#438)", {
  skip_on_cran() # full screens
  # #438's reproduction: the score looked the column up by its backquoted
  # label, reported it "not found in `data`", and skipped the strongest
  # covariate, while wald entered it first.
  withr::local_seed(11L)
  n <- 400L
  d <- data.frame(age = stats::rnorm(n))
  d[["_X1"]] <- stats::rnorm(n)
  d$time <- stats::rexp(n) * exp(-1.0 * d[["_X1"]])
  d$status <- 1L
  base <- suppressWarnings(do.call(hazard, list(
    formula = survival::Surv(time, status) ~ 1, data = d, dist = "weibull",
    theta = c(0.5, 1), fit = TRUE
  )))
  for (sc in list(~ age + `_X1`, c("age", "`_X1`"))) {
    w <- character()
    sw <- withCallingHandlers(
      hzr_stepwise(base, scope = sc, data = d, direction = "forward",
                   criterion = "score", slentry = 0.5, trace = FALSE),
      warning = function(cnd) {
        w <<- c(w, conditionMessage(cnd))
        invokeRestart("muffleWarning")
      }
    )
    expect_false(any(grepl("not found", w, fixed = TRUE)))
    expect_identical(sw$steps$variable[1L], "`_X1`")
    expect_true("`_X1`" %in% rrt_final_terms(sw))
  }
})

test_that("score does not read a literal `age:mal` column for the interaction (#449)", {
  skip_on_cran() # full screens
  # The interaction candidate is spelled like the literal column, and the
  # score read that column for it. An interaction is not a column the score
  # can test, so it is declined as it is when no such column exists.
  withr::local_seed(1L)
  d0 <- rrt_avc()
  d <- data.frame(d0, `age:mal` = d0$mal * 2 + stats::rnorm(nrow(d0), sd = 0.5),
                  check.names = FALSE)
  run <- function(dd) {
    base <- rrt_fit(dd, "age + mal", c(0.1, 1, 0, 0))
    suppressWarnings(hzr_stepwise(base, data = dd, scope = ~ age:mal,
                                  direction = "forward", criterion = "score",
                                  slentry = 0.99, trace = FALSE))
  }
  without <- run(d0) # the known behaviour with no look-alike column
  with <- run(d)
  expect_identical(nrow(with$steps), nrow(without$steps))
  expect_identical(with$criteria$uncomputable_reasons,
                   without$criteria$uncomputable_reasons)
})

test_that("hzr_bootstrap() replicates enter the column `age ` (#449)", {
  skip_on_cran() # a bootstrap
  d <- rrt_lookalike()
  base <- rrt_fit(d, "1", c(0.1, 1))
  bs <- suppressWarnings(suppressMessages(hzr_bootstrap(
    base, n_boot = 3L, seed = 1L, scope = c("age", "age "),
    criterion = "wald", slentry = 0.99
  )))
  expect_identical(bs$n_success, 3L)
  # Both columns selected in every replicate: before, `age ` refit as `age`,
  # added no column, and was never selected.
  pct <- stats::setNames(bs$summary$pct, bs$summary$parameter)
  expect_identical(unname(pct[c("age", "`age `")]), c(100, 100))
  # Each replicate's `age ` coefficient is that column's. The two columns'
  # effects are far apart in this fixture (`age ` is built from `mal`), so a
  # replicate that fitted `age`'s data under the `age ` name would show it.
  est <- bs$replicates
  a1 <- est$estimate[est$parameter == "`age `"]
  a0 <- est$estimate[est$parameter == "age"]
  expect_length(a1, 3L)
  expect_length(a0, 3L)
  direct <- rrt_fit(d, "age + `age `", c(0.1, 1, 0, 0))
  expect_true(all(a1 > 0.05))       # direct fit on the full data: see below
  expect_true(all(abs(a0) < 0.05))
  expect_gt(direct$fit$theta[4L], 0.05)
  expect_lt(abs(direct$fit$theta[3L]), 0.05)
})

test_that("an entry spelled `_X1` and its drop are one variable to max_move (#441)", {
  skip_on_cran() # a two-way screen that oscillates
  # Once a bare "_X1" can enter, the entry carries the spelling and the drop
  # the label `` `_X1` ``. Counted under both, the variable took twice
  # max_move moves to freeze. A noise column enters (p < slentry) and fails
  # slstay, so it oscillates: frozen at step 6, dropped at step 7, stopped.
  withr::local_seed(2L)
  d0 <- rrt_avc()
  d <- data.frame(d0, `_X1` = stats::rnorm(nrow(d0)), check.names = FALSE)
  fit <- rrt_fit(d, "age + mal", c(0.1, 1, 0, 0))
  sw <- suppressWarnings(hzr_stepwise(fit, data = d, scope = "_X1",
                                      direction = "both", criterion = "wald",
                                      slentry = 0.99, slstay = 0.2,
                                      max_steps = 20L, trace = FALSE))
  expect_identical(sum(sw$steps$action == "enter"), 3L)
  expect_identical(sw$steps$action[6L], "frozen")
  expect_identical(nrow(sw$steps), 7L)
  expect_false(sw$criteria$hit_max_steps)
  expect_identical(sw$scope$frozen, "`_X1`")
})

test_that("wald, aic and score agree beside the age:mal interaction (#442, #449)", {
  skip_on_cran() # nine full screens
  # The #442 fixture: the interaction age:mal is in the model and a data
  # column is literally named `age:mal`. That column is a variable of its
  # own, so every criterion must reach the same final model. At 17a5b009
  # wald and aic entered it while score compared its column NAME with the
  # model's term labels, matched the interaction, and declined it as
  # not_expandable. Run over a weak (nearly collinear with age), a noise and
  # a strong column, because a weak enough column CAN hide the divergence:
  # whether it does depends on how weak it is.
  d0 <- rrt_avc()
  make <- list(
    weak   = function() {
      as.numeric(scale(d0$age)) * 0.5 + 1 + stats::rnorm(nrow(d0), sd = 0.01)
    },
    noise  = function() stats::rnorm(nrow(d0)),
    strong = function() d0$mal * 2 + stats::rnorm(nrow(d0), sd = 0.5)
  )
  for (k in names(make)) {
    withr::local_seed(1L)
    d <- data.frame(d0, `age:mal` = make[[k]](), check.names = FALSE)
    fit <- rrt_fit(d, "age * mal", c(0.1, 1, 0, 0, 0))
    final <- lapply(c(wald = "wald", aic = "aic", score = "score"),
                    function(crit) rrt_screen(fit, d, criterion = crit))
    terms_by <- lapply(final, rrt_final_terms)
    expect_setequal(terms_by$wald, c("age", "mal", "age:mal", "`age:mal`"))
    expect_identical(terms_by$aic, terms_by$wald, info = k)
    expect_identical(terms_by$score, terms_by$wald, info = k)
    for (crit in names(final)) {
      expect_identical(length(final[[crit]]$criteria$uncomputable_reasons),
                       0L, info = paste(k, crit))
    }
  }
})

test_that("an untested entry is reported by the name a drop would use (#441)", {
  skip_on_cran() # full screens
  # $criteria$wald_untested_entries is public. An entry and a removal of the
  # same variable must be written the same way there, the resolved label,
  # whatever spelling the scope used, or a reader matching the two lists, or
  # matching either against $steps, misses it. The entry's Wald test is
  # masked to NA, as test-stepwise-backward-uncomputable.R does for x3.
  withr::local_seed(2L)
  d0 <- rrt_avc()
  d <- data.frame(d0, `_X1` = stats::rnorm(nrow(d0)), check.names = FALSE)
  fit <- rrt_fit(d, "age + mal", c(0.1, 1, 0, 0))
  orig <- .hzr_candidate_score
  local_mocked_bindings(.hzr_candidate_score = function(...) {
    a <- list(...)
    s <- orig(...)
    if (identical(a$mode, "entry")) {
      cols <- colnames(a$candidate$data$x)
      var <- cols[match(a$names, paste0("beta", seq_along(cols)))]
      if (identical(var, "`_X1`")) {
        s$score <- NA_real_
        s$p_value <- NA_real_
        s$stat <- NA_real_
      }
    }
    s
  })
  for (sc in list("_X1", "`_X1`", ~ `_X1`)) {
    sw <- suppressWarnings(hzr_stepwise(fit, data = d, scope = sc,
                                        direction = "forward",
                                        criterion = "wald", trace = FALSE))
    expect_identical(sw$criteria$wald_untested_entries, "`_X1`")
  }
})

test_that("a column no formula can name is skipped, saying so (#449)", {
  skip_on_cran() # full screens
  # A column literally named "." has no term label: terms() reads `.` as
  # "every other column". Its placeholder identity was pasted into formula
  # text, so a multiphase default scope stopped with the parser's
  # "unexpected '<'", and a single-distribution one failed three refits the
  # same way. It is not a candidate; the screen says why, once, and goes on.
  withr::local_seed(1L)
  d <- rrt_avc()
  d[["."]] <- stats::rnorm(nrow(d))
  mp <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(constant = hzr_phase("constant")),
    fit = TRUE, control = list(n_starts = 1L)
  ))
  sd1 <- rrt_fit(d, "1", c(0.1, 1))
  for (base in list(mp, sd1)) {
    w <- character()
    sw <- NULL
    expect_no_error(sw <- withCallingHandlers(
      hzr_stepwise(base, data = d, direction = "forward", criterion = "wald",
                   slentry = 0.99, trace = FALSE,
                   control = list(n_starts = 1L)),
      warning = function(cnd) {
        w <<- c(w, conditionMessage(cnd))
        invokeRestart("muffleWarning")
      }
    ))
    expect_identical(sum(grepl("\".\"", w, fixed = TRUE) &
                           grepl("cannot be a stepwise candidate", w,
                                 fixed = TRUE)), 1L)
    expect_false(any(grepl("unexpected", w, fixed = TRUE)))
    expect_identical(sw$criteria$n_refit_failures, 0L)
  }
})

test_that("score says why it declines an interaction beside a literal column (#449)", {
  skip_on_cran() # a full screen
  # The interaction is a term, not a column; that a column of the same
  # spelling exists must not make the decline read "not found in `data`".
  withr::local_seed(1L)
  d0 <- rrt_avc()
  d <- data.frame(d0, `age:mal` = d0$mal * 2 + stats::rnorm(nrow(d0), sd = 0.5),
                  check.names = FALSE)
  base <- rrt_fit(d, "age + mal", c(0.1, 1, 0, 0))
  w <- character()
  sw <- withCallingHandlers(
    hzr_stepwise(base, data = d, scope = ~ age:mal, direction = "forward",
                 criterion = "score", slentry = 0.99, trace = FALSE),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(any(grepl("not found in `data`", w, fixed = TRUE)))
  expect_true(any(grepl("is not a single column of `data`", w, fixed = TRUE)))
  expect_true(any(grepl("different variable", w, fixed = TRUE)))
  expect_identical(nrow(sw$steps), 0L)
})

test_that("$steps$variable is the term label on entry and drop alike (#449)", {
  skip_on_cran() # two-way screens
  # One variable, one name in the step table, whatever form the scope took:
  # an entry used to carry the scope's spelling (`_X1`) and a drop the label
  # (`` `_X1` ``), so filtering the table by either missed half the rows.
  withr::local_seed(2L)
  d0 <- rrt_avc()
  d <- data.frame(d0, `_X1` = stats::rnorm(nrow(d0)), check.names = FALSE)
  fit <- rrt_fit(d, "age + mal", c(0.1, 1, 0, 0))
  for (sc in list("_X1", NULL, "`_X1`", ~ `_X1`)) {
    sw <- suppressWarnings(hzr_stepwise(fit, data = d, scope = sc,
                                        direction = "both", criterion = "wald",
                                        slentry = 0.99, slstay = 0.2,
                                        max_steps = 4L, trace = FALSE))
    acts <- sw$steps$action
    expect_true(all(c("enter", "drop") %in% acts)) # known positive
    what <- paste(deparse(sc), collapse = "")
    # The entry row and the drop row for the same column, in one screen.
    expect_identical(unique(sw$steps$variable[acts == "enter"]), "`_X1`",
                     info = what)
    expect_identical(unique(sw$steps$variable[acts == "drop"]), "`_X1`",
                     info = what)
  }
})

test_that("score records an interaction it declines as not_single_column (#449)", {
  skip_on_cran() # full screens
  # An interaction is a term, not a column the score can test. Its reason
  # was `non_numeric`, which describes a column, with or without a literal
  # column of the same spelling beside it.
  withr::local_seed(1L)
  d0 <- rrt_avc()
  d <- data.frame(d0, `age:mal` = d0$mal * 2 + stats::rnorm(nrow(d0), sd = 0.5),
                  check.names = FALSE)
  for (dd in list(d0, d)) {
    base <- rrt_fit(dd, "age + mal", c(0.1, 1, 0, 0))
    sw <- suppressWarnings(hzr_stepwise(base, data = dd, scope = ~ age:mal,
                                        direction = "forward",
                                        criterion = "score", slentry = 0.99,
                                        trace = FALSE))
    expect_identical(sw$criteria$uncomputable_reasons,
                     c(not_single_column = 1L))
  }
  expect_match(.hzr_score_reason_text("not_single_column"),
               "not a single column", fixed = TRUE)
})
