# tests/testthat/test-stepwise-dot-formula.R
# hzr_stepwise() read a `~ .` base formula with terms() and no data. That
# errors, the error was turned into "no terms", and the screen reported zero
# steps as if it had finished (#279). It now stops and asks for the terms.

avc_279 <- na.omit(avc[, c("int_dead", "dead", "age", "mal")])
# Deterministic, and unrelated to survival: the screen should drop it.
avc_279$noise <- sin(seq_len(nrow(avc_279)))

test_that("a `~ .` base fit is refused before the screen prints anything", {
  fit <- hazard(survival::Surv(int_dead, dead) ~ ., data = avc_279,
                dist = "weibull", theta = c(0.3, 1, 0, 0, 0), fit = TRUE)
  err <- NULL
  out <- utils::capture.output(
    err <- tryCatch(
      hzr_stepwise(fit, data = avc_279, direction = "backward",
                   criterion = "wald", slstay = 0.2, trace = TRUE),
      error = conditionMessage
    )
  )
  expect_match(err, "write the terms out")
  expect_identical(out, character(0))
})

test_that("the term reader refuses `.` instead of returning no terms", {
  expect_error(
    .hzr_formula_rhs_terms(survival::Surv(int_dead, dead) ~ .),
    "write the terms out"
  )
  expect_identical(.hzr_formula_rhs_terms(~ age + mal), c("age", "mal"))
})

test_that("a `scope` of `~ .` is refused rather than read as empty", {
  # It used to give a screen with no candidates, and a zero-step result.
  fit <- hazard(survival::Surv(int_dead, dead) ~ age, data = avc_279,
                dist = "weibull", theta = c(0.3, 1, 0), fit = TRUE)
  expect_error(
    hzr_stepwise(fit, data = avc_279, scope = ~ ., direction = "forward",
                 criterion = "wald", trace = FALSE),
    "a `scope` can list its variables"
  )
})

mp_279 <- function(early_formula = NULL, constant_formula = NULL) {
  suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ .,
    data = avc_279[, c("int_dead", "dead", "age")], dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.15, nu = 1.4, m = 1,
                           fixed = "m", formula = early_formula),
      constant = hzr_phase("constant", formula = constant_formula)
    ),
    fit = TRUE, control = list(n_starts = 1L, conserve = FALSE)
  ))
}

test_that("a multiphase global `~ .` is refused when a phase inherits it", {
  # A phase with no formula inherits the global design, and every candidate
  # refit re-expands the global `.` against the screen's data: entering
  # `mal` into `early` also put it in `constant`, and nothing said so.
  expect_error(
    hzr_stepwise(mp_279(), data = avc_279[, c("int_dead", "dead", "age",
                                              "mal")],
                 scope = list(early = ~ mal), direction = "forward",
                 criterion = "wald", slentry = 0.99, trace = FALSE),
    "write the terms out"
  )
})

test_that("a multiphase global `~ .` runs when every phase has a formula", {
  # Then no phase reads the global design, so its `.` changes nothing.
  sw <- suppressWarnings(hzr_stepwise(
    mp_279(~ age, ~ age),
    data = avc_279[, c("int_dead", "dead", "age", "mal")],
    scope = list(early = ~ mal), direction = "forward", criterion = "wald",
    slentry = 0.99, trace = FALSE
  ))
  entered <- sw$steps[sw$steps$action == "enter", ]
  expect_identical(entered$variable, "mal")
  expect_identical(entered$phase, "early")
  expect_true("early.mal" %in% names(coef(sw)))
  expect_false("constant.mal" %in% names(coef(sw)))
})

test_that("the same screen with its terms written out drops noise", {
  fit <- hazard(survival::Surv(int_dead, dead) ~ age + mal + noise,
                data = avc_279, dist = "weibull",
                theta = c(0.3, 1, 0, 0, 0), fit = TRUE)
  sw <- hzr_stepwise(fit, data = avc_279, direction = "backward",
                     criterion = "wald", slstay = 0.2, trace = FALSE)
  expect_identical(sw$steps$action, "drop")
  expect_identical(sw$steps$variable, "noise")
  expect_identical(colnames(sw$data$x), c("age", "mal"))
})
