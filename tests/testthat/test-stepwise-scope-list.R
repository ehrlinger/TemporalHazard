# A multiphase `scope` is a named list of one-sided formulas. A list whose
# element was a character vector stopped with "formula must be a `formula`
# object", which names neither `scope` nor the phase (#328, item 4).

scope_list_fit <- function(d) {
  hazard(survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
         phases = list(early = hzr_phase("cdf", t_half = 0.15, nu = 1.4,
                                         m = 1, fixed = "m"),
                       late = hzr_phase("constant")),
         fit = TRUE, control = list(n_starts = 1L, conserve = FALSE))
}

test_that("a multiphase scope element that is not a formula names the phase (#328)", {
  skip_on_cran() # a multiphase fit
  d <- stats::na.omit(avc[, c("int_dead", "dead", "age", "mal")])
  fit <- scope_list_fit(d)
  msg <- tryCatch(
    hzr_stepwise(fit, scope = list(late = c("age", "mal")), data = d,
                 trace = FALSE),
    error = conditionMessage
  )
  expect_identical(
    msg,
    paste0("`scope$late` must be a one-sided formula such as `~ age + mal`, ",
           "or NULL, not a character vector. A multiphase `scope` is a ",
           "named list of formulas keyed by phase.")
  )
})
