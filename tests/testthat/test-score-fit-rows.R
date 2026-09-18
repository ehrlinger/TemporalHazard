# The score path's row guard counts the rows the FIT used (#372).
#
# A multiphase fit drops every row with a missing value in a phase covariate,
# but keeps the caller's full `time` in `$data`. .hzr_score_q() guarded row
# alignment against `length(current$data$time)` -- the full count -- so it
# passed, every candidate then failed to line up with the fit's design, and a
# score screen stopped with zero steps, labelling each candidate
# `not_expandable`: the candidate blamed for a base-data problem. The guard's
# own comment names this case and says to fail loudly; it now can.

avc_372 <- function() avc[, c("int_dead", "dead", "age", "com_iv", "opmos")]

fit_372 <- function(d) {
  suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes",
                        formula = ~ age),
      constant = hzr_phase("constant", formula = ~ age)
    ),
    fit = TRUE, control = list(n_starts = 1L, maxit = 500L)
  ))
}

screen_372 <- function(base, d) {
  suppressWarnings(hzr_stepwise(
    base, scope = list(early = ~ com_iv + opmos, constant = ~ com_iv + opmos),
    data = d, direction = "forward", criterion = "score", slentry = 0.99,
    max_steps = 1L, trace = FALSE
  ))
}

na_rows_372 <- function(k) {
  d <- avc_372()
  set.seed(1)
  idx <- sample(nrow(d), k)
  d$age[idx] <- NA
  list(d = d, idx = idx)
}

test_that("a base fit that dropped rows stops the score screen, naming the cause (#372)", {
  skip_on_cran() # multiphase fits
  x <- na_rows_372(15L)
  base <- fit_372(x$d)
  # The fixture's mechanism: the fit used fewer rows than `time` records.
  expect_identical(NROW(base$fit$x_list$early), 295L)
  expect_length(base$data$time, 310L)
  expect_error(screen_372(base, x$d), "dropped 15 rows", fixed = TRUE)
  expect_error(screen_372(base, x$d), "Refit the base model on complete cases",
               fixed = TRUE)
})

test_that("passing complete-case data without refitting is refused the same way (#372)", {
  skip_on_cran() # multiphase fits
  # The obvious workaround: the data now has the fit's 295 rows, but the fit
  # still carries 310 in `time`, so the candidate cannot line up either.
  x <- na_rows_372(15L)
  base <- fit_372(x$d)
  expect_error(screen_372(base, x$d[-x$idx, ]), "dropped 15 rows",
               fixed = TRUE)
})

test_that("the exact-multiple case is refused too, not scored from misaligned rows (#372)", {
  skip_on_cran() # multiphase fits
  # 155 of 310: the lengths are exact multiples, so R's recycling warning
  # never fires. The guard must not depend on that warning.
  x <- na_rows_372(155L)
  base <- fit_372(x$d)
  expect_error(screen_372(base, x$d), "dropped 155 rows", fixed = TRUE)
})

test_that("refitting on complete cases scores normally, so the fixture is live (#372)", {
  skip_on_cran() # multiphase fits
  # Control: the remedy the error names works, and the screen it unblocks
  # makes a real selection.
  x <- na_rows_372(15L)
  cc <- x$d[-x$idx, ]
  sw <- screen_372(fit_372(cc), cc)
  expect_identical(sw$steps$variable, "com_iv")
  expect_identical(sw$steps$phase, "early")
  expect_equal(sw$criteria$n_uncomputable_scores %||% 0L, 0L)
})
