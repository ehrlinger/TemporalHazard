# A logical column is an ordinary 0/1 predictor -- .hzr_modellable_vars() says
# so explicitly, and offers logicals as candidates under `scope = NULL`. Both
# criteria then refused them, so a screen died before step 1 on a column the
# package had chosen for itself. The assertion that matters is not that the
# screen runs, but that a logical column and its numeric twin are the same
# model: same selection, same coefficient, to the optimizer's tolerance.

make_pair <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  d_num <- avc[, c("int_dead", "dead", "age", "mal")]
  d_lgl <- d_num
  d_lgl$mal <- as.logical(d_num$mal)
  list(num = d_num, lgl = d_lgl)
}

base_fit <- function(d) {
  hazard(survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
         phases = list(constant = hzr_phase("constant")), fit = TRUE)
}

for (crit in c("score", "wald")) {
  test_that(paste0("a logical candidate screens as its numeric twin (", crit, ")"), {
    skip_on_cran()
    p <- make_pair()

    # The fixture must actually exercise the logical path, or this proves nothing.
    expect_type(p$num$mal, "integer")
    expect_type(p$lgl$mal, "logical")
    expect_true("mal" %in% TemporalHazard:::.hzr_modellable_vars(p$lgl, "mal"))

    r_num <- hzr_stepwise(base_fit(p$num), scope = NULL, data = p$num,
                          direction = "forward", criterion = crit)
    r_lgl <- hzr_stepwise(base_fit(p$lgl), scope = NULL, data = p$lgl,
                          direction = "forward", criterion = crit)

    # The logical column must actually enter, or comparing the two runs
    # compares two screens that both ignored it and proves nothing.
    expect_true("mal" %in% r_num$steps$variable)
    expect_true("mal" %in% r_lgl$steps$variable)

    # Same variables entered, in the same order, on the same evidence.
    expect_equal(r_lgl$steps$variable, r_num$steps$variable)
    expect_equal(r_lgl$steps$action, r_num$steps$action)
    expect_equal(r_lgl$steps$p_value, r_num$steps$p_value, tolerance = 1e-6)

    # And the same fit underneath, not merely the same bookkeeping. The
    # coefficient NAMES differ by design (model.matrix() calls the logical's
    # column malTRUE), so compare the values.
    expect_equal(unname(stats::coef(r_lgl)), unname(stats::coef(r_num)),
                 tolerance = 1e-6)
  })
}

test_that("the score criterion's refusal names an escape hatch that works", {
  skip_on_cran()

  # A two-level factor expands to one design-matrix column, so the Wald path
  # refits and tests it. The score path cannot expand it and refuses -- and the
  # refusal used to tell the caller that `criterion = "wald"` "rejects it too,
  # so switching criterion will not help". That stopped being true once the
  # coefficient-name lookup returned the expansion it had found, so the message
  # was steering people away from the one thing that works.

  data(avc, package = "TemporalHazard", envir = environment())
  d <- avc[, c("int_dead", "dead", "age", "mal")]
  d$mal <- factor(ifelse(d$mal > 0, "b", "a"))
  sc <- stats::setNames(list(~ age + mal), "constant")
  base <- function() {
    hazard(survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
           phases = list(constant = hzr_phase("constant")), fit = TRUE)
  }

  err <- tryCatch(
    hzr_stepwise(base(), scope = sc, data = d, direction = "forward",
                 criterion = "score"),
    error = conditionMessage
  )
  expect_match(err, "not numeric", fixed = TRUE)
  expect_false(grepl("switching criterion will not help", err, fixed = TRUE))
  expect_match(err, 'criterion = "wald"', fixed = TRUE)

  # Naming an escape hatch that does not work would be worse than naming none.
  r <- suppressWarnings(
    hzr_stepwise(base(), scope = sc, data = d, direction = "forward",
                 criterion = "wald")
  )
  expect_true("mal" %in% r$steps$variable)
})
