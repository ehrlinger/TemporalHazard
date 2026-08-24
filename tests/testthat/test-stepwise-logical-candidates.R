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
