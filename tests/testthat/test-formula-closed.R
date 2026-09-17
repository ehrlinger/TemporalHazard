# Pins the closed-formula check's refusals before its two walkers are folded
# into one (#271). `:::` reaches a namespace's internals and is never on the
# safe list; today that refusal is incidental, so a refactor that treats
# `:::` as `::` would pass it silently.

closed_cases <- list(
  list(f = ~ base:::log(age), closed = FALSE),
  list(f = ~ splines:::ns(age, df = 3), closed = FALSE),
  list(f = ~ stats:::log(age), closed = FALSE),
  list(f = ~ survival:::log(age), closed = FALSE),
  # `::` into a namespace with no list.
  list(f = ~ survival::strata(grp), closed = FALSE),
  # `log` is listed under base, not stats.
  list(f = ~ stats::log(age), closed = FALSE),
  # A namespace or function written as a string, and a qualified name that
  # is not called: the two walkers disagreed on these before the fold.
  list(f = ~ "base"::log(age), closed = FALSE),
  list(f = ~ base::"log"(age), closed = FALSE),
  list(f = ~ I(base::log), closed = FALSE),
  # Controls.
  list(f = ~ splines::ns(age, df = 3), closed = TRUE),
  list(f = ~ base::log(age), closed = TRUE)
)

closed_frame <- data.frame(age = c(40, 60, 80), grp = c("a", "b", "a"))

test_that("the global closed-formula check refuses ::: and off-list qualified calls", {
  got <- vapply(closed_cases, function(k) {
    f <- stats::as.formula(paste("y", paste(deparse(k$f), collapse = "")))
    environment(f) <- globalenv()
    .hzr_formula_closed(f, c(names(closed_frame), "."),
                        .hzr_rebuild_functions)
  }, logical(1))
  expect_identical(got, vapply(closed_cases, `[[`, logical(1), "closed"))
})

test_that("the phase closed-formula check refuses ::: and off-list qualified calls", {
  got <- vapply(closed_cases, function(k) {
    f <- k$f
    environment(f) <- globalenv()
    .hzr_formula_closed(f, names(closed_frame),
                        .hzr_rebuild_functions)
  }, logical(1))
  expect_identical(got, vapply(closed_cases, `[[`, logical(1), "closed"))
})

test_that("both legacy rebuild routes check against the one function list", {
  # Two identical lists must be kept identical by hand; one list cannot
  # drift. This fails if either route is given its own list again.
  ns <- asNamespace("TemporalHazard")
  lists <- grep("rebuild_functions$", ls(ns, all.names = TRUE), value = TRUE)
  expect_identical(lists, ".hzr_rebuild_functions")

  seen <- list()
  local_mocked_bindings(.hzr_formula_closed = function(formula, columns,
                                                       functions) {
    seen[[length(seen) + 1L]] <<- functions
    FALSE
  })
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)
  g <- hazard(survival::Surv(int_dead, dead) ~ age, data = d,
              dist = "weibull", theta = c(mu = 0.01, nu = 0.5, 0.1))
  g$data$x_design <- NULL
  .hzr_recover_x_design(g)
  m <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes",
                        formula = ~ age),
      constant = hzr_phase("constant")
    ),
    fit = FALSE
  ))
  .hzr_phase_design_from_frame(m, "early", m$spec$phases$early)
  expect_length(seen, 2L)
  expect_identical(seen[[1L]], .hzr_rebuild_functions)
  expect_identical(seen[[2L]], .hzr_rebuild_functions)
})
