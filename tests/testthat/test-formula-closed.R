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
  # Controls.
  list(f = ~ splines::ns(age, df = 3), closed = TRUE),
  list(f = ~ base::log(age), closed = TRUE)
)

closed_frame <- data.frame(age = c(40, 60, 80), grp = c("a", "b", "a"))

test_that("the global closed-formula check refuses ::: and off-list qualified calls", {
  got <- vapply(closed_cases, function(k) {
    f <- stats::as.formula(paste("y", paste(deparse(k$f), collapse = "")))
    environment(f) <- globalenv()
    .hzr_formula_closed(f, closed_frame)
  }, logical(1))
  expect_identical(got, vapply(closed_cases, `[[`, logical(1), "closed"))
})

test_that("the phase closed-formula check refuses ::: and off-list qualified calls", {
  got <- vapply(closed_cases, function(k) {
    f <- k$f
    environment(f) <- globalenv()
    .hzr_phase_formula_closed(f, names(closed_frame),
                              .hzr_phase_rebuild_functions)
  }, logical(1))
  expect_identical(got, vapply(closed_cases, `[[`, logical(1), "closed"))
})
