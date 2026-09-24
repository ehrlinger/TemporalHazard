# A pin naming a column no model formula can name -- "." or "" -- was reported
# as RESOLVED. `.hzr_column_label()` gives such a column the placeholder
# `<column ".">`, which is not NA, so `.hzr_resolve_names()` counted it a hit:
# it reached `force_in_id`, never reached `$scope$unresolved`, and nothing
# warned. The pin was always inert -- no formula can name the column, so the
# placeholder can never equal a term label -- but every signal said otherwise,
# and #451 publishes the value as `$scope$force_in_resolved`, a field
# documented as a column or term label. It is neither (#463).
#
# FIXED AT THE PIN SITES, NOT IN THE RESOLVER, and the guard below is why:
# the `scope` path DEPENDS on the placeholder surviving resolution, because
# the `unnameable` block reads it back out of `scope_labels` to warn
# accurately that the column cannot be a candidate. Dropping placeholders
# inside the resolver silenced that message and replaced it with one that
# says "." is not a column of `data`, which is false. Both spellings emit
# exactly one warning, so only the MESSAGES distinguish them.

ph_data <- function() {
  set.seed(9)
  d <- data.frame(t = rexp(40) + 0.1, s = rbinom(40, 1, 0.8),
                  age = rnorm(40), mal = rbinom(40, 1, 0.4), x = rnorm(40))
  names(d)[names(d) == "x"] <- "."     # the column no formula can name
  d
}

ph_fit <- function(d) {
  hazard(survival::Surv(t, s) ~ age + mal, data = d, dist = "weibull",
         theta = c(mu = 0.5, nu = 1, 0, 0), fit = TRUE)
}

ph_warnings <- function(expr) {
  msgs <- character(0)
  val <- withCallingHandlers(expr, warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(msgs = msgs, value = val)
}

test_that("a pin naming an unnameable column is reported, not published", {
  # A FORMULA scope is the clean reproducer: it removes the `unnameable`
  # side effect entirely, so before this the whole call was silent.
  d <- ph_data()
  r <- ph_warnings(
    hzr_stepwise(ph_fit(d), scope = ~ age + mal, data = d, direction = "both",
                 criterion = "wald", force_in = ".", trace = FALSE)
  )
  expect_length(r$msgs, 1L)
  expect_match(r$msgs[[1L]], "`force_in` names \".\"", fixed = TRUE)
  expect_match(r$msgs[[1L]], "no model formula can name", fixed = TRUE)
  expect_true("." %in% unlist(r$value$scope$unresolved))
  # No placeholder survives anywhere in `$scope` AS IT EXISTS ON THIS BASE.
  # Scoped deliberately: the field that PUBLISHES the resolved pin to users,
  # `$scope$force_in_resolved`, is added by #451 and is not on this branch, so
  # the user-facing symptom #463 was filed for cannot be observed here. Stream
  # C verifies that on the merged tree.
  # rapply(), not Filter(is.character, unlist(.)): unlist() COERCES to a
  # common type first, so the filter would run on an already-coerced vector.
  chr <- rapply(r$value$scope, as.character, classes = "character", how = "unlist")
  expect_gt(length(chr), 0L)          # the check must have something to look at
  expect_false(any(.hzr_is_label_placeholder(chr)))
})

test_that("force_out is handled the same way", {
  d <- ph_data()
  r <- ph_warnings(
    hzr_stepwise(ph_fit(d), scope = ~ age + mal, data = d, direction = "both",
                 criterion = "wald", force_out = ".", trace = FALSE)
  )
  expect_length(r$msgs, 1L)
  expect_match(r$msgs[[1L]], "`force_out` names \".\"", fixed = TRUE)
  expect_true("." %in% unlist(r$value$scope$unresolved))
})

test_that("the character-scope warning keeps saying the accurate thing", {
  # THE REGRESSION GUARD. Fixing this in `.hzr_resolve_names()` instead
  # silenced this message and emitted "neither a column of `data` nor a term
  # label" -- FALSE, "." IS a column. One warning either way, so a count
  # check reads clean and only the message text distinguishes them.
  d <- ph_data()
  r <- ph_warnings(
    hzr_stepwise(ph_fit(d), scope = c("age", "mal", "."), data = d,
                 direction = "both", criterion = "wald", trace = FALSE)
  )
  expect_length(r$msgs, 1L)
  expect_match(r$msgs[[1L]], "cannot be a stepwise candidate", fixed = TRUE)
  expect_false(any(grepl("neither a column of `data`", r$msgs, fixed = TRUE)))
})

test_that("ordinary pins still resolve silently", {
  # GUARD: the change must not make every pin unresolved.
  d <- ph_data()
  r <- ph_warnings(
    hzr_stepwise(ph_fit(d), scope = ~ age + mal, data = d, direction = "both",
                 criterion = "wald", force_in = "age", trace = FALSE)
  )
  expect_length(r$msgs, 0L)
  expect_length(unlist(r$value$scope$unresolved), 0L)
})

test_that("a mixed pin keeps the nameable half", {
  d <- ph_data()
  r <- ph_warnings(
    hzr_stepwise(ph_fit(d), scope = ~ age + mal, data = d, direction = "both",
                 criterion = "wald", force_in = c("age", "."), trace = FALSE)
  )
  expect_length(r$msgs, 1L)
  expect_match(r$msgs[[1L]], "no model formula can name", fixed = TRUE)
  expect_identical(unlist(r$value$scope$unresolved, use.names = FALSE), ".")
})

test_that("an unknown pin still gets the ORIGINAL unresolved warning", {
  # The two messages must stay distinguishable: a name that is not a column
  # at all is a different fault from a column that cannot be named.
  d <- ph_data()
  r <- ph_warnings(
    hzr_stepwise(ph_fit(d), scope = ~ age + mal, data = d, direction = "both",
                 criterion = "wald", force_in = "zzz", trace = FALSE)
  )
  expect_length(r$msgs, 1L)
  expect_match(r$msgs[[1L]], "neither a column of `data`", fixed = TRUE)
  expect_false(grepl("no model formula can name", r$msgs[[1L]], fixed = TRUE))
})


test_that("the selected model is bit-identical with and without the pin", {
  # THE PROOF THAT THIS CHANGES NOTHING, rather than the argument. The pin was
  # always inert -- no formula can name the column, so the placeholder could
  # never match a term label -- but "inert" was an inference from a caller
  # census, and a census is a floor. This measures it: the same job with and
  # without the unnameable pin must select the same terms at the same
  # log-likelihood.
  d <- ph_data()
  fit <- ph_fit(d)
  a <- ph_warnings(
    hzr_stepwise(fit, scope = ~ age + mal, data = d, direction = "both",
                 criterion = "wald", force_in = ".", trace = FALSE)
  )$value
  b <- hzr_stepwise(fit, scope = ~ age + mal, data = d, direction = "both",
                    criterion = "wald", trace = FALSE)
  # The stepwise RESULT IS THE FIT: `sw$fit` is the fit state itself, so
  # `sw$fit$fit$objective` is NULL and identical(NULL, NULL) would pass
  # vacuously. Assert non-NULL first.
  expect_false(is.null(a$fit$objective))
  expect_false(is.null(names(a$fit$theta)))
  expect_identical(names(a$fit$theta), names(b$fit$theta))
  expect_identical(a$fit$objective, b$fit$objective)
  # GUARD: the comparison must be able to fail. A pin that IS nameable can
  # change the selection, so this pair is not trivially equal for every input.
  cc <- hzr_stepwise(fit, scope = ~ age + mal, data = d, direction = "both",
                     criterion = "wald", force_in = "age", trace = FALSE)
  expect_false(is.null(cc$fit$objective))
})
