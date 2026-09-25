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
  # No placeholder survives anywhere in `$scope`, which includes
  # `force_in_resolved` (#451), the field that PUBLISHES resolved pins to users.
  # A placeholder there is the symptom #463 was filed for.
  # rapply(), not Filter(is.character, unlist(.)): unlist() COERCES to a
  # common type first, so the filter would run on an already-coerced vector.
  chr <- rapply(r$value$scope, as.character, classes = "character", how = "unlist")
  expect_gt(length(chr), 0L)          # the check must have something to look at
  expect_false(any(.hzr_is_label_placeholder(chr)))
  # The scan cannot see two failures, so each gets its own line. If the field
  # VANISHED, the scan would have one field fewer to look at and would still
  # pass, and so would `expect_length(NULL, 0L)`; hence the `is.null` check.
  # If a NON-placeholder value such as the spelling "." were published, the
  # scan would pass it, because "." is not in placeholder form and `$scope`
  # already holds "." legitimately in `force_in`; hence the length check.
  expect_false(is.null(r$value$scope$force_in_resolved))
  expect_length(r$value$scope$force_in_resolved, 0L)
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
  # GUARD, and it has to COMPARE something. An earlier version built `cc`
  # with a nameable pin and then only asserted `cc$fit$objective` was
  # non-NULL -- it never compared it, so deleting the `force_in = "age"`
  # argument left the test passing and the guard could not see the mutation
  # it was written for.
  cc <- hzr_stepwise(fit, scope = ~ age + mal, data = d, direction = "both",
                     criterion = "wald", force_in = "age", trace = FALSE)
  expect_false(is.null(cc$fit$objective))
  expect_false(identical(cc$fit$objective, b$fit$objective))
})


test_that("the default scope emits BOTH warnings, and neither may be dropped", {
  # Two guards now answer two different questions about the same column --
  # can it be a CANDIDATE, and can it be PINNED -- and on the default scope
  # both fire. Pinned together deliberately: dropping either leaves a warning
  # count of 1, which reads clean, so only asserting both messages can catch
  # it. That is this change's own cascade lesson turned on its own code
  # (stream C, #463 review).
  d <- ph_data()
  r <- ph_warnings(
    hzr_stepwise(ph_fit(d), data = d, direction = "both",
                 criterion = "wald", force_in = ".", trace = FALSE)
  )
  expect_length(r$msgs, 2L)
  expect_true(any(grepl("cannot be a stepwise candidate", r$msgs, fixed = TRUE)))
  expect_true(any(grepl("cannot be pinned", r$msgs, fixed = TRUE)))
})

test_that("the literal placeholder text does not advise renaming nothing", {
  # `known` in .hzr_resolve_names() holds every column's label, placeholders
  # included, so a user typing `<column ".">` verbatim reaches the pin path.
  # There is no column of that name, so "Rename that column" would point at
  # nothing; the advice is omitted rather than made up.
  d <- ph_data()
  r <- ph_warnings(
    hzr_stepwise(ph_fit(d), scope = ~ age + mal, data = d, direction = "both",
                 criterion = "wald", force_in = '<column ".">', trace = FALSE)
  )
  expect_length(r$msgs, 1L)
  expect_match(r$msgs[[1L]], "cannot be pinned", fixed = TRUE)
  expect_false(grepl("Rename", r$msgs[[1L]], fixed = TRUE))
  # and a real column name still gets the advice
  r2 <- ph_warnings(
    hzr_stepwise(ph_fit(d), scope = ~ age + mal, data = d, direction = "both",
                 criterion = "wald", force_in = ".", trace = FALSE)
  )
  expect_match(r2$msgs[[1L]], "Rename that column", fixed = TRUE)
})

# UNRESOLVED PINS KEEP THE USER'S ORDER (#465 review). `.hzr_resolve_names()`
# keeps its input's order, but the pin filter above appended the unnameable
# pins AFTER the names the resolver had already left unresolved, so
# `force_in = c(".", "nope")` came back as "nope", ".". The set was right;
# the order was not the one the user wrote, in the field or in print().
# The repeated-name case is here because the easy fixes -- union(), unique()
# -- would restore the order and quietly drop the repeat.
test_that("unresolved pins keep the order they were written in", {
  d <- ph_data()
  f <- ph_fit(d)
  sw_pins <- function(...) {
    suppressWarnings(
      hzr_stepwise(f, scope = ~ age + mal, data = d, direction = "both",
                   criterion = "wald", trace = FALSE, ...)
    )
  }

  sw <- sw_pins(force_in = c(".", "nope"))
  expect_identical(sw$scope$unresolved$force_in, c(".", "nope"))
  # Where the user reads it, not only where it is produced.
  expect_true(any(grepl("(unresolved `force_in`, ignored: \".\", \"nope\")",
                        capture.output(print(sw)), fixed = TRUE)))

  # The other order stays as written too, so a fix that merely moved the
  # placeholders to the FRONT would fail here.
  expect_identical(sw_pins(force_in = c("nope", "."))$scope$unresolved$force_in,
                   c("nope", "."))

  # A repeated name is kept, in place.
  expect_identical(
    sw_pins(force_in = c(".", "nope", "."))$scope$unresolved$force_in,
    c(".", "nope", ".")
  )

  # Both pin sites share the filter; `force_out` is checked separately.
  expect_identical(sw_pins(force_out = c(".", "nope"))$scope$unresolved$force_out,
                   c(".", "nope"))
})

# A COLUMN NAMED "" (#465 review). NEWS promises that a pin on a column
# called "" warns and is listed as unresolved, and until now only "." had a
# test. A "" column reaches `hzr_stepwise()` only when its `data` differs
# from the fit's: `hazard()` itself stops on any data frame carrying a ""
# column, even an unused one. So the model is fitted on clean data here.
# `direction = "forward"` with nothing left to add means no refit, so the
# pin's own warning is the only one -- a refit on this data would hit that
# separate `hazard()` failure and add warnings that are not this test's.
test_that("a pin on a column named \"\" warns and is listed as unresolved", {
  d <- ph_data()
  clean <- d[, names(d) != "."]
  names(d)[names(d) == "."] <- ""
  f <- ph_fit(clean)
  r <- ph_warnings(
    hzr_stepwise(f, scope = ~ age + mal, data = d, direction = "forward",
                 criterion = "wald", force_in = "", trace = FALSE)
  )
  expect_length(r$msgs, 1L)
  expect_match(r$msgs[[1L]], "`force_in` names \"\", which no model formula",
               fixed = TRUE)
  expect_identical(r$value$scope$unresolved$force_in, "")
  # The pin is not published as resolved. The `is.null` line comes first
  # because `expect_length(NULL, 0L)` passes. An earlier draft had only the
  # length check, on a branch where the field did not exist yet, so it could
  # not fail (#469 review).
  expect_false(is.null(r$value$scope$force_in_resolved))
  expect_length(r$value$scope$force_in_resolved, 0L)
})
