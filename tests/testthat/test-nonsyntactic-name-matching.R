# A non-syntactic variable name is matched by NAME, not by the spelling
# `terms()` gives its label (#437).
#
# terms() backquotes a label whose variable is not a syntactic name: `_X1`
# for the column `_X1`. `force_in`, `force_out` and a character `scope` are
# documented as variables, and the translator emits bare names, so the two
# spellings never met. The pin was ignored and the variable was dropped,
# with nothing said -- a wrong model and no message.
#
# SAS names reach R this way routinely: a leading underscore, a dot, a
# reserved word.

# `age` is the column a free backward step drops here, so IT is the one
# renamed: a pin on `_X1` then has to REDIRECT the drop onto `mal`. Renaming
# the other column instead gives a fixture where an ignored pin and an
# honoured pin produce the same answer, and the test cannot fail (#437).
ns_data_437 <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc[, c("int_dead", "dead", "age", "mal")])
  names(d)[names(d) == "age"] <- "_X1"
  d
}

ns_fit_437 <- function(d) {
  suppressWarnings(hazard(survival::Surv(int_dead, dead) ~ `_X1` + mal,
                          data = d, dist = "weibull",
                          theta = c(0.1, 1, 0, 0), fit = TRUE))
}

ns_drop_437 <- function(fit, d, ...) {
  suppressWarnings(.hzr_stepwise_backward_step(fit, data = d,
                                               criterion = "wald",
                                               slstay = 1e-9, ...))
}

test_that("the fixture is live: unpinned, the NON-SYNTACTIC one is dropped (#437)", {
  skip_on_cran() # a backward step
  # Known positive, and it fixes which answer means "the pin was ignored".
  # Without it, a test that the pin holds could pass because nothing was
  # going to be dropped at all.
  d <- ns_data_437()
  st <- ns_drop_437(ns_fit_437(d), d)
  expect_true(st$accepted)
  expect_identical(st$variable, "`_X1`")
})

test_that("pinning the non-syntactic one redirects the drop (#437)", {
  skip_on_cran() # a backward step
  # The sharp form. `_X1` is what an unpinned step drops, so an IGNORED pin
  # drops `_X1` and an HONOURED one drops `mal`: the two answers differ, and
  # the unfixed code gives the first.
  d <- ns_data_437()
  st <- ns_drop_437(ns_fit_437(d), d, force_in = "_X1")
  expect_true(st$accepted)
  expect_identical(st$variable, "mal")
})

test_that("force_in pins a non-syntactic variable, written bare (#437)", {
  skip_on_cran() # a backward step
  # The documented contract, and what the translator emits.
  d <- ns_data_437()
  st <- ns_drop_437(ns_fit_437(d), d, force_in = c("mal", "_X1"))
  expect_false(st$accepted)
  expect_gt(nrow(st$all_scores), 0L) # or `all()` is vacuously true
  expect_true(all(st$all_scores$force_in))
})

test_that("force_in also accepts the label spelling (#437)", {
  skip_on_cran() # a backward step
  # A regression guard, not a demonstration: the label spelling matched
  # before the fix too. It pins the decision that the key map is
  # many-to-one, so someone reading the label off $steps names the same
  # variable as someone writing it bare.
  d <- ns_data_437()
  st <- ns_drop_437(ns_fit_437(d), d, force_in = c("mal", "`_X1`"))
  expect_false(st$accepted)
  expect_gt(nrow(st$all_scores), 0L) # or `all()` is vacuously true
  expect_true(all(st$all_scores$force_in))
})

test_that("force_out keeps a non-syntactic candidate out (#437)", {
  skip_on_cran() # a candidate enumeration
  # The mirror of force_in, through a formula scope, whose terms are labels.
  d <- ns_data_437()
  base <- suppressWarnings(hazard(survival::Surv(int_dead, dead) ~ 1,
                                  data = d, dist = "weibull",
                                  theta = c(0.1, 1), fit = TRUE))
  offered <- function(...) {
    vapply(.hzr_stepwise_candidates(base, data = d, ...),
           function(c) c$var, character(1))
  }
  # Known positive: it is offered when not excluded.
  expect_true(any(offered(scope = ~ `_X1` + mal) %in% c("_X1", "`_X1`")))
  expect_false(any(offered(scope = ~ `_X1` + mal, force_out = "_X1") %in%
                     c("_X1", "`_X1`")))
})

test_that("a character scope does not re-offer a variable in the model (#437)", {
  skip_on_cran() # a candidate enumeration
  # A character scope carries bare names; the model's own terms are labels.
  # Unmatched, `_X1` was offered as a candidate although the model already
  # had it, and the refit then failed on a duplicate column.
  d <- ns_data_437()
  fit <- ns_fit_437(d)
  offered <- vapply(.hzr_stepwise_candidates(fit, scope = c("mal", "_X1"),
                                             data = d),
                    function(c) c$var, character(1))
  expect_length(offered, 0L)
})

test_that("a syntactic name is unaffected (#437)", {
  skip_on_cran() # a backward step
  # The control: the change must not alter matching where label and name
  # already agree.
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc[, c("int_dead", "dead", "age", "mal")])
  fit <- suppressWarnings(hazard(survival::Surv(int_dead, dead) ~ age + mal,
                                 data = d, dist = "weibull",
                                 theta = c(0.1, 1, 0, 0), fit = TRUE))
  st_free <- ns_drop_437(fit, d)
  expect_true(st_free$accepted)
  st_pinned <- ns_drop_437(fit, d, force_in = c("age", "mal"))
  expect_false(st_pinned$accepted)
})

test_that("a scope naming a variable twice offers it once (#437)", {
  skip_on_cran() # a candidate enumeration
  # `setdiff()` de-duplicates as well as subtracting, and replacing it with a
  # key comparison silently dropped that: the variable was scored, refit and
  # reported TWICE, with nothing said. Found by review of this fix, not by
  # the fix's own tests.
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc[, c("int_dead", "dead", "age", "mal")])
  base <- suppressWarnings(hazard(survival::Surv(int_dead, dead) ~ 1,
                                  data = d, dist = "weibull",
                                  theta = c(0.1, 1), fit = TRUE))
  vars <- function(...) {
    vapply(.hzr_stepwise_candidates(base, data = d, ...),
           function(c) c$var, character(1))
  }
  expect_identical(vars(scope = c("age", "age", "mal")), c("age", "mal"))
  # The same variable written two ways is still one variable.
  expect_identical(vars(scope = c("age", "`age`")), "age")
})

test_that("the key is parsed, not stripped of backticks (#437)", {
  skip_on_cran() # a backward step
  # Stripping backticks is not the inverse of quoting. A column whose name
  # CONTAINS one is labelled `a\`b`, which strips to "a\\b" -- a name that
  # does not exist -- while parsing gives back "a`b". Measured on a live
  # fit, not reasoned from the grammar.
  nm <- "a`b"
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc[, c("int_dead", "dead", "age", "mal")])
  names(d)[names(d) == "mal"] <- nm
  fit <- suppressWarnings(hazard(
    stats::as.formula("survival::Surv(int_dead, dead) ~ age + `a\\`b`"),
    data = d, dist = "weibull", theta = c(0.1, 1, 0, 0), fit = TRUE
  ))
  labels <- attr(stats::terms(stats::formula(fit$call$formula)), "term.labels")
  expect_identical(.hzr_var_key(labels), c("age", nm))
  expect_false(identical(gsub("`", "", labels, fixed = TRUE)[2L], nm))
  # And it holds end to end: the pin is honoured.
  st <- ns_drop_437(fit, d, force_in = c("age", nm))
  expect_false(st$accepted)
})

test_that("a term that is not a single symbol keeps its label (#437)", {
  # An interaction or a function call has no single variable to name, so it
  # must be returned unchanged and go on matching by label as before.
  expect_identical(.hzr_var_key(c("z:f", "I(age > 50)", "log(age)")),
                   c("z:f", "I(age > 50)", "log(age)"))
  expect_identical(.hzr_var_key(character()), character())
})
