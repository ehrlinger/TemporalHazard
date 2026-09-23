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
  offered <- function(fit, ...) {
    vapply(.hzr_stepwise_candidates(fit, data = d, ...),
           function(c) c$var, character(1))
  }
  # Known positive FIRST: without it, `expect_length(0)` is satisfied by an
  # enumeration that returned nothing at all, which is indistinguishable
  # from one that matched correctly and had nothing to offer.
  base <- suppressWarnings(hazard(survival::Surv(int_dead, dead) ~ 1,
                                  data = d, dist = "weibull",
                                  theta = c(0.1, 1), fit = TRUE))
  expect_setequal(offered(base, scope = c("mal", "_X1")), c("mal", "_X1"))
  expect_length(offered(fit, scope = c("mal", "_X1")), 0L)
})

test_that("entry needs the QUOTED spelling; a bare one fails loudly (#437)", {
  skip_on_cran() # three forward screens
  # The distinction the NEWS bullet has to make, pinned so the prose cannot
  # drift from it. Matching a non-syntactic variable is not the same as
  # being able to ADD one: the refit pastes the candidate into a formula, so
  # it parses only when the candidate is already quoted. A formula scope
  # carries `terms()` labels and therefore does; a character scope of bare
  # names does not, and that failure is LOUD (#441).
  d <- ns_data_437()
  base <- suppressWarnings(hazard(survival::Surv(int_dead, dead) ~ 1,
                                  data = d, dist = "weibull",
                                  theta = c(0.1, 1), fit = TRUE))
  fwd <- function(...) {
    suppressWarnings(hzr_stepwise(base, data = d, direction = "forward",
                                  criterion = "aic", trace = FALSE, ...))
  }
  entered <- function(sw) sw$steps$variable[sw$steps$action == "enter"]

  # A formula scope: the label is already quoted, so the add succeeds.
  sw_f <- fwd(scope = ~ `_X1` + mal)
  expect_true("`_X1`" %in% entered(sw_f))
  expect_identical(sw_f$criteria$n_refit_failures, 0L)

  # A character scope written bare: it cannot be added, and says so.
  sw_c <- fwd(scope = c("_X1", "mal"))
  expect_false(any(entered(sw_c) %in% c("_X1", "`_X1`")))
  expect_gt(sw_c$criteria$n_refit_failures, 0L)

  # Written with backticks, the same character scope does add it: the
  # barrier is the spelling of the candidate, not the shape of the scope.
  expect_true("`_X1`" %in% entered(fwd(scope = c("`_X1`", "mal"))))
})

test_that("multiphase candidate enumeration matches by variable too (#437)", {
  skip_on_cran() # a multiphase fit
  # The two multiphase sites in .hzr_stepwise_candidates() are changed by
  # this fix and every other test here is single-distribution, so reverting
  # either would go unnoticed.
  d <- ns_data_437()
  mp <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE, control = list(n_starts = 2L, maxit = 500L)
  ))
  vars <- function(...) {
    vapply(.hzr_stepwise_candidates(mp, data = d, ...),
           function(c) c$var, character(1))
  }
  sc <- list(early = ~ `_X1` + mal, constant = NULL)
  # Known positive: both are offered when nothing is excluded.
  expect_setequal(vars(scope = sc), c("`_X1`", "mal"))
  # force_out, written bare, must reach the label.
  expect_identical(vars(scope = sc, force_out = "_X1"), "mal")
  # And the default-scope site, on a syntactic frame so the pre-existing
  # bare-name paste in that branch is not what is under test.
  data(avc, package = "TemporalHazard", envir = environment())
  dd <- stats::na.omit(avc[, c("int_dead", "dead", "age", "mal")])
  mp2 <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = dd, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE, control = list(n_starts = 2L, maxit = 500L)
  ))
  all_vars <- vapply(.hzr_stepwise_candidates(mp2, data = dd),
                     function(c) c$var, character(1))
  expect_true("age" %in% all_vars) # known positive
  kept <- vapply(.hzr_stepwise_candidates(mp2, data = dd, force_out = "age"),
                 function(c) c$var, character(1))
  expect_false("age" %in% kept)
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
  # Assert the RELATION, not the key's spelling: the label must key the same
  # as the bare name, and differently from the backtick-stripped text.
  expect_identical(.hzr_var_key(labels), .hzr_var_key(c("age", nm)))
  stripped <- gsub("`", "", labels, fixed = TRUE)[2L]
  expect_false(identical(stripped, nm))
  expect_false(identical(.hzr_var_key(labels)[2L], .hzr_var_key(stripped)))
  # And it holds end to end: the pin is honoured.
  st <- ns_drop_437(fit, d, force_in = c("age", nm))
  expect_false(st$accepted)
})

test_that("a term that is not a single symbol keeps its label (#437)", {
  # An interaction or a function call has no single variable to name, so it
  # must be returned unchanged and go on matching by label as before.
  # Two labels with the same text key the same...
  expect_identical(.hzr_var_key("z:f"), .hzr_var_key("z:f"))
  # ...but an expression NEVER keys as the literal column of the same text,
  # or the column disappears from the candidates (#442).
  for (e in c("z:f", "I(age > 50)", "log(age)")) {
    expect_false(identical(.hzr_var_key(e),
                           .hzr_var_key(paste0("`", e, "`"))),
                 info = e)
  }
  expect_identical(.hzr_var_key(character()), character())
})

# --- Codex review of #442: expression labels must stay distinct from literal
# column names (discussion_r4082123105). -------------------------------------

test_that("a literal column named like an interaction is still offered (#442)", {
  skip_on_cran() # a fit plus candidate enumeration
  # `terms()` backquotes the literal column `age:mal`, and parsing that label
  # gives the symbol `age:mal`. The INTERACTION label age:mal is a `:` call,
  # not a symbol. Keyed on the name alone the two collided, and the literal
  # column silently vanished from the candidates while the screen finished
  # normally. main offers it; the first version of this fix did not.
  data(avc, package = "TemporalHazard", envir = environment())
  d0 <- stats::na.omit(avc[, c("int_dead", "dead", "age", "mal")])
  d <- data.frame(d0, `age:mal` = as.numeric(scale(d0$age)) * 0.5 + 1,
                  check.names = FALSE)
  # The column is genuinely not the product, or the test proves nothing.
  expect_false(isTRUE(all.equal(d[["age:mal"]], d$age * d$mal)))
  fit <- suppressWarnings(hazard(survival::Surv(int_dead, dead) ~ age * mal,
                                 data = d, dist = "weibull",
                                 theta = c(0.1, 1, 0, 0, 0), fit = TRUE))
  expect_true("age:mal" %in%
                attr(stats::terms(stats::formula(fit$call$formula)),
                     "term.labels")) # the interaction IS in the model
  offered <- vapply(.hzr_stepwise_candidates(fit, scope = ~ `age:mal`, data = d),
                    function(c) c$var, character(1))
  expect_identical(offered, "`age:mal`")
})

test_that("the key preserves main's distinctions except where #437 merges them", {
  # A DIFFERENTIAL test, added because two silent regressions in this PR came
  # from generalising a comparison without enumerating what the old one
  # DISTINGUISHED. main compared raw strings. For every pair of label shapes,
  # the new key must agree with that string comparison, except for an
  # explicit allow-set of the pairs #437 deliberately merges.
  shapes <- c("age", "`age`", "_X1", "`_X1`", "a b", "`a b`",
              "age:mal", "`age:mal`", "log(age)", "`log(age)`",
              "I(age > 50)", "`I(age > 50)`",
              # Literals, which parse to something that is neither a symbol
              # nor a failure. "NULL" is the one that needs a sentinel rather
              # than is.null() to classify, since str2lang("NULL") RETURNS
              # NULL; a mutation to is.null() survived every other test here.
              "NULL", "`NULL`", "TRUE", "`TRUE`", "1", "`1`")
  # Each entry: a label spelling and the same VARIABLE written bare. Nothing
  # else may merge.
  allow <- list(c("age", "`age`"), c("_X1", "`_X1`"), c("a b", "`a b`"))
  allowed <- function(a, b) {
    any(vapply(allow, function(p) {
      setequal(c(a, b), p)
    }, logical(1)))
  }

  deviations <- character()
  merged_by_allow <- 0L
  for (i in seq_along(shapes)) {
    for (j in seq_along(shapes)) {
      if (j <= i) next
      a <- shapes[i]
      b <- shapes[j]
      main_equal <- identical(a, b)                       # main's comparison
      new_equal  <- identical(.hzr_var_key(a), .hzr_var_key(b))
      if (allowed(a, b)) {
        # An allow-set entry must EARN its place: it must merge now, and it
        # must not have merged before, or it is padding.
        expect_true(new_equal, info = paste("allow-set pair not merged:", a, b))
        expect_false(main_equal, info = paste("allow-set pair was already equal:", a, b))
        merged_by_allow <- merged_by_allow + 1L
      } else if (!identical(new_equal, main_equal)) {
        deviations <- c(deviations, paste0(a, " <> ", b,
                                           " (main=", main_equal,
                                           ", new=", new_equal, ")"))
      }
    }
  }
  expect_identical(merged_by_allow, length(allow)) # every listed pair was seen
  expect_identical(deviations, character())
})
