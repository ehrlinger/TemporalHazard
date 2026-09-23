# How `force_in`, `force_out` and a character `scope` name a variable
# (#437, #442).
#
# terms() backquotes a label whose variable is not a syntactic name: `_X1`
# for the column `_X1`. The three arguments are documented as variables, and
# the translator emits bare names, so a pin on "_X1" never met its own
# candidate and was silently ignored (#437).
#
# The first fix PARSED each string and guessed what it named. Four review
# rounds each found a new silent wrong answer in that guess (#442): a
# candidate scored twice; a literal column `age:mal` swallowed by the
# interaction; "age", "age " and "age # x" merged; a column literally named
# `x` (with backticks) merged with `x`; and a reserved word such as "TRUE"
# parsed to a constant, so a pin on the column TRUE was ignored. The same
# text is a raw column NAME at some sites and a term LABEL at others, and no
# function of the string alone is right at both.
#
# So a name is resolved by LOOKUP, once, when hzr_stepwise() is called:
#   1. exactly a column of `data`      -> that column;
#   2. else exactly a term label of the model or `scope` -> that term;
#   3. else a warning naming it, and it is ignored.
#
# The tests read what the user sees: which terms a screen dropped, and which
# candidates hzr_stepwise() OFFERED to its first forward step. The second is
# captured by wrapping the real enumeration, so the resolution under test is
# the one hzr_stepwise() performs, not a helper called in isolation.

nsn_avc <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  stats::na.omit(avc[, c("int_dead", "dead", "age", "mal")])
}

# Rename `age`. It is the term a free backward screen drops first here, so a
# pin on the renamed column has to REDIRECT the drop onto `mal`: an ignored
# pin and an honoured one give different answers (#437).
nsn_rename_age <- function(nm) {
  d <- nsn_avc()
  names(d)[names(d) == "age"] <- nm
  d
}

nsn_fit <- function(d, rhs = "1", theta = c(0.1, 1)) {
  suppressWarnings(hazard(
    stats::as.formula(paste("survival::Surv(int_dead, dead) ~", rhs)),
    data = d, dist = "weibull", theta = theta, fit = TRUE
  ))
}

# The candidates hzr_stepwise() offers its first forward step, as spelled in
# the candidate list. The real enumeration runs; the step is then handed an
# empty list, so nothing is scored or refit and the screen stops.
nsn_offered <- function(fit, data, ...) {
  seen <- NULL
  real <- .hzr_stepwise_candidates
  local_mocked_bindings(.hzr_stepwise_candidates = function(...) {
    out <- real(...)
    if (is.null(seen)) seen <<- vapply(out, function(c) c$var, character(1))
    list()
  })
  hzr_stepwise(fit, data = data, direction = "forward", criterion = "wald",
               trace = FALSE, ...)
  # NULL would mean the enumeration never ran, which must not read as "offered
  # nothing".
  expect_false(is.null(seen))
  seen
}

# The terms a backward screen dropped, in order. slstay is tiny, so every
# term that is not pinned is dropped.
nsn_dropped <- function(fit, data, ...) {
  sw <- suppressWarnings(hzr_stepwise(fit, data = data, direction = "backward",
                                      criterion = "wald", slstay = 1e-9,
                                      trace = FALSE, ...))
  sw$steps$variable[sw$steps$action == "drop"]
}

nsn_agemal <- function() {
  d0 <- nsn_avc()
  d <- data.frame(d0, `age:mal` = as.numeric(scale(d0$age)) * 0.5 + 1,
                  check.names = FALSE)
  # The column is genuinely not the product, or the tests prove nothing.
  expect_false(isTRUE(all.equal(d[["age:mal"]], d$age * d$mal)))
  d
}

# The terms of the model a screen FINISHED with. Read off the fitted object
# the screen returns, so every step, refit included, has run.
nsn_final_terms <- function(sw) {
  attr(stats::terms(stats::formula(sw$call$formula)), "term.labels")
}

nsn_screen <- function(fit, data, ...) {
  suppressWarnings(hzr_stepwise(fit, data = data, criterion = "wald",
                                trace = FALSE, ...))
}

# --- real screens, run to completion ---------------------------------------
#
# nsn_offered() below stops each screen before anything is scored or refit,
# which is right for the resolution layer and blind to everything after it.
# These run the whole screen and read the final model.

test_that("a real screen keeps a pinned `_X1` in the final model (#437)", {
  skip_on_cran() # full screens
  d <- nsn_rename_age("_X1")
  fit <- nsn_fit(d, "`_X1` + mal", c(0.1, 1, 0, 0))
  sw_free <- nsn_screen(fit, d, direction = "backward", slstay = 1e-9)
  expect_identical(nsn_final_terms(sw_free), character()) # known positive
  sw_pin <- nsn_screen(fit, d, direction = "backward", slstay = 1e-9,
                       force_in = "_X1")
  expect_identical(nsn_final_terms(sw_pin), "`_X1`")
})

test_that("a real screen never enters a force_out column (#437, #442)", {
  skip_on_cran() # full screens
  for (nm in c("_X1", "TRUE")) {
    d <- nsn_rename_age(nm)
    base <- nsn_fit(d)
    sc <- stats::as.formula(paste0("~ `", nm, "` + mal"))
    lab <- paste0("`", nm, "`")
    sw_free <- nsn_screen(base, d, direction = "forward", slentry = 0.99,
                          scope = sc)
    expect_setequal(nsn_final_terms(sw_free), c(lab, "mal")) # known positive
    sw_out <- nsn_screen(base, d, direction = "forward", slentry = 0.99,
                         scope = sc, force_out = nm)
    expect_identical(nsn_final_terms(sw_out), "mal", info = nm)
  }
})

test_that("a screen over a literal `age:mal` column completes (#442)", {
  skip_on_cran() # full screens
  # Regression found by review of 36214f79: the candidate was compared by
  # the COLUMN's label while the refit entered the term its spelling pastes
  # to, so the next iteration re-offered it and the no-op refit stopped the
  # screen with an error. 71277ff8 completed with one entry. WHICH term the
  # refit enters for such a column is a known limitation of the refit
  # (#449), so this asserts only what does not depend on it: the screen
  # finishes, enters once, fails no refit, and grows the model by one term.
  d <- nsn_agemal()
  fit <- nsn_fit(d, "age + mal", c(0.1, 1, 0, 0))
  for (sc in list(NULL, "age:mal")) {
    sw <- NULL
    expect_no_error(
      sw <- nsn_screen(fit, d, direction = "forward", slentry = 0.99,
                       scope = sc)
    )
    expect_identical(sum(sw$steps$action == "enter"), 1L)
    expect_identical(sw$criteria$n_refit_failures, 0L)
    expect_length(nsn_final_terms(sw), 3L)
  }
})

# --- #437: a bare non-syntactic name ---------------------------------------

test_that("force_in = '_X1' pins the column `_X1` (#437)", {
  skip_on_cran() # backward screens
  d <- nsn_rename_age("_X1")
  fit <- nsn_fit(d, "`_X1` + mal", c(0.1, 1, 0, 0))
  # Known positive: unpinned, `_X1` is the first term dropped.
  expect_identical(nsn_dropped(fit, d)[1L], "`_X1`")
  # Pinned bare, as documented and as the translator emits it.
  expect_identical(nsn_dropped(fit, d, force_in = "_X1"), "mal")
  # The label spelling names the same term.
  expect_identical(nsn_dropped(fit, d, force_in = "`_X1`"), "mal")
})

test_that("a character scope does not re-offer a variable in the model (#437)", {
  skip_on_cran() # fits
  d <- nsn_rename_age("_X1")
  base <- nsn_fit(d)
  fit <- nsn_fit(d, "`_X1` + mal", c(0.1, 1, 0, 0))
  expect_setequal(nsn_offered(base, d, scope = c("mal", "_X1")),
                  c("mal", "_X1")) # known positive
  expect_identical(nsn_offered(fit, d, scope = c("mal", "_X1")), character())
})

test_that("force_out keeps a non-syntactic column out, either spelling (#437)", {
  skip_on_cran() # a fit
  d <- nsn_rename_age("_X1")
  base <- nsn_fit(d)
  expect_true("_X1" %in% nsn_offered(base, d)) # known positive
  expect_identical(nsn_offered(base, d, force_out = "_X1"), "mal")
  expect_identical(nsn_offered(base, d, force_out = "`_X1`"), "mal")
  # Through a formula scope, whose candidates are labels.
  expect_setequal(nsn_offered(base, d, scope = ~ `_X1` + mal),
                  c("`_X1`", "mal"))
  expect_identical(nsn_offered(base, d, scope = ~ `_X1` + mal,
                               force_out = "_X1"), "mal")
})

test_that("entry needs the QUOTED spelling; a bare one fails loudly (#437)", {
  skip_on_cran() # three forward screens
  # Pins the NEWS paragraph that separates MATCHING from ENTRY. Matching a
  # non-syntactic variable is not being able to ADD one: the refit pastes
  # the candidate's spelling into the formula text. For `_X1` the quoted
  # spelling parses and the bare one does not, and that failure is LOUD
  # (#441). (A spelling that parses to a DIFFERENT term is #449.)
  d <- nsn_rename_age("_X1")
  base <- nsn_fit(d)
  fwd <- function(...) {
    suppressWarnings(hzr_stepwise(base, data = d, direction = "forward",
                                  criterion = "aic", trace = FALSE, ...))
  }
  entered <- function(sw) sw$steps$variable[sw$steps$action == "enter"]

  sw_f <- fwd(scope = ~ `_X1` + mal)
  expect_true("`_X1`" %in% entered(sw_f))
  expect_identical(sw_f$criteria$n_refit_failures, 0L)

  sw_c <- fwd(scope = c("_X1", "mal"))
  expect_false(any(entered(sw_c) %in% c("_X1", "`_X1`")))
  expect_gt(sw_c$criteria$n_refit_failures, 0L)

  expect_true("`_X1`" %in% entered(fwd(scope = c("`_X1`", "mal"))))
})

test_that("multiphase: force_out written bare reaches the label (#437)", {
  skip_on_cran() # a multiphase fit
  d <- nsn_rename_age("_X1")
  mp <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(
      early    = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                           fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE, control = list(n_starts = 2L, maxit = 500L)
  ))
  sc <- list(early = ~ `_X1` + mal, constant = NULL)
  expect_setequal(nsn_offered(mp, d, scope = sc), c("`_X1`", "mal"))
  expect_identical(nsn_offered(mp, d, scope = sc, force_out = "_X1"), "mal")
})

test_that("a syntactic name is unaffected (#437)", {
  skip_on_cran() # backward screens
  d <- nsn_avc()
  fit <- nsn_fit(d, "age + mal", c(0.1, 1, 0, 0))
  expect_setequal(nsn_dropped(fit, d), c("age", "mal"))
  expect_identical(nsn_dropped(fit, d, force_in = "age"), "mal")
  expect_identical(nsn_dropped(fit, d, force_in = c("age", "mal")),
                   character())
})

# --- #442 defect 1: de-duplication ----------------------------------------

test_that("a scope naming a variable twice offers it once (#442, 1)", {
  skip_on_cran() # a fit
  d <- nsn_avc()
  base <- nsn_fit(d)
  expect_identical(nsn_offered(base, d, scope = c("age", "age", "mal")),
                   c("age", "mal"))
  # A column and its own label are one variable.
  d2 <- nsn_rename_age("_X1")
  base2 <- nsn_fit(d2)
  expect_identical(nsn_offered(base2, d2, scope = c("_X1", "`_X1`", "mal")),
                   c("_X1", "mal"))
})

# --- #442 defect 2: a literal column named like an interaction ------------

test_that("a literal `age:mal` column is offered beside the interaction (#442, 2)", {
  skip_on_cran() # a fit
  d <- nsn_agemal()
  fit <- nsn_fit(d, "age * mal", c(0.1, 1, 0, 0, 0))
  # The interaction is in the model; the literal column is not.
  # scope = NULL offers columns: the literal one must be among them.
  expect_identical(nsn_offered(fit, d), "age:mal")
  # A formula scope carries its label.
  expect_identical(nsn_offered(fit, d, scope = ~ `age:mal`), "`age:mal`")
  # A character scope: the COLUMN wins, so the model's interaction does not
  # swallow it.
  expect_identical(nsn_offered(fit, d, scope = "age:mal"), "age:mal")
})

test_that("a string that is a column AND a term label names the column", {
  skip_on_cran() # a fit
  # The precedence, pinned where the two readings give different answers.
  d <- nsn_agemal()
  base <- nsn_fit(d)
  sc <- ~ `age:mal` + age:mal
  expect_setequal(nsn_offered(base, d, scope = sc),
                  c("`age:mal`", "age:mal")) # known positive
  # "age:mal" is the literal column; the interaction stays offered.
  expect_identical(nsn_offered(base, d, scope = sc, force_out = "age:mal"),
                   "age:mal")
  # The interaction is reachable only through a string that is not a
  # column: here, with no such column, it resolves to the term.
  d0 <- nsn_avc()
  base0 <- nsn_fit(d0)
  expect_identical(nsn_offered(base0, d0, scope = ~ age + age:mal,
                               force_out = "age:mal"), "age")
})

# --- #442 defect 3: whitespace and comments -------------------------------

test_that("'age', 'age ' and 'age # x' stay distinct (#442, 3)", {
  skip_on_cran() # a fit
  d0 <- nsn_avc()
  d <- data.frame(d0, `age ` = d0$age * 2, `age # x` = d0$age + 1,
                  check.names = FALSE)
  base <- nsn_fit(d)
  expect_setequal(nsn_offered(base, d), c("age", "mal", "age ", "age # x"))
  expect_setequal(nsn_offered(base, d, scope = c("age", "age ")),
                  c("age", "age "))
  # force_out excludes the column named, and only that one.
  expect_setequal(nsn_offered(base, d, force_out = "age"),
                  c("mal", "age ", "age # x"))
  expect_setequal(nsn_offered(base, d, force_out = "age "),
                  c("age", "mal", "age # x"))
})

test_that("a near-miss spelling warns rather than matching (#442, 3)", {
  skip_on_cran() # a fit
  d <- nsn_avc()
  base <- nsn_fit(d)
  # No column "age " exists here, so it names nothing, and says so.
  expect_warning(
    off <- nsn_offered(base, d, force_out = "age "),
    "\"age \"", fixed = TRUE
  )
  expect_setequal(off, c("age", "mal"))
})

# --- #442 defect 4: a column literally named `x` --------------------------

test_that("a column named `x` (with backticks) is not x (#442, 4)", {
  skip_on_cran() # a fit
  d <- nsn_rename_age("x")
  d <- data.frame(d, `\`x\`` = d$x * 2 + 1, check.names = FALSE)
  base <- nsn_fit(d)
  expect_setequal(nsn_offered(base, d), c("x", "mal", "`x`"))
  expect_setequal(nsn_offered(base, d, force_out = "x"), c("mal", "`x`"))
  expect_setequal(nsn_offered(base, d, force_out = "`x`"), c("x", "mal"))
})

# --- #442 defect 5: reserved words and numbers ----------------------------

test_that("force_in honours a column named like a constant (#442, 5)", {
  skip_on_cran() # backward screens, several fits
  for (nm in c("TRUE", "NULL", "NA", "Inf", "next", "1")) {
    d <- nsn_rename_age(nm)
    lab <- paste0("`", nm, "`")
    fit <- nsn_fit(d, paste(lab, "+ mal"), c(0.1, 1, 0, 0))
    # Known positive: unpinned, it is the first term dropped.
    expect_identical(nsn_dropped(fit, d)[1L], lab, info = nm)
    expect_identical(nsn_dropped(fit, d, force_in = nm), "mal", info = nm)
  }
})

test_that("force_out and scope honour a column named TRUE (#442, 5)", {
  skip_on_cran() # a fit
  d <- nsn_rename_age("TRUE")
  base <- nsn_fit(d)
  expect_setequal(nsn_offered(base, d), c("TRUE", "mal")) # known positive
  expect_identical(nsn_offered(base, d, force_out = "TRUE"), "mal")
  expect_identical(nsn_offered(base, d, scope = "TRUE"), "TRUE")
})

# --- an unresolved name is never silent -----------------------------------

test_that("an unresolved name warns, naming it, for each argument", {
  skip_on_cran() # a fit
  d <- nsn_avc()
  base <- nsn_fit(d)
  expect_warning(off <- nsn_offered(base, d, force_in = "nosuch_in"),
                 "`force_in`.*nosuch_in")
  expect_setequal(off, c("age", "mal"))
  expect_warning(off <- nsn_offered(base, d, force_out = "nosuch_out"),
                 "`force_out`.*nosuch_out")
  expect_setequal(off, c("age", "mal"))
  expect_warning(off <- nsn_offered(base, d, scope = c("age", "nosuch_sc")),
                 "`scope`.*nosuch_sc")
  expect_identical(off, "age")
  # A resolved name is quiet, and a term label over columns of `data`
  # resolves in a character scope although it is not a column.
  expect_no_warning(nsn_offered(base, d, force_out = "age"))
  expect_no_warning(off <- nsn_offered(base, d,
                                       scope = c("age", "log(age)", "age:mal")))
  expect_identical(off, c("age", "log(age)", "age:mal"))
})

test_that("hzr_bootstrap() passes the unresolved-name warning on", {
  skip_on_cran() # a bootstrap
  d <- nsn_avc()
  fit <- suppressWarnings(hazard(survival::Surv(int_dead, dead) ~ age,
                                 data = d, dist = "weibull",
                                 theta = c(0.1, 1, 0), fit = TRUE))
  expect_warning(
    suppressMessages(hzr_bootstrap(fit, n_boot = 2L, seed = 1L,
                                   scope = ~ age + mal, criterion = "wald",
                                   force_out = "nosuch_boot")),
    "nosuch_boot"
  )
})

# --- the resolution itself ------------------------------------------------

test_that("each column resolves to its own terms() label, and only it", {
  # Replaces #442's differential test, which compared a parsed key against a
  # string comparison and, through its allow-set, ASSERTED that "TRUE" and
  # "`TRUE`" were different variables. The property now is the design's:
  # a column name resolves to the label terms() gives that column, distinct
  # columns never merge, and a column's label names the column.
  nms <- c("age", "_X1", "a b", "a`b", "a\\b", "a\nb", "a\tb", ".x", "x",
           "`x`", "age ", "age # x", "age:mal", "log(age)",
           "TRUE", "FALSE", "NULL", "NA", "Inf", "NaN", "NA_integer_",
           "NA_real_", "NA_character_", "next", "break", "if", "1", "1e3")
  d <- as.data.frame(stats::setNames(
    lapply(seq_along(nms), function(i) as.numeric(i)), nms
  ), check.names = FALSE)
  expect_identical(names(d), nms) # the frame really carries every name
  labels <- attr(stats::terms(stats::reformulate(".", response = NULL),
                              data = d), "term.labels")
  expect_length(labels, length(nms))

  by_name <- .hzr_resolve_names(nms, d, arg = "`x`")
  expect_identical(by_name$id, labels)
  expect_identical(anyDuplicated(by_name$id), 0L)

  by_label <- .hzr_resolve_names(labels, d, arg = "`x`")
  expect_identical(by_label$id, labels)

  # A string that is neither a column nor a label is not guessed at.
  for (s in c("`age`", "age  ", "\"age\"", "age#x")) {
    expect_warning(r <- .hzr_resolve_names(s, d, arg = "`x`"),
                   encodeString(s, quote = "\""), fixed = TRUE)
    expect_length(r$id, 0L)
  }
})
