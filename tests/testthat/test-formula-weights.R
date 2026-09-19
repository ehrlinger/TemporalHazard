# On the formula interface, `weights = <name>` is looked up in `data` first,
# then the calling frame, as the vector interface has done since #151 (#392).
# The formula path used to skip `data`: a column-only name was not found, and
# a name bound both as a column and in the calling frame silently read the
# calling frame's vector. stats::lm() also looks in `data` first, but falls
# back to the formula's environment, not the calling frame.

# The weights column is named fw_wcol, a name no caller plausibly has: the
# column-only test asserts no ambiguity warning, which is correct only while
# no variable of that name is visible from the calling frame.
fw_data <- function() {
  set.seed(1)
  data.frame(t = rexp(40), s = rbinom(40, 1, 0.7), x = rnorm(40),
             fw_wcol = runif(40, 0.5, 1.5))
}

test_that("a weights name that is a column of data is read from data", {
  d <- fw_data()
  want <- hazard(survival::Surv(t, s) ~ x, data = d, weights = d$fw_wcol,
                 dist = "weibull", theta = c(1, 1, 0), fit = TRUE)
  expect_no_warning(
    got <- hazard(survival::Surv(t, s) ~ x, data = d, weights = fw_wcol,
                  dist = "weibull", theta = c(1, 1, 0), fit = TRUE)
  )
  expect_identical(got$fit$theta, want$fit$theta)
  expect_identical(got$fit$objective, want$fit$objective)
})

test_that("a name that is both a column and a calling-frame vector reads the column, and warns", {
  d <- fw_data()
  fw_wcol <- rep(1, 40) # disagrees with d$fw_wcol, so the reading is observable
  column <- hazard(survival::Surv(t, s) ~ x, data = d, weights = d$fw_wcol,
                   dist = "weibull", theta = c(1, 1, 0), fit = TRUE)
  frame <- hazard(survival::Surv(t, s) ~ x, data = d, weights = rep(1, 40),
                  dist = "weibull", theta = c(1, 1, 0), fit = TRUE)
  # The two candidates give different fits, so the check below can fail.
  expect_false(isTRUE(all.equal(column$fit$objective, frame$fit$objective)))
  expect_warning(
    got <- hazard(survival::Surv(t, s) ~ x, data = d, weights = fw_wcol,
                  dist = "weibull", theta = c(1, 1, 0), fit = TRUE),
    "'fw_wcol' \\(weights\\).*both a column of 'data'.*The column was used"
  )
  expect_identical(got$fit$objective, column$fit$objective)
  expect_identical(got$fit$theta, column$fit$theta)
})

test_that("a calling-frame name that is not a column still resolves there, without a warning", {
  d <- fw_data()
  w_local <- d$fw_wcol * 2
  want <- hazard(survival::Surv(t, s) ~ x, data = d, weights = d$fw_wcol * 2,
                 dist = "weibull", theta = c(1, 1, 0), fit = TRUE)
  expect_no_warning(
    got <- hazard(survival::Surv(t, s) ~ x, data = d, weights = w_local,
                  dist = "weibull", theta = c(1, 1, 0), fit = TRUE)
  )
  expect_identical(got$fit$objective, want$fit$objective)
})

test_that("both interfaces read the same vector and warn on the same ambiguous input", {
  # Interface parity is about behaviour, which is what the design decided:
  # the same column read and a warning raised, with one diagnosis. The
  # remedy clause differs on purpose: only the vector interface can drop
  # `data` to reach the calling frame.
  d <- fw_data()
  fw_wcol <- rep(1, 40)
  catch <- function(expr) {
    msg <- NULL
    val <- withCallingHandlers(expr, warning = function(w) {
      if (grepl("both a column of 'data'", conditionMessage(w))) {
        msg <<- conditionMessage(w)
      }
      invokeRestart("muffleWarning")
    })
    list(fit = val, msg = msg)
  }
  via_formula <- catch(hazard(survival::Surv(t, s) ~ x, data = d,
                              weights = fw_wcol, dist = "weibull",
                              theta = c(1, 1, 0), fit = TRUE))
  via_vector <- catch(hazard(data = d, time = t, status = s,
                             x = as.matrix(d["x"]), weights = fw_wcol,
                             dist = "weibull", theta = c(1, 1, 0), fit = TRUE))
  expect_false(is.null(via_formula$msg))
  expect_false(is.null(via_vector$msg))
  diagnosis <- function(m) sub("The column was used\\..*$", "", m)
  expect_identical(diagnosis(via_formula$msg), diagnosis(via_vector$msg))
  expect_identical(via_formula$fit$fit$objective, via_vector$fit$fit$objective)
  expect_match(via_vector$msg, "or omit 'data' to use the calling frame's value\\.$")
  expect_match(via_formula$msg, "not a column of 'data'\\.$")
})

test_that("a namespace-qualified call does not raise a false ambiguity warning", {
  # `::` names a namespace and an export; neither is looked up in data, so a
  # column that happens to share either name is not ambiguous (#401 review).
  # Found on the vector interface too, where it had shipped since #151.
  d <- fw_data()
  names(d)[names(d) == "fw_wcol"] <- "base"
  base <- "a caller variable named base"
  w_ok <- d$base
  want <- hazard(survival::Surv(t, s) ~ x, data = d, weights = w_ok,
                 dist = "weibull", theta = c(1, 1, 0), fit = TRUE)
  expect_no_warning(
    got <- hazard(survival::Surv(t, s) ~ x, data = d,
                  weights = base::abs(w_ok), dist = "weibull",
                  theta = c(1, 1, 0), fit = TRUE)
  )
  expect_identical(got$fit$objective, want$fit$objective)
  expect_no_warning(
    hazard(data = d, time = t, status = s, x = as.matrix(d["x"]),
           weights = base::abs(w_ok), dist = "weibull", theta = c(1, 1, 0),
           fit = TRUE)
  )
})

test_that("a named argument's tag is not looked up, so it cannot be ambiguous (#402)", {
  # `n` in stats::runif(n = 40) names runif()'s argument; it is never looked
  # up, in data or anywhere else. A column `n` beside a caller variable `n`
  # is therefore not ambiguous, on either interface. A collector that also
  # harvested argument tags would warn here and name a column the fit never
  # read (#401 review).
  d <- fw_data()
  d$n <- 2
  n <- 1
  set.seed(7)
  want_f <- hazard(survival::Surv(t, s) ~ x, data = d,
                   weights = stats::runif(40), dist = "weibull",
                   theta = c(1, 1, 0), fit = TRUE)
  set.seed(7)
  expect_no_warning(
    got_f <- hazard(survival::Surv(t, s) ~ x, data = d,
                    weights = stats::runif(n = 40), dist = "weibull",
                    theta = c(1, 1, 0), fit = TRUE)
  )
  expect_identical(got_f$fit$objective, want_f$fit$objective)
  expect_no_warning(
    hazard(data = d, time = t, status = s, x = as.matrix(d["x"]),
           weights = stats::runif(n = 40), dist = "weibull",
           theta = c(1, 1, 0), fit = FALSE)
  )
})

test_that("an ambiguous argument inside a namespace-qualified call still warns", {
  # The control: skipping the namespace and export names must not skip the
  # call's own arguments, which ARE looked up in data.
  d <- fw_data()
  fw_wcol <- rep(1, 40)
  base <- "a caller variable named base"
  d$base <- 1
  expect_warning(
    hazard(survival::Surv(t, s) ~ x, data = d, weights = base::abs(fw_wcol),
           dist = "weibull", theta = c(1, 1, 0), fit = TRUE),
    "'fw_wcol' \\(weights\\): the name is both a column"
  )
})

test_that("the ambiguity warning's advice works when followed, on both interfaces", {
  # It used to say "Write data$<name>", and in the caller's frame `data` is
  # usually utils::data(), so following it errored (#401 review). The advice
  # now names the caller's own data argument when that is a plain symbol.
  d <- fw_data()
  advice_value <- function(msg, env) {
    adv <- regmatches(msg, regexpr("Write [^ ]+ for the column", msg))
    expect_length(adv, 1L)
    expr <- sub("<name>", "fw_wcol", sub(" for the column$", "", sub("^Write ", "", adv)))
    eval(parse(text = expr), envir = env)
  }
  caller <- function(iface) {
    fw_wcol <- rep(1, 40)
    msg <- NULL
    withCallingHandlers(
      if (iface == "formula") {
        hazard(survival::Surv(t, s) ~ x, data = d, weights = fw_wcol,
               dist = "weibull", theta = c(1, 1, 0))
      } else {
        hazard(data = d, time = t, status = s, x = as.matrix(d["x"]),
               weights = fw_wcol, dist = "weibull", theta = c(1, 1, 0))
      },
      warning = function(w) {
        if (grepl("both a column", conditionMessage(w))) msg <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      }
    )
    advice_value(msg, environment())
  }
  # Following the advice reads the column, not the caller's variable.
  expect_identical(caller("formula"), d$fw_wcol)
  expect_identical(caller("vector"), d$fw_wcol)
})

test_that("the ambiguity warning falls back when data is not a plain name", {
  # An inline expression or magrittr's `.` cannot be written as a prefix, so
  # the advice says to use the data frame passed as `data` instead.
  d <- fw_data()
  fw_wcol <- rep(1, 40)
  msgs <- function(expr) {
    m <- character()
    withCallingHandlers(expr, warning = function(w) {
      m <<- c(m, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    m[grepl("both a column", m)]
  }
  inline <- msgs(hazard(survival::Surv(t, s) ~ x, data = transform(d),
                        weights = fw_wcol, dist = "weibull", theta = c(1, 1, 0)))
  . <- d
  dot <- msgs(hazard(survival::Surv(t, s) ~ x, data = ., weights = fw_wcol,
                     dist = "weibull", theta = c(1, 1, 0)))
  for (m in list(inline, dot)) {
    expect_length(m, 1L)
    expect_match(m, "the data frame passed as 'data'")
    expect_false(grepl("data$", m, fixed = TRUE))
  }
})

# The warning says "The column was used" only when the argument is the name
# alone. Inside a larger expression the name may never be evaluated (a lazy
# function argument), may be rebound first (a loop variable, an assignment)
# or may be evaluated elsewhere (with()), so the warning does not say which
# value was read (#401 review).
# Each case compares the fit's objective with the two candidate vectors,
# because a message test passes whether or not the column was used.
fw_fit <- function(d, wexpr, env = parent.frame()) {
  msgs <- character()
  fit <- withCallingHandlers(
    eval(bquote(hazard(survival::Surv(t, s) ~ x, data = .(d),
                       weights = .(wexpr), dist = "weibull",
                       theta = c(1, 1, 0), fit = TRUE)), env),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  amb <- msgs[grepl("both a column", msgs)]
  list(objective = fit$fit$objective, ambiguity = length(amb),
       definite = sum(grepl("The column was used", amb, fixed = TRUE)),
       unchecked = sum(grepl("which value it read, if any, is not checked",
                             amb, fixed = TRUE)))
}

test_that("a name inside a larger expression warns without saying which value was used", {
  d <- fw_data()
  fw_wcol <- rep(1, 40)
  column <- fw_fit(d, quote(d$fw_wcol))$objective
  frame <- fw_fit(d, quote(rep(1, 40)))$objective
  expect_false(isTRUE(all.equal(frame, column)))
  here <- environment()
  # The objective goes either way; the message must be true for both.
  shapes <- list(
    list(quote((function(a) rep(1, 40))(fw_wcol)), frame),
    list(quote((function(fw_wcol) fw_wcol)(rep(1, 40))), frame),
    list(quote({
      for (fw_wcol in list(rep(1, 40))) NULL
      fw_wcol
    }), frame),
    list(quote({
      fw_wcol <- rep(1, 40)
      fw_wcol
    }), frame),
    list(quote(with(list(fw_wcol = rep(1, 40)), fw_wcol)), frame),
    # `here` is bound before the call: environment() inside the argument
    # would return the data mask itself.
    list(quote(local(fw_wcol, envir = here)), frame),
    list(quote(local(fw_wcol)), column),
    list(quote((function() fw_wcol)()), column),
    list(quote(sapply(1:40, function(i, k = fw_wcol) k[i])), column),
    list(quote(fw_wcol * 1), column)
  )
  for (s in shapes) {
    got <- fw_fit(d, s[[1]])
    expect_identical(got$objective, s[[2]], label = deparse(s[[1]]))
    expect_identical(got$ambiguity, 1L, label = deparse(s[[1]]))
    expect_identical(got$definite, 0L, label = deparse(s[[1]]))
    expect_identical(got$unchecked, 1L, label = deparse(s[[1]]))
  }
  # The name alone is the one case where the column is known to be read.
  got <- fw_fit(d, quote(fw_wcol))
  expect_identical(got$objective, column)
  expect_identical(got$definite, 1L)
})

test_that("an empty index in a masked argument does not stop the fit", {
  # `m[, 2]` carries the empty symbol, which is no name to look up; it
  # errored "invalid first argument" (vector interface since 1.2.2).
  d <- fw_data()
  wm <- cbind(1, d$fw_wcol)
  want <- hazard(survival::Surv(t, s) ~ x, data = d, weights = d$fw_wcol,
                 dist = "weibull", theta = c(1, 1, 0), fit = TRUE)
  expect_no_warning(
    got <- hazard(survival::Surv(t, s) ~ x, data = d, weights = wm[, 2],
                  dist = "weibull", theta = c(1, 1, 0), fit = TRUE)
  )
  expect_identical(got$fit$objective, want$fit$objective)
  tm <- cbind(d$t, 0)
  want_v <- hazard(time = d$t, status = d$s, data = d, dist = "weibull",
                   theta = c(1, 1), fit = TRUE)
  expect_no_warning(
    got_v <- hazard(time = tm[, 1], status = d$s, data = d,
                    dist = "weibull", theta = c(1, 1), fit = TRUE)
  )
  expect_identical(got_v$fit$objective, want_v$fit$objective)
})

test_that("an expression the ambiguity check cannot walk is fitted unchecked", {
  # A generated expression can nest deeper than the walk can recurse. The
  # check is a diagnostic, so its failure must not stop a fit that main
  # made (#401 review). The walk's failure is simulated: the depth at which
  # recursion fails differs by platform.
  d <- fw_data()
  fw_wcol <- rep(1, 40)
  want <- hazard(survival::Surv(t, s) ~ x, data = d, weights = d$fw_wcol,
                 dist = "weibull", theta = c(1, 1, 0), fit = TRUE)
  local_mocked_bindings(
    .hzr_ambiguity_symbols = function(e) stop("evaluation nested too deeply")
  )
  expect_no_warning(
    got <- hazard(survival::Surv(t, s) ~ x, data = d, weights = fw_wcol,
                  dist = "weibull", theta = c(1, 1, 0), fit = TRUE)
  )
  expect_identical(got$fit$objective, want$fit$objective)
})
