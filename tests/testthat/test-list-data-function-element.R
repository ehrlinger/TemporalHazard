# `data` masks the calling frame while hazard() evaluates its arguments and
# the formula's response. R's function lookup walks past every binding that
# is not a function, so a FUNCTION-valued element of `data` was called in
# place of the function the expression names -- silently, changing the fit
# (#420, since 1.2.2). hazard() now refuses such an element, as stats::lm()
# does. A numeric element of the same name was never consulted, and a
# list-column is a list, which the lookup skips: both still fit.
#
# The nrow sweep below is not decoration. An earlier version of this fix
# guarded only the vector path, on the measured belief that the formula path
# was immune. It is immune only at nrow >= 2, and for an accidental reason:
# .hzr_numeric_frame_values() replicates each column to nrow and dies on a
# function while doing it. At nrow == 1 there is nothing to replicate, the
# function survives, and the formula path reached it. At n = 0 it dies again
# with a THIRD message, "object of type 'closure' is not subsettable". One
# apparent guard, three regimes, absent in the middle one -- so every case
# here runs at n = 0, 1, 2 and 40. A guard that works as a side effect of
# size is absent exactly where the side effect does not occur, and passing
# at both 0 and 40 would still have missed 1.

ldf_flag <- function() new.env(parent = emptyenv())

ldf_fn <- function(env, val) {
  function(...) {
    env$called <- TRUE
    val
  }
}

# A malformed-but-constructible frame: data.frame(), `$<-` and `[[<-` all
# refuse a function column, and structure() does not.
ldf_frame <- function(n, nm, fn) {
  cols <- list(tt = rep(2.5, n), ss = rep(1L, n), age = rep(55, n))
  cols[[nm]] <- fn
  structure(cols, class = "data.frame", row.names = seq_len(n))
}

ldf_list <- function(n, nm, fn) {
  cols <- list(tt = rep(2.5, n), ss = rep(1L, n), age = rep(55, n))
  cols[[nm]] <- fn
  cols
}

# The refusal's own words. The pre-existing upstream guard says "attempt to
# replicate an object of type ...", so every expectation below keys on THIS
# text: an assertion satisfied by the old guard would pass on main and pin
# nothing (it did, in the first version of this file).
ldf_msg <- "holds a function named"

test_that("a function-valued element is refused on BOTH interfaces, at every nrow", {
  for (n in c(0L, 1L, 2L, 40L)) {
    env <- ldf_flag()
    env$called <- FALSE
    fn <- ldf_fn(env, rep(7, n))
    for (shape in c("frame", "list")) {
      d <- if (shape == "frame") ldf_frame(n, "rep", fn) else ldf_list(n, "rep", fn)
      expect_error(
        hazard(time = tt, status = ss, data = d, weights = rep(1, n),
               dist = "weibull", theta = c(1, 1), fit = TRUE),
        ldf_msg,
        info = paste("vector interface,", shape, "n =", n)
      )
      # The formula path reaches `data` twice: the Surv() response and the
      # `weights` lookup. Both are behind the same refusal now.
      expect_error(
        hazard(survival::Surv(tt, ss) ~ age, data = d, weights = rep(1, n),
               dist = "weibull", theta = c(1, 1, 0), fit = TRUE),
        ldf_msg,
        info = paste("formula interface,", shape, "n =", n)
      )
    }
    # Refused before anything is evaluated, so the masking function never ran.
    expect_false(env$called)
  }
})

test_that("a masked function in the Surv() response is refused, not silently used", {
  # The worst shape: the RESPONSE is replaced. On 1.2.11 a 1-row frame
  # carrying `round` made Surv(round(tt), ss) read 99 where round(2.5) is 2.
  env <- ldf_flag()
  env$called <- FALSE
  d <- ldf_frame(1L, "round", ldf_fn(env, 99))
  expect_error(
    hazard(survival::Surv(round(tt), ss) ~ age, data = d, dist = "weibull",
           theta = c(1, 1, 0), fit = TRUE),
    ldf_msg
  )
  expect_false(env$called)
})

test_that("the refusal names every function element, once each, and says what to do", {
  d <- list(tt = rep(2.5, 3), ss = rep(1L, 3), rep = base::rep, seq = base::seq)
  msg <- tryCatch(
    hazard(time = tt, status = ss, data = d, dist = "weibull",
           theta = c(1, 1), fit = TRUE),
    error = conditionMessage
  )
  expect_match(msg, "'rep'", fixed = TRUE)
  expect_match(msg, "'seq'", fixed = TRUE)
  # A duplicated name is named once, not once per element.
  d2 <- list(tt = rep(2.5, 3), ss = rep(1L, 3), rep = base::rep, rep = base::seq)
  msg2 <- tryCatch(
    hazard(time = tt, status = ss, data = d2, dist = "weibull",
           theta = c(1, 1), fit = TRUE),
    error = conditionMessage
  )
  expect_equal(lengths(regmatches(msg2, gregexpr("'rep'", msg2, fixed = TRUE))), 1L)
})

test_that("an element no expression can look up is not refused", {
  # `eval()` cannot reach an element with no name, so it never had the
  # defect. Refusing it would be a false refusal, and the first version of
  # this fix made one.
  d <- list(rep(2.5, 3), rep(1L, 3), base::rep)
  fit <- hazard(time = c(2.5, 3.5, 4.5), status = c(1L, 1L, 0L), data = d,
                dist = "weibull", theta = c(1, 1), fit = TRUE)
  expect_true(is.finite(fit$fit$objective))
  # Partially named: the unnamed function element is still unreachable.
  d2 <- list(tt = c(2.5, 3.5, 4.5), ss = c(1L, 1L, 0L), base::rep)
  fit2 <- hazard(time = tt, status = ss, data = d2, dist = "weibull",
                 theta = c(1, 1), fit = TRUE)
  expect_true(is.finite(fit2$fit$objective))
})

test_that("a numeric element of the same name still fits, unchanged", {
  set.seed(1)
  n <- 40
  tt <- rexp(n)
  ss <- rbinom(n, 1, 0.7)
  ref <- hazard(time = tt, status = ss, weights = rep(1, n), dist = "weibull",
                theta = c(1, 1), fit = TRUE)
  num_list <- hazard(time = t, status = s,
                     data = list(t = tt, s = ss, rep = seq_len(n)),
                     weights = rep(1, n), dist = "weibull",
                     theta = c(1, 1), fit = TRUE)
  num_col <- hazard(time = t, status = s,
                    data = data.frame(t = tt, s = ss, rep = seq_len(n)),
                    weights = rep(1, n), dist = "weibull",
                    theta = c(1, 1), fit = TRUE)
  expect_equal(num_list$fit$objective, ref$fit$objective)
  expect_equal(num_col$fit$objective, ref$fit$objective)
})

test_that("a data frame with a list-column of functions still fits", {
  # A list-column is a LIST, so the function-call lookup walks past it. This
  # is the shape closest to the refusal that must NOT be refused.
  set.seed(1)
  n <- 40
  df <- data.frame(t = rexp(n), s = rbinom(n, 1, 0.7))
  df$rep <- replicate(n, function(...) stop("must not be called"),
                      simplify = FALSE)
  fit <- hazard(time = t, status = s, data = df, weights = rep(1, n),
                dist = "weibull", theta = c(1, 1), fit = TRUE)
  expect_true(is.finite(fit$fit$objective))
})

test_that("the deliberate helper in data now errors, with a workaround that works", {
  d <- list(tt = c(2.5, 3.5, 4.5), ss = c(1L, 1L, 0L), f = function(x) x * 2)
  msg <- tryCatch(
    hazard(time = f(tt), status = ss, data = d, dist = "weibull",
           theta = c(1, 1), fit = TRUE),
    error = conditionMessage
  )
  expect_match(msg, "'f'", fixed = TRUE)
  expect_match(msg, "calling environment", fixed = TRUE)
  helper <- function(x) x * 2
  fit <- hazard(time = helper(tt), status = ss,
                data = list(tt = c(2.5, 3.5, 4.5), ss = c(1L, 1L, 0L)),
                dist = "weibull", theta = c(1, 1), fit = TRUE)
  expect_true(is.finite(fit$fit$objective))
})


test_that("a `data` that is not a list is still refused for its SHAPE", {
  # The guard iterates `data`, so it must not be the first thing to touch an
  # object that is not list-like: `vapply()` would answer for it, with a
  # coercion error for an S4 object, or -- worse -- by finding a "function
  # element" in an environment and telling the user to remove it, which does
  # not make an environment acceptable `data`. Both were regressions against
  # 1.2.11's message, and neither is about functions at all.
  tt <- c(2.5, 3.5, 4.5)
  ss <- c(1L, 1L, 0L)
  shape <- "'data' must be a data frame or a list."
  setClass("LdfS4", representation(a = "numeric"))
  on.exit(removeClass("LdfS4"), add = TRUE)
  expect_error(
    hazard(time = tt, status = ss, data = new("LdfS4", a = 1),
           dist = "weibull", theta = c(1, 1), fit = TRUE),
    shape, fixed = TRUE
  )
  e <- new.env()
  e$g <- mean
  expect_error(
    hazard(time = tt, status = ss, data = e, dist = "weibull",
           theta = c(1, 1), fit = TRUE),
    shape, fixed = TRUE
  )
  # Unchanged since 1.2.11, and kept here so the three stay together.
  expect_error(
    hazard(time = tt, status = ss, data = matrix(1:4, 2), dist = "weibull",
           theta = c(1, 1), fit = TRUE),
    shape, fixed = TRUE
  )
})

test_that("a function nested inside a list element is not refused", {
  # `eval()` cannot reach it as a function, so it never had the defect.
  # The columns are NOT named after variables bound here: a local of the
  # same name raises the masked-ambiguity warning (#401), which would leak
  # into the suite's warning multiset for a reason unrelated to this test.
  fit <- hazard(time = a, status = b,
                data = list(a = c(2.5, 3.5, 4.5), b = c(1L, 1L, 0L),
                            inner = list(f = mean)),
                dist = "weibull", theta = c(1, 1), fit = TRUE)
  expect_true(is.finite(fit$fit$objective))
})
