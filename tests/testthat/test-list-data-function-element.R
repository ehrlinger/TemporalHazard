# On the vector interface `data` masks the calling frame, so a FUNCTION-valued
# element of a list `data` was called in place of the function the expression
# names -- silently, changing the fit (#420, since 1.2.2). R's function lookup
# skips bindings that are not functions, so a numeric element of the same name
# was always harmless, and a data frame cannot hold a function column at all.
# hazard() now refuses such an element, as stats::lm() does.

ldf_setup <- function() {
  set.seed(1)
  n <- 40
  list(t = rexp(n), s = rbinom(n, 1, 0.7))
}

# Records whether it ran, so a test asserts WHICH function was called rather
# than trusting a message.
ldf_flagged_rep <- function(env) {
  function(...) {
    env$called <- TRUE
    runif(40, 0.2, 3)
  }
}

test_that("a function-valued element of list data is refused, naming it", {
  d <- ldf_setup()
  env <- new.env()
  env$called <- FALSE
  data_list <- list(t = d$t, s = d$s, rep = ldf_flagged_rep(env))
  expect_error(
    hazard(time = t, status = s, data = data_list, weights = rep(1, 40),
           dist = "weibull", theta = c(1, 1), fit = TRUE),
    "rep"
  )
  # The refusal happens before any expression is evaluated, so the masking
  # function never runs. On 1.2.11 it ran and the fit silently changed.
  expect_false(env$called)
})

test_that("the refusal names every function element, and says what to do", {
  d <- ldf_setup()
  msg <- tryCatch(
    hazard(time = t, status = s,
           data = list(t = d$t, s = d$s, rep = base::rep, seq = base::seq),
           dist = "weibull", theta = c(1, 1), fit = TRUE),
    error = conditionMessage
  )
  expect_match(msg, "'rep'", fixed = TRUE)
  expect_match(msg, "'seq'", fixed = TRUE)
  expect_match(msg, "vector interface", fixed = TRUE)
})

test_that("a numeric element of the same name still fits, unchanged", {
  # The negative control the fix must not disturb: a non-function binding was
  # never consulted for a call, so these three and the no-`data` reference
  # must all give the same objective.
  d <- ldf_setup()
  env <- new.env()
  env$called <- FALSE
  ref <- hazard(time = d$t, status = d$s, weights = rep(1, 40),
                dist = "weibull", theta = c(1, 1), fit = TRUE)
  num_list <- hazard(time = t, status = s,
                     data = list(t = d$t, s = d$s, rep = seq_len(40)),
                     weights = rep(1, 40), dist = "weibull",
                     theta = c(1, 1), fit = TRUE)
  num_col <- hazard(time = t, status = s,
                    data = data.frame(t = d$t, s = d$s, rep = seq_len(40)),
                    weights = rep(1, 40), dist = "weibull",
                    theta = c(1, 1), fit = TRUE)
  expect_equal(num_list$fit$objective, ref$fit$objective)
  expect_equal(num_col$fit$objective, ref$fit$objective)
  expect_false(env$called)
})

test_that("a data frame with a list-column of functions still fits", {
  # A list-column is a LIST, not a function, so R's function lookup skips it
  # and the formula-free path was never exposed through a data frame. The
  # refusal must not start rejecting one.
  d <- ldf_setup()
  df <- data.frame(t = d$t, s = d$s)
  df$rep <- replicate(40, function(...) stop("must not be called"),
                      simplify = FALSE)
  fit <- hazard(time = t, status = s, data = df, weights = rep(1, 40),
                dist = "weibull", theta = c(1, 1), fit = TRUE)
  expect_true(is.finite(fit$fit$objective))
})

test_that("the deliberate helper in list data now errors, with a workaround", {
  # The cost of the refusal, pinned rather than left to be discovered: using
  # the mask to reach a helper is no longer possible. The message has to
  # carry the way out.
  d <- ldf_setup()
  msg <- tryCatch(
    hazard(time = f(t), status = s,
           data = list(t = d$t, s = d$s, f = function(x) x * 2),
           dist = "weibull", theta = c(1, 1), fit = TRUE),
    error = conditionMessage
  )
  expect_match(msg, "'f'", fixed = TRUE)
  expect_match(msg, "calling environment", fixed = TRUE)
  # The workaround the message prescribes must actually work.
  helper <- function(x) x * 2
  fit <- hazard(time = helper(t), status = s,
                data = list(t = d$t, s = d$s),
                dist = "weibull", theta = c(1, 1), fit = TRUE)
  expect_true(is.finite(fit$fit$objective))
})

test_that("a hand-built data frame carrying a function column does not fit", {
  # data.frame(), `$<-` and `[[<-` all refuse a function column, so a normally
  # built frame cannot hold one. structure(list(...), class = "data.frame")
  # can, and is.data.frame() is TRUE for it. That is why the check below is
  # not exempted for data frames -- an exemption would be written as though
  # frames were immune, which this object disproves.
  #
  # It does NOT reach the check: .hzr_numeric_frame_values() reads every
  # column first and dies on the function with "attempt to replicate an
  # object of type 'special'". So the shape is refused, opaquely, by an
  # earlier guard. What this test pins is the part that matters -- it does
  # not silently fit -- not the wording, which comes from elsewhere and is
  # outside this change's scope.
  d <- ldf_setup()
  df <- structure(list(t = d$t, s = d$s, rep = base::rep),
                  class = "data.frame", row.names = seq_len(40))
  expect_true(is.data.frame(df))
  expect_true(any(vapply(df, is.function, logical(1))))
  expect_error(
    hazard(time = t, status = s, data = df, weights = rep(1, 40),
           dist = "weibull", theta = c(1, 1), fit = TRUE)
  )
})
