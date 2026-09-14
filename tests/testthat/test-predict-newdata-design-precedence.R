# Which route builds the covariate design at newdata (#272).
#
# newdata may carry the formula's variables (grp = "old") or the fit's
# design columns (grpyoung = 1).  When it carried both and they disagreed,
# the design columns silently won: grp = "old" plus grpyoung = 1 gave the
# "young" value.  User newdata is now rebuilt from the formula's variables
# whenever they are all present; a design-named column is then an unused
# extra.  Design columns are used only when the variables are absent (a fit
# saved before the design was stored, or a caller passing design columns),
# and hzr_deciles() / hzr_gof(), which evaluate at the fitted design, say so.

.dp_avc <- local({
  data(avc, package = "TemporalHazard")
  d <- na.omit(avc)
  d$grp <- factor(ifelse(d$age > 100, "old", "young"))
  d
})

.dp_th <- c(mu = 0.01, nu = 0.5, b_age = 0.004, b_young = 0.7)

.dp_weibull <- function(eta, time) (0.01 * time)^0.5 * exp(eta)

test_that("a contradicting design column does not override the variable", {
  w <- hazard(survival::Surv(int_dead, dead) ~ age + grp, data = .dp_avc,
              dist = "weibull", theta = .dp_th)
  old <- .dp_weibull(0.004 * 60, 2)
  young <- .dp_weibull(0.004 * 60 + 0.7, 2)
  # The two answers differ, so the check below can fail.
  expect_gt(young / old, 1.5)

  for (type in c("cumulative_hazard", "survival")) {
    want <- if (type == "survival") exp(-old) else old
    got <- predict(w, type = type,
                   newdata = data.frame(time = 2, age = 60, grp = "old",
                                        grpyoung = 1))
    expect_equal(unname(got), want, tolerance = 1e-12, label = type)
  }
  expect_equal(
    predict(w, type = "linear_predictor",
            newdata = data.frame(age = 60, grp = "old", grpyoung = 1)),
    0.004 * 60, tolerance = 1e-12)
})

test_that("some variables beside all design columns is refused, not guessed", {
  # With `sex` missing the variables cannot be rebuilt, and taking the
  # design columns would silently ignore grp = "old" (it gave 0.362, the
  # "young" value).  Neither answer is safe, so it stops.
  d <- .dp_avc
  d$sex <- factor(ifelse(d$mal == 1, "M", "F"))
  w <- hazard(survival::Surv(int_dead, dead) ~ age + grp + sex, data = d,
              dist = "weibull", theta = c(mu = 0.01, nu = 0.5, 0.004, 0.7, 0.2))
  expect_error(
    predict(w, type = "cumulative_hazard",
            newdata = data.frame(time = 2, age = 60, grp = "old",
                                 grpyoung = 1, sexM = 0)),
    "gives the formula variable\\(s\\) 'grp'.*lacks 'sex'"
  )
  # A transform: log(age) is a design column, age the missing variable.
  w2 <- hazard(survival::Surv(int_dead, dead) ~ log(age) + grp, data = d,
               dist = "weibull", theta = c(mu = 0.01, nu = 0.5, 0.2, 0.7))
  nd2 <- data.frame(time = 2, check.names = FALSE, `log(age)` = log(60),
                    grp = "old", grpyoung = 1)
  expect_error(predict(w2, type = "cumulative_hazard", newdata = nd2),
               "gives the formula variable\\(s\\) 'grp'.*lacks 'age'")
})

test_that("a variable that other columns are built from is not left stale", {
  # Copy the design, change `age`: I(age^2) and age:grpyoung still hold the
  # old values.  Taking them as given gave 0.362 where age * grp with
  # grp "young" gives 0.660.  `age` is a formula variable feeding other
  # columns, so with `grp` missing this is a mix, and it stops.
  d <- .dp_avc
  w <- hazard(survival::Surv(int_dead, dead) ~ age * grp, data = d,
              dist = "weibull", theta = c(mu = 0.01, nu = 0.5, 0.004, 0.7, 0.01))
  nd <- data.frame(time = 2, age = 60, grpyoung = 1, `age:grpyoung` = 0,
                   check.names = FALSE)
  expect_error(predict(w, type = "cumulative_hazard", newdata = nd),
               "gives the formula variable\\(s\\) 'age'.*lacks 'grp'")
  w2 <- hazard(survival::Surv(int_dead, dead) ~ age + I(age^2) + grp, data = d,
               dist = "weibull", theta = c(mu = 0.01, nu = 0.5, 0.004, 0, 0.7))
  nd2 <- data.frame(time = 2, age = 60, `I(age^2)` = 0, grpyoung = 1,
                    check.names = FALSE)
  expect_error(predict(w2, type = "cumulative_hazard", newdata = nd2),
               "gives the formula variable\\(s\\) 'age'.*lacks 'grp'")
  # Given grp as well, the design is rebuilt and nothing is stale.
  nd3 <- data.frame(time = 2, age = 60, grp = "young")
  expect_equal(unname(predict(w, type = "cumulative_hazard", newdata = nd3)),
               .dp_weibull(0.004 * 60 + 0.7 + 0.01 * 60, 2), tolerance = 1e-12)
})

test_that("a fit saved before the design was stored rebuilds it (#301)", {
  # A formula fit from 1.2.10 or earlier has no x_design.  Taking its design
  # columns beside a contradicting grp swapped the answer, so #292 refused
  # any column but the design columns (main at 9ec83e7 refused both newdata
  # below).  The design is now rebuilt from the stored formula and data
  # frame, and used once it reproduces the fitted design, so the legacy fit
  # answers as a new one does.  Truth: (0.01 * 2)^0.5 for "old", times
  # exp(0.7) for "young"; the swap gave them the other way round.
  w <- hazard(survival::Surv(int_dead, dead) ~ grp, data = .dp_avc,
              dist = "weibull", theta = c(mu = 0.01, nu = 0.5, b_young = 0.7))
  leg <- w
  leg$data$x_design <- NULL
  truth <- c(0.141421356237, 0.284787639017)
  nd_var <- data.frame(time = 2, grp = c("old", "young"))
  nd_mix <- data.frame(time = 2, grp = c("old", "young"), grpyoung = c(1, 0))
  for (nd in list(nd_var, nd_mix)) {
    got <- unname(predict(leg, type = "cumulative_hazard", newdata = nd))
    expect_equal(got, truth, tolerance = 1e-10)
    expect_identical(
      got, unname(predict(w, type = "cumulative_hazard", newdata = nd))
    )
  }
  # An unused column beside the design columns is ignored, as for a new fit.
  got <- predict(leg, type = "cumulative_hazard",
                 newdata = data.frame(grpyoung = c(0, 1), time = 2, junk = 9))
  expect_equal(unname(got), truth, tolerance = 1e-10)
})

test_that("a legacy formula with a constant is not rebuilt, even one kept", {
  # `k` is copied with the call and has not moved here, but nothing in the
  # fitted design records it, so a constant is never trusted (#301
  # reviews).  Variables alone are refused, as on main at 9ec83e7; the
  # design columns are answered.
  k <- 50
  w <- hazard(survival::Surv(int_dead, dead) ~ age + I(age > k),
              data = .dp_avc, dist = "weibull",
              theta = c(mu = 0.01, nu = 0.5, 0.004, 0.3))
  leg <- w
  leg$data$x_design <- NULL
  expect_true(exists("k", envir = leg$call_env, inherits = FALSE))
  expect_error(
    predict(leg, type = "cumulative_hazard",
            newdata = data.frame(time = 2, age = c(40, 60), junk = 1)),
    "lacks the covariate column\\(s\\) 'I\\(age > k\\)TRUE'"
  )
  nd <- data.frame(time = 2, age = c(40, 60), `I(age > k)TRUE` = c(0, 1),
                   check.names = FALSE)
  expect_equal(unname(predict(leg, type = "cumulative_hazard", newdata = nd)),
               .dp_weibull(0.004 * c(40, 60) + 0.3 * c(0, 1), 2),
               tolerance = 1e-12)
})

test_that("a legacy formula passed by name is rebuilt when it is closed", {
  # hazard(f, data = ) stores the symbol `f`, looked up wherever `f` is
  # bound when predict() runs.  A later `f` with another cutoff `k` once
  # reproduced the fitted design and was used silently (#301 review); a
  # formula with a constant is no longer rebuilt at all.  A closed one
  # carries its constants in its column names, so an `f` that reproduces
  # the fitted design is that design, wherever it is bound.
  f <- survival::Surv(int_dead, dead) ~ grp
  w <- hazard(f, data = .dp_avc, dist = "weibull",
              theta = c(mu = 0.01, nu = 0.5, b_young = 0.7))
  leg <- w
  leg$data$x_design <- NULL
  nd <- data.frame(time = 2, grp = c("old", "young"))
  expect_true(exists("f", envir = leg$call_env, inherits = FALSE))
  expect_equal(unname(predict(leg, type = "cumulative_hazard", newdata = nd)),
               c(0.141421356237, 0.284787639017), tolerance = 1e-10)
  # `f` bound only beyond the copied bindings, as a top-level `f` is.
  leg$call_env <- new.env(parent = list2env(list(f = f)))
  expect_equal(unname(predict(leg, type = "cumulative_hazard", newdata = nd)),
               c(0.141421356237, 0.284787639017), tolerance = 1e-10)
  # Bound to anything but a formula, it is not used.
  leg$call_env <- new.env(parent = list2env(list(f = "grp")))
  expect_error(predict(leg, type = "cumulative_hazard", newdata = nd),
               "lacks the covariate column\\(s\\) 'grpyoung'")
})

test_that("a legacy formula's constant in a function frame is not trusted", {
  # A formula made in a function keeps that frame as its environment, and
  # the frame is saved as it stands when the object is saved, not when it
  # was fitted: here the fit used k = 50 and the frame then moved on to k2.
  # No fitted age separates them, so the rebuild reproduced the fitted
  # design and predicted with k2, silently (#301 third review).
  # No fitted age lies in (50, k2], so the fitted design cannot tell them apart.
  k2 <- (50 + min(.dp_avc$age[.dp_avc$age > 50])) / 2
  g <- function() {
    k <- 50
    f <- survival::Surv(int_dead, dead) ~ age + I(age > k)
    fit <- hazard(f, data = .dp_avc, dist = "weibull",
                  theta = c(mu = 0.01, nu = 0.5, 0.004, 0.3))
    k <- k2
    fit
  }
  leg <- g()
  leg$data$x_design <- NULL
  a <- c(40, (50 + k2) / 2)   # the second age lies between the cutoffs
  nd <- data.frame(time = 2, age = a, `I(age > k)TRUE` = c(0, 1),
                   check.names = FALSE)
  expect_equal(unname(predict(leg, type = "cumulative_hazard", newdata = nd)),
               .dp_weibull(0.004 * a + 0.3 * c(0, 1), 2), tolerance = 1e-12)
})

test_that("a closure copied with the call is not trusted as a value", {
  # thr, given by value to vapply(), is copied into call_env with the call,
  # but its body reads its own frame, which moved on after the fit (k = 50,
  # then k2).  Trusted, it rebuilt the design with k2 (#301 third review).
  k2 <- (50 + min(.dp_avc$age[.dp_avc$age > 50])) / 2
  g <- function() {
    k <- 50
    thr <- function(x) x > k
    fit <- hazard(survival::Surv(int_dead, dead) ~ age +
                    I(vapply(age, thr, logical(1))),
                  data = .dp_avc, dist = "weibull",
                  theta = c(mu = 0.01, nu = 0.5, 0.004, 0.3))
    k <- k2
    fit
  }
  leg <- g()
  leg$data$x_design <- NULL
  expect_true(is.function(get0("thr", envir = leg$call_env, inherits = FALSE)))
  col <- "I(vapply(age, thr, logical(1)))TRUE"
  expect_identical(colnames(leg$data$x), c("age", col))
  a <- c(40, (50 + k2) / 2)
  nd <- data.frame(time = 2, age = a, x = c(0, 1))
  names(nd)[3] <- col
  expect_equal(unname(predict(leg, type = "cumulative_hazard", newdata = nd)),
               .dp_weibull(0.004 * a + 0.3 * c(0, 1), 2), tolerance = 1e-12)
})

test_that("a design function shadowed by the user's own is not trusted", {
  # sqrt() is trusted only as base R's own.  Here the user's sqrt() reads a
  # `k` that moved after the fit, and no fitted age separates the two
  # cutoffs, so a rebuild matched the fitted design and predicted with the
  # new one (#301 fourth review).
  k2 <- (50 + min(.dp_avc$age[.dp_avc$age > 50])) / 2
  k <- 50
  sqrt <- function(x) x > k
  w <- hazard(survival::Surv(int_dead, dead) ~ age + sqrt(age),
              data = .dp_avc, dist = "weibull",
              theta = c(mu = 0.01, nu = 0.5, 0.004, 0.3))
  leg <- w
  leg$data$x_design <- NULL
  k <- k2
  expect_identical(colnames(leg$data$x), c("age", "sqrt(age)TRUE"))
  a <- c(40, (50 + k2) / 2)
  nd <- data.frame(time = 2, age = a, `sqrt(age)TRUE` = c(0, 1),
                   check.names = FALSE)
  expect_equal(unname(predict(leg, type = "cumulative_hazard", newdata = nd)),
               .dp_weibull(0.004 * a + 0.3 * c(0, 1), 2), tolerance = 1e-12)
})

test_that("a legacy formula calling the user's own function is not rebuilt", {
  # thr() in the workspace may have been redefined since the fit; moved to
  # a cutoff no fitted age separates from the old one, the rebuild matched
  # the fitted design and predicted with the new cutoff (#301 third review).
  k2 <- (50 + min(.dp_avc$age[.dp_avc$age > 50])) / 2
  had <- exists("thr", envir = globalenv(), inherits = FALSE)
  old <- if (had) get("thr", envir = globalenv())
  withr::defer(
    if (had) assign("thr", old, envir = globalenv())
    else rm("thr", envir = globalenv())
  )
  assign("thr", function(x) x > 50, envir = globalenv())
  w <- hazard(survival::Surv(int_dead, dead) ~ age + thr(age),
              data = .dp_avc, dist = "weibull",
              theta = c(mu = 0.01, nu = 0.5, 0.004, 0.3))
  leg <- w
  leg$data$x_design <- NULL
  assign("thr", function(x) x > k2, envir = globalenv())
  a <- c(40, (50 + k2) / 2)
  nd <- data.frame(time = 2, age = a, `thr(age)TRUE` = c(0, 1),
                   check.names = FALSE)
  expect_equal(unname(predict(leg, type = "cumulative_hazard", newdata = nd)),
               .dp_weibull(0.004 * a + 0.3 * c(0, 1), 2), tolerance = 1e-12)
})

test_that("a legacy formula of R's design functions and literals is rebuilt", {
  # splines::ns(), factor(), c(), I() and `>` are R's own, and the knots,
  # levels and cutoff are written into the formula, so the design is rebuilt
  # and the legacy fit answers as the same fit with its design (a pkg::fn
  # term was refused by review 2's rule, and c() by review 4's).
  w <- hazard(survival::Surv(int_dead, dead) ~
                splines::ns(age, knots = c(40, 60)) +
                factor(grp, levels = c("young", "old")) + I(age > 50),
              data = .dp_avc, dist = "weibull",
              theta = c(mu = 0.01, nu = 0.5, 0.1, 0.2, 0.3, 0.4, 0.5))
  leg <- w
  leg$data$x_design <- NULL
  nd <- data.frame(time = 2, age = c(40, 60), grp = c("young", "old"),
                   junk = 1)
  got <- unname(predict(leg, type = "cumulative_hazard", newdata = nd))
  expect_identical(
    got, unname(predict(w, type = "cumulative_hazard", newdata = nd))
  )
  expect_gt(abs(got[2] / got[1] - 1), 0.01)   # the rows differ
})

test_that("a legacy formula constant read from the workspace is not trusted", {
  # A top-level fit keeps no binding for `k` in I(age > k): at predict()
  # time it is whatever the workspace holds.  Moved to a value no fitted age
  # separates from the old one, it reproduced the fitted design and was used
  # silently, design-column newdata then being rebuilt from `age` with the
  # new k (#301 second review).  Such a fit is not rebuilt, so the columns
  # are taken as given, as on main at 9ec83e7.
  had <- exists("k", envir = globalenv(), inherits = FALSE)
  old <- if (had) get("k", envir = globalenv())
  withr::defer(
    if (had) assign("k", old, envir = globalenv())
    else rm("k", envir = globalenv())
  )
  assign("k", 50, envir = globalenv())
  w <- hazard(survival::Surv(int_dead, dead) ~ age + I(age > k),
              data = .dp_avc, dist = "weibull",
              theta = c(mu = 0.01, nu = 0.5, 0.004, 0.3))
  expect_false(exists("k", envir = w$call_env, inherits = FALSE))
  leg <- w
  leg$data$x_design <- NULL
  # No fitted age lies in (50, k2], so the fitted design cannot tell them apart.
  k2 <- (50 + min(.dp_avc$age[.dp_avc$age > 50])) / 2
  assign("k", k2, envir = globalenv())
  a <- c(40, (50 + k2) / 2)   # the second age lies between the cutoffs
  nd <- data.frame(time = 2, age = a, `I(age > k)TRUE` = c(0, 1),
                   check.names = FALSE)
  expect_equal(unname(predict(leg, type = "cumulative_hazard", newdata = nd)),
               .dp_weibull(0.004 * a + 0.3 * c(0, 1), 2), tolerance = 1e-12)
})

test_that("a function object pasted into a legacy formula is not trusted", {
  # bquote() can put a function itself, not its name, into a formula.  The
  # column name shows the function's code but not what its body reads, so
  # a rebuild matched the fitted design after `k` moved and predicted with
  # the new cutoff (#301 fifth review).
  k2 <- (50 + min(.dp_avc$age[.dp_avc$age > 50])) / 2
  k <- 50
  thr <- function(x) x > k
  f <- eval(bquote(survival::Surv(int_dead, dead) ~ age + I(.(thr)(age))))
  w <- hazard(f, data = .dp_avc, dist = "weibull",
              theta = c(mu = 0.01, nu = 0.5, 0.004, 0.3))
  leg <- w
  leg$data$x_design <- NULL
  k <- k2
  col <- colnames(leg$data$x)[2L]
  a <- c(40, (50 + k2) / 2)   # the second age lies between the cutoffs
  nd <- data.frame(time = 2, age = a, x = c(0, 1))
  names(nd)[3] <- col
  expect_equal(unname(predict(leg, type = "cumulative_hazard", newdata = nd)),
               .dp_weibull(0.004 * a + 0.3 * c(0, 1), 2), tolerance = 1e-12)
})

test_that("a legacy formula that looks a value up by string is not rebuilt", {
  # get("k") names no symbol a walker can see, and get() is not a design
  # function, so the formula is not closed.  Trusted as base R's own, it
  # rebuilt the design with the moved k (#301 fourth review).
  had <- exists("k", envir = globalenv(), inherits = FALSE)
  old <- if (had) get("k", envir = globalenv())
  withr::defer(
    if (had) assign("k", old, envir = globalenv())
    else rm("k", envir = globalenv())
  )
  assign("k", 50, envir = globalenv())
  w <- hazard(survival::Surv(int_dead, dead) ~ age + I(age > get("k")),
              data = .dp_avc, dist = "weibull",
              theta = c(mu = 0.01, nu = 0.5, 0.004, 0.3))
  leg <- w
  leg$data$x_design <- NULL
  col <- "I(age > get(\"k\"))TRUE"
  expect_identical(colnames(leg$data$x), c("age", col))
  # No fitted age lies in (50, k2], so the fitted design cannot tell them apart.
  k2 <- (50 + min(.dp_avc$age[.dp_avc$age > 50])) / 2
  assign("k", k2, envir = globalenv())
  a <- c(40, (50 + k2) / 2)
  nd <- data.frame(time = 2, age = a, x = c(0, 1))
  names(nd)[3] <- col
  expect_equal(unname(predict(leg, type = "cumulative_hazard", newdata = nd)),
               .dp_weibull(0.004 * a + 0.3 * c(0, 1), 2), tolerance = 1e-12)
})

test_that("a rebuilt legacy design carries the fit's predvars to unseen rows", {
  # poly(), scale(), ns() and bs() take coefficients, a centre and scale,
  # or knots from the fitting data, and predict() reuses them at new rows
  # through the recorded predvars.  Those live in no column name, so the
  # rebuild, from the same frame, must record the same ones: at ages not
  # in the fit, where a basis recomputed from newdata would differ, the
  # legacy fit answers as the same fit with its design.
  a <- c(33.3, 47.7, 61.1, 88.8)
  expect_false(any(a %in% .dp_avc$age))
  terms <- list(poly = quote(poly(age, 2)), scale = quote(scale(age)),
                ns = quote(splines::ns(age, df = 3)),
                bs = quote(splines::bs(age, df = 3)))
  p <- c(poly = 2, scale = 1, ns = 3, bs = 3)
  nd <- data.frame(time = 2, age = a)
  for (nm in names(terms)) {
    f <- eval(bquote(survival::Surv(int_dead, dead) ~ .(terms[[nm]])))
    w <- hazard(f, data = .dp_avc, dist = "weibull",
                theta = c(mu = 0.01, nu = 0.5,
                          seq(0.1, 0.3, length.out = p[[nm]])))
    leg <- w
    leg$data$x_design <- NULL
    got <- unname(predict(leg, type = "cumulative_hazard", newdata = nd))
    expect_identical(
      got, unname(predict(w, type = "cumulative_hazard", newdata = nd)),
      label = nm
    )
    expect_gt(max(got) / min(got), 1.01, label = nm)   # the rows differ
  }
})

test_that("a legacy formula computed in the call is not re-run", {
  # A computed `formula =` (as.formula(), reformulate(), maybe with
  # sample()) would run again, with its side effects and its current
  # inputs, on every predict().  It is not rebuilt; the refusal stands.
  w <- hazard(stats::as.formula("survival::Surv(int_dead, dead) ~ grp"),
              data = .dp_avc, dist = "weibull",
              theta = c(mu = 0.01, nu = 0.5, b_young = 0.7))
  w$data$x_design <- NULL
  expect_error(
    predict(w, type = "cumulative_hazard",
            newdata = data.frame(time = 2, grp = c("old", "young"),
                                 grpyoung = c(1, 0))),
    "earlier version.*also has 'grp'.*pass only the design columns"
  )
})

test_that("a legacy fit with a row-level vector outside `data` keeps its route", {
  # ~ age + wv, `wv` a vector beside the data frame, not a column: the
  # design is not rebuilt, and design-column newdata is answered as on main
  # at 9ec83e7 (a rebuild gave a raw "variable lengths differ" error; #301
  # review).
  wv <- seq_len(nrow(.dp_avc)) / nrow(.dp_avc)
  w <- hazard(survival::Surv(int_dead, dead) ~ age + wv, data = .dp_avc,
              dist = "weibull", theta = c(mu = 0.01, nu = 0.5, 0.004, 0.3))
  w$data$x_design <- NULL
  got <- predict(w, type = "cumulative_hazard",
                 newdata = data.frame(time = 2, age = 60, wv = 1))
  expect_equal(unname(got), .dp_weibull(0.004 * 60 + 0.3, 2), tolerance = 1e-12)
})

test_that("legacy design-column means follow the #272 rule, as a new fit's do", {
  # age and I(age^2) at the design-column means: `age` is a formula
  # variable, so the design is rebuilt from it and I(age^2) becomes
  # mean(age)^2, as for the same fit with its design.  Main at 9ec83e7 took
  # the columns as given (mean(age^2)); NEWS says so.
  w <- hazard(survival::Surv(int_dead, dead) ~ age + I(age^2), data = .dp_avc,
              dist = "weibull", theta = c(mu = 0.01, nu = 0.5, 0.004, 1e-4))
  leg <- w
  leg$data$x_design <- NULL
  nd <- as.data.frame(t(colMeans(w$data$x)), check.names = FALSE)
  nd$time <- 2
  got <- unname(predict(leg, type = "cumulative_hazard", newdata = nd))
  expect_identical(got,
                   unname(predict(w, type = "cumulative_hazard", newdata = nd)))
  m <- mean(.dp_avc$age)
  expect_equal(got, .dp_weibull(0.004 * m + 1e-4 * m^2, 2), tolerance = 1e-12)
})

test_that("a legacy design that is not rebuilt exactly keeps the refusal", {
  # The rebuilt design is used only if it reproduces the fitted one, same
  # columns and same values.  Anything missing or different leaves #292's
  # refusal in place rather than guessing (#301).
  w <- hazard(survival::Surv(int_dead, dead) ~ grp, data = .dp_avc,
              dist = "weibull", theta = c(mu = 0.01, nu = 0.5, b_young = 0.7))
  w$data$x_design <- NULL
  nd_mix <- data.frame(time = 2, grp = c("old", "young"), grpyoung = c(1, 0))
  broken <- list(
    # Same columns, one value differs: the frame is not the fitted data.
    value = function(f) {
      f$data$x[1L, "grpyoung"] <- 1 - f$data$x[1L, "grpyoung"]
      f
    },
    # Same values under another column name.
    names = function(f) {
      f$call$formula <- quote(survival::Surv(int_dead, dead) ~
                                I(grp == "young"))
      f
    },
    # The formula now builds other columns.
    columns = function(f) {
      f$call$formula <- quote(survival::Surv(int_dead, dead) ~ age + grp)
      f
    },
    # The formula no longer evaluates against the frame.
    error = function(f) {
      f$call$formula <- quote(survival::Surv(int_dead, dead) ~ log(grp))
      f
    },
    # The frame has lost a row since the fit.
    rows = function(f) {
      f$data$frame <- f$data$frame[-1L, ]
      f
    },
    # A 1.0.3-era fit kept no data frame.
    frame = function(f) {
      f$data$frame <- NULL
      f
    },
    # A fit saved before call_env kept no bindings for its call.
    call_env = function(f) {
      f$call_env <- NULL
      f
    }
  )
  for (nm in names(broken)) {
    b <- broken[[nm]](w)
    expect_error(
      predict(b, type = "cumulative_hazard", newdata = nd_mix),
      "earlier version.*also has 'grp'.*pass only the design columns",
      label = nm
    )
    # Design columns alone are still answered, by position-free name match.
    got <- predict(b, type = "cumulative_hazard",
                   newdata = data.frame(grpyoung = c(0, 1), time = 2))
    expect_equal(unname(got), c(0.141421356237, 0.284787639017),
                 tolerance = 1e-10, label = nm)
  }
})

test_that("a vector-interface fit has no formula to rebuild (#301)", {
  # Given data = (the masking path) it stores the frame, but no formula, so
  # nothing is rebuilt and an extra column stays unused, as before.
  d <- .dp_avc
  v <- hazard(time = int_dead, status = dead,
              x = cbind(age = d$age, mal = d$mal), data = d,
              dist = "weibull", theta = c(mu = 0.01, nu = 0.5, 0.004, 0.3))
  expect_true(is.data.frame(v$data$frame))
  expect_null(v$call$formula)
  expect_null(v$data$x_design)
  got <- predict(v, type = "cumulative_hazard",
                 newdata = data.frame(time = 2, age = 60, mal = 1,
                                      grp = "old", junk = 9))
  expect_equal(unname(got), .dp_weibull(0.004 * 60 + 0.3, 2), tolerance = 1e-12)
})

test_that("hzr_gof() and hzr_deciles() still run on a fit saved before x_design", {
  # Their newdata is design-level and marked, so the legacy refusal must
  # not touch it.  Values pinned from b0aa955 (before the refusal) on the
  # same legacy object; a fitted model, so 1e-4 (relative: all >= 0.02).
  w <- hazard(survival::Surv(int_dead, dead) ~ age + I(age^2) + grp,
              data = .dp_avc, dist = "weibull",
              theta = c(mu = 0.01, nu = 0.5, 0, 0, 0), fit = TRUE)
  w$data$x_design <- NULL
  g <- hzr_gof(w)
  expect_length(g$par_cumhaz, 270L)
  expect_equal(unname(g$par_cumhaz[c(1, 135, 270)]),
               c(0.02289393213, 0.2218081341, 0.3079777839), tolerance = 1e-4)
  expect_equal(sum(hzr_deciles(w, time = 60)$expected), 67.99999988,
               tolerance = 1e-4)
})

test_that("a vector-interface fit still ignores a column it does not use", {
  # No formula, so no variable that could contradict: an extra column is
  # just unused, as for a fit with a stored design.
  d <- .dp_avc
  v <- hazard(time = d$int_dead, status = d$dead,
              x = cbind(age = d$age, mal = d$mal),
              dist = "weibull", theta = c(mu = 0.01, nu = 0.5, 0.004, 0.3))
  got <- predict(v, type = "cumulative_hazard",
                 newdata = data.frame(time = 2, age = 60, mal = 1, junk = 9))
  expect_equal(unname(got), .dp_weibull(0.004 * 60 + 0.3, 2), tolerance = 1e-12)
})

test_that("design columns alone are still taken by name", {
  w <- hazard(survival::Surv(int_dead, dead) ~ age + grp, data = .dp_avc,
              dist = "weibull", theta = .dp_th)
  got <- predict(w, type = "cumulative_hazard",
                 newdata = data.frame(grpyoung = c(0, 1), time = 2, age = 60))
  expect_equal(unname(got),
               .dp_weibull(0.004 * 60 + 0.7 * c(0, 1), 2), tolerance = 1e-12)
  # Pinned from c4678f1, before #272: this route is unchanged.
  expect_equal(unname(got), c(0.179781779, 0.3620360441), tolerance = 1e-9)
})

test_that("hzr_gof() still evaluates at the design-column means", {
  # For I(age^2) the mean of the design column is mean(age^2), not
  # mean(age)^2: the rebuild route would give the latter.
  w <- hazard(survival::Surv(int_dead, dead) ~ age + I(age^2), data = .dp_avc,
              dist = "weibull", theta = c(mu = 0.01, nu = 0.5, 0, 0),
              fit = TRUE)
  th <- w$fit$theta
  x_bar <- colMeans(w$data$x)
  expect_gt(abs(x_bar[[2]] / x_bar[[1]]^2 - 1), 0.1)   # the routes differ
  gof <- hzr_gof(w)
  want <- (th[[1]] * gof$time)^th[[2]] * exp(sum(x_bar * th[3:4]))
  expect_equal(unname(gof$par_cumhaz), want, tolerance = 1e-10)
  # Pinned from c4678f1, before #272. A fitted, poorly scaled model, so
  # 1e-4 to absorb cross-platform optimizer noise; the rebuild route would
  # move it by about 31%.
  expect_length(gof$par_cumhaz, 270L)
  expect_equal(unname(gof$par_cumhaz[c(1, 135, 270)]),
               c(0.02303918854, 0.2220404487, 0.3080652665),
               tolerance = 1e-4)
})

test_that("hzr_deciles() still evaluates each subject at its own design row", {
  w <- hazard(survival::Surv(int_dead, dead) ~ age + I(age^2), data = .dp_avc,
              dist = "weibull", theta = c(mu = 0.01, nu = 0.5, 0, 0),
              fit = TRUE)
  th <- w$fit$theta
  dec <- hzr_deciles(w, time = 60)
  # Expected events are the sum of each subject's H at their own follow-up.
  h <- (th[[1]] * w$data$time)^th[[2]] *
    exp(as.numeric(w$data$x %*% th[3:4]))
  expect_equal(sum(dec$expected), sum(h), tolerance = 1e-8)
  # Pinned from c4678f1, before #272. 1e-4 (relative: the value is ~68)
  # absorbs cross-platform optimizer noise on a poorly scaled fit. This
  # pins the value; the check against the design rows above pins the route,
  # since here the two routes happen to agree.
  expect_equal(sum(dec$expected), 67.99859561, tolerance = 1e-4)
})

test_that("hzr_deciles() uses the fitted rows even if a formula constant changed", {
  # ~ age + I(age > k): its design frame carries `age`, a formula variable,
  # so without the design-level marker it would be rebuilt, and I(age > k)
  # re-evaluated with the current k.  The fitted rows must win.
  d <- .dp_avc
  k <- 50
  w <- hazard(survival::Surv(int_dead, dead) ~ age + I(age > k), data = d,
              dist = "weibull", theta = c(mu = 0.01, nu = 0.5, 0, 0),
              fit = TRUE)
  th <- w$fit$theta
  want <- sum((th[[1]] * w$data$time)^th[[2]] *
                exp(as.numeric(w$data$x %*% th[3:4])))
  k <- 500   # now every I(age > k) would be FALSE on a rebuild
  expect_gt(sum(w$data$x[, 2]), 0)            # the fitted column is not all 0
  expect_equal(sum(hzr_deciles(w, time = 60)$expected), want,
               tolerance = 1e-8)
})

test_that("a multiphase fit saved before the design was stored rebuilds it", {
  skip_on_cran()  # a multiphase fit
  set.seed(1)
  m <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ age + grp, data = .dp_avc,
    dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes"),
      constant = hzr_phase("constant")),
    fit = TRUE))
  leg <- m
  leg$data$x_design <- NULL
  nd_mix <- data.frame(time = 2, age = 60, grp = c("old", "young"),
                       grpyoung = c(1, 0))
  # The global design is rebuilt (#301), so grp wins as it does for the
  # same fit with its design; main at 9ec83e7 refused this newdata.
  got <- unname(predict(leg, type = "cumulative_hazard", newdata = nd_mix))
  expect_identical(
    got, unname(predict(m, type = "cumulative_hazard", newdata = nd_mix))
  )
  # Pinned from main at 4b68020 on the same object (a fitted model: 1e-4),
  # as design columns (grpyoung 0 then 1: "old" then "young").
  expect_equal(got, c(0.0140056, 0.320817), tolerance = 1e-4)
  # Without its data frame (a 1.0.3-era fit) nothing can be rebuilt, and
  # the refusal stands.
  leg$data$frame <- NULL
  expect_error(
    predict(leg, type = "cumulative_hazard", newdata = nd_mix),
    "earlier version.*also has 'grp'.*pass only the design columns"
  )
})

test_that("multiphase time_windows: the rebuilt global design is expanded", {
  skip_on_cran()  # a multiphase fit
  # The global design a formula-less phase inherits is window-expanded at
  # fit time (age_w1, age_w2).  Rebuilt unexpanded at newdata, two rows met
  # two per-window coefficients and %*% returned four values, silently.
  # Reference: predict() at the fitted data, which uses the stored design.
  d <- stats::na.omit(get("avc", envir = asNamespace("TemporalHazard")))
  set.seed(1)
  m <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ age, data = d, dist = "multiphase",
    time_windows = 12,
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes"),
      constant = hzr_phase("constant")),
    fit = TRUE))
  # One subject on each side of the cut, so both windows are exercised.
  rows <- c(which(d$int_dead <= 12)[1], which(d$int_dead > 12)[1])
  want <- predict(m, type = "cumulative_hazard")[rows]
  nd <- data.frame(time = d$int_dead[rows], age = d$age[rows])
  got <- predict(m, newdata = nd, type = "cumulative_hazard")
  expect_length(got, 2L)
  expect_equal(unname(got), unname(want), tolerance = 1e-10)
})

test_that("a phase-formula factor with reordered levels codes as the fit did", {
  skip_on_cran()  # a multiphase fit
  # Global ~ age plus a phase formula ~ grp.  #292 makes this phase path
  # reachable at newdata; before #290 kept the fit's levels, reversed
  # levels came back swapped (0.213738/0.21485) with no error.  The answer
  # must be the structural one, whatever order or type grp arrives in.
  set.seed(1)
  pf <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ age, data = .dp_avc, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes",
                        formula = ~ grp),
      constant = hzr_phase("constant")),
    fit = TRUE))
  th <- pf$fit$theta
  tt <- c(2, 20)
  base <- predict(pf, newdata = data.frame(time = tt), type = "cumulative_hazard",
                  decompose = TRUE)
  # Row 1 is "old" (early: no grp effect), row 2 "young"; the constant
  # phase inherits the global age effect at age 60.
  ref <- base$early * exp(th[["early.grpyoung"]] * c(0, 1)) +
    base$constant * exp(sum(th[grep("^constant\\.age", names(th))]) * 60)
  expect_gt(abs(th[["early.grpyoung"]]), 0.5)          # levels matter here
  for (g in list(factor(c("old", "young"), levels = c("young", "old")),
                 c("old", "young"),
                 factor(c("old", "young"), levels = c("old", "young")))) {
    got <- predict(pf, newdata = data.frame(time = tt, age = 60, grp = g),
                   type = "cumulative_hazard")
    expect_equal(unname(got), ref, tolerance = 1e-10)
  }
  # The structural values, pinned (a fitted model: 1e-4, relative).
  expect_equal(ref, c(0.05608905, 0.3599925), tolerance = 1e-4)
})

test_that("a phase formula the fit did not use is not rebuilt at newdata", {
  skip_on_cran()  # multiphase fits
  # The fit uses a phase's own formula only on the formula interface; a
  # vector-interface fit ignores hzr_phase(formula = ), the phase inherits
  # the global x, and attr(fit$x_list, "from_formula") records FALSE.
  # predict() routed on !is.null(ph$formula) alone and rebuilt the unused
  # formula: log(age) times a coefficient fitted on age, and under
  # time_windows eight values for four rows -- on main 7e50e2b as well.
  # Reference: predict() at the fitted rows (the stored design).
  d <- stats::na.omit(get("avc", envir = asNamespace("TemporalHazard")))
  rows <- c(which(d$int_dead <= 12)[1:2], which(d$int_dead > 12)[1:2])
  nd <- data.frame(time = d$int_dead[rows], age = d$age[rows])
  for (tw in list(NULL, 12)) {
    set.seed(1)
    f <- suppressWarnings(hazard(
      time = d$int_dead, status = d$dead, x = cbind(age = d$age),
      dist = "multiphase", time_windows = tw,
      phases = list(
        early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                          fixed = "shapes", formula = ~ log(age)),
        constant = hzr_phase("constant")),
      fit = TRUE))
    # The fit recorded that it ignored the formula (not merely no record).
    expect_identical(unname(attr(f$fit$x_list, "from_formula")),
                     c(FALSE, FALSE))
    want <- predict(f, type = "cumulative_hazard")[rows]
    got <- predict(f, newdata = nd, type = "cumulative_hazard")
    label <- if (is.null(tw)) "no windows" else "time_windows = 12"
    expect_length(got, length(rows))
    expect_equal(unname(got), unname(want), tolerance = 1e-10, label = label)
  }
})

test_that("a fit without the from_formula record routes by formula and columns", {
  skip_on_cran()  # multiphase fits
  # Every fit saved before the record existed (it is in no tag through
  # v1.2.9) is routed by the fallback hzr_gof() shares: a phase uses its
  # formula iff it has one AND its stored columns are not the inherited
  # ones.  The vector fit still takes the global route, a formula fit
  # still takes its phase formula.
  d <- stats::na.omit(get("avc", envir = asNamespace("TemporalHazard")))
  d$grp <- factor(ifelse(d$age > 100, "old", "young"))
  rows <- c(1, 50, 150)
  set.seed(1)
  vf <- suppressWarnings(hazard(
    time = d$int_dead, status = d$dead, x = cbind(age = d$age),
    dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                        fixed = "shapes", formula = ~ log(age)),
      constant = hzr_phase("constant")),
    fit = TRUE))
  attr(vf$fit$x_list, "from_formula") <- NULL
  expect_equal(
    unname(predict(vf, type = "cumulative_hazard",
                   newdata = data.frame(time = d$int_dead[rows],
                                        age = d$age[rows]))),
    unname(predict(vf, type = "cumulative_hazard")[rows]), tolerance = 1e-10)
  set.seed(1)
  ff <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ age, data = d, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                        fixed = "shapes", formula = ~ grp),
      constant = hzr_phase("constant")),
    fit = TRUE))
  attr(ff$fit$x_list, "from_formula") <- NULL
  expect_equal(
    unname(predict(ff, type = "cumulative_hazard",
                   newdata = data.frame(time = d$int_dead[rows],
                                        age = d$age[rows],
                                        grp = d$grp[rows]))),
    unname(predict(ff, type = "cumulative_hazard")[rows]), tolerance = 1e-10)
})

test_that("a 1.0.3-era fit (no record, no frame, no design) keeps its phase formula", {
  skip_on_cran()  # a multiphase fit
  # Fitted on its phase formula ~ log(age), but saved before x_design,
  # data$frame (1.1.0) and the from_formula record existed.  A frame-based
  # fallback sent it down the global route (age times a log(age)
  # coefficient); the column test keeps it on its phase formula, as main
  # 7e50e2b did.  Values pinned from main on the same stripped object.
  d <- stats::na.omit(get("avc", envir = asNamespace("TemporalHazard")))
  set.seed(1)
  f <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ age, data = d, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                        fixed = "shapes", formula = ~ log(age)),
      constant = hzr_phase("constant")),
    fit = TRUE))
  f$data$x_design <- NULL
  f$fit$x_design <- NULL
  f$data$frame <- NULL
  attr(f$fit$x_list, "from_formula") <- NULL
  rows <- c(which(d$int_dead <= 12)[1:2], which(d$int_dead > 12)[1:2])
  got <- predict(f, type = "cumulative_hazard",
                 newdata = data.frame(time = d$int_dead[rows], age = d$age[rows]))
  expect_equal(unname(got), unname(predict(f, type = "cumulative_hazard")[rows]),
               tolerance = 1e-10)
  expect_equal(unname(got),
               c(0.01628361325, 0.08219660741, 0.1292144135, 0.1674520008),
               tolerance = 1e-4)
})

test_that("time_windows: an unused phase formula without the record goes global", {
  skip_on_cran()  # a multiphase fit
  # The phase inherited the window-expanded global x (age_w1, age_w2).
  # The fallback compares with those window names, not with data$x's age,
  # so it keeps the phase on the global route: four values, matching.
  d <- stats::na.omit(get("avc", envir = asNamespace("TemporalHazard")))
  set.seed(1)
  f <- suppressWarnings(hazard(
    time = d$int_dead, status = d$dead, x = cbind(age = d$age),
    dist = "multiphase", time_windows = 12,
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                        fixed = "shapes", formula = ~ log(age)),
      constant = hzr_phase("constant")),
    fit = TRUE))
  attr(f$fit$x_list, "from_formula") <- NULL
  rows <- c(which(d$int_dead <= 12)[1:2], which(d$int_dead > 12)[1:2])
  got <- predict(f, type = "cumulative_hazard",
                 newdata = data.frame(time = d$int_dead[rows], age = d$age[rows]))
  expect_length(got, 4L)
  expect_equal(unname(got), unname(predict(f, type = "cumulative_hazard")[rows]),
               tolerance = 1e-10)
})

test_that("coinciding global and phase formulas give the same answer either way", {
  skip_on_cran()  # a multiphase fit
  # Global ~ age with phase ~ age: without the record the phase's columns
  # equal the inherited ones, so the fallback takes the global route; with
  # the record it takes the phase route.  Both must be the fitted answer.
  d <- stats::na.omit(get("avc", envir = asNamespace("TemporalHazard")))
  set.seed(1)
  f <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ age, data = d, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                        fixed = "shapes", formula = ~ age),
      constant = hzr_phase("constant")),
    fit = TRUE))
  rows <- c(1, 50, 150, 250)
  nd <- data.frame(time = d$int_dead[rows], age = d$age[rows])
  want <- unname(predict(f, type = "cumulative_hazard")[rows])
  via_phase <- unname(predict(f, newdata = nd, type = "cumulative_hazard"))
  g <- f
  attr(g$fit$x_list, "from_formula") <- NULL
  via_global <- unname(predict(g, newdata = nd, type = "cumulative_hazard"))
  expect_equal(via_phase, want, tolerance = 1e-10)
  expect_equal(via_global, want, tolerance = 1e-10)
})

test_that("time_windows with a factor and several covariates matches the fit", {
  skip_on_cran()  # a multiphase fit
  d <- stats::na.omit(get("avc", envir = asNamespace("TemporalHazard")))
  d$grp <- factor(ifelse(d$age > 100, "old", "young"))
  set.seed(1)
  m <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ age + grp, data = d,
    dist = "multiphase", time_windows = 12,
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes"),
      constant = hzr_phase("constant")),
    fit = TRUE))
  rows <- c(which(d$int_dead <= 12)[1:2], which(d$int_dead > 12)[1:2])
  # grp arrives with its levels reversed, so ignoring the fit's levels
  # would code it wrong; both levels are present among the rows.
  nd <- data.frame(time = d$int_dead[rows], age = d$age[rows],
                   grp = factor(as.character(d$grp[rows]),
                                levels = c("young", "old")))
  expect_setequal(as.character(nd$grp), c("old", "young"))
  got <- predict(m, newdata = nd, type = "cumulative_hazard")
  expect_length(got, length(rows))
  expect_equal(unname(got),
               unname(predict(m, type = "cumulative_hazard")[rows]),
               tolerance = 1e-10)
})

test_that("an extra newdata column cannot mask a formula constant", {
  # cutoff is a constant of the formula's environment, not a fitting-data
  # variable, so a same-named newdata column is an unused extra.  Passing
  # all of newdata to model.frame() let it replace the constant: an extra
  # cutoff = 0 turned I(30 > 50) into I(30 > 0), silently (0.425 for 0.191).
  cutoff <- 50
  w <- hazard(survival::Surv(int_dead, dead) ~ I(age > cutoff) + mal,
              data = .dp_avc, dist = "weibull",
              theta = c(mu = 0.01, nu = 0.5, b_gt = 0.8, b_mal = 0.3))
  want_eta <- 0.8 * (30 > 50) + 0.3
  nd <- data.frame(time = 2, age = 30, mal = 1, cutoff = 0)
  expect_equal(unname(predict(w, newdata = nd, type = "cumulative_hazard")),
               .dp_weibull(want_eta, 2), tolerance = 1e-12)
  expect_equal(predict(w, newdata = nd[, -1], type = "linear_predictor"),
               want_eta, tolerance = 1e-12)
  # The extra column is ignored, as documented: same as without it.
  expect_equal(predict(w, newdata = nd, type = "cumulative_hazard"),
               predict(w, newdata = nd[, c("time", "age", "mal")],
                       type = "cumulative_hazard"), tolerance = 1e-12)
})

# ---- formula constants, outside-data row covariates, and the row backstop ---
# One rule in both rebuild helpers.  A formula variable that is not a
# fitting-data column is classified by the length of its value in the
# formula's environment: one value per fitting row makes it a covariate taken
# from outside `data`, which must then come from newdata; any other length
# makes it a constant, which a same-named newdata column must never mask.
# A rebuilt design with a row count other than newdata's is refused.

.oc_avc <- function() {
  stats::na.omit(get("avc", envir = asNamespace("TemporalHazard")))
}

test_that("an extra newdata column cannot mask a phase-formula constant", {
  skip_on_cran()  # a multiphase fit
  d <- .oc_avc()
  cutoff <- 50
  set.seed(1)
  m <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes",
                        formula = ~ I(age > cutoff)),
      constant = hzr_phase("constant")),
    fit = TRUE))
  # I(30 > 50) is FALSE, so the answer is the covariate-free baseline.
  base <- predict(m, newdata = data.frame(time = 2), type = "cumulative_hazard")
  masked <- predict(m, type = "cumulative_hazard",
                    newdata = data.frame(time = 2, age = 30, cutoff = 0))
  expect_equal(unname(masked), unname(base), tolerance = 1e-10)
  # The indicator matters, so a masked constant would have shown.
  expect_gt(abs(m$fit$theta[[grep("^early\\.I", names(m$fit$theta))]]), 0.1)
})

test_that("a term with row-level values from outside data is refused", {
  # predict(newdata =) takes only the columns of the model's `data`. A term
  # built from an object kept outside `data` cannot be rebuilt for new rows,
  # even when newdata carries a column of the same name, so it is refused.
  d <- .oc_avc()
  n <- nrow(d)
  set.seed(3)
  zz <- rnorm(n)
  M <- cbind(zz, zz^2)
  Lz <- list(z = zz, k = 3)
  en <- new.env()
  en$z <- zz
  ext <- data.frame(z = zz)
  rv <- rev(seq_len(n))
  nd <- d[rv, ]
  nd$zz <- zz[rv]
  nd$M <- M[rv, ]
  msg <- "uses row-level values taken from outside `data`"
  fit <- function(rhs, beta) {
    hazard(stats::as.formula(paste("survival::Surv(int_dead, dead) ~", rhs)),
           data = d, dist = "weibull", theta = c(mu = 0.01, nu = 0.5, beta))
  }
  lp <- function(w) predict(w, newdata = nd, type = "linear_predictor")
  for (rhs in c("zz", "M[, 1]", "Lz$z", "en$z", "ext$z")) {
    expect_error(lp(fit(rhs, 0.3)), msg)
  }
  expect_error(lp(fit("M", c(0.8, 0.3))), msg)
  expect_error(lp(fit("poly(M[, 1], 2)", c(0.3, 0.1))), msg)
  expect_error(lp(fit("age + zz", c(0.01, 0.3))), "term 'zz' of the model")
})

test_that("a phase term with row-level values from outside data is refused", {
  skip_on_cran()  # multiphase fits
  d <- .oc_avc()
  set.seed(7)
  zz <- d$age[sample(nrow(d))] / 100   # a row covariate, not a column of d
  Lz <- list(z = zz, k = 3)
  rv <- rev(seq_len(nrow(d)))
  nd <- data.frame(time = d$int_dead[rv], zz = zz[rv])
  mp <- function(glob, form) {
    set.seed(1)
    suppressWarnings(hazard(
      glob, data = d, dist = "multiphase",
      phases = list(
        early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1,
                          fixed = "shapes", formula = form),
        constant = hzr_phase("constant")),
      fit = TRUE))
  }
  f_zz <- mp(survival::Surv(int_dead, dead) ~ 1, ~ zz)
  f_lz <- mp(survival::Surv(int_dead, dead) ~ 1, ~ Lz$z)
  g_zz <- mp(survival::Surv(int_dead, dead) ~ zz, NULL)
  ch <- function(f, x) predict(f, newdata = x, type = "cumulative_hazard")
  expect_error(ch(f_zz, nd), "term 'zz' of phase 'early' uses row-level")
  expect_error(ch(f_lz, nd), "term 'Lz\\$z' of phase 'early' uses row-level")
  expect_error(ch(g_zz, nd), "term 'zz' of the model uses row-level")

  # A fit saved without the phase design but with its data knows zz is not a
  # data column, and refuses the same way.
  f18 <- f_zz
  f18$fit$x_design <- NULL
  expect_error(ch(f18, nd), "term 'zz' of phase 'early' uses row-level")
  # A 1.0.3-era fit (no design, frame or record) cannot tell zz from a data
  # column, so it takes zz from newdata, as it did then.
  f19 <- f18
  f19$data$frame <- NULL
  attr(f19$fit$x_list, "from_formula") <- NULL
  rows <- c(1, 50, 150)
  expect_equal(unname(ch(f19, nd[match(rows, rv), ])),
               unname(predict(f_zz, type = "cumulative_hazard")[rows]),
               tolerance = 1e-10)
})

test_that("a rebuilt design never has more rows than newdata", {
  skip_on_cran()  # a multiphase fit
  # A 1.0.3-era fit (no stored design, frame or record) cannot classify zz,
  # so its rebuild would take the 305-row fitting vector for a one-row
  # newdata.  The backstop refuses instead of returning 305 values.
  d <- .oc_avc()
  set.seed(7)
  zz <- d$age[sample(nrow(d))] / 100
  set.seed(1)
  f <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ 1, data = d, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes",
                        formula = ~ zz),
      constant = hzr_phase("constant")),
    fit = TRUE))
  f$fit$x_design <- NULL
  f$data$frame <- NULL
  attr(f$fit$x_list, "from_formula") <- NULL
  expect_error(predict(f, newdata = data.frame(time = 2, age = 60),
                       type = "cumulative_hazard"),
               "rows for 1 row\\(s\\) of 'newdata': a term uses row-level")
})

test_that("a vector formula constant (spline knots) still predicts", {
  d <- .oc_avc()
  k <- c(40, 120)
  w <- hazard(survival::Surv(int_dead, dead) ~ splines::ns(age, knots = k),
              data = d, dist = "weibull",
              theta = c(mu = 0.01, nu = 0.5, 0.3, -0.2, 0.1))
  rows <- c(1, 50, 150)
  want <- unname(predict(w, type = "cumulative_hazard")[rows])
  nd <- data.frame(time = d$int_dead[rows], age = d$age[rows])
  expect_equal(unname(predict(w, newdata = nd, type = "cumulative_hazard")),
               want, tolerance = 1e-10)
  # A newdata column named like the knots is an unused extra.
  nd$k <- 999
  expect_equal(unname(predict(w, newdata = nd, type = "cumulative_hazard")),
               want, tolerance = 1e-10)
})

test_that("the outside-data refusal at one row, duplicate rows, one fit row", {
  d <- .oc_avc()
  n <- nrow(d)
  set.seed(3)
  zz <- rnorm(n)
  Lz <- list(z = zz, k = 3)
  cutoff <- 50
  beta <- c(0.8, 0.3)
  wb <- function(f, data = d, b = beta) {
    hazard(f, data = data, dist = "weibull", theta = c(mu = 0.01, nu = 0.5, b))
  }
  w_zz <- wb(survival::Surv(int_dead, dead) ~ zz, b = 0.3)
  w_lz <- wb(survival::Surv(int_dead, dead) ~ Lz$z, b = 0.3)
  w_ct <- wb(survival::Surv(int_dead, dead) ~ I(age > cutoff) + mal)
  truth <- function(x) drop(cbind(x$age > cutoff, x$mal) %*% beta)
  lp <- function(w, x) unname(predict(w, newdata = x, type = "linear_predictor"))

  # One row cannot be shifted, so the row-count backstop refuses.
  one <- d[5, ]
  one$zz <- zz[5]
  backstop <- paste0("has ", n, " rows for 1 row")
  expect_error(lp(w_zz, one), backstop)
  expect_error(lp(w_lz, d[5, ]), backstop)
  expect_equal(lp(w_ct, d[5, ]), truth(d[5, ]), tolerance = 1e-12)

  # Duplicate rows: a design built from newdata's columns moves with its
  # rows, so it passes; an outside term does not, however alike the rows.
  dup <- d[c(1, 1, 2, 2, 3), ]
  expect_equal(lp(w_ct, dup), truth(dup), tolerance = 1e-12)
  expect_error(lp(w_lz, d[rep(1, n), ]), "uses row-level values")

  # One fitting row: a scalar constant is still a constant.
  w1 <- wb(survival::Surv(int_dead, dead) ~ I(age > cutoff) + mal,
           data = d[1, , drop = FALSE])
  expect_equal(lp(w1, d[1:3, ]), truth(d[1:3, ]), tolerance = 1e-12)
})

test_that("the shift check scales its tolerance to each column's spread", {
  # A tolerance floored at an absolute sqrt(eps) let an outside covariate
  # that varies by less than that pass, and its fitting-order values were
  # used silently (3.47 off here). Data-only designs at extreme scales must
  # still predict.
  d <- .oc_avc()
  n <- nrow(d)
  rv <- rev(seq_len(n))
  set.seed(1)
  zt <- 0.5 + rnorm(n) * 1e-9
  w <- hazard(survival::Surv(int_dead, dead) ~ age + zt, data = d,
              dist = "weibull", theta = c(mu = 0.01, nu = 0.5, 0.01, 1e9))
  nd <- d[rv, ]
  nd$zt <- zt[rv]
  expect_error(predict(w, newdata = nd, type = "linear_predictor"),
               "term 'zt' of the model uses row-level")
  for (f in list(survival::Surv(int_dead, dead) ~ scale(age),
                 survival::Surv(int_dead, dead) ~ I(age * 1e-12),
                 survival::Surv(int_dead, dead) ~ I(age * 1e8))) {
    wf <- hazard(f, data = d, dist = "weibull",
                 theta = c(mu = 0.01, nu = 0.5, 0.1))
    expect_equal(
      unname(predict(wf, newdata = d[rv, ], type = "linear_predictor")),
      unname(predict(wf, type = "linear_predictor"))[rv], tolerance = 1e-10)
  }
})

test_that("data_vars skips `$` names and is unchanged for ordinary formulas", {
  # data_vars feeds the route rule, so ordinary formulas must record exactly
  # what they recorded at e61d602.  `cfg$time` names a list element, not the
  # data column `time`, which must not become a required variable.
  d <- .oc_avc()
  d$grp <- factor(ifelse(d$age > 100, "old", "young"))
  cutoff <- 50
  dv <- function(f, dd = d) hazard(f, data = dd, dist = "weibull")$data$x_design$data_vars
  expect_identical(dv(survival::Surv(int_dead, dead) ~ age + grp), c("age", "grp"))
  expect_identical(dv(survival::Surv(int_dead, dead) ~ I(age > cutoff) + mal),
                   c("age", "mal"))
  expect_identical(dv(survival::Surv(int_dead, dead) ~ poly(age, 2)), "age")
  d2 <- d
  d2$time <- d2$int_dead
  cfg <- list(time = 50)
  expect_identical(dv(survival::Surv(int_dead, dead) ~ I(age > cfg$time) + mal, d2),
                   c("age", "mal"))
  w <- hazard(survival::Surv(int_dead, dead) ~ I(age > cfg$time) + mal,
              data = d2, dist = "weibull",
              theta = c(mu = 0.01, nu = 0.5, 0.8, 0.3))
  expect_equal(predict(w, newdata = data.frame(age = 30, mal = 1),
                       type = "linear_predictor"), 0.3, tolerance = 1e-12)
})

test_that("the time check ignores a phase formula the fit did not use", {
  skip_on_cran()  # a multiphase fit
  d <- .oc_avc()
  set.seed(1)
  f <- suppressWarnings(hazard(
    time = d$int_dead, status = d$dead, x = cbind(age = d$age),
    dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes",
                        formula = ~ time),
      constant = hzr_phase("constant")),
    fit = TRUE))
  rows <- c(1, 50, 150)
  got <- predict(f, type = "cumulative_hazard",
                 newdata = data.frame(time = d$int_dead[rows], age = d$age[rows]))
  expect_equal(unname(got), unname(predict(f, type = "cumulative_hazard")[rows]),
               tolerance = 1e-10)
})

test_that("a multiphase global design takes the variable over its column", {
  skip_on_cran()  # a multiphase fit
  fit <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ grp, data = .dp_avc,
    dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes"),
      constant = hzr_phase("constant")),
    fit = TRUE))
  as_old <- predict(fit, type = "cumulative_hazard",
                    newdata = data.frame(time = c(1, 4), grp = "old"))
  both <- predict(fit, type = "cumulative_hazard",
                  newdata = data.frame(time = c(1, 4), grp = "old",
                                       grpyoung = 1))
  as_young <- predict(fit, type = "cumulative_hazard",
                      newdata = data.frame(time = c(1, 4), grp = "young"))
  expect_gt(max(abs(as_young / as_old - 1)), 0.05)
  expect_equal(both, as_old, tolerance = 1e-12)
})
