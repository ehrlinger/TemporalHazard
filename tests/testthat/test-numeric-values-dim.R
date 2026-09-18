# A classed numeric column whose stored doubles are not its values is read as
# its values whether or not it has a dim (#371). A column with a dim is kept
# as a matrix, because a genuine matrix column (I(cbind(p, q)), a Surv) is a
# legitimate part of a model and flattening it would change that model.

# The class is local, so this needs no bit64 dependency: its as.double()
# returns values the stored doubles do not hold.
registerS3method("as.double", "hzr_dim_wrapped",
                 function(x, ...) attr(x, "values"))

dim_wrap <- function(v) {
  structure(matrix(9e-300, length(v), 1L), values = v,
            class = "hzr_dim_wrapped")
}

dim_fit_data <- function() {
  data(avc, package = "TemporalHazard", envir = environment())
  d <- stats::na.omit(avc)[1:120, ]
  data.frame(int_dead = d$int_dead, dead = d$dead, age = d$age)
}

dim_fit <- function(dd) {
  suppressWarnings(hazard(survival::Surv(int_dead, dead) ~ age, data = dd,
                          dist = "weibull",
                          theta = c(mu = 0.01, nu = 0.5, 0.004), fit = TRUE))
}

test_that("a classed numeric matrix column fits on its values, not its bits", {
  plain <- dim_fit_data()
  wrapped <- plain
  wrapped$age <- dim_wrap(plain$age)
  raw <- plain
  raw$age <- as.vector(unclass(dim_wrap(plain$age)))

  want <- dim_fit(plain)
  got <- dim_fit(wrapped)
  expect_equal(unname(got$fit$theta), unname(want$fit$theta),
               tolerance = 1e-12)
  expect_equal(got$fit$objective, want$fit$objective, tolerance = 1e-12)
  expect_equal(unname(got$data$x), unname(want$data$x), tolerance = 1e-12)

  # Control: reading the stored doubles gives a different fit, so the
  # assertions above cannot pass over an unconverted column.
  bits <- dim_fit(raw)
  expect_false(isTRUE(all.equal(unname(bits$fit$theta),
                                unname(want$fit$theta))))
})

test_that("a classed numeric matrix column predicts on its values", {
  plain <- dim_fit_data()
  fit <- dim_fit(plain)
  m <- suppressWarnings(hazard(
    survival::Surv(int_dead, dead) ~ age, data = plain, dist = "multiphase",
    phases = list(
      early = hzr_phase("cdf", t_half = 0.5, nu = 1, m = 1, fixed = "shapes"),
      constant = hzr_phase("constant")
    ),
    fit = TRUE
  ))
  nd <- data.frame(time = c(1, 2, 5), age = c(60, 90, 30))
  wrapped <- nd
  wrapped$age <- dim_wrap(nd$age)
  raw <- nd
  raw$age <- as.vector(unclass(dim_wrap(nd$age)))

  for (object in list(fit, m)) {
    for (type in c("survival", "cumulative_hazard")) {
      want <- predict(object, newdata = nd, type = type)
      expect_identical(predict(object, newdata = wrapped, type = type), want)
      # Control: the stored doubles give a different answer.
      expect_false(isTRUE(all.equal(
        predict(object, newdata = raw, type = type), want
      )))
    }
  }
})

test_that("a genuine matrix column keeps its shape and values", {
  # What keep_dim = TRUE exists for: flattening these would change the model.
  p <- c(1, 2, 3)
  q <- c(4, 5, 6)
  df <- data.frame(id = 1:3)
  df$m <- cbind(p = p, q = q)
  df$sv <- survival::Surv(c(5, 3, 8), c(1, 0, 1))
  out <- .hzr_numeric_frame_values(df)

  expect_identical(dim(out$m), c(3L, 2L))
  expect_identical(out$m, df$m)
  # A Surv column is a classed numeric with a dim whose stored doubles ARE
  # its values; it must keep its class, or it is no longer a response.
  expect_s3_class(out$sv, "Surv")
  expect_identical(out$sv, df$sv)
})

test_that("a classed numeric matrix column keeps its shape when converted", {
  v <- c(60, 90, 30)
  df <- data.frame(id = 1:3)
  df$w <- dim_wrap(v)
  out <- .hzr_numeric_frame_values(df)

  expect_identical(dim(out$w), c(3L, 1L))
  expect_equal(as.vector(out$w), v, tolerance = 1e-12)
  expect_false(inherits(out$w, "hzr_dim_wrapped"))
})

test_that("an integer AsIs matrix column keeps its class and storage", {
  # Its values equal its storage, only in another storage mode; comparing
  # with identical() on the raw storage called that a difference and
  # rewrote the column as a bare double matrix (#381 review).
  df <- data.frame(id = 1:3)
  df$m <- I(cbind(1:3, 4:6))
  out <- .hzr_numeric_frame_values(df)

  expect_s3_class(out$m, "AsIs")
  expect_identical(storage.mode(out$m), "integer")
  expect_identical(out$m, df$m)
})
