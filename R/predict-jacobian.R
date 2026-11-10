#' Derivatives behind a prediction's standard errors
#'
#' **Experimental.** Returns the Jacobian that [predict.hazard()] used for
#' the standard errors of a `se.fit = TRUE` result, row for row with the
#' result as you hold it now: after selecting, reordering or repeating its
#' rows. With it, any weighted sum of the predictions has a delta-method
#' standard error, such as the expected number of events in a cohort, which
#' is the sum of the cumulative hazards at each patient's own time.
#'
#' @details
#' The derivatives are stored on the result when [predict.hazard()] returns
#' it. R keeps that stored copy unchanged when you select or bind rows, so it
#' can describe rows the result no longer has, or miss rows it now has. This
#' function matches each row of `object` to the row of the prediction it came
#' from, by its row name (a repeated row, named `"2.1"` by `[`, matches row
#' `"2"`), and checks that its `fit` and `se.fit` are exactly the values that
#' row was given. It stops, naming the first row it cannot match, rather than
#' return derivatives for the wrong rows. Results bound together from more
#' than one call to `predict()`, or whose `fit` or `se.fit` columns were
#' changed or rounded, are refused: compute their derivatives from each
#' original result and combine the pieces by `names`.
#'
#' The columns are the parameters the standard errors used. A parameter held
#' fixed in the fit without a variance, or one the prediction does not depend
#' on (such as a coefficient whose covariate is 0 in every row), is left out,
#' so two predictions can have different columns. A Weibull fit's first
#' column is `log(mu)`, not `mu`.
#'
#' @param object A data frame returned by [predict.hazard()] with
#'   `se.fit = TRUE`, possibly with its rows selected, reordered or repeated.
#' @return A list with three elements: `J`, a matrix with one row per row of
#'   `object` and one column per parameter, holding the derivative of that
#'   row's `fit` with respect to the parameter; `vcov`, the covariance matrix
#'   of those parameters; and `names`, the parameters' names as
#'   [summary.hazard()] names its coefficient rows. For weights `w`, the
#'   standard error of `sum(w * object$fit)` is
#'   `sqrt(t(g) %*% vcov %*% g)` with `g = colSums(w * J)`.
#' @seealso [predict.hazard()]
#' @export
#' @examples
#' set.seed(1)
#' dat <- data.frame(time = rexp(60, 0.3), status = rbinom(60, 1, 0.7),
#'                   age = rnorm(60))
#' fit <- hazard(survival::Surv(time, status) ~ age, data = dat,
#'               dist = "weibull", theta = c(0.3, 1, 0), fit = TRUE)
#' # Each patient's cumulative hazard at their own follow-up time.
#' res <- predict(fit, newdata = dat, type = "cumulative_hazard",
#'                se.fit = TRUE)
#' j <- hzr_jacobian(res)
#' g <- colSums(j$J)
#' c(expected = sum(res$fit), se = sqrt(drop(t(g) %*% j$vcov %*% g)))
hzr_jacobian <- function(object) {
  jac <- attr(object, "jacobian", exact = TRUE)
  if (!is.data.frame(object) || is.null(jac)) {
    stop("'object' carries no derivatives: pass a result of predict() with ",
         "se.fit = TRUE whose standard errors could be computed.",
         call. = FALSE)
  }
  if (is.null(jac$fit) || is.null(jac$se.fit)) {
    stop("'object' carries derivatives in an older form that cannot be ",
         "matched to its rows; predict() again.", call. = FALSE)
  }
  if (!all(c("fit", "se.fit") %in% names(object))) {
    stop("'object' has no 'fit' and 'se.fit' columns, so its rows cannot ",
         "be matched to the derivatives; keep both columns.", call. = FALSE)
  }
  rn <- rownames(object)
  jr <- rownames(jac$J)
  idx <- match(rn, jr)
  # `[` names a repeated row "2.1", "2.2"; it came from row "2".
  rep_row <- is.na(idx)
  idx[rep_row] <- match(sub("\\.[0-9]+$", "", rn[rep_row]), jr)
  same <- function(a, b) (a == b & !is.na(a) & !is.na(b)) |
    (is.na(a) & is.na(b))
  ok <- !is.na(idx)
  ok[ok] <- same(object$fit[ok], jac$fit[idx[ok]]) &
    same(object$se.fit[ok], jac$se.fit[idx[ok]])
  if (!all(ok)) {
    bad <- which(!ok)[1L]
    stop("Row ", bad, " of 'object' (row name \"", rn[bad], "\") does not ",
         "match a row of the prediction its derivatives came from: rows ",
         "were bound from another call to predict(), renamed, or their ",
         "'fit' or 'se.fit' changed. ", sum(!ok), " of ", length(ok),
         " rows do not match. Take the derivatives from each original ",
         "result.", call. = FALSE)
  }
  J <- jac$J[idx, , drop = FALSE]
  rownames(J) <- rn
  list(J = J, vcov = jac$vcov, names = jac$names)
}
