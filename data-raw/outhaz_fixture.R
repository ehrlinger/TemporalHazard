# NOTE: This is a maintainer-only script, NOT part of the installed package.
# It lives in data-raw/ (which is .Rbuildignore'd) so it is never shipped to
# CRAN, checked by R CMD check, or reachable by users. It regenerates
# inst/extdata/outhaz-fixture.rds, the fixture behind test-read-outhaz.R.

# outhaz_fixture.R - a synthetic SAS `outhaz` dataset
#
# PURPOSE
# -------
# `PROC HAZARD`'s `outhaz=` dataset is the full-precision parity reference: it
# stores converged estimates and the asymptotic variance-covariance matrix at
# double precision where the printed .lst carries about seven significant
# figures. `hzr_read_outhaz()` parses it, so its test needs a file in exactly
# that layout.
#
# WHY THE NUMBERS ARE INVENTED
# ----------------------------
# The layout below was mapped from a real `outhaz` dataset (a four-parameter
# early+constant fit). That file contains no PHI -- it is model parameters
# only -- but it belongs to an unpublished study, and this package is public.
# So the STRUCTURE is reproduced faithfully and the VALUES are invented. Every
# free parameter and every vcov entry is derived from a reciprocal of an
# irrational or a prime, so that each carries a full double mantissa and no
# reader can mistake it for a fit.
#
# THE LAYOUT, and what `hzr_read_outhaz()` relies on
# --------------------------------------------------
#   * 17 rows: six model-structure flags, then eleven parameters.
#   * Columns: _NAME_, _EST_, _STATUS_, then one column per parameter. The
#     parameter columns hold the vcov, addressed by name in both directions.
#   * Flag rows carry NA in _STATUS_ and in every parameter column. That NA is
#     how the reader separates flags from parameters -- there is no type
#     column to ask.
#   * _STATUS_ is 1 for a free parameter and 0 for a fixed one, and is stored
#     as a double, not an integer.
#   * Fixed parameters occupy a row and a column of the vcov block, filled
#     with zeros. Only the free-by-free submatrix is meaningful.
#   * The two triangles of the vcov are written independently, so the block is
#     symmetric to about 1e-18 but NOT bit-identical across the diagonal. That
#     asymmetry is reproduced deliberately: it is why the test compares with
#     expect_equal() and not identical().

make_outhaz_fixture <- function(path = "inst/extdata/outhaz-fixture.rds") {
  flags <- c(G1FLAG = 2, FIXDEL0 = 1, FIXMNU1 = 0,
             G3FLAG = 0, FIXGE2 = 0, FIXGAE2 = 0)

  params <- c("DELTA", "THALF", "NU", "M", "TAU", "GAMMA",
              "ALPHA", "ETA", "E0", "C0", "L0")
  free <- c("THALF", "NU", "E0", "C0")

  # Invented estimates. Reciprocals of irrationals: full mantissas, obviously
  # not a fit. Fixed parameters are zero, as SAS writes them.
  est <- stats::setNames(numeric(length(params)), params)
  est[["THALF"]] <- pi / 100
  est[["NU"]]    <- sqrt(2)
  est[["E0"]]    <- -exp(1)
  est[["C0"]]    <- -log(10) / log(2)

  status <- stats::setNames(rep(0, length(params)), params)
  status[free] <- 1

  # A positive-definite vcov over the free parameters, as L %*% t(L) so it
  # cannot accidentally be indefinite.
  lt <- matrix(0, 4, 4, dimnames = list(free, free))
  diag(lt) <- c(1 / pi, 1 / exp(1), 1 / sqrt(7), 1 / sqrt(11))
  lt[2, 1] <- 1 / 13
  lt[3, 1] <- -1 / 17
  lt[3, 2] <- 1 / 19
  lt[4, 1] <- 1 / 23
  lt[4, 2] <- -1 / 29
  lt[4, 3] <- 1 / 31
  vc <- lt %*% t(lt)

  # Reproduce SAS's independently-written triangles (see header).
  vc["THALF", "NU"] <- vc["NU", "THALF"] + 3.469446951953614e-18

  # Assemble in the column order SAS emits.
  block <- matrix(0, length(params), length(params),
                  dimnames = list(params, params))
  block[free, free] <- vc

  d <- data.frame(
    `_NAME_`   = c(names(flags), params),
    `_EST_`    = c(unname(flags), unname(est)),
    `_STATUS_` = c(rep(NA_real_, length(flags)), unname(status)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  for (p in params) {
    d[[p]] <- c(rep(NA_real_, length(flags)), unname(block[, p]))
  }

  stopifnot(nrow(d) == 17, ncol(d) == 3 + length(params))
  saveRDS(d, path)
  invisible(d)
}
