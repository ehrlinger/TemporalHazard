# Evaluate a translated job's emitted chunks the way rendering the .qmd
# would: in a clean environment holding only the input data. Asserting the
# SHAPE of emitted calls cannot distinguish a document that runs from one
# that errors -- every defect in #151/#152 was invisible to shape assertions
# and visible to this in one line.

# The render environment's parent. globalenv() is not one: an unrelated
# global `fit`, or any leftover symbol from another test file, satisfies a
# broken translation and this helper reports "ok" -- an oracle that ambient
# state can satisfy is not an oracle. The attached package environment sits
# AFTER globalenv() on the search path, so parenting on it reaches the
# package's exports, stats and base, and nothing of the user's.
.render_parent <- function() {
  if ("package:TemporalHazard" %in% search()) {
    return(as.environment("package:TemporalHazard"))
  }
  # Belt and braces for a load path that does not attach: rebuild the same
  # layer by hand rather than silently fall back to globalenv().
  e <- new.env(parent = baseenv())
  for (pkg in c("stats", "TemporalHazard")) {
    ns <- asNamespace(pkg)
    for (s in getNamespaceExports(pkg)) assign(s, get(s, envir = ns), envir = e)
  }
  e
}

render_sim <- function(job, data = list()) {
  env <- new.env(parent = .render_parent())
  for (nm in names(data)) assign(nm, data[[nm]], envir = env)
  res <- character(0)
  for (nm in names(job$calls)) {
    # job$calls holds quoted calls emitted by this package's own
    # hzr_translate_sas(), not arbitrary external input -- evaluating them
    # is exactly what rendering the translated .qmd does.
    out <- tryCatch({
      eval(job$calls[[nm]], env)
      "ok"
    }, error = function(e) paste0("ERROR: ", conditionMessage(e)))
    res[[nm]] <- out
  }
  # all(character(0) == "ok") is TRUE, so an empty `res` would report ok
  # having evaluated nothing -- require at least one call ran.
  list(ok = length(res) > 0 && all(res == "ok"), results = res, env = env)
}

# --- synthetic data for a translated SAS job ---------------------------
#
# render_sim(job) with nothing bound cannot exercise anything: every job
# whose PROC HAZARD read a SAS DATA= dataset opens with a guard chunk that
# stop()s when that dataset is absent, so the guard fails first and for
# every job, and the fit/grid/pred chunks are never really tested. The
# helpers below read the column names back off the emitted calls and build a
# small data frame per dataset so the guard passes and the rest actually
# runs. The fit is meaningless -- the point is that the chunks execute.

.sas_strip_assign <- function(cl) {
  if (is.call(cl) && identical(cl[[1L]], as.name("<-"))) cl[[3L]] else cl
}

# The called function of a chunk, looking through the `x <- ...` binding, so
# that classification does not depend on whether the chunk assigns. That
# matters: whether the fit chunk binds its result is exactly the defect this
# test has to be able to catch, so it must not also decide what gets tested.
.sas_head <- function(cl) {
  cl <- .sas_strip_assign(cl)
  if (is.call(cl) && is.name(cl[[1L]])) as.character(cl[[1L]]) else ""
}

# Every symbol in argument position; the function position is dropped, so
# `ifelse`, `hazard`, `~` and friends do not come back as column names.
.sas_syms <- function(x) {
  if (is.name(x)) return(as.character(x))
  if (is.call(x)) return(unlist(lapply(as.list(x)[-1L], .sas_syms)))
  character(0)
}

#' Classify a translated job's chunks.
#'
#' `eligible` says whether synthetic data can make the whole document run:
#' it needs at least one chunk that calls `hazard()`, and every `predict()`
#' must name one of those chunks. `refusals` are the chunks the translator
#' deliberately emits as `stop()` (SELECTION, LCENSOR + ICENSOR).
sas_job_shape <- function(job) {
  nms <- names(job$calls)
  heads <- vapply(job$calls, .sas_head, character(1))
  fit_chunks <- nms[heads %in% c("hazard", "hzr_read_outhaz")]
  pred_i <- which(heads == "predict")
  refs <- vapply(pred_i, function(i) as.character(.sas_strip_assign(job$calls[[i]])[[2L]]),
                 character(1))
  dangling <- setdiff(unique(refs), fit_chunks)
  reason <- if (!length(job$calls)) {
    "no calls emitted"
  } else if (!length(fit_chunks)) {
    "no PROC HAZARD block that translates to a fit"
  } else if (length(dangling)) {
    paste0("predict() references a chunk that holds no fit: ", paste(dangling, collapse = ", "))
  } else {
    ""
  }
  list(eligible = !nzchar(reason), reason = reason, heads = heads,
       fit_chunks = fit_chunks, hazard_chunks = nms[heads == "hazard"],
       refusals = nms[heads == "stop"], preds = nms[pred_i])
}

#' Build a synthetic data frame per dataset the job's emitted calls name.
#'
#' Column roles are read off the argument they appear in -- `time=`,
#' `status=`, `time_lower=`, `weights=`, and the phase formulas -- so the
#' values are at least of the right kind (positive times, 0/1 status, an
#' entry time below the event time).
#'
#' Only the job's SAS *input* datasets are built. A `newdata=` prediction
#' grid is deliberately NOT: the document is supposed to build its own grid,
#' so fabricating one hands the job the very object it fails to create and
#' counts it as rendered when it cannot run. That is how this oracle went
#' vacuous the last time -- a job whose grid the translator refused rendered
#' anyway, on data no chunk of it produced.
sas_synth_data <- function(job, n = 24L) {
  nms <- names(job$calls)
  heads <- vapply(job$calls, .sas_head, character(1))

  roles <- list()
  add <- function(d, slot, v) {
    v <- setdiff(unique(v), c(".hzr_status", d))
    if (!length(v)) return(invisible(NULL))
    r <- roles[[d]]
    if (is.null(r)) r <- list()
    r[[slot]] <- unique(c(r[[slot]], v))
    roles[[d]] <<- r
  }

  for (k in seq_along(job$calls)) {
    rhs <- .sas_strip_assign(job$calls[[k]])
    if (grepl("^status(_[0-9]+)?$", nms[[k]])) {
      # ICENSOR hoists status into its own chunk, derived into the dataset
      # with transform(<data>, .hzr_status = .) so hazard()'s data mask
      # cannot shadow it; with no DATA= it is a bare local binding.
      d <- if (identical(heads[[k]], "transform")) as.character(rhs[[2L]]) else ""
      add(d, "status", .sas_syms(rhs))
    }
    if (identical(heads[[k]], "hazard")) {
      a <- as.list(rhs)[-1L]
      d <- if (!is.null(a$data)) as.character(a$data) else ""
      add(d, "time", .sas_syms(a$time))
      add(d, "status", .sas_syms(a$status))
      add(d, "lower", .sas_syms(a$time_lower))
      add(d, "wt", .sas_syms(a$weights))
      add(d, "cov", .sas_syms(a$phases))
    }
  }

  time_col <- seq_len(n) * (60 / n)
  covs <- function(i) sin(seq_len(n) * i)
  mk <- function(r) {
    cols <- list()
    i <- 0L
    put <- function(v, gen) {
      for (nm in v) {
        if (is.null(cols[[nm]])) {
          i <<- i + 1L
          cols[[nm]] <<- gen(i)
        }
      }
    }
    put(r$time, function(i) time_col)
    # One residue class per count variable, so no two of them ever fire on
    # the same row and one class is left over where none does. EVENT,
    # ICENSOR's C3 and RCENSOR's C2 are SEPARATE counts: a row where two
    # fire is two observations at once, which the translator refuses at fit
    # time (#157, #162). The patterns here used to be `%% (i %% 3L + 2L)`,
    # which coincide at every common multiple -- that generates data the
    # refusal is right to reject, and the rejection then reads as a
    # rendering failure, testing this generator rather than the translation.
    # Disjoint classes still reach every branch of the emitted ifelse(),
    # which is what the old comment was after.
    n_status <- length(r$status)
    status_i <- 0L
    put(r$status, function(i) {
      status_i <<- status_i + 1L
      as.numeric(seq_len(n) %% (n_status + 1L) == status_i - 1L)
    })
    put(r$lower, function(i) time_col * 0.25)
    put(r$wt, function(i) rep(1, n))
    put(r$cov, covs)
    cols
  }

  data <- list()
  free <- list()
  for (d in names(roles)) {
    cols <- mk(roles[[d]])
    if (!nzchar(d)) {
      free <- c(free, cols)
    } else {
      if (is.null(cols[["time"]])) cols[["time"]] <- time_col
      data[[d]] <- as.data.frame(cols)
    }
  }

  c(data, free)
}
