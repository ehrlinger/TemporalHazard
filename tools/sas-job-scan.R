#!/usr/bin/env Rscript
# sas-job-scan.R -- find SAS HAZARD jobs that widen TemporalHazard's coverage.
#
# WHY
#   tools/sas-hazard-profile.R measures a job tree in aggregate and deliberately
#   names no file, so it cannot say which jobs to run. This scan gives one row
#   per job, says whether the job can be compared with SAS at all, and picks a
#   small set that covers each feature.
#
# WHAT IT READS AND WRITES
#   Reads *.sas, the job's .lst and .log, and the existence of its datasets. It
#   opens no dataset. It writes only under --out, which must be outside the tree
#   scanned. It follows no symbolic link below a ROOT, so a job reached only
#   through one is not scanned; the summary counts the links passed over.
#
#   --out/jobs.tsv and --out/selected.tsv hold paths and stay on the server.
#   --out/summary.txt and standard output hold counts and feature names only.
#   --out/refusals.tsv tallies the translator's reasons with quoted text, paths and
#   capitalised names removed. A lower-case name can survive: read it before sharing.
#
# USAGE (on the server, never over a desktop mount)
#   Rscript tools/sas-job-scan.R --out DIR [--per-feature 2]
#     [--allow-old-listings] ROOT [ROOT ...]
#   --allow-old-listings is for a tree whose file times mean nothing, such as a
#   git checkout, where a listing can look older than the job that wrote it.
#
# READING THE RESULT
#   A feature at zero means no job showed it to this scan, not that none uses
#   it: a job that reaches HAZARD only through its own macro is counted under
#   macro_only, and its features are not seen. Scan a tree whose features are
#   known first (the hazard repository's tests/ and dist/examples/).
#
# EXIT STATUS
#   0 report written   1 usage error   2 no job with a HAZARD block found

.scan_rule <- function(feature, matrix, pattern) {
  data.frame(feature = feature, matrix = matrix, pattern = pattern)
}

# One rule per observable feature, read from the normalised text of the job's
# HAZARD blocks. `matrix` names the row of tools/feature-matrix-map.csv.
.scan_rules <- function() {
  r <- .scan_rule
  # A phase is active when its MU is written and positive, the rule the
  # translator follows (R/sas-parse-parms.R, mu_active). A statement only names
  # covariates: an intercept-only phase has none, and MUE=0 with EARLY X; has no
  # phase. A macro value (&M) cannot be read, so it counts as written.
  mu <- function(name) paste0("\\b", name, " *= *[+]?([0-9.]*[1-9]|&)")
  rbind(
    r("phase_early", "Multiphase likelihood and optimizer", mu("MUE")),
    r("phase_constant", "Multiphase likelihood and optimizer", mu("MUC")),
    r("phase_late", "Multiphase likelihood and optimizer", mu("MUL")),
    r("covar_early", "Multiphase likelihood and optimizer", "; *EARLY[ ;]"),
    r("covar_constant", "Multiphase likelihood and optimizer", "; *CONSTANT[ ;]"),
    r("covar_late", "Multiphase likelihood and optimizer", "; *LATE[ ;]"),
    r("shape_late_tau", "Phase specification and shapes (hzr_phase)", "\\bTAU *="),
    r("shape_late_gamma", "Phase specification and shapes (hzr_phase)", "\\bGAMMA *="),
    r("shape_fixed", "Phase specification and shapes (hzr_phase)",
      "\\bFIX(M|NU|DELTA|TAU|GAMMA|ALPHA|ETA|THALF|MUE|MUC|MUL)\\b"),
    r("interval_censoring", "Interval and left censoring", "; *ICENSOR[ =]"),
    r("right_censor_var", "Right censoring and data handling", "; *RCENSOR[ =]"),
    r("left_truncation", "Entry times (left truncation)", "; *LCENSOR[ =]"),
    r("weights", "Observation weights", "; *WEIGHT[ =]"),
    r("selection", "Stepwise selection (hzr_stepwise)", "; *SELECTION\\b"),
    r("selection_backward", "Stepwise selection (hzr_stepwise)",
      "; *SELECTION[^;]*\\bBACKWARD\\b"),
    r("selection_fast", "Stepwise selection (hzr_stepwise)",
      "; *SELECTION[^;]*\\bFAST\\b"),
    r("selection_move", "Stepwise selection (hzr_stepwise)",
      "; *SELECTION[^;]*\\bMOVE\\b"),
    r("restrict", "Stepwise selection (hzr_stepwise)", "; *RESTRICT\\b"),
    r("conserve", "Conservation of Events", "\\bCONSERVE\\b"),
    r("nocov", "Standard errors and the Hessian", "\\bNOCOV\\b"),
    r("robust", "Standard errors and the Hessian", "\\b(SEMI)?ROBUST\\b"),
    r("outhaz", "SAS output readers (hzr_read_outhaz and listings)", "\\bOUTHAZ *="),
    r("steepest", "Convergence and degraded-fit reporting", "\\bSTEEPEST\\b"),
    r("quasinewton", "Convergence and degraded-fit reporting",
      "\\bQUASI(NEWTON)?\\b")
  )
}

# Features seen in the whole job, not one block: the other procedures and the
# package's macros.
.scan_job_rules <- function() {
  r <- .scan_rule
  rbind(
    r("hazpred", "Prediction (predict.hazard)", "(PROC |%)HAZPRED\\b"),
    r("inhaz", "SAS output readers (hzr_read_outhaz and listings)", "\\bINHAZ *="),
    r("repeated_events", "Repeated events", "%REPEAT *[(]"),
    r("bootstrap", "Bootstrap (hzr_bootstrap)", "%(HAZBOOT|BOOTSTRAP)[A-Z_.]* *[(]"),
    r("deciles", "Goodness of fit and deciles", "%DECILES[A-Z_.]* *[(]"),
    r("chisqgf", "Goodness of fit and deciles", "%CHISQGF *[(]"),
    r("hazplot", "Prediction (predict.hazard)", "%HAZPLOT *[(]")
  )
}

.scan_features <- function() c(.scan_rules()$feature, .scan_job_rules()$feature)

# The text of the job with comments gone, and its HAZARD blocks' text. The body
# of a %MACRO HAZ... or %MACRO BOOT... definition is the library, not the job, so
# it is set aside: what is left says whether the file also runs anything.
.scan_text <- function(lines) {
  txt <- TemporalHazard:::.hzr_sas_normalise(lines)
  definition <- "%MACRO +(HAZ|BOOT).*?%MEND[^;]*;"
  defines <- grepl(definition, txt, perl = TRUE)
  if (defines) txt <- gsub(definition, " ", txt, perl = TRUE)
  blocks <- TemporalHazard:::.hzr_sas_blocks(txt)
  procs <- vapply(blocks, function(b) b$proc, character(1L))
  list(txt = txt, defines = defines, n_hazard = sum(procs == "HAZARD"),
       n_hazpred = sum(procs == "HAZPRED"), n_repeat = sum(procs == "REPEAT"),
       hazard = vapply(blocks[procs == "HAZARD"], function(b) b$text,
                       character(1L)))
}

# Whether the data of every HAZARD block in the job can be found. FALSE when any
# block's dataset is known to be gone, because the job's features are credited
# together and that block cannot be re-run. Otherwise NA when this scan cannot
# tell for any block, and TRUE only when all were found. A job with no HAZARD
# block has told this scan nothing, so it is NA.
.scan_data_found <- function(hazard, lines, job) {
  if (!length(hazard)) return(NA)
  # Comments off but case kept: the normalised text is upper case, a Unix path not.
  src <- paste(TemporalHazard:::.hzr_sas_strip_comments(lines), collapse = " ")
  found <- vapply(hazard, .scan_block_data, logical(1L), src, job,
                  USE.NAMES = FALSE)
  if (any(!found, na.rm = TRUE)) return(FALSE)
  if (anyNA(found)) NA else TRUE
}

# One block: TRUE when its dataset exists, FALSE when its library resolves and
# the dataset is not there, NA when this scan cannot tell. NA is not a pass: only
# the job's author knows what !MACRO or a %LET path means.
.scan_block_data <- function(block, src, job) {
  d <- regmatches(block, regexpr("\\bDATA *= *[A-Z0-9_]+[.][A-Z0-9_]+", block))
  if (!length(d)) return(NA)
  d <- strsplit(sub("^DATA *= *", "", d[[1L]]), ".", fixed = TRUE)[[1L]]
  m <- regmatches(src, gregexec(
    paste0("(?i)libname +", d[1L], "( +[a-z0-9]+)? *['\"]([^'\"]+)['\"]"),
    src, perl = TRUE
  ))[[1L]]
  # Two declarations of one library depend on the order the job runs them.
  if (!length(m) || length(unique(m[3L, ])) != 1L) return(NA)
  path <- m[3L, 1L]
  if (grepl("[!&%$]", path)) return(NA)
  # A job's own library can sit on a desktop mount, and a stale one hangs.
  if (startsWith(path, "/Volumes/")) return(NA)
  if (!grepl("^/", path)) path <- file.path(dirname(job), path)
  if (!dir.exists(path)) return(FALSE)
  length(list.files(path, pattern = paste0("^", d[2L], "[.]sas7bdat$"),
                    ignore.case = TRUE)) > 0L
}

# The translator's refusals in a job's calls, as reasons named by chunk. A refusal
# is a chunk that is itself stop(...), bound or bare. Searching the text cannot
# tell it from the guards a translated job carries, which hold a stop() inside an
# if().
.scan_stops <- function(calls) {
  strings <- function(x) {
    if (is.character(x)) x else if (is.call(x)) unlist(lapply(as.list(x)[-1L], strings))
  }
  rhs <- lapply(calls, function(cl) {
    if (is.call(cl) && identical(cl[[1L]], as.name("<-"))) cl[[3L]] else cl
  })
  stops <- Filter(function(r) is.call(r) && identical(r[[1L]], as.name("stop")), rhs)
  vapply(stops, function(r) paste(strings(r), collapse = " "), character(1L))
}

# "ok": a fit and nothing refused. This says the translation emits a fit, not that
# the fit runs or matches the job: n_untranslated counts what it left out. A fit
# with a prediction that refuses, or that reads a fit the translator could not
# place, is "fit_only", because the fits can still be compared. "refused": a fit or
# a step it reads refuses, or there is no fit at all.
.scan_fit_status <- function(calls, inhaz_resolved = TRUE) {
  stops <- .scan_stops(calls)
  beyond_pred <- !startsWith(as.character(names(stops)), "pred")
  if (any(beyond_pred) || !any(startsWith(as.character(names(calls)), "fit"))) {
    return("refused")
  }
  if (length(stops) || !isTRUE(inhaz_resolved)) "fit_only" else "ok"
}

# Why a job has that status. A stop that decides "refused" comes before a refused
# prediction, so a fit's refusal is not tallied under a prediction's reason.
.scan_status_reason <- function(calls, status) {
  if (identical(status, "ok")) return("")
  stops <- .scan_stops(calls)
  beyond_pred <- stops[!startsWith(as.character(names(stops)), "pred")]
  if (length(beyond_pred)) return(beyond_pred[[1L]])
  if (length(stops)) return(stops[[1L]])
  if (identical(status, "refused")) "no fit in the translation" else "INHAZ not placed"
}

# The translator prints SAS's words and the job's own names alike in capitals. Of
# the words in capitals only these are kept, so a variable or dataset name cannot
# leave, and one reason tallies as one reason whichever job gave it. A word missing
# here costs clarity, never privacy.
.scan_sas_words <- c(
  "A", "CONSTANT", "DATA", "EARLY", "EVENT", "FAST", "HAZARD", "HAZPRED", "ICENSOR",
  "INHAZ", "LATE", "MAXVARS", "MOVE", "MU", "NU", "OUTHAZ", "PARMS", "PROC",
  "REPEAT", "RESTRICT", "SAS", "SELECTION", "SLE", "SLS", "THALF", "TIME"
)

# The job's own path and file name go first, whole, because they can hold spaces
# no pattern below would cross. Then quoted text and other paths.
.scan_reason <- function(m, path = NULL) {
  for (known in c(path, if (!is.null(path)) basename(path))) {
    m <- gsub(known, "<path>", m, fixed = TRUE)
  }
  m <- gsub("`[^`]*`|'[^']*'|\"[^\"]*\"", "<x>", m)
  m <- gsub("/[^[:space:]]*", "<path>", m)
  m <- gsub("[^[:space:]]+[.]sas\\b", "<file>", m, ignore.case = TRUE)
  at <- gregexpr("(?<![A-Za-z0-9_])[A-Z0-9_]*[A-Z_][A-Z0-9_]*(?![A-Za-z0-9_])", m, perl = TRUE)
  regmatches(m, at) <- lapply(regmatches(m, at), function(w) {
    ifelse(w %in% .scan_sas_words, w, "<name>")
  })
  substr(gsub("[[:space:]]+", " ", m), 1L, 100L)
}

# Without out_dir, hzr_translate_sas() parses only and writes nothing.
.scan_translate <- function(path) {
  on.exit(setTimeLimit(), add = TRUE)
  tryCatch({
    setTimeLimit(elapsed = 60, transient = TRUE)
    j <- suppressWarnings(suppressMessages(
      TemporalHazard::hzr_translate_sas(path)
    ))
    status <- .scan_fit_status(j$calls, j$inhaz_resolved)
    list(status = status,
         reason = .scan_reason(.scan_status_reason(j$calls, status), path),
         n_untranslated = NROW(j$untranslated))
  }, error = function(e) {
    list(status = "error", reason = .scan_reason(conditionMessage(e), path),
         n_untranslated = NA_integer_)
  })
}

# One row for one file. `lst_after_job` is the check that the listing can be
# from this version of the job; `log_errors` that it finished.
.scan_job <- function(path, max_bytes = 5e6) {
  empty <- data.frame(path = path, kind = "no_block", md5 = NA_character_,
                      n_hazard = 0L, n_hazpred = 0L)
  if (isTRUE(file.size(path) > max_bytes)) {
    empty$kind <- "too_large"
    return(empty)
  }
  lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) NULL,
                    warning = function(w) NULL)
  if (is.null(lines)) {
    empty$kind <- "unreadable"
    return(empty)
  }
  s <- tryCatch(.scan_text(lines), error = function(e) NULL)
  if (is.null(s)) {
    empty$kind <- "lexer_error"
    return(empty)
  }
  empty$md5 <- unname(tools::md5sum(path))
  # A job can reach the package's macros with no block of its own, such as a
  # bootstrap driver, and it still counts.
  by_macro <- any(vapply(.scan_job_rules()$pattern, grepl, logical(1L), s$txt))
  if (s$n_hazard == 0L && s$n_hazpred == 0L && s$n_repeat == 0L && !by_macro) {
    if (s$defines) {
      empty$kind <- "macro_definition"
      return(empty)
    }
    # %HAZARD( lexes as a block, so a macro call that reaches here is the job's
    # own wrapper, whose features this scan cannot see.
    if (grepl("%[A-Z_][A-Z0-9_.]*(HAZ|BOOT)[A-Z0-9_.]* *[(]", s$txt)) {
      empty$kind <- "macro_only"
    }
    return(empty)
  }

  companion <- function(ext) {
    hit <- paste0(sub("[.][^.]*$", "", path), c(ext, toupper(ext)))
    hit <- hit[file.exists(hit)]
    if (length(hit)) hit[1L] else NA_character_
  }
  lst <- companion(".lst")
  log <- companion(".log")
  log_lines <- if (!is.na(log)) {
    tryCatch(readLines(log, warn = FALSE), error = function(e) NULL)
  }

  blocks <- paste(s$hazard, collapse = " ; ")
  flag <- function(rules, text) {
    stats::setNames(as.list(vapply(rules$pattern, grepl, logical(1L), text,
                                   USE.NAMES = FALSE)), rules$feature)
  }
  tr <- .scan_translate(path)
  cbind(
    data.frame(
      path = path, kind = "job", md5 = empty$md5, n_hazard = s$n_hazard,
      n_hazpred = s$n_hazpred,
      has_lst = !is.na(lst),
      lst_after_job = !is.na(lst) && file.mtime(lst) >= file.mtime(path),
      has_log = !is.na(log),
      log_errors = if (is.null(log_lines)) NA_integer_
                   else sum(grepl("^ERROR", log_lines)),
      data_found = .scan_data_found(s$hazard, lines, path),
      translate = tr$status, reason = tr$reason, n_untranslated = tr$n_untranslated
    ),
    as.data.frame(flag(.scan_rules(), blocks)),
    as.data.frame(flag(.scan_job_rules(), s$txt))
  )
}

# A job can be compared with SAS when a listing not older than the job exists,
# no log reports an error, and its data are not known to be gone. A missing log
# is counted in the summary, not held against the job, because many jobs keep
# only the listing. A job the translator cannot read still counts: that is a
# finding about the job.
.scan_eligible <- function(jobs, allow_old = FALSE) {
  jobs$has_lst & (jobs$lst_after_job | allow_old) &
    (is.na(jobs$log_errors) | jobs$log_errors == 0L) &
    (is.na(jobs$data_found) | jobs$data_found)
}

# A comparison needs a fit to compare, so only these jobs can be selected. A job
# without one is still in jobs.tsv, and its reason is tallied in refusals.tsv.
.scan_selectable <- function(jobs) {
  jobs$eligible & jobs$translate %in% c("ok", "fit_only")
}

# What a job's translation can show. A fit_only job's predictions refuse, so it
# proves nothing about prediction, or about reading the fit it predicts from,
# whatever its text uses. Other macros (bootstrap, deciles, chisqgf) stay credited:
# the translator never emits them, so a selected job has a fit to start from and
# the macro's step is run by hand.
.scan_credit <- function(jobs) {
  jobs[jobs$translate == "fit_only", c("hazpred", "inhaz", "hazplot")] <- FALSE
  jobs
}

# Greedy cover: repeatedly take the job that adds the most still-needed
# features, each worth more the rarer it is. Ties go to the smaller job, then
# the path, so two runs over one tree choose the same set.
.scan_select <- function(jobs, features, per_feature = 2L) {
  if (!nrow(jobs)) return(integer())
  m <- as.matrix(jobs[features]) * 1
  rarity <- 1 / pmax(colSums(m), 1)
  need <- pmin(colSums(m), per_feature)
  picked <- integer()
  open <- seq_len(nrow(jobs))
  ord <- order(jobs$n_hazard + jobs$n_hazpred, jobs$path)
  while (any(need > 0) && length(open)) {
    gain <- as.vector(m[open, need > 0, drop = FALSE] %*% rarity[need > 0])
    best <- order(-gain, match(open, ord))[1L]
    if (gain[best] <= 0) break
    i <- open[best]
    picked <- c(picked, i)
    need <- pmax(need - m[i, ], 0)
    open <- setdiff(open, i)
  }
  picked
}

.scan_write_summary <- function(jobs, n_jobs, rows, selected, path, walk) {
  features <- .scan_features()
  rule_rows <- rbind(.scan_rules(), .scan_job_rules())
  credit <- .scan_credit(jobs)
  feat <- data.frame(
    feature = features,
    jobs = vapply(features, function(f) sum(jobs[[f]]), integer(1L)),
    comparable = vapply(features, function(f) sum(jobs[[f]] & jobs$eligible),
                        integer(1L)),
    fits = vapply(features, function(f) sum(credit[[f]] & .scan_selectable(jobs)),
                  integer(1L)),
    selected = vapply(features, function(f) sum(credit[[f]][selected]), integer(1L)),
    matrix_row = rule_rows$matrix[match(features, rule_rows$feature)],
    row.names = NULL
  )
  kinds <- table(rows$kind)
  con <- file(path, "w")
  on.exit(close(con))
  w <- function(...) writeLines(paste0(...), con)
  w("TemporalHazard ", as.character(utils::packageVersion("TemporalHazard")),
    "; scanned ", as.character(Sys.Date()))
  w("symbolic links not followed: ", walk$links)
  w("directories that could not be read: ", walk$unreadable)
  w("files: ", nrow(rows))
  for (k in names(kinds)) w("  ", k, ": ", kinds[[k]])
  w("jobs: ", n_jobs, " (", nrow(jobs), " after removing identical copies)")
  w("comparable with SAS: ", sum(jobs$eligible))
  w("  no listing: ", sum(!jobs$has_lst))
  w("  listing older than job: ", sum(jobs$has_lst & !jobs$lst_after_job))
  w("  no log: ", sum(!jobs$has_log))
  w("  log with errors: ", sum(!is.na(jobs$log_errors) & jobs$log_errors > 0L))
  w("  dataset missing: ", sum(!is.na(jobs$data_found) & !jobs$data_found))
  w("  dataset unknown (kept): ", sum(is.na(jobs$data_found)))
  w("translator: ", sum(jobs$translate == "ok"), " ok, ",
    sum(jobs$translate == "fit_only"), " fit only, ",
    sum(jobs$translate == "refused"), " refused, ",
    sum(jobs$translate == "error"), " error")
  w("comparable and fits: ", sum(.scan_selectable(jobs)))
  w("  of those, with untranslated rows: ",
    sum(.scan_selectable(jobs) & !is.na(jobs$n_untranslated) & jobs$n_untranslated > 0L))
  w("translator reasons: ", length(unique(jobs$reason[jobs$translate != "ok"])),
    " (see refusals.tsv)")
  w("selected: ", length(selected))
  w("")
  w(paste(names(feat), collapse = "\t"))
  for (i in seq_len(nrow(feat))) w(paste(feat[i, ], collapse = "\t"))
  w("")
  none <- feat$feature[feat$jobs == 0L]
  gap <- feat$feature[feat$jobs > 0L & feat$comparable == 0L]
  w("seen in no job (unproven, not absent): ",
    if (length(none)) paste(none, collapse = ", ") else "none")
  w("seen, but in no comparable job: ",
    if (length(gap)) paste(gap, collapse = ", ") else "none")
  nofit <- feat$feature[feat$comparable > 0L & feat$fits == 0L]
  w("seen, but in no job that fits: ",
    if (length(nofit)) paste(nofit, collapse = ", ") else "none")
  invisible(feat)
}

# A path that may not exist yet, with its existing part's symlinks resolved.
.scan_resolve <- function(path) {
  path <- path.expand(path)
  rest <- character()
  while (!file.exists(path) && dirname(path) != path) {
    rest <- c(basename(path), rest)
    path <- dirname(path)
  }
  do.call(file.path, as.list(c(normalizePath(path), rest)))
}

# The *.sas files under the roots. list.files(recursive = TRUE) follows a link to
# a directory, which leaves the tree asked for and can hang on a stale mount, so
# this walks one directory at a time and reads a link without touching its
# target. It also says what it could not read, which list.files() reports as an
# empty directory.
.scan_walk <- function(roots) {
  files <- list()
  links <- 0L
  unreadable <- 0L
  todo <- roots
  while (length(todo)) {
    d <- todo[length(todo)]
    todo <- todo[-length(todo)]
    if (file.access(d, 5L) != 0L) {
      unreadable <- unreadable + 1L
      next
    }
    entries <- list.files(d, full.names = TRUE)
    # NA is a platform with no links, which has nothing to follow either.
    is_link <- !Sys.readlink(entries) %in% c("", NA_character_)
    links <- links + sum(is_link)
    entries <- entries[!is_link]
    is_dir <- dir.exists(entries)
    todo <- c(todo, entries[is_dir])
    files[[length(files) + 1L]] <-
      grep("[.]sas$", entries[!is_dir], ignore.case = TRUE, value = TRUE)
  }
  list(files = sort(unlist(files)), links = links, unreadable = unreadable)
}

# Returns the exit status, so a test can run it without quitting.
scan_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  fail <- function(status, ...) {
    cat(..., "\n", sep = "", file = stderr())
    status
  }
  usage <- paste("usage: Rscript sas-job-scan.R --out DIR [--per-feature N]",
                 "[--allow-old-listings] ROOT...")
  allow_old <- "--allow-old-listings" %in% args
  args <- args[args != "--allow-old-listings"]
  take <- function(flag) {
    i <- which(args == flag)
    if (!length(i) || i == length(args)) return(NULL)
    v <- args[i + 1L]
    args <<- args[-c(i, i + 1L)]
    v
  }
  out <- take("--out")
  per_arg <- take("--per-feature")
  per <- if (is.null(per_arg)) 2L else suppressWarnings(as.integer(per_arg))
  if (is.null(out) || !length(args)) return(fail(1L, usage))
  if (is.na(per) || per < 1L) return(fail(1L, "--per-feature must be at least 1"))
  # Before anything touches the path: a stale mount can hang dir.exists().
  if (any(startsWith(c(args, out), "/Volumes/"))) {
    return(fail(1L, "refusing a desktop mount: run this on the server"))
  }
  if (any(!dir.exists(args))) return(fail(1L, "a ROOT is not a directory"))
  roots <- normalizePath(args)
  if (any(startsWith(roots, "/Volumes/"))) {
    return(fail(1L, "refusing a desktop mount: run this on the server"))
  }
  # Decided before anything is created: refusing --out ROOT/new must not leave
  # the new directory in the tree it was told to read only.
  if (any(startsWith(paste0(.scan_resolve(out), "/"), paste0(roots, "/")))) {
    return(fail(1L, "--out must be outside the tree scanned"))
  }
  dir.create(out, recursive = TRUE, showWarnings = FALSE, mode = "0700")
  out <- normalizePath(out)
  Sys.chmod(out, "0700")

  walk <- .scan_walk(roots)
  files <- walk$files
  if (!length(files)) {
    return(fail(2L, "FAIL LOUD: no *.sas files under the ROOT given (",
                walk$links, " symbolic links not followed, ", walk$unreadable,
                " directories unreadable). ",
                "This is a failed scan, not an empty result."))
  }
  # Internals this scan reads. Without them every job would land in lexer_error,
  # which reads as a tree with no jobs.
  needed <- c(".hzr_sas_normalise", ".hzr_sas_blocks", ".hzr_sas_strip_comments")
  have <- vapply(needed, exists, logical(1L),
                 envir = asNamespace("TemporalHazard"), inherits = FALSE)
  if (!all(have)) {
    return(fail(1L, "this TemporalHazard is too old: it has no ",
                paste(needed[!have], collapse = ", "), ". Install a current one."))
  }
  rows <- vector("list", length(files))
  for (i in seq_along(files)) {
    rows[[i]] <- .scan_job(files[i])
    if (i %% 500L == 0L) cat(i, "of", length(files), "files\n", file = stderr())
  }
  brief <- do.call(rbind, lapply(rows, `[`, c("path", "kind")))
  all <- do.call(rbind, Filter(function(r) identical(r$kind, "job"), rows))
  if (is.null(all)) {
    return(fail(2L, "FAIL LOUD: ", length(files),
                " .sas file(s), none with a HAZARD block."))
  }
  all$eligible <- .scan_eligible(all, allow_old)
  # A copy adds no coverage. Keep a comparable copy if any copy is comparable.
  all <- all[order(!all$eligible, all$path), ]
  jobs <- all[!duplicated(all$md5), ]
  features <- .scan_features()
  pool <- which(.scan_selectable(jobs))
  selected <- pool[.scan_select(.scan_credit(jobs[pool, ]), features, per)]
  utils::write.table(jobs, file.path(out, "jobs.tsv"), sep = "\t",
                     row.names = FALSE, quote = FALSE)
  utils::write.table(.scan_credit(jobs)[selected, c("path", features, "translate", "n_untranslated")],
                     file.path(out, "selected.tsv"), sep = "\t",
                     row.names = FALSE, quote = FALSE)
  no_fit <- jobs[jobs$translate != "ok", c("translate", "reason")]
  key <- paste(no_fit$translate, no_fit$reason)
  n_key <- table(key)
  tally <- no_fit[match(names(n_key), key), ]
  tally$jobs <- as.integer(n_key)
  utils::write.table(tally[order(-tally$jobs, tally$reason), c("jobs", "translate", "reason")],
                     file.path(out, "refusals.tsv"), sep = "\t",
                     row.names = FALSE, quote = FALSE)
  .scan_write_summary(jobs, nrow(all), brief, selected,
                      file.path(out, "summary.txt"), walk)
  writeLines(readLines(file.path(out, "summary.txt")))
  0L
}

if (sys.nframe() == 0L) {
  if (!requireNamespace("TemporalHazard", quietly = TRUE)) {
    cat("TemporalHazard is not installed for this R.\n", file = stderr())
    quit(status = 1L)
  }
  quit(status = scan_main())
}
