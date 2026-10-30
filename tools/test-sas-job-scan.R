#!/usr/bin/env Rscript
# test-sas-job-scan.R -- tools/sas-job-scan.R on planted jobs.
#
# Every feature rule gets a job that uses it, and the same job with the feature
# only in a comment. The positive is what stops a zero in the real scan being
# read as absence; the negative is what stops a rule firing on comment prose.
#
# USAGE (from the package root, with the package installed)
#   Rscript tools/test-sas-job-scan.R
# Prints one line per case and "ALL n CASES PASSED"; exits 1 on any failure.

if (!exists("%||%", baseenv())) `%||%` <- function(a, b) if (is.null(a)) b else a
sys.source(file.path("tools", "sas-job-scan.R"), envir = environment())

n <- 0L
ok <- function(label, cond) {
  n <<- n + 1L
  cat(if (isTRUE(cond)) "pass " else "FAIL ", label, "\n", sep = "")
  if (!isTRUE(cond)) quit(status = 1L)
}
tree <- function() {
  d <- tempfile("scan-tree-")
  dir.create(d)
  d
}
job <- function(body = "", opts = "", stmts = "") {
  c("libname study '/nonexistent/place';",
    "%HAZARD(",
    paste("PROC HAZARD DATA=WORK.PLANTED", opts, ";"),
    "TIME T; EVENT E;",
    "PARMS MUE=0.1 THALF=1 NU=1 M=1 MUC=0.1;",
    stmts,
    ");",
    body)
}
plant <- function(dir, name, lines, lst = TRUE, log = "NOTE: done") {
  p <- file.path(dir, paste0(name, ".sas"))
  writeLines(lines, p)
  if (lst) {
    writeLines("listing", sub("[.]sas$", ".lst", p))
    Sys.setFileTime(sub("[.]sas$", ".lst", p), Sys.time() + 60)
  }
  if (!is.null(log)) writeLines(log, sub("[.]sas$", ".log", p))
  invisible(p)
}

# Each rule: what a job writes to show the feature, and whether it is an option
# on the PROC, a statement in the block, or text in the job.
uses <- list(
  phase_early = list(stmts = "EARLY X;"),
  phase_constant = list(stmts = "CONSTANT X;"),
  phase_late = list(stmts = "LATE X;"),
  shape_late_tau = list(stmts = "PARMS TAU=1;"),
  shape_late_gamma = list(stmts = "PARMS GAMMA=1;"),
  shape_fixed = list(stmts = "PARMS FIXM;"),
  interval_censoring = list(stmts = "ICENSOR C;"),
  right_censor_var = list(stmts = "RCENSOR C;"),
  left_truncation = list(stmts = "LCENSOR C;"),
  weights = list(stmts = "WEIGHT W;"),
  selection = list(stmts = "SELECTION SLE=0.1;"),
  selection_backward = list(stmts = "SELECTION BACKWARD;"),
  selection_fast = list(stmts = "SELECTION FAST;"),
  selection_move = list(stmts = "SELECTION MOVE=2;"),
  restrict = list(stmts = "SELECTION; RESTRICT X;"),
  conserve = list(opts = "CONSERVE"),
  nocov = list(opts = "NOCOV"),
  robust = list(stmts = "SELECTION ROBUST;"),
  outhaz = list(opts = "OUTHAZ=OUTEST"),
  steepest = list(opts = "STEEPEST"),
  quasinewton = list(opts = "QUASI"),
  hazpred = list(body = "%HAZPRED(PROC HAZPRED DATA=G INHAZ=OUTEST; TIME T;);"),
  inhaz = list(body = "%HAZPRED(PROC HAZPRED DATA=G INHAZ=OUTEST; TIME T;);"),
  repeated_events = list(body = "%REPEAT(IN=A, OUT=B);"),
  bootstrap = list(body = "%HAZBOOT(DATA=A);"),
  deciles = list(body = "%DECILES(DATA=A);"),
  chisqgf = list(body = "%CHISQGF(DATA=A);"),
  hazplot = list(body = "%HAZPLOT(DATA=A);")
)
ok("every rule has a planted use", setequal(names(uses), .scan_features()))

for (f in names(uses)) {
  u <- uses[[f]]
  d <- tree()
  r <- .scan_job(plant(d, "positive", job(u$body %||% "", u$opts %||% "", u$stmts %||% "")))
  ok(paste("rule fires:", f), identical(r$kind, "job") && isTRUE(r[[f]]))
  said <- paste("/*", u$opts %||% "", u$stmts %||% "", u$body %||% "", "*/")
  r <- .scan_job(plant(d, "negative", job(said)))
  ok(paste("comment does not fire:", f), identical(r$kind, "job") && isFALSE(r[[f]]))
}

# Spaces around `=` are valid SAS, and the normaliser keeps them.
d <- tree()
r <- .scan_job(plant(d, "spaced", job(
  body = "%HAZPRED(PROC HAZPRED DATA=G INHAZ = OUTEST; TIME T;);",
  opts = "OUTHAZ = OUTEST", stmts = "PARMS TAU = 1 GAMMA = 2;"
)))
ok("spaced assignments are seen",
   all(unlist(r[c("shape_late_tau", "shape_late_gamma", "outhaz", "inhaz")])))
ok("spaced DATA= is read",
   isFALSE(.scan_data_found("PROC HAZARD DATA = STUDY.X;",
                            "libname study '/nonexistent/place';",
                            file.path(d, "j.sas"))))

# Which library a job means.
lib <- tree()
file.create(file.path(lib, "x.sas7bdat"))
ok("a commented libname is not the job's",
   isTRUE(.scan_data_found("PROC HAZARD DATA=STUDY.X;",
                           c("/* libname study '/old/missing'; */",
                             paste0("libname study '", lib, "';")),
                           file.path(d, "j.sas"))))
ok("a macro reference inside a path is unknown",
   is.na(.scan_data_found("PROC HAZARD DATA=STUDY.X;",
                          "libname study \"/data/&study\";", file.path(d, "j.sas"))))
ok("a library on a desktop mount is unknown, and is not touched",
   is.na(.scan_data_found("PROC HAZARD DATA=STUDY.X;",
                          "libname study '/Volumes/stale/data';",
                          file.path(d, "j.sas"))))
ok("two declarations of one library are unknown",
   is.na(.scan_data_found("PROC HAZARD DATA=STUDY.X;",
                          c("libname study '/a';", "libname study '/b';"),
                          file.path(d, "j.sas"))))

# Kinds a job can fail into.
d <- tree()
ok("no block", identical(.scan_job(plant(d, "a", "data x; set y; run;"))$kind, "no_block"))
ok("own wrapper is macro_only",
   identical(.scan_job(plant(d, "b", "%MYHAZARD(DATA=A);"))$kind, "macro_only"))
ok("a macro library file is not a job",
   identical(.scan_job(plant(d, "c", c("%macro hazboot(data=);", "%mend;")))$kind,
             "macro_definition"))
ok("a file defining a macro and running a job is a job",
   identical(.scan_job(plant(d, "f", c("%macro hazboot(data=);", "%mend;", job())))$kind,
             "job"))
ok("a commented macro definition does not hide a job",
   identical(.scan_job(plant(d, "g", c("/* %macro hazboot(data=); */", job())))$kind,
             "job"))
ok("a bootstrap driver with no block is a job",
   identical(.scan_job(plant(d, "e", "%HAZBOOT(DATA=A);"))$kind, "job"))

# Eligibility, each condition on its own against one accepted row.
d <- tree()
r <- .scan_job(plant(d, "good", job()))
ok("known good job is comparable", isTRUE(.scan_eligible(r)))
ok("a missing dataset library fails the data check",
   isFALSE(.scan_data_found("PROC HAZARD DATA=STUDY.X;", "libname study '/nonexistent/place';",
                            file.path(d, "j.sas"))))
ok("a logical path is unknown, not missing",
   is.na(.scan_data_found("PROC HAZARD DATA=STUDY.X;", "libname study '!STUDY/data';",
                          file.path(d, "j.sas"))))
ok("no listing", isFALSE(.scan_eligible(.scan_job(plant(d, "n1", job(), lst = FALSE)))))
ok("log with an error",
   isFALSE(.scan_eligible(.scan_job(plant(d, "n2", job(), log = "ERROR: bad")))))
ok("no log is allowed", isTRUE(.scan_eligible(.scan_job(plant(d, "n3", job(), log = NULL)))))
ok("no listing, even with old listings allowed",
   isFALSE(.scan_eligible(.scan_job(plant(d, "n5", job(), lst = FALSE)), allow_old = TRUE)))
gone <- sub("WORK.PLANTED", "STUDY.GONE", job(), fixed = TRUE)
ok("a dataset known to be gone",
   isFALSE(.scan_eligible(.scan_job(plant(d, "n6", gone)))))
old <- plant(d, "n4", job())
Sys.setFileTime(sub("[.]sas$", ".lst", old), Sys.time() - 86400)
r <- .scan_job(old)
ok("old listing refused by default", isFALSE(.scan_eligible(r)))
ok("old listing allowed when asked", isTRUE(.scan_eligible(r, allow_old = TRUE)))

# Selection.
m <- data.frame(path = c("a", "b", "c", "d"), n_hazard = 1L, n_hazpred = 0L,
                f1 = c(TRUE, TRUE, FALSE, FALSE), f2 = c(FALSE, TRUE, FALSE, TRUE),
                f3 = c(FALSE, FALSE, TRUE, FALSE))
s1 <- .scan_select(m, c("f1", "f2", "f3"), 1L)
ok("selection covers every feature", all(colSums(as.matrix(m[s1, c("f1", "f2", "f3")])) >= 1L))
ok("selection is small", length(s1) == 2L)
ok("selection is stable", identical(s1, .scan_select(m, c("f1", "f2", "f3"), 1L)))
s2 <- .scan_select(m, c("f1", "f2", "f3"), 2L)
ok("per_feature is met where it can be",
   identical(unname(colSums(as.matrix(m[s2, c("f1", "f2", "f3")]))), c(2, 2, 1)))

# The whole run: read only, and nothing identifying in what may leave the server.
root <- tree()
out <- tree()
secret <- file.path(root, "STUDYSECRET")
dir.create(secret)
plant(secret, "jobsecret", job(stmts = "EARLY VARSECRET; SELECTION SLE=0.1;"))
plant(secret, "copy", job(stmts = "EARLY VARSECRET; SELECTION SLE=0.1;"))
state <- function(p) {
  f <- sort(list.files(p, recursive = TRUE, full.names = TRUE))
  list(f, unname(tools::md5sum(f)), file.mtime(f))
}
before <- state(root)
sink(tempfile())
status <- scan_main(c("--out", out, root))
sink()
ok("run exits 0", identical(status, 0L))
ok("the tree scanned is untouched", identical(before, state(root)))
ok("selected names the job", any(grepl("jobsecret|copy", readLines(file.path(out, "selected.tsv")))))
ok("an identical copy is not two jobs", nrow(utils::read.delim(file.path(out, "jobs.tsv"))) == 1L)
summ <- paste(readLines(file.path(out, "summary.txt")), collapse = "\n")
ok("summary reports the job", grepl("jobs: 2 (1 after", summ, fixed = TRUE))
ok("summary has no path, job or variable name", !grepl("SECRET|secret|/", gsub("Stepwise|hzr_[a-z_]+", "", summ)))
ok("output directory is private", identical(format(file.mode(out)), "700"))

# Refusals.
q <- function(...) {
  con <- file(tempfile(), "w")
  sink(con, type = "message")
  on.exit({
    sink(type = "message")
    close(con)
  })
  scan_main(c(...))
}
ok("--out inside the tree", identical(q("--out", file.path(root, "o"), root), 1L))
inside <- file.path(root, "new", "deeper")
ok("--out in a new directory inside the tree", identical(q("--out", inside, root), 1L))
ok("that refusal created nothing", !file.exists(file.path(root, "new")))
ok("--out on a desktop mount", identical(q("--out", "/Volumes/not-a-real-mount/o", root), 1L))
ok("--out required", identical(q(root), 1L))
ok("bad --per-feature", identical(q("--out", out, "--per-feature", "0", root), 1L))
ok("a tree with no .sas is a failed scan", identical(q("--out", tree(), tree()), 2L))
blocks <- tree()
plant(blocks, "x", "data x; run;")
ok("a tree with no job is a failed scan", identical(q("--out", tree(), blocks), 2L))
# A path that need not exist: the refusal must come before the path is touched.
ok("a desktop mount", identical(q("--out", tree(), "/Volumes/not-a-real-mount"), 1L))

cat("ALL ", n, " CASES PASSED\n", sep = "")
