# Measure hzr_translate_sas()'s PROC HAZPRED grids on the HAZPRED binary (#494).
#
# PROC HAZPRED predicts at every row of its DATA= dataset: the time is the
# column its TIME statement names (hazpred/timeprc.c:16-25) and each covariate
# is the column of the same name, which must be there
# (hazpred/hazpred.c:153-173). The translator emits R that builds that dataset
# from the job's own DATA steps. This script runs the emitted grid chunk for
# hp.death.AVC.hm1 and hm2, hands the result to the SAS/C HAZPRED binary with
# the job's own INHAZ= estimates, and compares the predictions with the ones
# SAS printed in the job's .lst. It then runs the three defects #494 names:
#
#   (a) the grid without its SET DESIGN covariates;
#   (b) TIME YEARS where YEARS = MONTHS/12;
#   (c) hm1's first definition of PREDICT, not its last.
#
# Run from the package root, with TemporalHazard loaded from this checkout:
#
#   Rscript -e 'devtools::load_all(quiet = TRUE); source("data-raw/hazpred-grid-oracle.R")'
#
# The binary reads its input from $TMPDIR/hzp.J<id>.X<ix>.{sas,dta,haz} and
# writes .out, as the %HAZPRED macro arranges (hazpred.sas:74-119).

repo <- path.expand(Sys.getenv("HAZARD_REPO", "~/Documents/GitHub/hazard"))
bin <- file.path(repo, "dist", "bin", "hazpred")
inhaz_file <- file.path(repo, "examples", "sasest", "hmdeath.sas7bdat")
stopifnot(file.exists(bin), file.exists(inhaz_file))
inhaz <- as.data.frame(haven::read_sas(inhaz_file))
work <- tempfile("hazpred-grid-oracle-")
dir.create(work)

run_hazpred <- function(grid, stmts, tag) {
  pre <- file.path(work, sprintf("hzp.J494.X%s", tag))
  haven::write_xpt(grid, paste0(pre, ".dta"), version = 5, name = "HZRCALL")
  haven::write_xpt(inhaz, paste0(pre, ".haz"), version = 5, name = "HZRCALL")
  writeLines(c(
    paste0("(PROC HAZPRED DATA=PREDICT INHAZ=EXAMPLES.HMDEATH OUT=PREDICT; ",
           stmts, " )"),
    sprintf("OBSCOUNT %d ;", nrow(grid)), sprintf("HAZCOUNT %d ;", nrow(inhaz)),
    "JOBID J494;", sprintf("JOBIX X%s;", tag)), paste0(pre, ".sas"))
  system(paste0("TMPDIR=", shQuote(work), " ", shQuote(bin), " < ",
                shQuote(paste0(pre, ".sas")), " > ", shQuote(paste0(pre, ".lst")),
                " 2>&1"))
  log <- readLines(paste0(pre, ".lst"), warn = FALSE)
  out <- tryCatch(as.data.frame(haven::read_xpt(paste0(pre, ".out"))),
                  error = function(e) NULL)
  list(log = log, out = out)
}

# The digital nomogram SAS printed: MONTHS and _SURVIV, in listing order.
listed <- function(lst) {
  x <- iconv(readLines(lst, warn = FALSE), "latin1", "ASCII", sub = " ")
  x <- x[seq(grep("Obs +MONTHS +YEARS +OPYEAR", x)[[1L]], length(x))]
  f <- strsplit(trimws(x), " +")
  f <- f[lengths(f) == 17L & grepl("^[0-9]+$", vapply(f, `[`, "", 1L))]
  data.frame(MONTHS = as.numeric(vapply(f, `[`, "", 2L)),
             SURVIV = as.numeric(vapply(f, `[`, "", 12L)))
}

emitted_grid <- function(sas, n = 1L) {
  job <- suppressWarnings(hzr_translate_sas(sas))
  cl <- job$calls[[grep("^grid", names(job$calls), value = TRUE)[[n]]]]
  env <- new.env(parent = baseenv())
  eval(cl, env)
  g <- get(as.character(cl[[2L]]), envir = env)
  g$time <- NULL # the R-side copy of TIME; HAZPRED reads the named column
  list(job = job, grid = g)
}

for (nm in c("hm1", "hm2")) {
  sas <- file.path(repo, "tests", sprintf("hp.death.AVC.%s.sas", nm))
  e <- emitted_grid(sas)
  r <- run_hazpred(e$grid, "TIME MONTHS;", toupper(nm))
  want <- listed(sub("[.]sas$", ".lst", sas))
  got <- r$out[r$out$DIGITAL == 1, ]
  stopifnot(nrow(want) > 0L, nrow(got) == nrow(want))
  cat(sprintf(paste("%s: emitted grid %d rows (%d digital); HAZPRED on it",
                    "vs SAS listing: max |MONTHS diff| %.2g, max |_SURVIV diff| %.2g\n"),
              nm, nrow(e$grid), nrow(got), max(abs(got$MONTHS - want$MONTHS)),
              max(abs(got[["_SURVIV"]] - want$SURVIV))))
}

# (a) The grid without the covariates SET DESIGN brings in, as main emitted it.
sas1 <- file.path(repo, "tests", "hp.death.AVC.hm1.sas")
g1 <- emitted_grid(sas1)$grid
ra <- run_hazpred(g1[c("MONTHS", "YEARS", "DIGITAL")], "TIME MONTHS;", "A")
cat("(a) no SET DESIGN columns:",
    paste(grep("ERROR|EXITING", ra$log, value = TRUE), collapse = "\n    "), "\n")

# (b) TIME YEARS predicts at the YEARS column, not at MONTHS.
rm_ <- run_hazpred(g1, "TIME MONTHS;", "BM")
ry <- run_hazpred(g1, "TIME YEARS;", "BY")
i <- which(rm_$out$DIGITAL == 1 & rm_$out$MAL == 1 & rm_$out$MONTHS == 12)
cat(sprintf("(b) MAL=1, MONTHS=12 (YEARS=1): _SURVIV %.5f under TIME MONTHS, %.5f under TIME YEARS\n",
            rm_$out[["_SURVIV"]][i], ry$out[["_SURVIV"]][i]))

# (c) hm1's first PREDICT (the log grid alone) vs its last (log + DIGITAL).
cat(sprintf("(c) hm1 grid rows: %d with the DIGITAL rows SAS's last PREDICT adds; %d in the first definition\n",
            nrow(g1), sum(g1$DIGITAL == 0)))
