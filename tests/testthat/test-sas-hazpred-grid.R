# #494: a PROC HAZPRED grid is the job's DATA steps, not the first DO loop.
#
# PROC HAZPRED predicts at every row of its DATA= dataset, reading the time
# from the TIME variable (hazpred/timeprc.c:16-25) and each covariate from
# its column (hazpred/hazpred.c:153-173). Every test below runs the emitted
# grid chunk and asserts the rows, covariates and times it actually holds.

translate_494 <- function(lines, env = parent.frame()) {
  f <- withr::local_tempfile(fileext = ".sas", .local_envir = env)
  writeLines(lines, f)
  suppressWarnings(hzr_translate_sas(f))
}

# Run a job's n-th grid chunk the way the rendered document does, and return
# the data frame it binds.
grid_494 <- function(job, n = 1L) {
  slots <- grep("^grid", names(job$calls), value = TRUE)
  cl <- job$calls[[slots[[n]]]]
  env <- new.env(parent = baseenv())
  eval(cl, env)
  get(as.character(cl[[2L]]), envir = env)
}

hazpred_494 <- function(data, time = "MONTHS", out = data) {
  sprintf("%%HAZPRED( PROC HAZPRED DATA=%s INHAZ=E.H OUT=%s; TIME %s; );",
          data, out, time)
}

test_that("SET DESIGN carries the design covariates into the grid (#494a)", {
  job <- translate_494(c(
    "DATA DESIGN; AGE=6; OPMOS=180; OP_AGE=OPMOS*AGE; DO MAL=0,1; OUTPUT; END;",
    "DATA PREDICT; SET DESIGN; DIGITAL=0; DO MONTHS=1,6,12; OUTPUT; END;",
    hazpred_494("PREDICT")
  ))
  g <- grid_494(job)
  # SAS writes each DESIGN row's loop together: MAL=0 at 1, 6, 12, then MAL=1.
  expect_equal(nrow(g), 6L)
  expect_equal(g$MAL, c(0, 0, 0, 1, 1, 1))
  expect_equal(g$AGE, rep(6, 6))
  expect_equal(g$OP_AGE, rep(1080, 6))
  expect_equal(g$DIGITAL, rep(0, 6))
  expect_equal(g$time, rep(c(1, 6, 12), 2))
  expect_equal(nrow(job$untranslated), 0L)
})

test_that("predict() runs at the SET DESIGN covariates, not at zero (#494a)", {
  skip_on_cran()
  set.seed(494)
  AVCS <- data.frame(INT_DEAD = stats::rexp(300, 0.2),
                     DEAD = rep(c(1, 1, 0), length.out = 300),
                     MAL = rep(c(0, 1), each = 150))
  job <- translate_494(c(
    "%HAZARD( PROC HAZARD DATA=AVCS CONDITION=14 OUTHAZ=E.H;",
    "EVENT DEAD; TIME INT_DEAD; PARMS MUE=0.2 THALF=0.15 NU=1 MUC=0.0005;",
    "EARLY MAL; );",
    "DATA DESIGN; MAL=1;",
    "DATA PREDICT; SET DESIGN; DO MONTHS=1,6,12; OUTPUT; END;",
    hazpred_494("PREDICT")
  ))
  sim <- suppressWarnings(render_sim(job, list(AVCS = AVCS)))
  expect_true(sim$ok, info = paste(names(sim$results), sim$results, collapse = " | "))
  # This fit masks three free variances as non-positive, so the HAZPRED
  # chunk's limits are withheld, with a warning (#586); its estimates stand.
  expect_warning(pred <- eval(job$calls$pred, sim$env),
                 "no variance for an estimated parameter", fixed = TRUE)
  got <- pred$fit
  fit <- get("fit", envir = sim$env)
  at_mal1 <- predict(fit, newdata = data.frame(time = c(1, 6, 12), MAL = 1),
                     type = "survival")
  at_zero <- predict(fit, newdata = data.frame(time = c(1, 6, 12), MAL = 0),
                     type = "survival")
  # The known positive: MAL moves the prediction, so the comparison can fail.
  expect_gt(max(abs(at_mal1 - at_zero)), 1e-3)
  expect_equal(got, at_mal1)
})

test_that("TIME names the grid variable predict() reads as time (#494b)", {
  job <- translate_494(c(
    "DATA PREDICT; DO MONTHS=6,12,24; YEARS=MONTHS/12; OUTPUT; END;",
    hazpred_494("PREDICT", time = "YEARS")
  ))
  g <- grid_494(job)
  expect_equal(g$MONTHS, c(6, 12, 24))
  expect_equal(g$time, c(0.5, 1, 2))
})

test_that("a time derived in a later step is the one predicted at (#494b)", {
  job <- translate_494(c(
    "DATA PREDICT; DO MONTHS=6,12,24; OUTPUT; END;",
    "DATA PREDICT; SET PREDICT; YEARS=MONTHS/12;",
    hazpred_494("PREDICT", time = "YEARS")
  ))
  expect_equal(grid_494(job)$time, c(0.5, 1, 2))
})

test_that("a grid is the last definition before the HAZPRED, as in SAS (#494c)", {
  job <- translate_494(c(
    "DATA PREDICT; DO MONTHS=1,2; OUTPUT; END;",
    "DATA DIGITAL; DO MONTHS=5,6,7; OUTPUT; END;",
    "DATA PREDICT; SET PREDICT DIGITAL; YEARS=MONTHS/12;",
    hazpred_494("PREDICT"),
    # Written after the HAZPRED above, so only the second HAZPRED reads it.
    "DATA PREDICT; DO MONTHS=9; OUTPUT; END;",
    hazpred_494("PREDICT")
  ))
  g1 <- grid_494(job, 1L)
  expect_equal(g1$time, c(1, 2, 5, 6, 7))
  expect_equal(g1$YEARS, c(1, 2, 5, 6, 7) / 12)
  expect_equal(grid_494(job, 2L)$time, 9)
})

test_that("SET A B stacks A over B, missing where only one has a variable (#494c)", {
  job <- translate_494(c(
    "DATA A; DIGITAL=0; DO MONTHS=1,2; OUTPUT; END;",
    "DATA B; DO MONTHS=3; OUTPUT; END;",
    "DATA PREDICT; SET A B;",
    hazpred_494("PREDICT")
  ))
  g <- grid_494(job)
  expect_equal(g$time, c(1, 2, 3))
  expect_equal(g$DIGITAL, c(0, 0, NA))
})

test_that("a HAZPRED OUT= keeps its DATA= rows for a later grid (#494c)", {
  job <- translate_494(c(
    "DATA DESIGN; DO MAL=0,1; OUTPUT; END;",
    "DATA PREDICT; SET DESIGN; DO MONTHS=1,2; OUTPUT; END;",
    hazpred_494("PREDICT", out = "RESULT"),
    "DATA NEXT; SET RESULT; S1=_SURVIV;",
    hazpred_494("NEXT")
  ))
  g1 <- grid_494(job, 1L)
  g2 <- suppressWarnings(grid_494(job, 2L))
  expect_equal(nrow(g2), 4L)
  expect_equal(g2$MAL, g1$MAL)
  expect_equal(g2$time, g1$time)
  # The first HAZPRED's _SURVIV is not carried: S1 is NA and recorded.
  expect_equal(g2$S1, rep(NA_real_, 4))
  expect_true(any(job$untranslated$construct == "S1=_SURVIV"))
})

test_that("PROC SORT reorders the grid it sorts (#494)", {
  job <- translate_494(c(
    "DATA PREDICT; DO MONTHS=3,1; OUTPUT; END;",
    "DATA DIGITAL; DO MONTHS=2; OUTPUT; END;",
    "DATA PREDICT; SET PREDICT DIGITAL;",
    "PROC SORT DATA=PREDICT; BY MONTHS;",
    hazpred_494("PREDICT")
  ))
  expect_equal(grid_494(job)$time, c(1, 2, 3))
})

test_that("PROC SORT puts a missing key first, as SAS does (#494)", {
  # SAS orders a numeric missing value below every number, so an ascending
  # sort puts it first; R's order() puts NA last by default.
  job <- translate_494(c(
    "DATA A; DO MONTHS=2,1; K=MONTHS; OUTPUT; END;",
    "DATA B; MONTHS=3;",
    "DATA PREDICT; SET A B;",
    "PROC SORT DATA=PREDICT; BY K;",
    hazpred_494("PREDICT")
  ))
  g <- grid_494(job)
  expect_equal(g$K, c(NA, 1, 2))
  expect_equal(g$time, c(3, 1, 2))
})

test_that("two OUTPUTs write each input row twice, in SAS's order (#494)", {
  job <- translate_494(c(
    "DATA BASE; DO B=1,2; OUTPUT; END;",
    "DATA PREDICT; SET BASE; MONTHS=6; A=1; OUTPUT; A=2; B2=5; OUTPUT;",
    hazpred_494("PREDICT")
  ))
  g <- grid_494(job)
  expect_equal(g$B, c(1, 1, 2, 2))
  expect_equal(g$A, c(1, 2, 1, 2))
  # B2 is first set after the first OUTPUT, so that row has it missing.
  expect_equal(g$B2, c(NA, 5, NA, 5))
})

test_that("an untranslated DATA-step statement is recorded, warned and NA (#494)", {
  job <- translate_494(c(
    "DATA DESIGN; AGE=6; RISK=0; IF AGE>5 THEN RISK=1; X=ROUND(AGE);",
    "DATA PREDICT; SET DESIGN; DO MONTHS=1,2; OUTPUT; END;",
    hazpred_494("PREDICT")
  ))
  # The grid is still emitted and predicted over.
  expect_identical(as.character(job$calls$pred[[1L]]), "predict")
  u <- job$untranslated
  expect_true(any(u$construct == "IF AGE>5 THEN RISK=1" & grepl("RISK, NA", u$reason)))
  expect_true(any(u$construct == "X=ROUND(AGE)" & grepl("calls ROUND()", u$reason, fixed = TRUE)))
  expect_warning(g <- grid_494(job), "RISK")
  expect_equal(g$RISK, c(NA_real_, NA_real_))
  expect_equal(g$X, c(NA_real_, NA_real_))
  # The translated neighbours are unaffected.
  expect_equal(g$AGE, c(6, 6))
})

test_that("a value carried over from the previous loop pass is NA (#494)", {
  # SAS gives PREV = 0, 1, 2. Run once per row, the assignment would read
  # LAST = 0 on every row, a wrong value rather than a missing one.
  job <- translate_494(c(
    "DATA PREDICT; LAST=0; DO MONTHS=1,2,3; PREV=LAST; LAST=MONTHS; OUTPUT; END;",
    hazpred_494("PREDICT")
  ))
  g <- suppressWarnings(grid_494(job))
  expect_equal(g$PREV, rep(NA_real_, 3))
  expect_equal(g$LAST, c(1, 2, 3))
  expect_true(any(grepl("previous loop pass", job$untranslated$reason)))
})

refused_494 <- function(job) {
  heads <- vapply(job$calls, function(x) as.character(x[[1L]])[[1L]], "")
  expect_false("predict" %in% heads)
  expect_true("stop" %in% heads)
  job$untranslated$reason[grepl("^DATA=", job$untranslated$construct)]
}

test_that("a TIME variable the grid does not carry refuses it (#494b)", {
  job <- translate_494(c(
    "DATA PREDICT; DO MONTHS=1,2; OUTPUT; END;",
    hazpred_494("PREDICT", time = "YEARS")
  ))
  expect_match(refused_494(job), "TIME variable YEARS", fixed = TRUE)
})

test_that("a HAZPRED with no TIME statement refuses its grid (#494b)", {
  job <- translate_494(c(
    "DATA PREDICT; DO MONTHS=1,2; OUTPUT; END;",
    "%HAZPRED( PROC HAZPRED DATA=PREDICT INHAZ=E.H OUT=P; );"
  ))
  expect_match(refused_494(job), "no TIME statement", fixed = TRUE)
})

test_that("a grid another procedure rewrote is refused, not read from before it (#494c)", {
  job <- translate_494(c(
    "DATA PREDICT; DO MONTHS=1,2; OUTPUT; END;",
    "PROC MEANS DATA=PREDICT NOPRINT; OUTPUT OUT=PREDICT MEAN=;",
    hazpred_494("PREDICT")
  ))
  expect_match(refused_494(job), "PREDICT is written by PROC MEANS", fixed = TRUE)
})

test_that("a loop that changes a variable after its OUTPUT is refused (#494)", {
  # X is 1 on every row but the first in SAS; no per-row statement says so.
  job <- translate_494(c(
    "DATA PREDICT; X=0; DO MONTHS=1,2; OUTPUT; X=1; END;",
    hazpred_494("PREDICT")
  ))
  expect_match(refused_494(job), "OUTPUT or DELETE inside a loop", fixed = TRUE)
})

test_that("a macro call inside a grid's DATA step is refused (#494)", {
  job <- translate_494(c(
    "DATA PREDICT; DO MONTHS=1,2; OUTPUT; END; %ADDCOLS;",
    hazpred_494("PREDICT")
  ))
  expect_match(refused_494(job), "`%ADDCOLS`", fixed = TRUE)
})

test_that("a DO WHILE list is refused, not run to its TO bound (#494)", {
  job <- translate_494(c(
    "DATA PREDICT; DO MONTHS=1 TO 12 WHILE(MONTHS<6); OUTPUT; END;",
    hazpred_494("PREDICT")
  ))
  expect_match(refused_494(job), "calls WHILE()", fixed = TRUE)
})

test_that("a loop bound read from SET data is refused (#494)", {
  # Each row would run its own loop; the first row's bound is not all rows'.
  # T is also base R's TRUE, so folding it without checking that it is a
  # constant of the step would read T*1 as 1.
  job <- translate_494(c(
    "DATA DESIGN; DO T=2,3; OUTPUT; END;",
    "DATA PREDICT; SET DESIGN; DO MONTHS=1 TO T*1; OUTPUT; END;",
    hazpred_494("PREDICT")
  ))
  expect_match(refused_494(job), "not a constant set earlier", fixed = TRUE)
  # The same bound through an assignment: U is data, not a constant.
  job <- translate_494(c(
    "DATA DESIGN; DO T=2,3; OUTPUT; END;",
    "DATA PREDICT; SET DESIGN; U=T*1; DO MONTHS=1 TO U; OUTPUT; END;",
    hazpred_494("PREDICT")
  ))
  expect_match(refused_494(job), "not a constant set earlier", fixed = TRUE)
})

test_that("a range SAS runs no times is refused; a descending one runs (#494)", {
  job <- translate_494(c(
    "DATA PREDICT; DO MONTHS=5 TO 1; OUTPUT; END;",
    hazpred_494("PREDICT")
  ))
  expect_match(refused_494(job), "runs no times", fixed = TRUE)
  down <- translate_494(c(
    "DATA PREDICT; DO MONTHS=5 TO 1 BY -2; OUTPUT; END;",
    hazpred_494("PREDICT")
  ))
  expect_equal(grid_494(down)$time, c(5, 3, 1))
})

test_that("statements that decide rows this cannot read are refused (#494)", {
  cases <- list(
    # An unconditional DELETE writes no rows at all.
    c("DATA P1; DO MONTHS=1,2; OUTPUT; END;", "DATA PREDICT; SET P1; DELETE;"),
    # Two loops each write rows; SAS gives 1, 2, 3.
    c("DATA PREDICT; DO MONTHS=1,2; OUTPUT; END; DO MONTHS=3; OUTPUT; END;"),
    # A conditional OUTPUT decides which rows exist.
    c("DATA P1; DO MONTHS=1,2; OUTPUT; END;", "DATA PREDICT; SET P1; IF MONTHS=1 THEN OUTPUT;"),
    # A DROP list this cannot read may drop a covariate.
    c("DATA PREDICT; A1=1; A2=2; DO MONTHS=1,2; OUTPUT; END; DROP A1-A2;"),
    # A DATA statement writing two datasets.
    c("DATA PREDICT; DO MONTHS=1; OUTPUT; END;", "DATA PREDICT OTHER; DO MONTHS=2; OUTPUT; END;"),
    # PROC SORT NODUPKEY deletes rows as well as sorting them.
    c("DATA PREDICT; DO MONTHS=1,1,2; OUTPUT; END;", "PROC SORT DATA=PREDICT NODUPKEY; BY MONTHS;")
  )
  for (k in seq_along(cases)) {
    job <- translate_494(c(cases[[k]], hazpred_494("PREDICT")))
    heads <- vapply(job$calls, function(x) as.character(x[[1L]])[[1L]], "")
    expect_false("predict" %in% heads, info = paste("case", k))
    expect_true(any(grepl("^DATA=", job$untranslated$construct)), info = paste("case", k))
  }
})

test_that("a subsetting IF or IF ... THEN DELETE refuses the grid (#494)", {
  # Both decide which rows exist. Keeping them all predicted at rows SAS
  # deleted: hp.dthip.PAIVS.time gave 24504 rows where SAS wrote 23483.
  for (del in c("IF MONTHS<3;", "IF MONTHS=3 THEN DELETE;",
                "IF MONTHS<3 THEN X=1; ELSE DELETE;")) {
    job <- translate_494(c(
      "DATA PREDICT; DO MONTHS=1,2,3; OUTPUT; END;",
      paste("DATA PREDICT; SET PREDICT;", del),
      hazpred_494("PREDICT")
    ))
    expect_match(refused_494(job), "decides which rows exist", fixed = TRUE, info = del)
  }
})

test_that("a variable read but never set is a missing column, as in SAS (#494)", {
  # SAS puts every variable a step names in the dataset, missing where never
  # set. hp.death.COMPARISON reads _PLAD50 without setting it.
  job <- translate_494(c(
    "DATA PREDICT; DO MONTHS=1,2; X=Y*2; OUTPUT; END;",
    "DATA NEXT; SET PREDICT; IF Z>1 THEN W=1; V=U*0; U=0;",
    hazpred_494("NEXT")
  ))
  g <- suppressWarnings(grid_494(job))
  # V reads U before U is set, so V is missing; U is then set.
  for (v in c("X", "Y", "Z", "W", "V")) {
    expect_equal(g[[v]], c(NA_real_, NA_real_), info = v)
  }
  expect_equal(g$U, c(0, 0))
})

test_that("PROC SQL INSERT, DELETE, UPDATE and ALTER rewrite the grid, and refuse it (#494)", {
  for (dml in c("INSERT INTO PREDICT SET MONTHS=9;", "DELETE FROM PREDICT WHERE MONTHS=1;",
                "UPDATE PREDICT SET MONTHS=9;", "ALTER TABLE PREDICT DROP MONTHS;",
                "insert into predict set months=9;")) {
    job <- translate_494(c(
      "DATA PREDICT; DO MONTHS=1,2; OUTPUT; END;",
      paste("PROC SQL;", dml, "QUIT;"),
      hazpred_494("PREDICT")
    ))
    expect_match(refused_494(job), "PREDICT is written by PROC SQL", fixed = TRUE, info = dml)
  }
})

test_that("a WORK. libref names the same dataset; another libref does not (#494)", {
  # SAS's one-level names live in WORK, so WORK.PREDICT is PREDICT.
  job <- translate_494(c(
    "DATA WORK.PREDICT; DO MONTHS=1,2; OUTPUT; END;",
    "DATA PREDICT; SET WORK.PREDICT; YEARS=MONTHS/12;",
    hazpred_494("WORK.PREDICT")
  ))
  expect_equal(grid_494(job)$time, c(1, 2))
  expect_equal(grid_494(job)$YEARS, c(1, 2) / 12)
  # A HAZPRED reading WORK.PREDICT passes its rows on through OUT=.
  job <- translate_494(c(
    "DATA PREDICT; DO MONTHS=1,2; OUTPUT; END;",
    hazpred_494("WORK.PREDICT", out = "RESULT"),
    "DATA NEXT; SET RESULT;",
    hazpred_494("NEXT")
  ))
  expect_equal(grid_494(job, 2L)$time, c(1, 2))
  # A later write through WORK. is the last definition, and refuses the grid.
  for (w in c("PROC APPEND BASE=WORK.PREDICT DATA=X; RUN;",
              "PROC MEANS DATA=PREDICT; OUTPUT OUT=WORK.PREDICT; RUN;",
              "PROC SQL; CREATE TABLE WORK.PREDICT AS SELECT * FROM PREDICT; QUIT;")) {
    job <- translate_494(c(
      "DATA PREDICT; DO MONTHS=1,2; OUTPUT; END;", "DATA X; MONTHS=9;", w,
      hazpred_494("PREDICT")
    ))
    expect_match(refused_494(job), "PREDICT is written by", fixed = TRUE, info = w)
  }
  # SASUSER.PREDICT is another dataset: writing it leaves PREDICT alone.
  job <- translate_494(c(
    "DATA PREDICT; DO MONTHS=1,2; OUTPUT; END;",
    "PROC APPEND BASE=SASUSER.PREDICT DATA=PREDICT; RUN;",
    hazpred_494("PREDICT")
  ))
  expect_equal(grid_494(job)$time, c(1, 2))
})

test_that("a value an IF changes inside a loop is NA for the whole loop (#494)", {
  # SAS gives Y = 0, 0, 5: C carries the IF's 5 into the next pass.
  job <- translate_494(c(
    "DATA PREDICT; C=0; DO MONTHS=1,2,3; Y=C; IF MONTHS>1 THEN C=5; OUTPUT; END;",
    hazpred_494("PREDICT")
  ))
  g <- suppressWarnings(grid_494(job))
  expect_equal(g$Y, rep(NA_real_, 3))
  expect_true(any(job$untranslated$construct == "Y=C" &
                    grepl("previous loop pass", job$untranslated$reason)))
})

test_that("LOG of 0, EXP overflow and division by zero are missing, as in SAS (#494)", {
  job <- translate_494(c(
    "DATA PREDICT; DO MONTHS=0,1; L=LOG(MONTHS); E=EXP(1000*MONTHS); R=1/MONTHS; OUTPUT; END;",
    hazpred_494("PREDICT")
  ))
  g <- grid_494(job)
  expect_equal(g$L, c(NA, 0))
  expect_equal(g$E, c(1, NA))
  expect_equal(g$R, c(NA, 1))
})

test_that("PROC SQL CREATE TABLE and PROC APPEND rewrite the grid, and refuse it (#494)", {
  cases <- list(
    c("DATA PREDICT; DO MONTHS=1,2,3; OUTPUT; END;",
      "PROC SQL; CREATE TABLE PREDICT AS SELECT * FROM PREDICT WHERE MONTHS<2; QUIT;"),
    c("DATA PREDICT; DO MONTHS=1,2; OUTPUT; END;", "DATA X; MONTHS=9;",
      "PROC APPEND BASE=PREDICT DATA=X; RUN;"),
    c("DATA PREDICT; DO MONTHS=1,2; OUTPUT; END;", "PROC DATASETS; DELETE PREDICT; RUN;",
      "DATA OLD; MONTHS=9;", "PROC DATASETS; CHANGE OLD=PREDICT; RUN;")
  )
  for (k in seq_along(cases)) {
    job <- translate_494(c(cases[[k]], hazpred_494("PREDICT")))
    expect_match(refused_494(job), "PREDICT is written by", fixed = TRUE, info = paste("case", k))
  }
})

test_that("statements after a step's only OUTPUT do not reach the grid (#494)", {
  job <- translate_494(c(
    "DATA PREDICT; MONTHS=1; OUTPUT; MONTHS=2;",
    hazpred_494("PREDICT")
  ))
  expect_equal(grid_494(job)$time, 1)
})

test_that("a grid SET from a dataset no DATA step builds is refused (#494a)", {
  job <- translate_494(c(
    "DATA PREDICT; SET COHORT; DO MONTHS=1,2; OUTPUT; END;",
    hazpred_494("PREDICT")
  ))
  expect_match(refused_494(job), "reads COHORT", fixed = TRUE)
})

# The public corpus, against SAS's own listing of the same grid: the
# digital-nomogram rows each job prints after PROC HAZPRED ran.
lst_rows_494 <- function(lst) {
  # The printer plots use box-drawing bytes that are not valid UTF-8.
  x <- iconv(readLines(lst, warn = FALSE), "latin1", "ASCII", sub = " ")
  # The nomogram is the PROC PRINT whose header begins with these columns;
  # the DESIGN listing before it prints numbers too, with another header.
  x <- x[seq(grep("Obs +MONTHS +YEARS +OPYEAR", x)[[1L]], length(x))]
  f <- strsplit(trimws(x), " +")
  f <- f[lengths(f) == 17L & grepl("^[0-9]+$", vapply(f, `[`, "", 1L))]
  s <- do.call(rbind, lapply(f, function(r) r[1:11]))
  colnames(s) <- c("OBS", "MONTHS", "YEARS", "OPYEAR", "OPMOS", "AGE", "COM_IV",
                   "MAL", "INC_SURG", "ORIFICE", "STATUS")
  s
}

# Each listed value is exact to half a unit in its last printed place, and
# pages print to different places (0.16427 on one, 0.164 on the next).
expect_listed_494 <- function(actual, printed, info) {
  places <- nchar(sub("^[^.]*[.]?", "", printed))
  gap <- abs(actual - as.numeric(printed))
  expect_true(all(gap <= 0.5 * 10^-places + 1e-9),
              info = paste(info, "worst", format(max(gap - 0.5 * 10^-places))))
}

test_that("hp.death.AVC.hm1 and hm2 grids match SAS's listing (#494)", {
  skip_on_cran()
  repo <- path.expand(Sys.getenv("HAZARD_REPO", "~/Documents/GitHub/hazard"))
  for (job_name in c("hm1", "hm2")) {
    sas <- file.path(repo, "tests", paste0("hp.death.AVC.", job_name, ".sas"))
    lst <- sub("[.]sas$", ".lst", sas)
    skip_if_not(file.exists(sas) && file.exists(lst), "hazard checkout not available")
    job <- suppressWarnings(hzr_translate_sas(sas))
    expect_equal(nrow(job$untranslated), 0L, info = job_name)
    g <- grid_494(job)
    want <- lst_rows_494(lst)
    # hm1 lists 29 DIGITAL rows for each MAL, hm2 16 for each COM_IV.
    expect_equal(nrow(want), if (job_name == "hm1") 58L else 32L, info = job_name)
    dig <- g[g$DIGITAL == 1, ]
    expect_equal(nrow(dig), nrow(want), info = job_name)
    for (v in setdiff(colnames(want), "OBS")) {
      expect_listed_494(dig[[v]], want[, v], info = paste(job_name, v))
    }
    # Every row, not only the digital ones, is predicted at TIME MONTHS.
    expect_equal(g$time, g$MONTHS, info = job_name)
    # The log grid is there too: 101 rows for each design row.
    expect_equal(sum(g$DIGITAL == 0), 202L, info = job_name)
  }
})

# Macro scope (1.2.13 release review). A %MACRO ... %MEND body is a
# definition: its DATA steps run where the macro is called, not where they
# stand. And a call this cannot read, between the grid's DATA step and the
# PROC HAZPRED, may rewrite the grid.
fit_macro_grid <- c(
  "%HAZARD( PROC HAZARD DATA=AVC OUTHAZ=OUTEST; TIME INT_DEAD; EVENT DEAD;",
  "PARMS MUE=0.3504743 THALF=0.1905077 NU=1.437416 M=1 FIXM MUC=4.391673E-07;",
  "EARLY AGE; );",
  "DATA PRED; AGE=50; DO INT_DEAD=1,2,3; OUTPUT; END; RUN;"
)
hazpred_macro_grid <- "%HAZPRED( PROC HAZPRED DATA=PRED INHAZ=OUTEST OUT=OUT; TIME INT_DEAD; );"
macro_rows <- function(job) job$untranslated[grepl("^%", job$untranslated$construct), ]

test_that("a DATA step inside a macro that is never called is not the grid", {
  job <- translate_494(c(
    fit_macro_grid,
    "%MACRO ALTGRID; DATA PRED; AGE=70; DO INT_DEAD=1,2,3; OUTPUT; END; RUN;",
    "%MEND ALTGRID;",
    hazpred_macro_grid
  ))
  g <- grid_494(job)
  # SAS never runs ALTGRID, so it predicts at AGE = 50 (main emitted 70).
  expect_equal(g$AGE, c(50, 50, 50))
  expect_equal(g$time, c(1, 2, 3))
  expect_equal(nrow(macro_rows(job)), 0L)
})

test_that("calling a macro whose body writes the grid refuses the grid", {
  job <- translate_494(c(
    fit_macro_grid,
    "%MACRO ALTGRID; DATA PRED; AGE=70; DO INT_DEAD=1,2,3; OUTPUT; END; RUN;",
    "%MEND ALTGRID;",
    "%ALTGRID;",
    hazpred_macro_grid
  ))
  expect_match(refused_494(job), "PRED is written by %ALTGRID", fixed = TRUE)
})

test_that("a call before the grid's step is recorded; macro statements add no row", {
  # Approved change (2026-10-09): the calls before the grid's step each get a
  # row, since they may set session state the grid reads. %LET, %PUT and a
  # macro whose body writes nothing add none.
  job <- translate_494(c(
    "%INCLUDE 'lib.sas'; %MAKEGRID;",
    fit_macro_grid,
    "%LET X=1; %PUT &X;",
    "%MACRO QUIET; %PUT HELLO; %MEND QUIET;",
    "%QUIET;",
    hazpred_macro_grid
  ))
  expect_equal(macro_rows(job)$construct, c("%INCLUDE 'LIB.SAS'", "%MAKEGRID"))
  expect_equal(grid_494(job)$AGE, c(50, 50, 50))
})

test_that("a job with no macros records no macro row", {
  job <- translate_494(c(fit_macro_grid, hazpred_macro_grid))
  expect_equal(nrow(macro_rows(job)), 0L)
  expect_equal(grid_494(job)$AGE, c(50, 50, 50))
})

test_that("a PROC HAZPRED inside a macro definition is refused", {
  # SAS runs the procedure where the macro is called, so its grid is
  # whatever PRED holds at the call. Here the call follows a second DATA
  # PRED with AGE = 70; main read the grid where the text stands and
  # emitted AGE = 50 with no row. The translation does not follow calls, so
  # it refuses. Approved change to #621's test, which asserted the grid
  # built inside the macro, and emitted it for a macro never called.
  job <- translate_494(c(
    fit_macro_grid,
    "%MACRO P;",
    hazpred_macro_grid,
    "%MEND P;",
    "DATA PRED; AGE=70; DO INT_DEAD=1,2,3; OUTPUT; END; RUN;",
    "%P;"
  ))
  expect_match(refused_494(job), "inside a %MACRO definition", fixed = TRUE)
  expect_false(any(grepl("^grid", names(job$calls))))

  # The same when the definition builds the grid itself, and when the macro
  # is never called: the refusal does not depend on what the body holds.
  for (call in c("%ALTPRED;", "")) {
    job <- translate_494(c(
      fit_macro_grid,
      "%MACRO ALTPRED; DATA PRED; AGE=70; DO INT_DEAD=1,2,3; OUTPUT; END; RUN;",
      hazpred_macro_grid,
      "%MEND ALTPRED;",
      call
    ))
    expect_match(refused_494(job), "inside a %MACRO definition", fixed = TRUE)
  }
})

test_that("a PROC HAZPRED outside a definition still emits the grid (control)", {
  # The definition above holds no PROC HAZPRED, so nothing is refused for
  # scope, and the uncalled macro's DATA step is not the grid.
  job <- translate_494(c(
    fit_macro_grid,
    "%MACRO P; %PUT HELLO; %MEND P;",
    "%P;",
    hazpred_macro_grid
  ))
  expect_equal(grid_494(job)$AGE, c(50, 50, 50))
  expect_false(any(grepl("%MACRO", job$untranslated$reason, fixed = TRUE)))
})

# A macro that sorts the grid in place rewrites it, as a DATA step in the
# body does: PROC SORT NODUPKEY drops the duplicate INT_DEAD = 1 row, and SAS
# predicts on 1, 2, 3. Main counted only DATA names and OUT= in a body, so
# the grid was emitted with all four rows and no row.
fit_sort_grid <- c(fit_macro_grid[1:3],
                   "DATA PRED; AGE=50; DO INT_DEAD=1,1,2,3; OUTPUT; END; RUN;")

test_that("calling a macro that sorts the grid in place refuses the grid", {
  job <- translate_494(c(
    fit_sort_grid,
    "%MACRO S; PROC SORT DATA=PRED NODUPKEY; BY INT_DEAD; RUN; %MEND S;",
    "%S;",
    hazpred_macro_grid
  ))
  expect_match(refused_494(job), "PRED is written by %S", fixed = TRUE)
  # A plain sort refuses too: the call is opaque, as a DATA step in it is.
  job <- translate_494(c(
    fit_sort_grid,
    "%MACRO S; PROC SORT DATA=WORK.PRED; BY DESCENDING INT_DEAD; RUN; %MEND S;",
    "%S;",
    hazpred_macro_grid
  ))
  expect_match(refused_494(job), "PRED is written by %S", fixed = TRUE)
})

test_that("a macro sort that writes elsewhere, or is never called, leaves the grid", {
  # OUT= names the dataset written; DATA= is then only read.
  job <- translate_494(c(
    fit_sort_grid,
    "%MACRO S; PROC SORT DATA=PRED OUT=SORTED NODUPKEY; BY INT_DEAD; RUN; %MEND S;",
    "%S;",
    hazpred_macro_grid
  ))
  expect_equal(grid_494(job)$INT_DEAD, c(1, 1, 2, 3))
  expect_equal(nrow(macro_rows(job)), 0L)
  # An uncalled macro sorts nothing.
  job <- translate_494(c(
    fit_sort_grid,
    "%MACRO S; PROC SORT DATA=PRED NODUPKEY; BY INT_DEAD; RUN; %MEND S;",
    hazpred_macro_grid
  ))
  expect_equal(grid_494(job)$INT_DEAD, c(1, 1, 2, 3))
  expect_equal(nrow(macro_rows(job)), 0L)
  # A sort of another dataset is no write of the grid.
  job <- translate_494(c(
    fit_sort_grid,
    "%MACRO S; PROC SORT DATA=OTHER NODUPKEY; BY INT_DEAD; RUN; %MEND S;",
    "%S;",
    hazpred_macro_grid
  ))
  expect_equal(grid_494(job)$INT_DEAD, c(1, 1, 2, 3))
  expect_equal(nrow(macro_rows(job)), 0L)
})

test_that("a macro sort with no DATA= refuses the grid", {
  # PROC SORT without DATA= sorts the most recent dataset, which depends on
  # where the macro is called. This cannot say which, so the call may rewrite
  # the grid, as an %INCLUDE may. Approved change to #621's test, which
  # recorded the call and emitted the grid.
  job <- translate_494(c(
    fit_sort_grid,
    "%MACRO S; PROC SORT NODUPKEY; BY INT_DEAD; RUN; %MEND S;",
    "%S;",
    hazpred_macro_grid
  ))
  expect_match(refused_494(job), "`%S` calls a macro whose body", fixed = TRUE)
  expect_match(refused_494(job), "may rewrite PRED", fixed = TRUE)
})

test_that("an open-code sort of the grid is unchanged by the macro rule (control)", {
  # Outside a macro, a plain BY sort is translated and NODUPKEY is refused,
  # as before.
  job <- translate_494(c(
    fit_sort_grid,
    "PROC SORT DATA=PRED; BY AGE INT_DEAD; RUN;",
    hazpred_macro_grid
  ))
  expect_equal(grid_494(job)$INT_DEAD, c(1, 1, 2, 3))
  job <- translate_494(c(
    fit_sort_grid,
    "PROC SORT DATA=PRED NODUPKEY; BY INT_DEAD; RUN;",
    hazpred_macro_grid
  ))
  expect_match(refused_494(job), "PRED is written by PROC SORT", fixed = TRUE)
})

# The grid allow-list (1.2.13 release review). The grid is emitted only when
# every statement that can change it is one the translation reads: three
# reviews each found a construct the old deny-list let through, and every
# one emitted a wrong grid with no $untranslated row. Each case below did
# that on main.
grid_rows_allow <- function(job) {
  job$untranslated[grepl("^DATA=", job$untranslated$construct), ]
}

test_that("a macro call with no semicolon refuses the grid it swallowed (T1)", {
  # `%SETUP` then `DATA PRED;` is one statement, `%SETUP DATA PRED`. SAS runs
  # the macro and then the DATA step, so it predicts at AGE = 70. Main lost
  # the DATA statement and emitted the earlier grid, AGE = 50, with no row.
  job <- translate_494(c(
    fit_macro_grid,
    "%MACRO SETUP; OPTIONS LS=80; %MEND SETUP;",
    "%SETUP",
    "DATA PRED; AGE=70; DO INT_DEAD=1,2,3; OUTPUT; END; RUN;",
    hazpred_macro_grid
  ))
  expect_match(refused_494(job), "%SETUP DATA PRED", fixed = TRUE)
  expect_false(any(grepl("^grid", names(job$calls))))
  # The same for a sort it swallows: SAS sorts PRED descending.
  job <- translate_494(c(
    fit_macro_grid,
    "%MACRO SETUP; OPTIONS LS=80; %MEND SETUP;",
    "%SETUP",
    "PROC SORT DATA=PRED; BY DESCENDING INT_DEAD; RUN;",
    hazpred_macro_grid
  ))
  expect_match(refused_494(job), "%SETUP PROC SORT DATA=PRED", fixed = TRUE)
})

test_that("a merged statement before the grid's steps refuses it too (T1)", {
  # The statement the parser could not read stands before PRED is built, but
  # the parser has lost its place in the job, so nothing after it is trusted.
  job <- translate_494(c(
    "%MACRO SETUP; OPTIONS LS=80; %MEND SETUP;",
    "%SETUP",
    "DATA OTHER; X=1; RUN;",
    fit_macro_grid,
    hazpred_macro_grid
  ))
  expect_match(refused_494(job), "%SETUP DATA OTHER", fixed = TRUE)
})

test_that("a WHERE statement in a PROC SORT of the grid refuses it (T2)", {
  # SAS keeps the rows with INT_DEAD > 1 (2 rows); main read only BY and
  # emitted all 3.
  for (sort in c("PROC SORT DATA=PRED0 OUT=PRED; WHERE INT_DEAD > 1; BY INT_DEAD; RUN;",
                 "PROC SORT DATA=PRED0 OUT=PRED; BY INT_DEAD; WHERE INT_DEAD > 1; RUN;")) {
    job <- translate_494(c(
      fit_macro_grid[1:3],
      "DATA PRED0; AGE=50; DO INT_DEAD=3,1,2; OUTPUT; END; RUN;",
      sort,
      hazpred_macro_grid
    ))
    expect_match(refused_494(job), "WHERE INT_DEAD > 1", fixed = TRUE, info = sort)
  }
})

test_that("a macro body's PROC DATASETS that replaces the grid refuses it (T3)", {
  # SAS deletes PRED and renames PRED2 to PRED, so it predicts at AGE = 70.
  # Main saw no name the body writes and emitted AGE = 50 with no row.
  job <- translate_494(c(
    fit_macro_grid,
    "DATA PRED2; AGE=70; DO INT_DEAD=1,2,3; OUTPUT; END; RUN;",
    "%MACRO SW; PROC DATASETS LIB=WORK NOLIST; DELETE PRED; CHANGE PRED2=PRED; RUN; QUIT;",
    "%MEND SW;",
    "%SW;",
    hazpred_macro_grid
  ))
  expect_match(refused_494(job), "PRED is written by %SW", fixed = TRUE)
})

test_that("a grid built inside an open-code %IF or %DO is refused (T6)", {
  # %IF 0 never runs, so SAS predicts at AGE = 50; main applied the step
  # unconditionally and emitted AGE = 60.
  for (cond in c("%IF 0 %THEN %DO;", "%DO I = 1 %TO 0;")) {
    job <- translate_494(c(
      fit_macro_grid,
      cond,
      "DATA PRED; AGE=60; DO INT_DEAD=1,2,3; OUTPUT; END; RUN;",
      "%END;",
      hazpred_macro_grid
    ))
    expect_match(refused_494(job), "PRED is built inside an open-code %IF/%DO block",
                 fixed = TRUE, info = cond)
  }
  # With the PROC HAZPRED inside the block too, no macro statement stands
  # between the grid's step and the PROC HAZPRED: the step's place in the
  # block is what refuses it.
  job <- translate_494(c(
    fit_macro_grid,
    "%IF &RUN %THEN %DO;",
    "DATA PRED; AGE=60; DO INT_DEAD=1,2,3; OUTPUT; END; RUN;",
    hazpred_macro_grid,
    "%END;"
  ))
  expect_match(refused_494(job), "PRED is built inside an open-code %IF/%DO block",
               fixed = TRUE)
  # A %IF block that only stands between the grid's step and the PROC
  # HAZPRED refuses too: whatever it runs may rewrite the grid.
  job <- translate_494(c(
    fit_macro_grid,
    "%IF 1 %THEN %DO; %PUT HELLO; %END;",
    hazpred_macro_grid
  ))
  expect_match(refused_494(job), "%IF 1 %THEN %DO", fixed = TRUE)
})

test_that("an %INCLUDE or an unread macro call after the grid refuses it", {
  # Approved change to #621's test, which recorded each call and still
  # emitted the grid the job's DATA steps show. Each call below may rewrite
  # PRED, and the maintainer prefers a loud refusal to a silent wrong grid.
  for (call in c("%INCLUDE 'makegrid.sas';", "%MAKEGRID;", "%MAKEGRID2(DS=PRED, AGE=70);")) {
    job <- translate_494(c(fit_macro_grid, call, hazpred_macro_grid))
    expect_match(refused_494(job), "may rewrite PRED", fixed = TRUE, info = call)
    expect_false(any(grepl("^grid", names(job$calls))), info = call)
  }
})

test_that("a procedure this cannot read, after the grid, refuses it", {
  # PROC COPY writes the members it copies, which no option names.
  job <- translate_494(c(
    fit_macro_grid,
    "PROC COPY IN=SAVED OUT=WORK; SELECT PRED; RUN;",
    hazpred_macro_grid
  ))
  expect_match(refused_494(job), "PROC COPY", fixed = TRUE)
})

test_that("a DATA step that runs other code, after the grid, refuses it", {
  for (stmt in c("CALL EXECUTE('DATA PRED; AGE=70; RUN;')",
                 "RC = DOSUBL('DATA PRED; AGE=70; RUN;')")) {
    job <- translate_494(c(
      fit_macro_grid,
      paste0("DATA _NULL_; ", stmt, "; RUN;"),
      hazpred_macro_grid
    ))
    expect_match(refused_494(job), "runs other SAS code", fixed = TRUE, info = stmt)
  }
})

test_that("a statement outside every step refuses the grid, wherever it stands", {
  job <- translate_494(c(
    "AGE=70;",
    fit_macro_grid,
    hazpred_macro_grid
  ))
  expect_match(refused_494(job), "`AGE=70`", fixed = TRUE)
})

test_that("OPTIONS OBS= or FIRSTOBS= before the PROC HAZPRED refuses the grid", {
  for (opt in c("OPTIONS OBS=2;", "OPTIONS NODATE FIRSTOBS=2;")) {
    job <- translate_494(c(opt, fit_macro_grid, hazpred_macro_grid))
    expect_match(refused_494(job), "OBS", fixed = TRUE, info = opt)
  }
  # OBS=MAX is the default, and reads every row.
  job <- translate_494(c("OPTIONS OBS=MAX;", fit_macro_grid, hazpred_macro_grid))
  expect_equal(grid_494(job)$AGE, c(50, 50, 50))
})

test_that("a SET that names no dataset refuses the grid", {
  # `SET;` reads the most recent dataset; main read it as no SET at all.
  job <- translate_494(c(
    fit_macro_grid[1:3],
    "DATA BASE; AGE=60; RUN;",
    "DATA PRED; SET; DO INT_DEAD=1,2,3; OUTPUT; END; RUN;",
    hazpred_macro_grid
  ))
  expect_match(refused_494(job), "SET", fixed = TRUE)
})

test_that("statements the allow-list reads leave the grid as it was (control)", {
  # Global statements, %LET and %PUT, a procedure that writes nothing, a
  # DATA step that writes another dataset, a macro whose body writes nothing,
  # and a call before the grid's step: none of them can change PRED.
  job <- translate_494(c(
    "%INCLUDE 'lib.sas';",
    "%PLOT( ID L=\"x\", END; LABELX L=\"y\", END; );",
    fit_macro_grid,
    "TITLE1 'A grid'; FOOTNOTE 'n'; OPTIONS LS=132 PS=60 NODATE;",
    "LIBNAME EX ('!HZEXAMPLES/sasest'); FILENAME F ('x');",
    "%LET X=1; %PUT &X;",
    "PROC PRINT DATA=PRED; VAR AGE; RUN;",
    "DATA OTHER; SET PRED; Y=AGE*2; RUN;",
    "%MACRO QUIET; %PUT HELLO; TITLE2 'q'; %MEND QUIET;",
    "%QUIET;",
    hazpred_macro_grid
  ))
  expect_equal(grid_494(job)$AGE, c(50, 50, 50))
  expect_equal(grid_494(job)$time, c(1, 2, 3))
  expect_equal(nrow(grid_rows_allow(job)), 0L)
  # Only the two calls before the grid's step are recorded (2026-10-09).
  expect_equal(macro_rows(job)$construct,
               c("%INCLUDE 'LIB.SAS'", "%PLOT( ID L=\"X\", END; LABELX L=\"Y\", END; )"))
})

test_that("a numeric LENGTH below 8 in the grid's step refuses it", {
  # SAS stores 0.1 in 3 bytes as 0.0999755859375; the grid would hold 0.1.
  for (len in c("LENGTH INT_DEAD 3;", "LENGTH DEFAULT=4;", "ATTRIB INT_DEAD LENGTH=3;")) {
    job <- translate_494(c(
      fit_macro_grid[1:3],
      paste("DATA PRED;", len, "AGE=50; DO INT_DEAD=0.1,0.2; OUTPUT; END; RUN;"),
      hazpred_macro_grid
    ))
    expect_match(refused_494(job), "numeric length below 8", fixed = TRUE, info = len)
  }
  # A character length, or a numeric 8, changes no value.
  job <- translate_494(c(
    fit_macro_grid[1:3],
    "DATA PRED; LENGTH LBL $ 3 INT_DEAD 8; AGE=50; DO INT_DEAD=0.1,0.2; OUTPUT; END; RUN;",
    hazpred_macro_grid
  ))
  expect_equal(grid_494(job)$INT_DEAD, c(0.1, 0.2))
})

test_that("OPTIONS OBS= set by a macro refuses the grid, wherever it is called", {
  # The option holds for every later step, so a call before the grid's step
  # changes the rows PROC HAZPRED reads as surely as one after it.
  job <- translate_494(c(
    "%MACRO LIM; OPTIONS OBS=2; %MEND LIM;",
    "%LIM;",
    fit_macro_grid,
    hazpred_macro_grid
  ))
  expect_match(refused_494(job), "`%LIM` calls a macro whose body holds", fixed = TRUE)
})

test_that("PROC CONTENTS OUT2= writes the dataset it names", {
  job <- translate_494(c(
    fit_macro_grid,
    "PROC CONTENTS DATA=X OUT2=PRED; RUN;",
    hazpred_macro_grid
  ))
  expect_match(refused_494(job), "PRED is written by PROC CONTENTS", fixed = TRUE)
})

test_that("option and programming statements in a listed procedure leave the grid", {
  # PROC IMPORT's GETNAMES= and PROC PHREG's programming statements are
  # assignments inside a procedure, not a lost DATA step.
  job <- translate_494(c(
    "PROC IMPORT DATAFILE='a.csv' OUT=X DBMS=CSV REPLACE; GETNAMES=YES; RUN;",
    "PROC PHREG DATA=X; MODEL T*D(0)=A2; A2=A*A; RUN;",
    fit_macro_grid,
    hazpred_macro_grid
  ))
  expect_equal(grid_494(job)$AGE, c(50, 50, 50))
  expect_equal(nrow(grid_rows_allow(job)), 0L)
})

test_that("a call before the grid's first step keeps the grid and records a row", {
  # Maintainer decision (2026-10-09). The call cannot write a dataset the
  # grid reads, because the grid's steps build them afresh, but it may set
  # session state they read (OPTIONS OBS=). The grid is emitted, and the row
  # says it is SAS's grid only if the call leaves that state alone.
  for (call in c("%INCLUDE 'setup.sas';", "%SETUP;")) {
    job <- translate_494(c(call, fit_macro_grid, hazpred_macro_grid))
    expect_equal(grid_494(job)$AGE, c(50, 50, 50), info = call)
    rows <- macro_rows(job)
    expect_equal(nrow(rows), 1L, info = call)
    expect_equal(rows$construct, sub(";$", "", toupper(call)), info = call)
    expect_match(rows$reason, "session state", fixed = TRUE, info = call)
    expect_match(rows$reason, "OPTIONS OBS=", fixed = TRUE, info = call)
    expect_match(rows$reason, "DATA=PRED", fixed = TRUE, info = call)
  }
  # Control: no call before the grid, no row.
  job <- translate_494(c(fit_macro_grid, hazpred_macro_grid))
  expect_equal(nrow(job$untranslated), 0L)
})
