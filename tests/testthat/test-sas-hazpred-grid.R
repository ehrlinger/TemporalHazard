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
  got <- eval(job$calls$pred, sim$env)$fit
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
