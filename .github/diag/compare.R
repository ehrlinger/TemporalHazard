out <- Sys.getenv("DIAG_OUT")
s <- read.csv(file.path(out, "serial.csv")); p <- read.csv(file.path(out, "parallel.csv"))
agg <- function(d) {
  d$key <- paste(d$file, d$test, sep = " :: ")
  d$n <- 1
  aggregate(cbind(passed, failed, skipped, warning, n) ~ key + file, data = d, FUN = sum)
}
m <- merge(agg(s), agg(p), by = c("key", "file"), all = TRUE, suffixes = c(".s", ".p"))
m[is.na(m)] <- 0
d <- m[m$passed.s != m$passed.p | m$failed.s != m$failed.p | m$skipped.s != m$skipped.p | m$n.s != m$n.p, ]
cat("rows serial", nrow(s), "parallel", nrow(p), "\n")
cat("differing (file, test) keys:", nrow(d), "\n")
byfile <- aggregate(cbind(dpass = passed.p - passed.s, dfail = failed.p - failed.s, dn = n.p - n.s) ~ file, data = d, FUN = sum)
print(byfile[order(-abs(byfile$dpass)), ], row.names = FALSE)
print(utils::head(d[order(-abs(d$passed.p - d$passed.s)), ], 60), row.names = FALSE)
write.csv(d, file.path(out, "diff.csv"), row.names = FALSE)
