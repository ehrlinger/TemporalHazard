# bh.dead SAS parity harness (development only)

`inst/dev/` is `.Rbuildignore`d: nothing here ships to CRAN.

**This repository is public.** Neither the row-level data nor the study's
aggregate results live here. The fixture is built from SAS listings on a
secure volume and written back to that volume. Committed code contains no
study values; the parity test derives every expectation from the fixture at
run time and skips when the volume is not mounted.

## Paths

| Environment variable | Purpose |
|----------------------|---------|
| `TEMPORALHAZARD_BHBLT` | Row-level SAS dataset (`bhblt.sas7bdat`) |
| `TEMPORALHAZARD_BHDEAD_FIXTURE` | Built fixture (`bhdead.rds`) |

Both are read from the environment only -- there is no default path baked
into this repository. Set them to point at the secure-volume locations
before rebuilding or loading the fixture.

## Rebuilding the fixture

With the volume mounted, pass the two listing paths explicitly:

```r
devtools::load_all()
source("inst/dev/bhdead-parity/parse-bhdead-lst.R")
.hzr_write_bhdead_fixture(
  hz_lst = "<path to the shape-fit listing, hz.dead.lst>",
  bh_lst = "<path to the bootstrap listing, bh.dead.lst>"
)
```

`hz_lst` and `bh_lst` are required arguments -- this script does not store or
default to a study directory. The output path defaults to
`Sys.getenv("TEMPORALHAZARD_BHDEAD_FIXTURE")`; if that variable is unset,
`.hzr_write_bhdead_fixture()` stops and asks you to set it.

This reads the shape-fit listing and the bootstrap listing, and writes
`bhdead.rds` to the configured output path.

## Running the parity test

```r
devtools::load_all()
testthat::test_file("tests/testthat/test-bhdead-sas-parity.R")
```

It skips unless the data and the fixture are both readable and `haven` is
installed.

## What the SAS side does

* Shape fit: `PROC HAZARD ... noconserve`, fixing only `M`.
* Bootstrap: `%hazboot(seed=-1, resampl=1000, sle=0.12, sls=0.1)`, base model
  fixing `nu` and `m`.

`seed=-1` means SAS seeds from the time of day, so its selection frequencies
are one random realisation and cannot be reproduced exactly. Bootstrap
assertions are therefore statistical, not exact.
