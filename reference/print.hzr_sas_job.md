# Print method for hzr_sas_job

Shows what
[`hzr_translate_sas()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_translate_sas.md)
made of one SAS job: the chunks it emitted, how many of the job's tokens
it mapped, and one line per construct it could not translate. Read the
untranslated lines before trusting the document – each one names
something the emitted R does differently from PROC HAZARD, or not at
all.

## Usage

``` r
# S3 method for class 'hzr_sas_job'
print(x, ...)
```

## Arguments

- x:

  An `hzr_sas_job` object, as returned by
  [`hzr_translate_sas()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_translate_sas.md).

- ...:

  Additional arguments (ignored).

## Value

The object `x`, invisibly.

## See also

[`hzr_translate_sas()`](https://ehrlinger.github.io/TemporalHazard/reference/hzr_translate_sas.md),
which builds the object.
