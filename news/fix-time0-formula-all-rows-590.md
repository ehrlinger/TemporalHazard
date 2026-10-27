* **A formula helper that reads an outside vector by position is now refused
  when rows at time 0 are dropped (#590).** `hazard()` drops a row at time 0
  before fitting, and builds the design on the rows that remain. A helper
  such as `function(a) a + g[seq_along(a)]`, with `g` aligned to every row
  given, then paired the remaining rows with the first values of `g`, with
  no warning. In one example all 48 of 48 rows were paired with the wrong
  value; in another, an alternating treatment arm read this way fitted its
  coefficient with the sign reversed. The existing check reads only the
  variables the formula names, so it could not see `g` inside the function.
  `hazard()` now also evaluates the formula's variables on every row in
  order and again with the remaining rows first, as the fit sees them, and
  refuses a variable that then gives a remaining row a value other than its
  own, naming it. Terms computed from the rows as a set, such as `scale()`,
  `poly()`, `ns()` and factors, pass and still use the remaining rows. A
  term that depends on the order of the rows, such as `cumsum(age)`, or is
  random, such as `jitter(age)`, is now refused when rows other than the
  last are dropped, where it was fitted before. A helper that depends on
  how many rows it is given, such as `tail(g, length(a))` or
  `a * mean(g[seq_along(a)])`, is not caught. Put such a vector in `data`,
  where it is dropped with its rows.
