* **`predict(newdata = )` no longer blames a term for a failure that was not
  its doing (#446).** When rebuilding the design failed, the refusal named
  whichever model term did not give one value per row of `newdata` and raised
  that *instead of* the failure. The naming was never verified to explain
  anything, so it could blame an innocent term and prescribe a remedy that
  could not be followed: a term returning a list failed on its type, and the
  message sent you to move a `zz` into `data` that was not a column of `data`
  at all. A condition raised inside one of your own terms was destroyed
  outright, class and all.

  What happens now depends on who raised the failure. **A condition carrying
  any class of your own passes through unchanged**: the same object, with its
  class, call and every field, so a `tryCatch()` on it fires and its fields
  are intact. *It gains no note*: nothing is appended to your condition,
  because that cannot be done reliably (`conditionMessage()` is a generic,
  and a class with its own method never reads the field that would be
  written). A plain error from the frame build keeps the term naming,
  **first**, with the original text quoted after it, so a message blaming a
  column you supplied correctly is no longer the first thing you read. And a
  term whose value is not a legal model-frame column at all (a list, say) is
  no longer diagnosed as a row-count problem; you get the type error by
  itself.

  One case is unchanged: a failure raised by a `model.frame()` call *inside*
  one of your terms, when that failure is a plain base error, is
  indistinguishable from our own frame assembly, so a row-mismatched term is
  still named alongside it. Your own text is quoted rather than discarded.
