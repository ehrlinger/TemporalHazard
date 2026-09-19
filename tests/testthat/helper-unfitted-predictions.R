# predict() on a model built with fit = FALSE warns that its numbers come
# from starting values rather than estimates (#144). Whole files here
# predict from such a model on purpose -- it is an intended, tested
# capability -- and one warning per call would bury the warnings that mean
# something. Switch it off for the suite by default; the tests that assert
# the warning switch it back on for their own scope, so the assertion still
# has something to catch.
options(TemporalHazard.warn_unfitted_prediction = FALSE)
