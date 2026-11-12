# OptSampling 0.1.0

- Package the supplied VNS and COS/ACOS R implementations, preserving their
  search engines, function names and fitted-object structures.
- Move the embedded C++ code to compiled package sources; register native
  routines through Rcpp. Loading the package no longer compiles code or writes
  functions into the global environment.
- Add vns_multistart() for supplied/random starts and per-start continuation
  budgets, retaining every attempt and choosing the lowest final design loss.
- Add systematic_support() and draw_systematic_sample() public helpers.
- Add a Spanish tutorial with synthetic examples, R help and invariant tests.
