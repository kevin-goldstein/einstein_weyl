# Reproducible core-validation repairs

This directory contains a reviewable patch, a generator, and focused numerical
checks. These are newly written repairs to the supplied source, not recovered
author files. The archival source in `reference/` is preserved. The production
source is now distributed in readable form; the complete calculation is checked
by `make verify`, while the focused tests here exercise its arithmetic helpers.

## Generating the production source

From the ancillary directory:

```sh
make regenerate-core
```

The generator defaults to the archival core source. An alternative input file
can be passed as the positional argument, including one whose header paths have
already been adapted. The generator makes the two include paths portable and
removes the obsolete development marker from the active copy. Each expected
replacement must occur exactly once; the
generator fails instead of silently applying a partial patch. It refuses to
overwrite its input. `core_validation.patch` records the algorithm changes;
`core_readability.patch` separately records the formatting, identifier changes,
comments and removal of four unused helpers. Regeneration uses clang-format;
ordinary `make verify` needs no formatter and never rewrites the source.
See `READABILITY.md` for the source map and token-consistency check.

## Mathematical reasons for the changes

1. **Sensitivity tube uses the actual Jacobian enclosure.** In
   `validated_parameter_substep` (formerly `qstep_centeredQ_TM`), `M` encloses the Taylor coefficients of the true
   center-orbit Jacobian and `Mp` contains selected scalar midpoints.
   Consequently, `si=qsenscoef(S0,M)` encloses coefficients of the actual
   variational solution initialized at the stored column center, while `sp`
   describes a distinct polynomial obtained using `Mp`. The polynomial `Pr`
   entering the sensitivity tube must use `si`. A Taylor remainder generated
   from the true equation does not account for the difference between the
   lower-order coefficients of `si` and `sp`. The patch retains `sp` for choosing
   the endpoint center, where the separately computed radius accounts for that
   choice. It changes only the polynomial used to establish the enclosure.

2. **Radius is measured about the represented center.** For an interval
   `[lo,hi]` and a stored binary128 center `c`, an enclosing radius is
   `max(up(c-lo),up(hi-c))`. Half the interval width bounds distances from the
   mathematical midpoint, which need not equal the rounded stored midpoint.
   The difference matters when an interval spans very few representable numbers.
   The new `qradius_about` rounds each subtraction upward. It is used in the
   initial sensitivity errors, local endpoint errors, and inflation of the
   sensitivity tube. The original midpoint choices remain unchanged.

3. **Tube inflation encloses its input hull.** Computing a rounded midpoint,
   then inflating half the width by 1.02, does not universally enclose a narrow
   interval about that midpoint. The inflation now starts with the endpoint
   distance radius. The multiplicative and additive inflation constants, and
   both subsequent arithmetic operations, are rounded upward. Forming the final
   lower and upper endpoints retains the existing outward rounding.

4. **Negative reciprocal endpoints have the same order as positive ones.** On
   either connected component of the nonzero real numbers, `x -> 1/x` is
   decreasing. The reciprocal of `[lo,hi]` is therefore
   `[down(1/hi),up(1/lo)]` for both signs. The zero-containing denominator check
   remains in place. The removed negative branch reversed the endpoints and
   applied the wrong rounding directions before multiplication. This is a
   generic interval-arithmetic repair; the intended physical core denominators
   are positive.

5. **Published positivity minima are lower endpoints.** `minA`, `minD`, and
   `minU` are parsed downstream as certified lower bounds. The explicit helper
   `mdownstr` uses MPFR's `%.*RDg` format, rounding decimal output toward minus
   infinity. It avoids depending on the replacement header's default scalar
   formatting policy. This change does not affect the internal integration.

6. **The sensitivity tube requires strict inclusion.** The predicate
   `strict_quad_subset(need, St[i][a])` requires both
   `need.lo > St[i][a].lo` and `need.hi < St[i][a].hi`. Equality at either
   endpoint is now a failed inclusion and causes further tube inflation or
   failure of the step. This implements the strict inequality used in Appendix C,
   equation label `eq:core-sensitivity-tube`. The focused tests call the actual
   production predicate and cover equality, escape and strictly interior cases.

## Focused verification

`check_core_repairs.cpp` includes the production source with its entry point
renamed, so every tested helper is the one used by the certificate. Build and run
it with:

```sh
make build/check_core_repairs
build/check_core_repairs
```

The checks demonstrate an actual failure of the old half-width formula for two
adjacent binary128 endpoints; verify containment with the repaired radius; check
negative reciprocal enclosures against 512-bit MPFR bounds for 99 intervals; and
verify directed decimal formatting for positive and negative values; and check
eight sensitivity-tube boundary cases. Full sensitivity propagation, the weight
schedule and the final Poincare--Miranda margins are checked by `make verify`.
