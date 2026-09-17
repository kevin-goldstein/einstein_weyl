# Reconstructed interval compatibility header

`mpfr_interval.hpp` was independently implemented from the interfaces visible in
the supplied Einstein–Weyl ancillary C++ callers on 2026-09-14. It is a compatible
replacement, not a recovery of the missing original source, and passing tests
does not authenticate the original archived output.

The implementation supports finite closed intervals. It throws for reversed
endpoints, invalid/nonfinite inputs, nonfinite results, and division by an
interval containing zero. Each interval arithmetic endpoint uses MPFR's directed
rounding. Multiplication evaluates all four endpoint products in both directions;
division reverses the reciprocal endpoints for either sign of a nonzero divisor.

`Real` stores a single MPFR dyadic scalar. Copy/move operations preserve its value
and precision. Ordinary scalar arithmetic uses the current `mpi::PREC`, which
defaults to 512 bits. Decimal `Real` construction rounds nearest unless supplied
an explicit rounding argument. Decimal `I` construction always encloses the exact
decimal using RNDD/RNDU. A caller requiring an exact-decimal bound must use the
latter or pass a direction explicitly.

`width(I)` rounds upward. `mid(I)` computes the exact dyadic midpoint using enough
temporary precision; this prevents the midpoint/half-width mismatch in the
callers' MPFR recentering routine without increasing the configured arithmetic
precision. The separate binary128 centering code is outside this header.

`str(I,n)` rounds the lower decimal downward and upper decimal upward using MPFR.
`str(Real,n,rnd)` accepts a rounding mode and defaults to nearest for diagnostic
scalars. Lower-bound outputs need `MPFR_RNDD`; upper-bound outputs need
`MPFR_RNDU`. The header cannot infer which role an arbitrary scalar plays.

Validation uses `tests/test_mpfr_interval.cpp`, with GMP exact rational endpoint
arithmetic as an independent oracle. It checks 1,400 randomly selected interval
pairs at 24, 53, 113, and 512 bits, all basic interval operations, negative
division, zero-divisor rejection, aliasing, copy/move across precision changes,
exact midpoints including widely separated exponents, malformed input rejection,
and decimal enclosure at 1, 2, 7, 30, and 50 significant digits. The test passed
with GCC 15, MPFR 4.2.2, and GMP 6.3.0 on this machine.

Example build from `work/reconstruction`:

```
/opt/local/bin/g++-mp-15 -std=c++17 -O2 -Wall -Wextra \
  -Iinclude -I/opt/local/include tests/test_mpfr_interval.cpp \
  -L/opt/local/lib -lmpfr -lgmpxx -lgmp -o tests/test_mpfr_interval
tests/test_mpfr_interval
```

No original manuscript or ancillary source was edited by this reconstruction.
