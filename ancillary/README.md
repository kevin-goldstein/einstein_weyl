# Einstein–Weyl existence certificate: revised sources and verification

This directory contains the complete numerical and symbolic support for the
revised manuscript. The interval header and positive weight schedule are newly
constructed implementations/data. They are not recovered or authenticated copies
of the two files absent from the earlier archive. Their provenance is recorded
in `include/MPFR_INTERVAL_PROVENANCE.md` and `weights/README.md`.

## Reproduce the complete calculation

Requirements: GNU GCC with IEEE binary128 and libquadmath, MPFR, GMP, GNU Make,
Python 3.9 or newer, and SymPy. Optional construction of the tail approximation
uses `mpmath`; it is not needed for the certificate checks.

From this directory run:

```sh
make verify
```

For installations that require explicit dependency locations:

```sh
make verify CXX=g++ MPFR_PREFIX=/usr PYTHON=python3
```

The verified environment is macOS arm64, GNU GCC 15.2.0, MPFR 4.2.2, GMP 6.3.0,
GNU Make 3.81, Python 3.14.7, and SymPy 1.14.0. Compilation uses GNU C++17,
`-O2 -fno-fast-math -ffp-contract=off`. Do not enable unsafe floating-point
optimizations. The supplied binary128/MPFR numeric bridge handles the conversion
entry points absent from the tested MPFR build. Its assumptions and tests are
included in the source.

`make verify` stops on failure and performs the following operations:

1. Check that the distributed readable core agrees with the reproducible repair
   and readability transformations, then build and run interval arithmetic,
   binary128 conversion and focused tests of the production core helpers.
2. Rebuild and run horizon positivity, horizon derivatives and the compact core.
3. Verify the initial weighted Hessian bound against the actual chosen schedule;
   run the horizon majorant, tail-jet, recurrence-derivation and initialization
   bounds; then run the centered tail and scalar-constraint checks.
4. Run exact-rational Poincare–Miranda matching on the newly generated core file.
5. Verify the constraint identities in SymPy and generate outward-rounded
   manuscript data for the parameter box, mass, temperature, literal-action
   Wald entropy and matching table.

The logs are written under `logs/`. The source-to-Python decimal transfer adds
one unit in the last printed place to every serialized core bound and center,
with the appropriate outward direction. This includes nearest formatting of
binary128 column centers, independently of the producer's formatting mode.
The table generator uses exact rational arithmetic and rational bounds for pi.

The paper fragments are generated in `../fragments/`, next to this ancillary
directory. To rerun just the exact identities or regenerate manuscript numbers
from an existing regenerated core output, use `make symbolic` or `make paper-data`.

## Production sources and supporting data

- `src/verified_core_taylor_model.repaired.cpp` is the production core source.
  Verification checks its executable tokens against the historical reference
  and the deterministic transformations in `patches/`, without rewriting the
  source or requiring a formatter. The maintenance target `make regenerate-core`
  rebuilds it using clang-format and the supplied `.clang-format` configuration.
- `include/mpfr_interval.hpp` supplies directed MPFR interval arithmetic and
  outward decimal endpoint formatting.
- `include/fixed_weight_schedule.hpp` contains 1591 positive rational weight
  vectors. Every logarithmic-norm and norm-change bound is recomputed by the core.
- `python/tail/tail_shape_data.py` is an explicit finite rational approximation;
  the centered-tail certificate checks its full residual independently.
- `python/symbolic/verify_constraints.py` verifies constraint propagation, its
  coordinate identification, the horizon constraint and the fixed-v derivative.
- `python/matching/generate_paper_data.py` generates the displayed outward bounds.
- `reference/verified_core_taylor_model.cpp` is the archived source read by the
  reproducible source-consistency check. The historical core output files in
  `reference/` are retained for diagnostic comparison; those old intervals are
  not inputs to the new final matching calculation.

The production core repairs interval sensitivity-tube polynomials, errors about
actual stored centers, tube inflation, negative reciprocal ordering and downward
printing of positivity minima. It also rejects equality at either boundary of
the sensitivity tube. Descriptive stage names and comments link the readable
production source to Appendix C. `patches/READABILITY.md` records the source map;
the algorithm and readability diffs are separate. Four unused development
helpers and the obsolete marker are omitted from the active generated source.

## Construction tools versus proof checks

`make weights` selects a new schedule using approximate backward integration and
positive comparison-matrix eigenvectors. Any positive weights are admissible
choices, but only a subsequent successful `make verify` establishes useful
bounds for them. Accuracy of the approximate trajectory is not a proof premise.
The recorded initial Hessian check depends on the selected initial weights.

The `mpmath` tail-shape generator and corrected-recurrence exploration scripts
construct candidate data. They do not supply rigorous bounds. The exact rational
centered-tail certificate checks the stored shape without trusting those
exploratory calculations. The numeric Python verification scripts use only the
standard library; the separate algebra checks use SymPy.

## Interpretation

A successful run establishes the program inequalities and identities for the
specified inputs, as part of the analytic argument in the revised manuscript.
It is not a bit-for-bit reproduction or authentication of the earlier archive.
See `VERIFICATION.md` and the current logs for measured results. The distributed
source, compiled paper and checksum manifest must be regenerated together after
any further edit.
