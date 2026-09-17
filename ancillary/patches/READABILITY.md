# Readable production core and source provenance

`src/verified_core_taylor_model.repaired.cpp` is the source compiled by the
certificate and the focused helper tests. The archived compressed source remains
in `reference/verified_core_taylor_model.cpp`. The transformation has two stages:

1. `patch_core_validation.py:patched` applies the numerical repairs listed in
   `README.md`. `core_validation.patch` records these changes before formatting.
   The second-round numerical change is strict sensitivity-tube containment:
   equality at either boundary fails the inclusion test.
2. `prepare_readable_core.py:readable` removes four unused helpers, renames
   identifiers by a fixed table, and inserts comments. It does not rewrite
   arithmetic expressions or reorder stages. The implicit successful return
   from `main` is made explicit as `return 0`, allowing helper tests to rename
   the entry point without a return-type warning. clang-format lays out the result;
   `core_readability.patch` records this separate stage.

The removed helpers are `qtruncate`, `qbinom`, `qshiftstate` and the interval
`qs(QI, int)` overload. The first and third were unused Taylor truncation/shift
experiments; `qbinom` was used only by the discarded shift helper. The scalar
`qs(__float128, int)` serializer remains in production.

## Computational source map

The source comments cite stable LaTeX labels in Appendix C so that equation
renumbering does not invalidate the references.

| Production function or type | Previous name | Purpose / proof label |
|---|---|---|
| `core_rhs` | `rhs` | Six-state ODE; main `eq:F`, `eq:first-order-rest` |
| `MpfrStateTaylor`, `mpfr_state_taylor_coefficients` | `MTS`, `mcoef` | State Taylor recurrence; `eq:core-order56-step` |
| `validated_center_step` | `mstep` | Strict Picard inclusion and state remainder; `eq:core-picard-inclusion` |
| `quad_jacobian_taylor_coefficients` | `qjacseries` | Jacobian series; `eq:core-first-variation` |
| `quad_sensitivity_taylor_coefficients` | `qsenscoef` | Variational recurrence; `eq:core-sensitivity-recurrence` |
| `validated_parameter_substep` | `qstep_centeredQ_TM` | Sensitivity tube, endpoint error, norm transfer and parameter bootstrap |
| `strict_quad_subset` | newly named predicate | Strict tube test; `eq:core-sensitivity-tube` |
| `StateHessianDual`, `interval_rhs_hessian` | `QAD2`, `qrhs_hessian` | Interval Hessian source; `eq:core-hessian-source` |
| `load_outward_weights`, `weight_index` | `qweights`, `qindex` | Directed conversion and schedule selection; `eq:core-error-update` |
| `initialize_horizon_first_variations` | `mhorizon` | Scaled parameter columns; `eq:core-initial-map`, `eq:initial-dual-padding` |
| `initialize_horizon_center` | `center_horizon` | Center initialization and analytic horizon tail |
| `main` | unchanged | Coarse/substep schedule and serialization; `eq:core-coarse-schedule` |

The parameter substep has comments at its existing phase boundaries: midpoint
predictor, strict interval inclusion, Hessian source, logarithmic norm, weight
transition and error update. The midpoint predictor still chooses only the stored
endpoint center; the inclusion polynomial uses the interval Jacobian series.

The paper-data generator reads the exact `b0` and `R0` decimal initializers with
a whitespace-tolerant parser, so line wrapping does not change this transfer.
Its regression tests cover compressed and multiline forms, and rejection of
missing or ambiguous initializers.

## Checking and regenerating

Normal verification never changes the distributed source:

```sh
make check-core-source
make verify
```

The consistency check repeats the deterministic repairs and named readability
transformations, then compares C++ lexical tokens with the distributed source.
Comments and whitespace are ignored; identifiers, literals and operators are
retained. A difference causes failure. This check requires Python but no
formatter. It also runs before compiling the production core and its focused
helper tests. The tests call the production radius, reciprocal, directed decimal
formatter and strict containment helpers directly.

Only source maintenance requires a formatter:

```sh
make regenerate-core
# Or select a formatter explicitly:
CLANG_FORMAT=/path/to/clang-format make regenerate-core
```

The recorded formatting run used Apple clang-format 16.0.0
(`clang-1600.0.26.6`) and the root `.clang-format` configuration. The generator
compares tokens before and after formatting and refuses output if they change.
Long exact numeric strings are deliberately not split. Formatter versions may
change whitespace without changing the checked program. After any maintenance
edit, rerun `make verify` and regenerate the release checksum manifest.
