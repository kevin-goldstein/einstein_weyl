# Verification of the revised certificate chain

The complete `make verify` pipeline passed on 14 September 2026 using macOS
arm64, GNU GCC 15.2.0, MPFR 4.2.2, GMP 6.3.0, GNU Make 3.81, Python 3.14.7,
and SymPy 1.14.0. The final verification transcript is `logs/full_make_verify.txt`.
The production programs were rebuilt from a clean build directory before this
run. The strict sensitivity-tube change and readable production source passed;
the core output and exact `paper_data.json` are byte-identical to the first
revision. The numerical enclosures are unchanged; timing and test/provenance
logs record the new run.

| Check | Result |
|---|---|
| MPFR interval exact-rational tests | Pass: 1400 interval pairs and edge cases |
| Binary128 numeric bridge | Pass: 10000 round trips, 3000 signed directed cases, boundary cases |
| Reproducible readable source | Pass: executable tokens agree with the repaired reference and fixed readability transformations |
| Focused production core helpers | Pass, including 8 strict tube-boundary cases |
| Exact decimal source transfer | Pass: line wrapping, missing and ambiguous initializer tests |
| Rebuilt horizon positivity | Pass |
| Rebuilt horizon derivatives | Pass |
| New initial weights versus `Binit=1e-50` | Pass; largest weighted Hessian bound below `2.327e-52` |
| Horizon majorant and tail-jet bounds | Pass |
| Direct and triangular horizon recurrences | Agree exactly in Q[b] through order 12 |
| All-index majorant constants / geometric-tail ratios | Pass |
| Horizon-to-core initialization amplification | Pass; `1e-70` dual padding suffices |
| Repaired compact core | Reaches t=40; 485 coarse steps and 1590 sensitivity substeps |
| Sensitivity bootstrap ratio | Approximately 0.00400696257435, below 1 |
| State bootstrap ratio | Approximately 0.00238332468681, below 1 |
| Revised centered tail | Pass; same certified ball radius and weights |
| Fixed-v scalar constraint | Pass |
| Exact-rational matching with decimal-transfer allowances | All five margins positive |
| Five exact symbolic constraint identities | Pass |
| Outward theorem, mass, temperature, literal-action entropy and face table | Generated and checked with exact rationals |

The following regenerated face-margin fractions are rounded downward:

```
b        0.540619770207
R        0.574495623641
A_inf    0.597534509130
beta     0.547484807780
eta      0.647565186013
```

The matching parser explicitly widens every serialized core endpoint, center,
error and remainder by a unit in its last displayed decimal place before
exact-rational post-processing. The usual one-percent error inflation is then
applied as documented. The minimum unrounded fraction exceeds
`0.5406197702071599`.

The regenerated center intervals have the same rational midpoints as the
supplied historical intervals and small outward differences in final decimal
places. This comparison is diagnostic. The historical intervals are not used
as proof inputs to the new match. The complete regeneration, including the new
positive schedule and repairs, is a new calculation rather than authentication
of the earlier missing-file archive.

The symbolic log is `logs/symbolic_constraints.txt`; the exact rational physical
and matching data are in `logs/paper_data.json`; the main Python transcript is
`logs/verification.txt`. The compact-core run with strict sensitivity-tube
inclusion took approximately 51.11 seconds on the tested host, as recorded in
`logs/core_reconstructed.log`.

These computational checks accompany the analytic arguments in the revised
paper. They do not establish local uniqueness of the global match, stability,
a complete cold branch, or regularity of the black-hole interior.
