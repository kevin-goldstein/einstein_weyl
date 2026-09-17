# Replacement weight schedule

`../include/fixed_weight_schedule.hpp` is a new replacement for a missing
dependency. It is **not** a recovered copy of the original author's schedule.

Run `python3 weights/generate_fixed_weights.py` from the reconstruction folder
to regenerate it. The generator uses only the Python standard library. It
reads the supplied central endpoint intervals, takes approximate midpoints,
integrates the center equations backwards with ordinary RK4, and chooses the
positive Perron right eigenvectors of the pointwise Metzler comparison matrix.
The matrix keeps the actual Jacobian diagonal and replaces its off-diagonal
entries with their absolute values. The last component is normalized to one.
The schedule has 1591 nodes on the caller's 1590 substeps, between 20/19 and 40.

The approximate integration and eigenvectors are optimization choices, not
validated data. The proof may use *any* positive weights: each listed decimal
string defines an exact positive rational number. For an interval Jacobian M,
the caller calculates the bound

    mu = max_i (upper(M_ii) + sum_{j != i} mag(M_ij) q_j/q_i)

with outward rounding. At a change of weights q to q_next it multiplies error
bounds by `max_i q_i/q_next_i`, also outward rounded. These independent checks,
not closeness to a Perron eigenvector or accuracy of the approximate orbit,
are what make a schedule admissible.

`weight_generation_diagnostics.json` records approximate orbit samples and
eigenvector residuals for inspection. It must not be used as a source of
rigorous state bounds. The core certificate and matching certificate must be
rerun with this replacement before claiming successful reconstruction.
