# Portable ancillary Python checks

These are revised portable verification programs. The stored exact-decimal
tail shape is unchanged. The tail estimates have been made explicit and repaired,
and the matching parser now widens every serialized core quantity by a decimal
unit. Additional programs verify the horizon derivations and initialization,
symbolic identities, and outward manuscript data.

From any directory, run:

```sh
python3 /path/to/ancillary/python/check_reconstruction.py
```

The default regenerated input is `../logs/core_reconstructed.txt`, relative to
this directory. An explicit different generated input can be selected with
`--core /path/to/core-output.txt`. The checker:

1. Compares the six regenerated center intervals with the supplied historical
   output in `../reference/core_tm_final2_t40.txt` using exact rational arithmetic.
2. Checks that the regenerated initial horizon Hessian bounds, divided by the
   newly chosen initial component weights, are smaller than the core's
   `Binit = 1e-50`. This uses the regenerated
   `../logs/horizon_second_derivative.txt` and accounts conservatively for
   binary128 weight conversion. It is an exact rational check.
3. Runs the independent rational horizon, recurrence-derivation, horizon
   initialization, tail, and constraint Python checks.
4. Invokes the Poincare–Miranda script with the regenerated core path explicitly.

The comparison with the supplied output is diagnostic. The supplied intervals
are not made rigorous by copying them, overlapping them, or reproducing their
decimal values. In particular, the independently reconstructed missing headers
are new implementations and a newly chosen weight schedule, not recovered
original files. The original missing header contents cannot be inferred uniquely
from the archive. Review and regeneration of the C++ core validation remain
necessary before treating its output as validated input to matching.

`--matching-only` skips the separate rational checks after they have already
passed, but still compares centers, checks the initial Hessian weights, and
passes the selected core file explicitly to matching. Using the supplied
reference itself as `--core` is only a smoke
test of the portable scripts, and the checker labels it accordingly.

The numeric certificate scripts use only the Python standard library. The
additional `symbolic/verify_constraints.py` requires SymPy and is run by the
full `make verify` pipeline. The auxiliary
tail-shape generation and recurrence exploration scripts require `mpmath`.
The stored exact-decimal shape can be certified directly without regenerating
it, because the tail certificate checks the residual of that supplied finite
approximation independently.
