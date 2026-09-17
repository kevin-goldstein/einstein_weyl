# Computational certificate for an Einstein--Weyl black hole

This repository contains the complete computational artifact supporting the existence proof in *A computer assisted existence proof for a non-Schwarzschild black hole in Einstein--Weyl gravity* by Kevin Goldstein and Vishnu Jejjala.

The manuscript is deliberately not included.
The repository contains the production sources, exact input data, interval arithmetic support, fixed weight schedule, tests, verification scripts, successful logs, and generated numerical fragments used by the paper.

## Reproduce the certificate

The required software and the verified platform are documented in [`ancillary/README.md`](ancillary/README.md).
In brief, the calculation requires GNU GCC with IEEE binary128 and libquadmath, MPFR, GMP, GNU Make, Python 3.9 or newer, and SymPy.

From the repository root, run:

```sh
cd ancillary
make verify
```

The command rebuilds the numerical programs and reruns the arithmetic tests, horizon and core enclosures, tail certificate, constraint identities, exact matching calculation, and generated numerical data.
It stops if any required check fails.
The final verified results and software versions are summarized in [`ancillary/VERIFICATION.md`](ancillary/VERIFICATION.md).

## Repository contents

- `ancillary/src/` contains the production C++ sources.
- `ancillary/include/` contains the interval arithmetic header and fixed weight schedule.
- `ancillary/python/` contains the exact and validated Python checks.
- `ancillary/tests/` contains focused regression tests.
- `ancillary/patches/` records the deterministic source repair and readability transformations.
- `ancillary/reference/` contains the archived source used by the source consistency check and historical diagnostic output.
- `ancillary/logs/` contains the successful verification output.
- `fragments/` contains the generated numerical macros and matching table used in the paper.

The active production core enforces strict containment of the sensitivity tube.
The complete verification passes with 485 coarse steps and 1,590 sensitivity substeps, and all five Poincaré--Miranda face margins are positive.

## Citation and public URL

Author metadata is supplied in [`CITATION.cff`](CITATION.cff).

See [`UPLOAD_INSTRUCTIONS.md`](UPLOAD_INSTRUCTIONS.md) for GitHub and Zenodo publication routes.

## License

No license has been selected in this bundle because the choice belongs to the authors.
Add the chosen license before making the repository public; without a license, ordinary copyright restrictions apply even though the files are visible.
# einstein_weyl
