# Environment and reproduction

The proof-critical checker, mutation suite and short-chain check use Python 3
standard-library modules only:

```sh
python3 -I verify.py --report /tmp/wmy5-verification.json
python3 -I test_mutations.py
python3 -I check_chain_short_identity.py
```

The monomial inventory and uniform Hessian-minor audit use SymPy (revision
execution used SymPy 1.14.0):

```sh
python3 discovery/grading_audit.py
python3 discovery/generalized_hessian_check.py
```

The checks were also run with `python3 -O -I -B` where applicable, so the main
verifier and all three added audits use explicit failures rather than removable
Python assertions, so success does not depend on active assertions. LaTeX uses the
standard article class, AMS packages and pdfLaTeX. `supplementary/` contains a
Wolfram Language source and recorded exploratory results; those calculations are
not part of the standalone certificate obligations.
