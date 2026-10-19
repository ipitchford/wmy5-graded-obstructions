# AI research index

## Entry and identity

Start with [README.md](README.md), then read [CLAIMS.json](CLAIMS.json),
[research_note.tex](research_note.tex), and [ASSURANCE.md](ASSURANCE.md). This is
version 1.0.0, creator Anonymous, status unrefereed candidate. Publication
metadata is in [CITATION.cff](CITATION.cff). Do not infer that the global
five-dimensional weak Markus--Yamabe problem has been solved.

## Claim-to-evidence map

- `WMY5-BASELINE`: Section 2; `verify.py`; exact six-dimensional characteristic
  identity and three exact zeros.
- `WMY5-NAIVE`: Proposition 2.1; `verify.py`; exact positive-eigenvalue witness.
- `WMY5-COVERAGE`: Proposition 3.1; `discovery/grading_audit.py`; exhaustive
  monomial inventories only for the two displayed bigradings.
- `WMY5-PARALLEL`: Theorem 4.1; `verify.py` and the cross certificates; exact
  exclusion of the complete displayed 17-parameter family.
- `WMY5-CHAIN`: Theorem 5.1; `verify.py`, the chain certificate, and
  `check_chain_short_identity.py`; exact exclusion of the complete displayed
  16-parameter family.
- `WMY5-HESSIAN`: Theorem 6.1 and Corollary 6.2;
  `discovery/generalized_hessian_check.py`; written triangular proof plus an
  exact symbolic nonzero minor for the specified uniform potential family.

## Reproduce and interpret

The standard-library checks are `python3 -I verify.py --report
/tmp/wmy5-verification.json`, `python3 -I test_mutations.py`, and `python3 -I
check_chain_short_identity.py`. The monomial and uniform-Hessian audits require
SymPy 1.14 or a compatible release. See [ENVIRONMENT.md](ENVIRONMENT.md).
The external report path keeps replay from mutating the manifest-bound package.

Exact replay establishes the stated polynomial identities and finite certificate
relations. It is not proof-assistant formalisation, independent human peer
review, or exhaustive priority certification. The supplementary Wolfram
Language searches are exploratory because small independent certificates are
not shipped for those larger systems.

## Reuse constraints and open directions

Preserve the fixed base terms, gradings, auxiliary architecture and nonzero
parameter hypotheses when reusing the family exclusions. Preserve the
affine-in-auxiliary restriction and affine source--target boundary when reusing
the Hessian theorem. The global five-dimensional construction/nonexistence
problem, nonlinear changes of variables and unrelated architectures remain open.

Read [SOURCES.md](SOURCES.md), [PROVENANCE.md](PROVENANCE.md), and
[LICENSES.md](LICENSES.md) before asserting novelty, authorship or reuse rights.
