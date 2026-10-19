# Five-variable weak Markus–Yamabe: graded compression and affine-Hessian obstructions

**Status: partial. No five-dimensional counterexample was constructed.**

The requested target was a real polynomial field X:R^5 -> R^5 with
`det(lambda I-DX)=(lambda+1)^5` identically and two distinct real zeros.
This package proves obstructions to two complete bigraded compression families
and to invertible affine normalisation of a uniform Hessian-gradient family.
It does not prove nonexistence of a five-dimensional counterexample.

## Replay

Python 3, standard library only:

```bash
python3 -I verify.py --report /tmp/wmy5-verification.json
python3 -I test_mutations.py
python3 -I check_chain_short_identity.py
python3 discovery/grading_audit.py
python3 discovery/generalized_hessian_check.py
```

The expected final messages are `ALL EXACT CHECKS PASSED` and
`ALL CORRUPTIONS REJECTED`.

The report path is deliberately outside the extracted package. Replaying the
checks must not overwrite the manifest-bound historical `verification_report.json`.

`verify.py` independently constructs and differentiates the polynomial maps,
expands their characteristic polynomials, checks all certificate inputs against
the resulting coefficients, and replays rational polynomial identities. It does
not rerun a Groebner search. It also reproduces the six-dimensional baseline,
checks all three stated zeros, and rejects the tempting five-dimensional map by
an exact positive-eigenvalue witness.

The original Hessian trace certificate uses a 15-by-15 integer minor reconstructed
from the nonconstant coefficients of the explicit potential's Hessian. The uniform
extension checks a symbolic minor with determinant
`-1304596316160*kappa^15*r^10*t^6`, nonzero when `kappa*r*t != 0`.

## Files

- `research_note.pdf`, `research_note.tex`: mathematical statements, assumptions,
  proofs and scope.
- `models.py`, `poly.py`, `verify.py`: standalone exact checker and explicit maps.
- `check_chain_short_identity.py`: standard-library replay in which eight
  affine-linear coefficient equations force the variables and the ninth has
  exact nonzero residual 16,367,616.
- `certificates/`: rational identity data and the Hessian coefficient-minor data.
- `verification_report.json`, `mutation_report.json`, `regeneration_report.json`: executed checks, including separate-copy certificate regeneration.
- `STATUS.json`: machine-readable result and exclusions.
- `discovery/`: optional SymPy certificate regeneration, exact monomial-inventory
  audit and symbolic Hessian-minor check. These require SymPy, unlike the
  standard-library verifier, and remain separate from it.
- `supplementary/`: broader Wolfram Language searches and their reported status.
  They are not verified by `verify.py`.
- `sources.json`: inspected primary sources, access date and GitHub blob hash.
- `MANIFEST.sha256`: SHA-256 hashes of the release files, excluding this manifest.

## Boundaries

No proof assistant was used. The written reductions and verifier have not had
external mathematical review. The literature search was targeted, not an
exhaustive priority audit. The new results concern the precise families stated in
the research note; neither polynomial degree alone nor ambient dimension five is
ruled out. The uniform Hessian result concerns affine source-target changes and
potentials affine in the three auxiliary variables below the specified quadratic
leading term; it is not a nonlinear-equivalence theorem.

## Persistent identifiers

- DOI: <https://doi.org/10.5281/zenodo.23477226>
- Repository: <https://github.com/ipitchford/wmy5-graded-obstructions>
- Evidence Press release: <https://evidencepress.org/releases/wmy5-graded-obstructions/>
