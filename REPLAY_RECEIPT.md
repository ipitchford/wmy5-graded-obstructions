# Replay receipt

Fresh revision replay was performed on 19 October 2026 from this candidate tree.

| Obligation | Result |
|---|---|
| Standard-library exact verifier, ordinary Python | Passed all five groups |
| Four deliberate certificate corruptions, ordinary Python | All rejected |
| Standard-library exact verifier, `-O -I -B` | Passed all five groups |
| Four deliberate certificate corruptions, `-O -I -B` | All rejected |
| Short chain contradiction | Passed; eight equations force the variables and the ninth has residual 16,367,616 |
| Parallel and chain monomial inventories | Passed; totals 27 and 28 |
| Uniform Hessian selected minor | Passed; determinant `-1304596316160*kappa^15*r^10*t^6` |
| pdfLaTeX build and cross references | Passed; nine-page PDF |

This receipt records producer-side execution. The supplied reviewer package also
contains a separately implemented SymPy reconstruction and byte-identical
certificate regeneration, but both remain within an AI-assisted coordinated
review workflow rather than external peer review.
