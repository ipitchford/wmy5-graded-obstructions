# Supplemental exact symbolic computations

Date: 14 October 2026. Wolfram Language kernel accessed through the Wolfram connector.
These are supplemental results, not obligations checked by `verify.py`.

The searches represented by `broader_search.wl` returned the unit Groebner basis
`{1}` in each of the following cases:

1. Parallel auxiliary weights (2,2), first component x^2*y+x^3*z, adjustable lower
   components, and the extra cross-coupling t*x*b.
2. Parallel weights (2,2), first component x^2*y, adjustable lower components,
   without the extra cross-coupling.
3. Parallel weights (2,2), first component x^3*z, adjustable lower components,
   without the extra cross-coupling.
4. Unequal weights (2,3), first component x^2*y+x^3*z, adjustable lower components
   and the listed cubic, quartic and auxiliary terms.

The supplied script reconstructs the equations rather than importing unexplained
coefficient data. The kernel's unit-ideal results are exact CAS outputs. No small
independent rational proof certificate for these broader systems is included.
They must not be described as covered by the standard-library verifier.

A syntax error in an initial multi-case Wolfram loop was corrected before the
successful evaluations; that failed call supplied no mathematical result.
Several local SymPy Groebner attempts on the broader systems exceeded the tool
execution cap. Such interruptions are not nonexistence certificates.

None of these families exhausts all polynomial maps in five variables.
