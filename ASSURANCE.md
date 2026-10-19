# Assurance and evidence boundaries

## Written mathematics

The principal warrant is the manuscript's exact argument. The two compression
theorems quantify only over their displayed complete bigraded families. The
uniform Hessian theorem quantifies over nonzero `kappa`, `r`, and `t`, arbitrary
polynomials in `(x,y)` as coefficients of terms affine in `(a,b,c)`, and
invertible affine source--target changes. It does not quantify over nonlinear
equivalences or all five-variable polynomial maps.

## Executed checks

On 19 October 2026 the standard-library verifier and four mutation tests passed
in ordinary and optimised Python. The short chain contradiction, both monomial
inventories and the symbolic Hessian minor also passed. The review package's
separate-copy regeneration and independent SymPy calculations were replayed
during revision. Exact arithmetic supports the displayed identities; it does
not convert the work into proof-assistant formalisation.

The reviewer-supplied audit and subsequent coordinated AI checks are internal
research assurance, not external specialist review or human peer review.

## Evidence classes

- Exact analytic proof: manuscript propositions and theorems.
- Exact finite replay: standard-library certificate checker, short linear
  contradiction, grading inventory and symbolic Hessian minor.
- Negative controls: four deliberate certificate corruptions, all rejected.
- Independent implementation within the same AI-assisted workflow: reviewer
  SymPy reconstruction and separate-copy certificate regeneration.
- Exploratory evidence only: broader Wolfram Language unit-ideal calculations
  without shipped small independent certificates.

## Novelty and limitations

A bounded source audit located close antecedents for chain shortening,
degree-filtered retargeting and Hessian coefficient-span reasoning, but no exact
predecessor for the two displayed bigraded exclusions or the uniform perturbation
statement. This is not exhaustive priority certification. Nonzero variation in
`r,t` is coordinate robustness, not a separate novelty claim. The requested
five-dimensional counterexample remains neither constructed nor ruled out.
