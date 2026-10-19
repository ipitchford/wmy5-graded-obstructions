#!/usr/bin/env python3
"""Exact symbolic check of the uniform Hessian-span minor."""
import json
from pathlib import Path

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


x, y, a, b, c, r, t, kappa = sp.symbols("x y a b c r t kappa")
variables = [x, y, a, b, c]
u = 1 + x * y
A = a * u**3 + r * x * b * u**2 - t * x**3 * c
phi = sp.expand(kappa * A**2)
pairs = [(i, j) for i in range(5) for j in range(i, 5)]
columns = [
    sp.expand((2 if i != j else 1) * sp.diff(phi, variables[i], variables[j]))
    for i, j in pairs
]
selectors = [b**2, a**2, x**2 * a**2, b, a, x * a, x**2 * a, x**3 * b, x**4 * a, x * y, x, x**3, x**2, x**4, x**6]
matrix = sp.zeros(15, 15)
rows = []
for row_index, monomial in enumerate(selectors):
    exponent = sp.Poly(monomial, *variables).monoms()[0]
    entries = []
    for column_index, polynomial in enumerate(columns):
        value = sp.factor(sp.Poly(polynomial, *variables).coeff_monomial(exponent))
        matrix[row_index, column_index] = value
        if value != 0:
            entries.append({"pair": list(pairs[column_index]), "coefficient": str(value)})
    rows.append({"selector": str(monomial), "entries": entries})

determinant = sp.factor(matrix.det())
expected = -1304596316160 * kappa**15 * r**10 * t**6
require(sp.expand(determinant - expected) == 0, "selected determinant factorisation mismatch")
require(abs(int(expected.subs({kappa: 1, r: 3, t: 1}))) == 77035107872931840, "fixed specialisation determinant mismatch")
report = {
    "status": "PASS",
    "arithmetic": "exact symbolic SymPy",
    "claim": "The selected Hessian-coefficient minor is nonzero whenever kappa, r and t are nonzero.",
    "A": "a*(1+x*y)^3 + r*x*b*(1+x*y)^2 - t*x^3*c",
    "selectors": [str(item) for item in selectors],
    "pairs": [list(pair) for pair in pairs],
    "rows": rows,
    "determinant": str(determinant),
    "specialisation_absolute_determinant": "77035107872931840",
    "boundary": "The manuscript's auxiliary-degree proof handles arbitrary affine-in-(a,b,c) perturbations; this script checks the load-bearing exact minor.",
}
Path(__file__).with_name("generalized_hessian_check_report.json").write_text(json.dumps(report, indent=2) + "\n")
print(f"PASS: selected determinant = {determinant}")
