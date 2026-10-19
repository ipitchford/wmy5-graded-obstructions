#!/usr/bin/env python3
"""Replay the nine-equation contradiction in the chain exclusion.

Usage: python3 -I check_chain_short_identity.py [certificate.json]
Only the Python standard library is required. Certificate indices are zero-based.
"""
import json
import sys
from fractions import Fraction
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


certificate_path = (
    Path(sys.argv[1])
    if len(sys.argv) > 1
    else Path(__file__).with_name("certificates") / "chain_final_certificate.json"
)
certificate = json.loads(certificate_path.read_text())
variables = certificate["variables"]
require(variables == ["p0", "p1", "p2", "p3", "v0", "v1", "v2", "v3", "v4", "w0", "w1", "w2"], "unexpected certificate variables")
inputs = certificate["inputs"]
require(len(inputs) == 27, "expected 27 final equations")
require(all(sum(monomial) <= 1 for polynomial in inputs for monomial, numerator, denominator in polynomial if numerator), "a final equation is not affine-linear")

weights = {
    0: Fraction(33177600),
    1: Fraction(5806080),
    2: Fraction(1019904),
    5: Fraction(11520),
    6: Fraction(-92897280, 31),
    7: Fraction(6336),
    8: Fraction(-556608),
    12: Fraction(-2073600, 31),
    25: Fraction(1),
}
combined = {}
for index, multiplier in weights.items():
    for monomial, numerator, denominator in inputs[index]:
        key = tuple(monomial)
        combined[key] = combined.get(key, Fraction(0)) + multiplier * Fraction(numerator, denominator)
combined = {key: value for key, value in combined.items() if value}
require(combined == {(0,) * len(variables): Fraction(16367616)}, "nine-input linear identity does not yield the required residual")

forced = {
    "w1": 0,
    "w2": 0,
    "w0": 0,
    "p1": 132,
    "p2": -744,
    "p3": 1392,
    "v3": -1920,
    "v4": 3264,
}
derived = {}
for index, target in [(6, "w1"), (8, "w2"), (12, "w0"), (0, "p1"), (1, "p2"), (2, "p3"), (5, "v3"), (7, "v4")]:
    coefficient = Fraction(0)
    constant = Fraction(0)
    for monomial, numerator, denominator in inputs[index]:
        value = Fraction(numerator, denominator)
        names = [name for name, degree in zip(variables, monomial) if degree]
        if not names:
            constant += value
        elif names == [target]:
            coefficient += value
        else:
            require(len(names) == 1 and names[0] in derived, f"equation {index} does not have the expected triangular form")
            constant += value * derived[names[0]]
    require(coefficient != 0, f"equation {index} does not determine {target}")
    derived[target] = -constant / coefficient
require(derived == forced, "forced values do not match the manuscript")


def evaluate(polynomial):
    result = Fraction(0)
    for monomial, numerator, denominator in polynomial:
        term = Fraction(numerator, denominator)
        for name, degree in zip(variables, monomial):
            if degree:
                term *= forced[name] ** degree
        result += term
    return result


require(all(evaluate(inputs[index]) == 0 for index in weights if index != 25), "a selected forcing equation does not vanish")
require(evaluate(inputs[25]) == 16367616, "the ninth equation does not yield the required residual")
require(25344 * 3264 + 34560 * (-1920) == 16367616, "displayed integer arithmetic mismatch")
print("PASS: all 27 final equations are affine-linear")
print("PASS: eight selected equations force the stated values")
print("PASS: input 25 then equals 16,367,616")
print("PASS: nine-input rational linear identity equals 16,367,616")
