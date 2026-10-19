#!/usr/bin/env python3
"""Enumerate the monomials permitted by the manuscript's two gradings."""
import json
from itertools import product
from pathlib import Path


VARIABLES = ("x", "y", "z", "a", "b")


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def describe(monomial):
    parts = []
    for name, exponent in zip(VARIABLES, monomial):
        if exponent:
            parts.append(name if exponent == 1 else f"{name}^{exponent}")
    return "*".join(parts) or "1"


def enumerate_family(weights, second):
    output = []
    for component in range(5):
        target = weights[component] + 1
        terms = []
        ranges = [range(target // weight + 1) for weight in weights[1:]]
        for exponents in product(*ranges):
            if sum(weight * exponent for weight, exponent in zip(weights[1:], exponents)) != target:
                continue
            x_power = second[component] - sum(degree * exponent for degree, exponent in zip(second[1:], exponents))
            if x_power < 0:
                continue
            monomial = (x_power, *exponents)
            require(sum(weight * exponent for weight, exponent in zip(weights, monomial)) == target, "first grading mismatch")
            require(sum(degree * exponent for degree, exponent in zip(second, monomial)) == second[component], "second grading mismatch")
            terms.append(monomial)
        output.append({
            "component": component + 1,
            "monomials": terms,
            "terms": [describe(monomial) for monomial in terms],
            "count": len(terms),
        })
    return output


def main():
    result = {
        "parallel": enumerate_family((0, 1, 1, 2, 2), (1, -1, -2, -1, -2)),
        "chain": enumerate_family((0, 1, 1, 2, 3), (1, -1, -2, -2, -2)),
    }
    parallel_total = sum(component["count"] for component in result["parallel"])
    chain_total = sum(component["count"] for component in result["chain"])
    require([component["count"] for component in result["parallel"]] == [2, 5, 4, 8, 8], "parallel component counts mismatch")
    require([component["count"] for component in result["chain"]] == [2, 4, 4, 7, 11], "chain component counts mismatch")
    require(parallel_total - 8 - 2 == 17, "parallel parameter count mismatch")
    require(chain_total - 8 - 2 - 2 == 16, "chain parameter count mismatch")
    result["parameter_counts"] = {
        "parallel_all_bigraded_coefficients": parallel_total,
        "parallel_after_fixed_base_and_normalised_feedthroughs": 17,
        "parallel_after_trace_and_second_coefficient_relations": 13,
        "chain_all_bigraded_coefficients": chain_total,
        "chain_after_fixed_base_normalised_feedthroughs_and_trace": 16,
        "chain_after_displayed_forced_zero_relations": 12,
    }
    output = Path(__file__).with_name("grading_audit_report.json")
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"parameter_counts": result["parameter_counts"]}, indent=2))


if __name__ == "__main__":
    main()
