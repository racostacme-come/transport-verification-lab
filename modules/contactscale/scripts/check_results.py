"""Recompute all recorded numerical outputs; allow platform-level roundoff."""

import argparse
import json
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose

from contactscale.study import run

ROOT = Path(__file__).resolve().parents[1]
ROUND_OFF_BUDGETS = {
    "dense_exponential_error_K": 1e-8,
    "heat_quadrature_error_J_m2": 1e-6,
    "max_energy_balance_residual_J_m2": 0.02,
    "domain_doubling_flux_relative_at_10s": 1e-8,
}


def compare(actual, expected, key=""):
    if key in ROUND_OFF_BUDGETS:
        assert abs(actual) < ROUND_OFF_BUDGETS[key], (key, actual)
    elif isinstance(expected, dict):
        assert actual.keys() == expected.keys()
        for name in expected:
            compare(actual[name], expected[name], name)
    elif isinstance(expected, list):
        assert len(actual) == len(expected)
        for a, e in zip(actual, expected, strict=True):
            compare(a, e, key)
    else:
        assert_allclose(actual, expected, rtol=2e-5, atol=2e-9, err_msg=key)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "build/reproduction")
    args = parser.parse_args()
    actual = run(args.output)
    expected = json.loads((ROOT / "results/metrics.json").read_text())
    compare(actual, expected)
    for name in ("history.csv", "refinement.csv", "smooth.csv"):
        a = np.genfromtxt(args.output / name, names=True, delimiter=",")
        e = np.genfromtxt(ROOT / "results" / name, names=True, delimiter=",")
        assert a.dtype.names == e.dtype.names
        for column in a.dtype.names:
            assert_allclose(a[column], e[column], rtol=2e-5, atol=2e-8, err_msg=f"{name}:{column}")
    assert min(actual["contact_flux_orders_at_1s"][-2:]) > 1.95
    assert min(actual["smooth_orders"]) > 1.98
    assert actual["refinement"][-1]["flux_error_at_1s_relative"] < 1e-3
    print("All JSON quantities and 514 CSV rows reproduced within documented tolerances.")


if __name__ == "__main__":
    main()
