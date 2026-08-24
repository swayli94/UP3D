"""Locks for the 2026-08-24 M1 re-spec (user ruling). Zero solves.

The ruling: (a) deleted as a criterion, (b) and (c) relaxed from 3 % to 20 %, and split
into two INDEPENDENT metrics M1b and M1c.

★ Why these locks exist here rather than in bench/: `bench/run_m1_gate.py` is in no test
collection, so a spec change there would run on no cadence -- the same G-CADENCE argument
as tests/F/test_F05_usability_anchor.py. `_criteria` takes rows, so the committed gate CSV can be
fed straight through it with no solve.
"""
import csv
import importlib.util
import os
import sys

import pytest

from tests.conftest import REPO_ROOT
ROOT = str(REPO_ROOT)   #: ★ 与目录深度无关（2026-08-24 重编号）
CSV = os.path.join(ROOT, "bench/gate_results/m1_gate_default.csv")


@pytest.fixture(scope="module")
def gate():
    spec = importlib.util.spec_from_file_location(
        "m1_gate", os.path.join(ROOT, "bench/run_m1_gate.py"))
    m = importlib.util.module_from_spec(spec)
    argv, sys.argv = sys.argv, ["m1_gate"]
    try:
        spec.loader.exec_module(m)
    finally:
        sys.argv = argv
    return m


@pytest.fixture(scope="module")
def rows():
    out = []
    for r in csv.DictReader(open(CSV)):
        out.append(dict(level=r["level"], C=float(r["C"]),
                        n_picard_seed=int(r["n_picard_seed"]),
                        cl_p=float(r["cl_p"]) if r["cl_p"] else None,
                        x_shock=float(r["x_shock"]) if r["x_shock"] else None,
                        converged=r["converged"] == "True"))
    return out


def test_thresholds_are_the_ruling(gate):
    """★ The ruling's numbers, so a silent drift back to 3 % (or on to something looser)
    is a red test rather than a quiet change."""
    assert gate.M1B_TOL == 0.20
    assert gate.M1C_TOL == 0.20


def test_shock_reference_is_the_committed_one(gate):
    """★ The script used to carry 0.61 +- 0.02, which appears in nine documents and in NO
    reference file. The committed reference is 0.62 +- 0.03
    (cases/reference_data/naca0012_m080/shock_reference.csv)."""
    assert (gate.SHOCK_REF, gate.SHOCK_TOL) == (0.62, 0.03)


def test_criteria_returns_two_independent_verdicts(gate, rows):
    """The split is the structural half of the ruling: two metrics, not one M1."""
    v = gate._criteria([r for r in rows if r["n_picard_seed"] == 0], 0)
    assert set(v) == {"M1b", "M1c"}


@pytest.mark.parametrize("seed,m1b,m1c", [(0, True, False), (5, True, True)])
def test_measured_verdicts_on_the_committed_data(gate, rows, seed, m1b, m1c):
    """★★ The re-spec's measured outcome, locked. M1b passes on both seeds; M1c passes on
    seed 5 and fails on seed 0 ONLY because medium has a single converged leg there -- a
    COVERAGE failure, which no tolerance can cure."""
    v = gate._criteria([r for r in rows if r["n_picard_seed"] == seed], seed)
    assert v["M1b"] is m1b and v["M1c"] is m1c, v


def test_the_margins_are_recorded(gate, rows):
    """★ M1c's medium margin at seed 5 is thin (17.86 % against 20 %); locking it means a
    regression that eats 2 pp becomes visible instead of flipping a verdict silently."""
    sub = [r for r in rows if r["n_picard_seed"] == 5 and r["level"] == "medium"
           and r["converged"] and r["cl_p"] is not None]
    v = [r["cl_p"] for r in sub]
    spread = (max(v) - min(v)) / min(v)
    assert abs(spread - 0.178563) < 1e-4, spread
    assert spread < gate.M1C_TOL, "margin gone"
    b = {r["level"]: r["cl_p"] for r in rows
         if r["n_picard_seed"] == 0 and r["C"] == 1.5 and r["cl_p"] is not None}
    d = abs((b["medium"] - b["coarse"]) / b["coarse"])
    assert abs(d - 0.16333) < 1e-4, d
