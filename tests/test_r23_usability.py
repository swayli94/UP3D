"""Locks for bench/usability.py -- the answer anchor.

★ G-CADENCE: the tool lives in bench/ and its locks live here, because CLAUDE.md records
that bench scripts are in no test collection and therefore run on no cadence. A usability
criterion without a tests lock would reproduce the problem it exists to fix.

No solve is needed: every case below is a committed number.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "bench"))
from usability import ANCHOR_RATIO_MAX, assess, assess_set   # noqa: E402

#: R22's six medium legs at alpha 1.25 / seed 0 (bench/studies/r22_c_interval/results)
R22 = [dict(C=1.0, seed=0, cl_p=0.857947, converged=False, n_limited=1135, n_floored=0),
       dict(C=1.1, seed=0, cl_p=0.040267, converged=True, n_limited=0, n_floored=0),
       dict(C=1.25, seed=0, cl_p=2.502974, converged=False, n_limited=5217, n_floored=1543),
       dict(C=1.5, seed=0, cl_p=0.341340, converged=True, n_limited=0, n_floored=0),
       dict(C=2.0, seed=0, cl_p=0.473876, converged=False, n_limited=0, n_floored=0),
       dict(C=2.5, seed=0, cl_p=0.554894, converged=False, n_limited=0, n_floored=0)]
#: the nine converged, zero-clamp legs of bench/gate_results/m1_gate_default.csv
GOOD = [dict(level="coarse", C=1.0, seed=0, cl_p=0.42473),
        dict(level="coarse", C=1.0, seed=5, cl_p=0.413008),
        dict(level="coarse", C=1.5, seed=0, cl_p=0.407978),
        dict(level="coarse", C=1.5, seed=5, cl_p=0.407108),
        dict(level="coarse", C=3.0, seed=0, cl_p=0.377355),
        dict(level="coarse", C=3.0, seed=5, cl_p=0.377442),
        dict(level="medium", C=1.5, seed=0, cl_p=0.34134),
        dict(level="medium", C=1.5, seed=5, cl_p=0.34134),
        dict(level="medium", C=3.0, seed=5, cl_p=0.289614)]
for g in GOOD:
    g.update(converged=True, n_limited=0, n_floored=0)


def test_catches_the_spurious_root():
    """U-CATCH: C=1.10 converged to 2.28e-13 with zero clamps and must still be rejected;
    C=1.5 must be accepted. This is the single thing the module exists to do."""
    res = assess_set(R22, axes=("seed",))
    by = {r["cl_p"]: r for r in res["legs"]}
    spur, good = by[0.040267], by[0.341340]
    assert not spur["usable"], spur
    assert "anchor_ratio" in spur["reason"] and "OUTLIER" in spur["reason"]
    assert spur["converged"] and not spur["clamped"], (
        "the point is that it passes the OLD checks")
    assert good["usable"], good


def test_rejects_none_of_the_committed_good_legs():
    """U-NOFALSE: a filter that rejects everything is trivially effective. Not one of the
    nine legs the project cites as good may be rejected."""
    res = assess_set(GOOD, axes=("level", "seed"))
    bad = [r for r in res["legs"] if not r["usable"]]
    assert not bad, bad
    assert res["n_usable"] == len(GOOD)


def test_single_leg_spread_is_undefined():
    """U-AXES: R13 -- a spread over fewer than two usable legs is UNDEFINED, never a
    number. Seed 0 at medium has exactly one usable C."""
    res = assess_set([l for l in R22], axes=("seed",))
    s = res["spreads"][("seed", 0)]
    assert s["n"] == 1 and s["spread"].startswith("UNDEFINED"), s


def test_two_usable_legs_give_a_spread_with_every_denominator():
    """R13: a spread must be quoted with its denominator, so all three are returned."""
    legs = [dict(seed=5, cl_p=0.341340, converged=True, n_limited=0, n_floored=0),
            dict(seed=5, cl_p=0.289614, converged=True, n_limited=0, n_floored=0)]
    s = assess_set(legs, axes=("seed",))["spreads"][("seed", 5)]
    assert s["n"] == 2
    assert set(s) >= {"rel_min", "rel_max", "rel_mean"}
    #: ★ compute the expectation from the two inputs rather than transcribing it --
    #: my first version copied 0.178563 off a rounded print and failed at 4e-05.
    exp = (0.341340 - 0.289614) / 0.289614
    assert abs(s["rel_min"] - exp) < 1e-12, (s, exp)


def test_clamped_ever_sees_transient_clamping():
    """U-HIST: R15 -- the scalars are final-step, so 0/0 does not mean never clamped."""
    r = assess(0.34, converged=True, n_limited=0, n_floored=0,
               clamp_history=[[0, 0], [7, 0], [0, 0]], consensus=0.34)
    assert r["usable"] and r["clamped_ever"], r


def test_missing_anchor_is_reported_not_silently_passed():
    """★ Passing an unanchored leg silently is exactly how C=1.10 got in."""
    r = assess(0.040267, converged=True, n_limited=0, n_floored=0, consensus=None)
    assert not r["usable"] and "NO ANCHOR" in r["reason"], r


def test_reason_is_specific():
    """U-REASON: never just False."""
    r = assess(0.34, converged=False, n_limited=12, n_floored=3, consensus=0.34)
    assert "not converged" in r["reason"] and "clamped 12/3" in r["reason"], r


def test_axes_are_mandatory():
    """★ R22's erratum: the omission that produced a wrong headline must be impossible
    to make silently."""
    with pytest.raises(ValueError, match="axes"):
        assess_set(GOOD, axes=())


def test_calibration_is_documented():
    """★ The threshold is a CALIBRATION; its derivation must stay in the docstring."""
    import usability
    d = usability.__doc__
    assert "CALIBRATION" in d and "1.47" in d and "8.48" in d
    assert 1.47 < ANCHOR_RATIO_MAX < 8.48, "the threshold must sit in the measured gap"
