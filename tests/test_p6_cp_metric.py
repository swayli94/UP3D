"""G6.1 surface-Cp smoothness metric (roadmap P6, design.md §3.1-3.2).

Correctness of `post.section_cut.cp_oscillation_metric` on synthetic curves --
no solve required. The metric must be ~0 on a smooth supersonic run, jump on a
2-cell sawtooth of known amplitude, shrink under h-refinement of a smooth curve,
and degrade gracefully when there is no resolved supersonic pocket.
"""

import numpy as np

from pyfp3d.post.section_cut import cp_oscillation_metric

# A representative sonic level; the supersonic run is cp < CP_STAR.
CP_STAR = -0.43


def _smooth_pocket(n, cp_min=-1.1):
    """A smooth suction pocket dipping below CP_STAR: cp(x) parabola in x/c,
    minimum cp_min at mid-pocket, rising back toward CP_STAR at both ends."""
    x = np.linspace(0.05, 0.6, n)
    xm = 0.5 * (x[0] + x[-1])
    span = 0.5 * (x[-1] - x[0])
    cp = cp_min + (CP_STAR - cp_min) * ((x - xm) / span) ** 2
    return x, cp


def test_smooth_curve_metric_near_zero():
    x, cp = _smooth_pocket(60)
    m = cp_oscillation_metric(x, cp, CP_STAR)
    assert m["n_super"] > 50
    # Monotone-derivative parabola: the slope reverses only once (at the
    # vertex), so the sign-gated sawtooth energy is ~0 -- the metric now
    # correctly ignores smooth curvature (unlike the raw second-difference RMS).
    assert m["n_reversals"] <= 2
    assert m["metric"] < 1e-3


def test_monotone_shock_not_counted_as_sawtooth():
    """A steep but MONOTONE compression (a shock-like ramp) must not register:
    no slope reversals -> ~0, even though its second difference is large."""
    x = np.linspace(0.05, 0.6, 40)
    cp = -1.1 + 0.9 / (1.0 + np.exp(-(x - 0.45) / 0.01))  # smooth monotone rise
    cp = np.minimum(cp, CP_STAR - 1e-6)                    # keep it supersonic
    m = cp_oscillation_metric(x, cp, CP_STAR)
    assert m["metric"] < 5e-3


def test_sawtooth_dominates_smooth():
    x, cp = _smooth_pocket(60)
    smooth = cp_oscillation_metric(x, cp, CP_STAR)["metric"]
    amp = 0.05
    saw = cp.copy()
    saw[1:-1:2] += amp   # alternating +/- on interior points -> 2-cell serration
    saw[2:-1:2] -= amp
    serrated = cp_oscillation_metric(x, saw, CP_STAR)["metric"]
    # The sawtooth must register far above the smooth baseline.
    assert serrated > 50 * smooth
    # A |d2| ~ 4*amp per point over a ~0.67 Cp range -> metric ~ 0.2-0.3.
    assert 0.1 < serrated < 0.5


def test_smooth_metric_near_zero_at_any_resolution():
    # A smooth curve has no sawtooth, so the sign-gated metric stays ~0 at both
    # coarse and fine sampling (it does NOT grow with resolution -- the failure
    # mode of a raw second-difference RMS on a strongly curved pocket).
    for n in (30, 60, 120):
        m = cp_oscillation_metric(*_smooth_pocket(n), CP_STAR)
        assert m["metric"] < 1e-3


def test_unordered_input_is_sorted():
    x, cp = _smooth_pocket(60)
    perm = np.random.default_rng(0).permutation(x.size)
    m_ref = cp_oscillation_metric(x, cp, CP_STAR)["metric"]
    m_shuf = cp_oscillation_metric(x[perm], cp[perm], CP_STAR)["metric"]
    assert np.isclose(m_ref, m_shuf, rtol=1e-12, atol=1e-15)


def test_no_supersonic_run_returns_nan():
    x = np.linspace(0.0, 1.0, 40)
    cp = np.full_like(x, 0.2)   # entirely subcritical (cp > cp_star)
    m = cp_oscillation_metric(x, cp, CP_STAR)
    assert m["n_super"] == 0
    assert np.isnan(m["metric"])


def test_too_short_run_returns_nan():
    x = np.array([0.1, 0.2, 0.3])
    cp = np.array([-0.6, -0.7, -0.6])  # only 3 supersonic points < min_points=5
    m = cp_oscillation_metric(x, cp, CP_STAR)
    assert m["n_super"] == 3
    assert np.isnan(m["metric"])
