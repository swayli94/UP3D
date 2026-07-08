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
    # Smooth parabola: only O(h^2 cp'') curvature, no odd-even content. A pure
    # parabola has constant nonzero cp'', so d2 = cp''*h^2 per point is a real
    # (small, refinement-shrinking) floor -- ~2.5e-3 here, ~100x below the
    # sawtooth of test_sawtooth_dominates_smooth.
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


def test_refinement_shrinks_smooth_metric():
    m_coarse = cp_oscillation_metric(*_smooth_pocket(30), CP_STAR)["metric"]
    m_fine = cp_oscillation_metric(*_smooth_pocket(120), CP_STAR)["metric"]
    # Smooth second difference ~ h^2: 4x points -> ~16x smaller.
    assert m_fine < 0.25 * m_coarse


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
