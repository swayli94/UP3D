"""
Track A / A1 -- the timing instrumentation contract (GA1.1 inertness +
faithfulness at test scale).

The A1 bottleneck study only means something if the `timings` dict every
driver now reports is (a) on the SAME schema across all four drivers, so
their phase breakdowns are comparable, and (b) FAITHFUL -- the accounted
phases explain almost all of the wall clock, so "the bottleneck is phase X"
is a real claim and not an artifact of unmeasured time. These are cheap
subsonic checks on NACA coarse; the expensive transonic faithfulness is
re-checked inside the demo (GA1.1 over the real case matrix).

The drivers are otherwise covered by their own P/B suites; here we only
assert the instrumentation, and that it did not perturb the answer (the
four methods still agree on gamma to the known <1%).
"""

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.solve.continuation import solve_transonic_lifting
from pyfp3d.solve.newton import solve_newton_lifting, solve_newton_transonic
from pyfp3d.solve.newton_ls import solve_multivalued_newton
from pyfp3d.solve.picard import solve_subsonic_lifting
from pyfp3d.solve.picard_ls import solve_multivalued_lifting
from pyfp3d.solve.timing import PHASES
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

MESH = "cases/meshes/naca0012_2.5d/coarse.msh"
M_INF = 0.5
ALPHA = 1.25


@pytest.fixture(scope="module")
def cut():
    mesh = read_mesh(MESH)
    mc, wc = cut_wake(mesh)
    return mc, wc


@pytest.fixture(scope="module")
def mvop():
    mesh = read_mesh(MESH)
    z = mesh.nodes[:, 2]
    wls = WakeLevelSet(
        np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
        direction=(1.0, 0.0, 0.0),
    )
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return mesh, MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)


def _assert_schema(t):
    """Canonical keys present, non-negative, and phases explain the wall."""
    for k in PHASES + ("wall", "other"):
        assert k in t, f"missing timings key {k!r}"
    for p in PHASES:
        assert t[p] >= -1e-9, f"negative phase {p}: {t[p]}"
    assert t["wall"] > 0.0
    # `other` = wall - sum(phases); faithfulness = it is small AND the phases
    # never over-count (other >= 0 up to timer jitter).
    assert t["other"] >= -1e-3, f"phases over-counted the wall: other={t['other']}"
    assert abs(t["wall"] - (sum(t[p] for p in PHASES) + t["other"])) < 1e-6
    # subsonic NACA coarse: >90% of the wall must be in the named phases
    assert t["other"] < 0.10 * t["wall"], (
        f"unaccounted time {t['other']:.3f}s is "
        f"{100 * t['other'] / t['wall']:.1f}% of {t['wall']:.3f}s wall")


def _assert_step_records(recs, n_iter):
    assert len(recs) == n_iter, f"{len(recs)} records vs {n_iter} iters"
    for r in recs:
        assert "residual" in r and "wall_cum_s" in r
        assert {"t_" + p for p in PHASES} <= set(r)
    # wall_cum_s is monotone non-decreasing
    cum = [r["wall_cum_s"] for r in recs]
    assert all(b >= a - 1e-6 for a, b in zip(cum, cum[1:]))


# --- single-level drivers: schema + step records --------------------------
#
# The LS drivers are run with farfield="vortex" throughout, to MATCH the
# conforming path's vortex far field -- the LS defaults differ (Picard
# vortex, Newton neumann) and neumann shifts gamma ~4% (a real O(Gamma/R)
# outer-boundary effect, not an instrumentation artifact), which would
# swamp the "the four methods agree" check below.

def _cp(cut):
    mc, wc = cut
    return solve_subsonic_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                                  n_picard_max=40)


def _cn(cut):
    mc, wc = cut
    return solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA)


def _lp(mvop):
    mesh, mv = mvop
    return solve_multivalued_lifting(mv, mesh, m_inf=M_INF, alpha_deg=ALPHA,
                                     n_outer_max=80, tol_residual=1e-8,
                                     farfield="vortex")


def _ln(mvop):
    mesh, mv = mvop
    return solve_multivalued_newton(mv, mesh, m_inf=M_INF, alpha_deg=ALPHA,
                                    n_seed=10, farfield="vortex")


def test_conforming_picard(cut):
    r = _cp(cut)
    assert r["converged"]
    _assert_schema(r["timings"])
    _assert_step_records(r["step_records"], r["n_picard"])
    # the Picard tell: the inner Kutta secant re-solves the frozen matrix,
    # so there is more than one linear solve per outer on average
    assert r["n_solves_total"] >= r["n_picard"]


def test_conforming_newton(cut):
    r = _cn(cut)
    assert r["converged"]
    t = r["timings"]
    _assert_schema(t)
    # legacy aliases preserved bit-for-bit for cases/demo/p8_newton
    assert t["jacobian"] == t["assembly"]
    assert t["amg_setup"] == t["precond"]
    assert t["gmres"] == t["linsolve"]
    _assert_step_records(r["step_records"], len(r["residual_history"]))
    assert r["n_refactor"] == 0            # amg path, no splu


def test_ls_picard(mvop):
    r = _lp(mvop)
    assert r["converged"]
    _assert_schema(r["timings"])
    _assert_step_records(r["step_records"], r["n_outer"])


def test_ls_newton(mvop):
    r = _ln(mvop)
    assert r["converged"]
    _assert_schema(r["timings"])
    _assert_step_records(r["step_records"], r["n_newton"])
    # A1 added the circulation trajectory the LS Newton never had
    assert all("gamma_mean" in s for s in r["step_records"])


def test_four_methods_agree(cut, mvop):
    """The comparison is only meaningful if the four drivers converge to the
    same circulation (GA1.3 at test scale): instrumentation did not perturb
    any of them, and conforming vs level-set agree to the known <1% under a
    matched far field."""
    gammas = np.array([
        float(_cp(cut)["gamma"][0]),
        float(_cn(cut)["gamma"][0]),
        float(_lp(mvop)["gamma"]),
        float(_ln(mvop)["gamma"]),
    ])
    spread = float(np.ptp(gammas)) / float(np.mean(gammas))
    assert spread < 0.01, f"gamma spread {spread:.4%} across methods: {gammas}"


# --- ramp wrappers: every level carries wall_s + timings, plus a total ----

def test_conforming_picard_ramp_levels(cut):
    """The conforming Picard ramp kept NO per-level record before A1."""
    mc, wc = cut
    r = solve_transonic_lifting(mc, wc, m_inf=0.5, alpha_deg=ALPHA,
                                m_start=0.45, dm=0.05, n_picard_seed=40,
                                n_picard_eval=60, max_gamma_evals=4)
    assert "level_results" in r and len(r["level_results"]) >= 1
    tot = r["timings_total"]
    _assert_schema(dict(tot, other=max(tot["other"], 0.0)))
    for lr in r["level_results"]:
        assert lr["wall_s"] > 0.0
        for k in PHASES:
            assert k in lr["timings"]
        assert "step_records" in lr
    # timings_total is the sum of the level timings (the anti-footgun)
    s = sum(lr["timings"]["wall"] for lr in r["level_results"])
    assert abs(tot["wall"] - s) < 1e-6


def test_conforming_newton_ramp_total(cut):
    mc, wc = cut
    r = solve_newton_transonic(mc, wc, m_inf=0.5, alpha_deg=ALPHA,
                               m_start=0.45, dm=0.05)
    assert "timings_total" in r and "level_results" in r
    for lr in r["level_results"]:
        assert "step_records" in lr and lr["wall_s"] > 0.0
    s = sum(lr["timings"]["wall"] for lr in r["level_results"])
    assert abs(r["timings_total"]["wall"] - s) < 1e-6
