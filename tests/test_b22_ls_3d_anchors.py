"""
Track B / B22 — the N3 gap closed: GATED anchor locks on the core 3-D
level-set numbers.

Why this file exists (the 2026-07-19 Kimi inspection, finding N3, confirmed
against the gated tier at its widest): the 3-D level-set numbers were locked
by NOTHING. B15's tests run entirely on the 2.5-D NACA mesh; the M6 ramp's
m_final / gamma / M_max / clamp counts were demo numbers no test asserted; the
b14 gated A/B compares two solver paths against EACH OTHER under the same
code, so a re-baseline that moves both arms together raises no alarm — and
that is exactly what happened twice in two days (B20: gamma 0.088338 ->
0.071909 with the ramp stalling at M0.6625; B21: restored to 0.088343 at
M0.84). The suite, with every gate enabled, was green through both.

These tests re-run the committed B15 M6 recipes and assert ABSOLUTE anchors
(band, not bitwise: the walk selection is discrete and thread-scheduling can
flip near-ties; the two B21 sweep variants agreed on gamma to 5.5e-6 relative,
so the 1e-4 bands carry ~20x margin while still catching any B20-sized move
by four orders of magnitude).

Anchor provenance (committed artifacts, 2026-07-19):
  medium — cases/analysis/c1_ls_jacobian_fd/results/n1_freeze_fix_sweep.csv
           (freeze_tol=1e-3 row = the committed recipe) and the refreshed
           cases/demo/b15_ls_newton_ramp/results/summary.csv;
  coarse — the refreshed cases/demo/b14_schur_precond/results/schur_ab.csv
           (part 4, trans/lagged arm = the committed recipe).

If a legitimate change moves these numbers: re-baseline them EXPLICITLY
(update the anchors, commit the regenerated evidence, and follow the
re-baseline erratum checklist in CLAUDE.md workflow step 5). That is the
point — the move must be a decision, not a silence.
"""

import os
from pathlib import Path

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.meshgen.wing3d import B_SEMI, x_te
from pyfp3d.solve.newton_ls import (
    B_NEWTON_M6_DEFAULTS,
    solve_multivalued_newton_transonic,
)
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

REPO_ROOT = Path(__file__).parent.parent
M6_DIR = REPO_ROOT / "cases" / "meshes" / "onera_m6_wakefree"
GATES = os.environ.get("PYFP3D_TRANSONIC_GATES", "0") == "1"
ALPHA = 3.06
RAMP = dict(m_target=0.84, alpha_deg=ALPHA, farfield="neumann",
            n_seed=40, n_newton_max=80, tol_residual=1e-10)

# The anchors. gamma/M_max as measured post-B21 (provenance in the docstring).
ANCHORS = {
    "coarse": dict(gamma=0.08493098, m_max=1.3684),
    "medium": dict(gamma=0.088343, m_max=2.4818),
}
GAMMA_RTOL, MMAX_RTOL = 1e-4, 1e-3


def _ramp(level):
    path = M6_DIR / f"{level}.msh"
    if not path.exists():
        pytest.skip(f"{path} not generated (gitignored)")
    mesh = read_mesh(path)
    a = np.radians(ALPHA)
    wls = WakeLevelSet(
        np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]]),
        direction=(np.cos(a), np.sin(a), 0.0))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)
    r = solve_multivalued_newton_transonic(mvop=mvop, mesh=mesh, **RAMP,
                                           **B_NEWTON_M6_DEFAULTS)
    return r


def _assert_anchors(r, level):
    a = ANCHORS[level]
    assert r["target_reached"], (
        f"M6 {level} ramp no longer reaches M0.84 (m_final={r['m_final']}) — "
        "a capability re-baseline; see this file's docstring")
    assert all(l["converged"] for l in r["levels"]), (
        f"levels {[int(l['converged']) for l in r['levels']]}")
    assert r["residual_history"][-1] < 1e-9
    assert r["n_limited"] == 0, f"n_limited={r['n_limited']}"
    assert r["n_floored"] <= 1, f"n_floored={r['n_floored']}"
    assert np.isclose(r["gamma"], a["gamma"], rtol=GAMMA_RTOL, atol=0.0), (
        f"gamma {r['gamma']:.8f} vs anchor {a['gamma']:.8f}")
    m_max = float(np.sqrt(r["mach2_max"]))
    assert np.isclose(m_max, a["m_max"], rtol=MMAX_RTOL, atol=0.0), (
        f"M_max {m_max:.4f} vs anchor {a['m_max']:.4f}")


@pytest.mark.skipif(not GATES, reason="heavy gated 3-D anchor (~35 s)")
def test_m6_coarse_ramp_anchor():
    """The committed M6 COARSE M0.84 ramp anchors (post-B20/B21 state)."""
    _assert_anchors(_ramp("coarse"), "coarse")


@pytest.mark.skipif(not GATES, reason="heavy gated 3-D anchor (~9 min)")
def test_m6_medium_ramp_anchor():
    """The committed M6 MEDIUM M0.84 ramp anchors — the number that moved
    silently through two re-baselines before this lock existed."""
    _assert_anchors(_ramp("medium"), "medium")
