"""M1a: the in-envelope mesh-convergence lock (phase two GS1.5a, decision D3).

What this locks and why it exists
---------------------------------
GS1.3b measured that the h-convergence of transonic LIFT degrades monotonically
with shock strength (docs/dev_phase_two/20260728-2133-s1-closure-stiffness-and-
envelope.md §4.3): over three mesh levels (h_wall 0.02 / 0.01 / 0.005 chord) the
successive cl changes and their ratio are

    M 0.60  (M_max 0.85)   +2.84 % -> +0.23 %   ratio 0.08
    M 0.72  (M_max 1.17)   +4.35 % -> +0.91 %   ratio 0.21
    M 0.75  (M_max 1.27)   +8.78 % -> +3.14 %   ratio 0.36
    M 0.7875(M_max 1.60)  +40.5 % -> +36.6 %    ratio 0.90   <- divergent

so the h-convergent envelope is M_max ~< 1.2. Product metric M1 (M0.80/alpha
1.25, M_max ~ 1.4) sits OUTSIDE it and is kept as a recorded FAIL; M1a locks
what the solver CAN do today, inside the envelope, so it cannot silently rot
while GS1.6 (entropy-fix full potential) works on widening the envelope.

Scope discipline
----------------
This is a MESH-CONVERGENCE lock, not an accuracy gate: the repo has no
committed transonic reference for M0.72/alpha1.25 (the reference data covers
M0.5 and the M0.80 shock position only). The Richardson extrapolation 0.2563 is
recorded for context and the accuracy anchor is an open item -- it needs an
external Euler/experimental reference, not a number of ours.

Assertions use relative tolerances (phase-two decision D1); the committed
values are drift references, NOT truth.
"""

import os

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.post.surface import wall_force_coefficients
from pyfp3d.solve.newton import solve_newton_lifting

M_INF, ALPHA = 0.72, 1.25

#: measured 2026-07-28 on DESKTOP-N6UP769 (16 threads), target-Mach direct
#: Newton (no Mach ramp -- GS1.2b), upwind_c 1.5, freeze_tol 1e-6.
#: ★ RE-ANCHORED 2026-07-31 (GS1b.11): the entropy correction became the DEFAULT, so
#: these are entropy-ON values. Measured at the runner-default 16 threads with this
#: file's own recipe (solve_newton_lifting direct at M0.72, upwind_c 1.5, m_crit 0.95,
#: freeze_tol 1e-6, precond direct).
#:
#: ANCHORED TO WHAT: this is a DRIFT LOCK, not a correctness claim -- M0.72 has no
#: external reference in this project. What carries external meaning is the
#: three-level CONSISTENCY criterion below (medium -> fine < 1 %), and it IMPROVED:
#: +0.91 % isentropic -> +0.63 % with the correction.
#:
#: SUPERSEDED isentropic values, kept per discipline #11 rather than overwritten:
#:     coarse cl 0.242797  x_shock 0.28742  m_max 1.1359
#:     medium cl 0.253351  x_shock 0.29150  m_max 1.1540
#:     fine   cl 0.255662  x_shock 0.28753  m_max 1.1658   RICHARDSON_CL 0.25631
#: The moves are +0.08 % / -0.04 % / -0.28 % in cl; the fine one exceeds CL_RTOL,
#: which is why this re-anchor was necessary (and why GS1b.6's "flipping the default
#: would not break the M1a lock" was imprecise -- it held for the consistency
#: criterion, not for the absolute locks).
LOCK = {
    "coarse": dict(cl=0.242997, x_shock=0.28706, m_max=1.1356),
    "medium": dict(cl=0.253237, x_shock=0.29075, m_max=1.1537),
    "fine": dict(cl=0.254823, x_shock=0.28358, m_max=1.1650),
}
CL_RTOL = 2.0e-3          # 0.2 %: well inside run-to-run scatter, far below
#                           the physics changes this is meant to catch
RICHARDSON_CL = 0.25511   # was 0.25631 isentropic; convergence ratio 6.45

run_fine = pytest.mark.skipif(
    os.environ.get("PYFP3D_TRANSONIC_GATES", "0") != "1",
    reason="the fine level is ~60 s; set PYFP3D_TRANSONIC_GATES=1")


def _solve(mesh_dir, level):
    path = mesh_dir / "naca0012_2.5d" / f"{level}.msh"
    if not path.exists():
        pytest.skip(f"{path.name} not generated "
                    "(cases/meshes/naca0012_2.5d/generate_naca0012.py)")
    mc, wc = cut_wake(read_mesh(path))
    r = solve_newton_lifting(mc, wc, m_inf=M_INF, alpha_deg=ALPHA,
                             upwind_c=1.5, m_crit=0.95, freeze_tol=1e-6,
                             freeze_refresh_max=8, precond="direct",
                             direct_refactor_every=4, n_newton_max=80)
    dz = float(np.ptp(mc.nodes[:, 2]))
    f = wall_force_coefficients(mc.nodes, mc.elements,
                                mc.boundary_faces["wall"], r["phi"],
                                alpha_deg=ALPHA, s_ref=dz, m_inf=M_INF)
    return r, float(f["cl"])


@pytest.fixture(scope="module")
def mesh_dir():
    from .conftest import REPO_ROOT
    return REPO_ROOT / "cases" / "meshes"


@pytest.mark.parametrize("level", ["coarse", "medium"])
def test_m1a_level_lock(mesh_dir, level):
    """Each in-envelope level reproduces its committed cl, converges to a
    genuine solution, and carries no clamped cells."""
    r, cl = _solve(mesh_dir, level)
    assert r["converged"], (
        f"{level} did not converge: |R| = {r['residual_history'][-1]:.2e}")
    assert r["n_limited"] == 0 and r["n_floored"] == 0, (
        f"{level} carries clamps ({r['n_limited']} limited / "
        f"{r['n_floored']} floored) -- a clamped state is not a solution")
    ref = LOCK[level]["cl"]
    assert abs(cl - ref) <= CL_RTOL * abs(ref), (
        f"{level} cl {cl:.6f} vs committed {ref:.6f} "
        f"({100 * (cl - ref) / ref:+.3f} %)")
    assert float(np.sqrt(r["mach2_max"])) < 1.25, (
        "M_max moved out of the in-envelope range this lock is about")


def test_m1a_coarse_to_medium_consistency(mesh_dir):
    """The cheap half of the envelope statement: coarse -> medium must move cl
    by less than 5 % (measured +4.35 %). Outside the envelope this number is
    +40 %, so it is a real discriminator even without the fine level."""
    _, cl_c = _solve(mesh_dir, "coarse")
    _, cl_m = _solve(mesh_dir, "medium")
    d = (cl_m - cl_c) / cl_c
    assert 0.0 < d < 0.05, (
        f"coarse->medium cl change {100 * d:+.2f} % outside (0, 5) % "
        f"(measured +4.35 % on 2026-07-28)")


@run_fine
def test_m1a_medium_to_fine_convergence(mesh_dir):
    """The binding M1a criterion: medium -> fine cl change < 1 %, and the
    convergence RATIO well below 1 (0.21 measured; >= 1 is divergence)."""
    _, cl_c = _solve(mesh_dir, "coarse")
    _, cl_m = _solve(mesh_dir, "medium")
    r_f, cl_f = _solve(mesh_dir, "fine")
    d1, d2 = cl_m - cl_c, cl_f - cl_m
    assert abs(d2 / cl_m) < 0.01, (
        f"medium->fine cl change {100 * d2 / cl_m:+.3f} % >= 1 %")
    ratio = d2 / d1
    assert 0.0 < ratio < 0.4, (
        f"convergence ratio {ratio:.3f} outside (0, 0.4) "
        "(0.21 measured; >= 1 means divergence)")
    rich = cl_f + d2 * ratio / (1.0 - ratio)
    assert abs(rich - RICHARDSON_CL) < 0.01 * RICHARDSON_CL, (
        f"Richardson extrapolation {rich:.5f} vs committed "
        f"{RICHARDSON_CL:.5f}")
    assert r_f["n_limited"] == 0 and r_f["n_floored"] == 0
