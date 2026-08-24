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
#: ★★ GATE RE-SPEC 2026-08-05, user ruling. Recorded in full in
#: phases/p2/docs/dev_phase_two/20260805-2330-m1a-respec.md, written BEFORE this edit as
#: roadmap sec 5 requires. Summary of what changed and what it COST:
#:
#: The ladder moved from (coarse, medium, fine) to (xcoarse, coarse, medium) because
#: `fine` was retired from all tests and demos on 2026-08-04. Measured consequence:
#: NONE of the three original assertions survives the move, because the two windows
#: are in DIFFERENT convergence regimes on this clean factor-2 ladder
#: (h_wall 0.040 / 0.020 / 0.010 / 0.005):
#:
#:     window                    d2/cl_m    ratio    observed p    Richardson
#:     coarse/medium/fine        +0.629 %   0.1556   2.68          0.255123
#:     xcoarse/coarse/medium     +4.214 %   0.5056   0.98          0.263707
#:
#: so the retired window's extrapolation matched the committed 0.25511 to five
#: digits while the new window overshoots by +3.4 %. The 1 %-level agreement and the
#: Richardson assertion are therefore RETIRED, not loosened -- the criterion below
#: measures different QUANTITIES rather than the same ones at a weaker threshold.
#:
#: The retired window's readings are kept as a historical MEASUREMENT (not a live
#: gate) so "this solver reached p ~ 2.7 over h 0.02 -> 0.005" is not lost with the
#: level: cl 0.242997 / 0.253237 / 0.254830 at coarse/medium/fine, seed 5.
LOCK = {
    "xcoarse": dict(cl=0.222742, m_max=1.0972),
    "coarse": dict(cl=0.242984, m_max=1.1356),
    "medium": dict(cl=0.253237, m_max=1.1537),
}
LADDER = ("xcoarse", "coarse", "medium")
CL_RTOL = 2.0e-3          # 0.2 %: well inside run-to-run scatter, far below
#                           the physics changes this is meant to catch
#: the re-spec's three criteria (sec 3 of the round file), with the out-of-envelope
#: readings that show each one DISCRIMINATES -- at M0.80 the same three-level sequence
#: gives d2 = -0.0666 (non-monotone), ratio -0.6003 and |d2/cl_m| 16.33 %.
#: ★ Note which criterion does the work: at M0.80 all three levels CONVERGE with zero
#: clamps (the 2026-08-05 cold-start seed fallback made that so), therefore
#: "everything converged" is NOT the discriminator -- MONOTONICITY is.
RATIO_MAX = 0.7           # measured 0.5065 in envelope, -0.6003 out
D2_REL_MAX = 0.05         # measured 4.22 % in envelope, 16.33 % out


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
    from tests.conftest import REPO_ROOT
    return REPO_ROOT / "cases" / "meshes"


@pytest.mark.parametrize("level", LADDER)
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


def test_m1a_three_level_convergence(mesh_dir):
    """★ The re-spec'd M1a criterion (2026-08-05): on (xcoarse, coarse, medium) the
    sequence must be MONOTONE, CONTRACTING, and settled to under 5 % on the last step.

    Each of the three discriminates out of envelope -- at M0.80 the same sequence gives
    d2 = -0.0666, ratio -0.6003 and 16.33 % -- while "all three converged" does NOT,
    since the cold-start seed fallback makes M0.80 converge too. So monotonicity is
    what carries the statement, and it is asserted first.

    What this criterion no longer claims, and why, is in the LOCK comment above and in
    phases/p2/docs/dev_phase_two/20260805-2330-m1a-respec.md sec 5.
    """
    cl = {}
    for level in LADDER:
        r, c = _solve(mesh_dir, level)
        assert r["converged"] and r["n_limited"] == 0 and r["n_floored"] == 0, (
            f"{level}: the sequence needs three genuine solutions")
        cl[level] = c
    d1 = cl["coarse"] - cl["xcoarse"]
    d2 = cl["medium"] - cl["coarse"]
    assert d1 > 0.0 and d2 > 0.0, (
        f"NOT monotone: d1 = {d1:+.6f}, d2 = {d2:+.6f}. This is the criterion that "
        f"separates in-envelope from out (M0.80 gives d2 = -0.0666)")
    ratio = d2 / d1
    assert 0.0 < ratio < RATIO_MAX, (
        f"convergence ratio {ratio:.4f} outside (0, {RATIO_MAX}) "
        f"(0.5065 measured in envelope; -0.6003 at M0.80)")
    rel = abs(d2 / cl["medium"])
    assert rel < D2_REL_MAX, (
        f"last-step change {100 * rel:.3f} % >= {100 * D2_REL_MAX:.0f} % "
        f"(4.22 % measured in envelope; 16.33 % at M0.80)")
