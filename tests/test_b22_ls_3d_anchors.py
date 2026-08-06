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

#: ★★ RE-SPEC 2026-08-06. Full record in docs/dev_phase_two/20260806-1200-b22-respec.md,
#: written before this edit as roadmap sec 5 requires. What happened and what it costs:
#:
#: `cases/meshes/onera_m6_wakefree/` was regenerated ROUND on 2026-08-04 while these
#: anchors were measured on the FLAT cap, and the 2026-08-06 gated run failed both.
#: Measured on the round meshes with the recipe otherwise untouched:
#:
#:   coarse   0.60/0.65/0.70/0.75/0.80/0.82 all converge with ZERO clamps; 0.84 fails
#:            twice at 0/2. Highest CLEAN level = 0.82 (was 0.84 on flat).
#:   medium   EVERY level is clamped, starting at M0.60 with 2/1, and M_max reads 3.31
#:            at M0.60 rising to 15.26 at M0.70 -- the round cap RESOLVES the P13 tip
#:            free-edge singularity. Highest converged = 0.6625 (4/1 clamps), and there
#:            is NO level with n_limited == 0 at all.
#:
#: ★ Independently corroborated: LE-15 measured the same ramp by a different route and
#: got the SAME per-level clamp counts -- 0.65 -> 3/1, 0.6625 -> 4/1, 0.675 -> 5/1,
#: 0.70 -> 9/2 -- and the same convergence boundary. Not a one-off.
#:
#: ★ And unlike G8.2, this CANNOT be bought back with the production taper:
#: solve_multivalued_newton_transonic has no tip_taper parameter at all (the taper
#: scales the CONFORMING path's Kutta rows), and B31 measured the LS-side tip cure as
#: closed-negative. So the round tip's cost on the level-set path is unrecoverable --
#: a capability fork between the two paths.
#:
#: SUPERSEDED flat-cap anchors, kept per discipline #11 rather than deleted:
#:     coarse gamma 0.08493098  M_max 1.3684   at M0.84
#:     medium gamma 0.088343    M_max 2.4818   at M0.84
ANCHORS = {
    #: the highest level that converges CLEAN (zero clamps), which is what GS1.4 lets
    #: us call a solution -- not merely the highest that "converged".
    "coarse": dict(m_clean=0.82, gamma=0.08138231, m_max=2.43637),
}
#: measured capability readings for medium, RECORDED not asserted (see the xfail below)
MEDIUM_RECORDED = dict(m_highest_converged=0.6625, gamma=0.07172634,
                       m_max=5.62615, n_limited=4, n_floored=1,
                       m_max_at_m060=3.31320, m_max_at_m070=15.25894)
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


def _clean_levels(r):
    """Levels that converged with ZERO clamps. GS1.4: a clamped state is not a
    solution, so "the highest level that converged" is the wrong thing to anchor --
    on the round medium mesh EVERY level converges-with-clamps and none is clean."""
    return [l for l in r["levels"]
            if l["converged"] and not l["n_limited"] and not l["n_floored"]]


@pytest.mark.skipif(not GATES, reason="heavy gated 3-D anchor (~3 min)")
def test_m6_coarse_ramp_anchor():
    """The M6 COARSE ramp's highest CLEAN level and its state (round-tip re-spec)."""
    r = _ramp("coarse")
    a = ANCHORS["coarse"]
    clean = _clean_levels(r)
    assert clean, "no level converged with zero clamps -- there is nothing to anchor"
    top = max(clean, key=lambda l: l["m_inf"])
    assert abs(top["m_inf"] - a["m_clean"]) < 1e-9, (
        f"highest CLEAN level {top['m_inf']:.4f} vs anchor {a['m_clean']} — a "
        f"capability re-baseline; see docs/dev_phase_two/20260806-1200-b22-respec.md")
    assert top["residual_norm"] < 1e-9, f"|R| = {top['residual_norm']:.2e}"
    assert np.isclose(top["gamma"], a["gamma"], rtol=GAMMA_RTOL, atol=0.0), (
        f"gamma {top['gamma']:.8f} vs anchor {a['gamma']:.8f}")
    assert np.isclose(top["mach_max"], a["m_max"], rtol=MMAX_RTOL, atol=0.0), (
        f"M_max {top['mach_max']:.5f} vs anchor {a['m_max']:.5f}")


@pytest.mark.xfail(strict=True, reason=(
    "RE-SPEC 2026-08-06, measured: on the ROUND-tip wakefree medium mesh this ramp has "
    "NO level with zero clamps -- M0.60 already carries 2 limited / 1 floored and M_max "
    "3.31, rising to 15.26 at M0.70, because the round cap resolves the P13 tip "
    "free-edge singularity. GS1.4 says a clamped state is not a solution, so there is "
    "nothing here to anchor, and anchoring a clamped state would void that contract "
    "inside this very test. STRICT on purpose: if a future change makes a clean level "
    "appear, this must go RED so the capability gets re-anchored rather than silently "
    "improving. Readings are in MEDIUM_RECORDED and in "
    "docs/dev_phase_two/20260806-1200-b22-respec.md."))
@pytest.mark.skipif(not GATES, reason="heavy gated 3-D anchor (~35 min)")
def test_m6_medium_ramp_anchor():
    """The M6 MEDIUM ramp — expected to have no clean level on the round tip."""
    r = _ramp("medium")
    assert _clean_levels(r), (
        "a clean level appeared on the round medium mesh — re-anchor this test")
