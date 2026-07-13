"""
Track B / B6 gate: transonic + Mach continuation on the level-set path
(docs/roadmap.md Track B B6; docs/design_track_b.md §5.2/D10).

What B6 adds over B3/B4 (subsonic lifting):

  * PER-SIDE artificial density on the cut elements -- a cut element has two
    velocity states (upper/lower DOF copies), so rho_tilde is evaluated once
    per side, and the upstream walk runs on a face-adjacency graph RESTRICTED
    TO THAT SIDE's own element set (the wake is a slip line: density
    information must not cross the tangential discontinuity).
    `MultivaluedOperator.element_rho_tilde` / `_side_upwind`.
  * LOCALIZED (supersonic-zone) damping. ★ The conforming P4 stabilizer
    (theta*diag on the whole system) does NOT transfer: diagonal damping is a
    Jacobi smoother, near-transparent to smooth global modes -- and on the B
    path the implicit Kutta makes the circulation a SOLUTION MODE, so global
    damping throttles Gamma itself (measured: Gamma crawls 0.0005 -> 0.017 in
    160 outers on a case that converges undamped in 35). B6 damps only the
    rows of nu > 0 elements (the pocket + shock), where the transonic
    instability actually lives (`damping_scope="supersonic"`).
  * `solve_multivalued_transonic`: the Mach ramp with NO Gamma secant -- a
    level is just a warm-started Picard solve (the Track B payoff: the P5
    st133-class per-station secant failure cannot occur structurally).

Gate anchors (roadmap B6, re-anchored 2026-07-11 P4-erratum aware): NACA
coarse M0.80 alpha 1.25 against the G8.1 Newton locks (shock 0.658 /
cl_p 0.459) as a Picard-QUALITY comparison; same-mesh same-recipe
Picard-vs-Picard A/B; dual-mesh rule (M0 embedded + M3 wake-free).

Cheap suite tests here run the SUBCRITICAL no-op clauses + a short M0.70
supercritical smoke on coarse; the M0.80 A/B numbers are asserted in the
gated test (PYFP3D_TRANSONIC_GATES=1) and recorded by the demo
(cases/demo/b6_transonic/).
"""

import os
from pathlib import Path

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.post.shock import shock_metrics
from pyfp3d.post.surface_ls import (
    cl_pressure_levelset,
    surface_curve_levelset,
    wall_cp_levelset,
)
from pyfp3d.solve.picard_ls import (
    solve_multivalued_lifting,
    solve_multivalued_transonic,
)
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

REPO_ROOT = Path(__file__).parent.parent
MESHES = REPO_ROOT / "cases" / "meshes"
M0_DIR = MESHES / "naca0012_2.5d"                 # wake-embedded
M3_DIR = MESHES / "naca0012_wakefree_2.5d"        # wake-free (workflow form)

GATES = os.environ.get("PYFP3D_TRANSONIC_GATES", "0") == "1"


def _load(directory, level):
    path = directory / f"{level}.msh"
    if not path.exists():
        pytest.skip(f"{path} not generated (gitignored)")
    return read_mesh(path)


def _build(mesh):
    z = mesh.nodes[:, 2]
    wls = WakeLevelSet(
        np.array([[1.0, 0.0, z.min()], [1.0, 0.0, z.max()]]),
        direction=(1.0, 0.0, 0.0),
    )
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    return MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)


# ---------------------------------------------------------------------------
# subcritical no-op: upwind machinery in the loop must not change B3/B4
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("mesh_dir", [M0_DIR, M3_DIR],
                         ids=["m0-embedded", "m3-wakefree"])
def test_subcritical_noop(mesh_dir):
    """upwind_c > 0 at a subcritical Mach reproduces the B3/B4 solve exactly
    in Gamma (nu == 0 everywhere -> rho_tilde == rho; the assembled matrices
    are identical, so any difference is pure solver noise ~1e-9)."""
    mesh = _load(mesh_dir, "coarse")
    mvop = _build(mesh)
    r_off = solve_multivalued_lifting(mvop, mesh, 0.5, alpha_deg=2.0)
    r_on = solve_multivalued_lifting(mvop, mesh, 0.5, alpha_deg=2.0,
                                     upwind_c=1.5, m_crit=0.95)
    assert r_on["nu_max"] == 0.0
    assert r_on["n_nu_active"] == 0
    assert r_on["n_limited"] == 0 and r_on["n_floored"] == 0
    assert abs(r_on["gamma"] - r_off["gamma"]) < 1e-7
    assert r_on["converged"]


def test_per_side_density_identity_subcritical():
    """At a fixed subcritical state the per-side ARTIFICIAL density equals the
    per-side isentropic density bitwise on each side's own element set
    (nu == 0 -> rho_tilde == rho exactly, the G4.2 property per side)."""
    mesh = _load(M0_DIR, "coarse")
    mvop = _build(mesh)
    r = solve_multivalued_lifting(mvop, mesh, 0.5, alpha_deg=2.0)
    phi = r["phi_ext"]
    ru, rl = mvop.element_densities(phi, 0.5)
    tu, tl = mvop.element_rho_tilde(phi, 0.5, 1.5, 0.95)
    in_upper, in_lower = mvop._side_element_sets()
    assert mvop.nu_max == 0.0
    np.testing.assert_array_equal(ru[in_upper], tu[in_upper])
    np.testing.assert_array_equal(rl[in_lower], tl[in_lower])


def test_side_upwind_graph_stays_on_side():
    """The same-side restriction (design_track_b.md §5.2/D10): the upper
    upwind graph never references a lower-only element and vice versa -- the
    wake is a slip line, density information does not cross it."""
    mesh = _load(M0_DIR, "coarse")
    mvop = _build(mesh)
    in_upper, in_lower = mvop._side_element_sets()
    upw_u, upw_l = mvop._side_upwind()
    for upw, keep in ((upw_u, in_upper), (upw_l, in_lower)):
        fn = upw.face_neighbors
        nb = fn[fn >= 0]
        assert keep[nb].all(), "upwind graph leaks across the wake"
        # rows of excluded elements must have no neighbors at all
        assert (fn[~keep] == -1).all()


def test_own_side_field_split_matches_assembly():
    """B6 fix: own_side_field's element split must be the assembly's
    (`_side_element_sets`), not a node-side test -- a te_lower element contains
    the eps-shifted '+' TE node, and the node test mislabeled it UPPER while
    the matrix weights it with rho_LOWER (measured M_max junk 2.56 vs 0.85)."""
    mesh = _load(M0_DIR, "coarse")
    mvop = _build(mesh)
    in_upper, _ = mvop._side_element_sets()
    up = np.ones(mvop.op.n_tets)
    lo = np.zeros(mvop.op.n_tets)
    own = mvop.own_side_field(up, lo)
    tel = np.asarray(mvop.cm.te_lower_elems, dtype=np.int64)
    assert (own[tel] == 0.0).all(), "te_lower elements must read the LOWER field"
    np.testing.assert_array_equal(own == 1.0, in_upper)


# ---------------------------------------------------------------------------
# supercritical smoke (cheap): M0.70 coarse -- a real shock-free supersonic
# pocket, converged, physical
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("mesh_dir", [M0_DIR, M3_DIR],
                         ids=["m0-embedded", "m3-wakefree"])
def test_supercritical_pocket_m070(mesh_dir):
    """M0.70 alpha1.25: supersonic pocket appears (M_max ~ 1.07), the per-side
    upwind + localized damping converge it cleanly (0 limited/floored), and
    the pocket sits on the wing -- not on the wake sheet."""
    mesh = _load(mesh_dir, "coarse")
    mvop = _build(mesh)
    r = solve_multivalued_transonic(mvop, mesh, 0.70, alpha_deg=1.25,
                                    n_outer_level=400)
    L = r["levels"][-1]
    assert L["converged"], f"M0.70 level did not converge: {L}"
    # element_mach2 default flipped side->main 2026-07-14: bit-identical on
    # the quasi-2D families (no beyond-tip cells touch aux dofs; re-read
    # cases/demo/b8_tip_taper_ls/run_b8_mmax_reread.py). Band unchanged.
    assert 1.0 < L["mach_max"] < 1.3
    assert L["n_limited"] == 0 and L["n_floored"] == 0
    # the supersonic pocket lives on the wing, not the cut strip
    ctr = mesh.nodes[mesh.elements].mean(axis=1)
    act = mvop.nu_active_elements()
    if len(act):
        is_cut = np.zeros(len(ctr), dtype=bool)
        is_cut[mvop.cm.cut_elems] = True
        assert not is_cut[act].any(), "nu active on the wake cut strip"


def test_global_damping_throttles_gamma_negative_result():
    """The recorded NEGATIVE result that shaped B6 (kept honest): GLOBAL fluid
    damping (the P4 stabilizer transplanted as-is) throttles the implicit-Kutta
    circulation mode -- after 60 outers Gamma is still < 25% of the value the
    undamped/localized path reaches, because diag damping is a Jacobi smoother
    and Gamma is a smooth global SOLUTION mode here (on the conforming path it
    is an outer secant unknown, outside the damped matrix)."""
    mesh = _load(M0_DIR, "coarse")
    mvop = _build(mesh)
    kw = dict(alpha_deg=1.25, upwind_c=1.5, n_outer_max=60)
    r_glob = solve_multivalued_lifting(mvop, mesh, 0.60, damping_theta=0.2,
                                       damping_scope="fluid", **kw)
    r_none = solve_multivalued_lifting(mvop, mesh, 0.60, **kw)
    assert r_none["gamma"] > 0.09          # undamped: essentially converged
    assert r_glob["gamma"] < 0.25 * r_none["gamma"]


# ---------------------------------------------------------------------------
# gated: the M0.80 B6 gate numbers (heavy, PYFP3D_TRANSONIC_GATES=1)
# ---------------------------------------------------------------------------

def test_live_vortex_feedback_runaway_negative_result():
    """The recorded B6 fold finding (kept honest, cheap version): with the
    LIVE Gamma -> far-field-vortex refresh (option a) at M0.80 coarse, Gamma
    climbs monotonically PAST the Newton solution (~0.2295) instead of
    settling -- the loop gain is > 1 near the fold, so the level never
    converges and the trajectory is monotone through the physical value.
    (The full runaway to M_max ~37 is the demo's job; here we only pin the
    monotone-past-the-solution character in a short budget.)"""
    mesh = _load(M0_DIR, "coarse")
    mvop = _build(mesh)
    r = solve_multivalued_transonic(mvop, mesh, 0.80, alpha_deg=1.25,
                                    n_outer_level=700)
    g = np.asarray(r["gamma_history"], dtype=np.float64)
    assert not r["levels"][-1]["converged"]
    assert g.max() > 0.26, "expected the runaway past the Newton value"


@pytest.mark.skipif(not GATES, reason="PYFP3D_TRANSONIC_GATES != 1")
@pytest.mark.parametrize("mesh_dir", [M0_DIR, M3_DIR],
                         ids=["m0-embedded", "m3-wakefree"])
def test_g_b6_m080_coarse_neumann(mesh_dir):
    """B6 gate, coarse M0.80 alpha1.25 (roadmap Track B B6): the transonic
    B-path recipe is farfield="neumann" (the Lopez outlet -- no Gamma
    feedback; the live option-a loop has gain > 1 near the fold and the
    lagged variant's outer map has no fixed point, design_track_b.md
    section 10.3). Clauses: physical field (0 limited/floored, M_max < 1.6),
    shock inside the Picard-quality band of the G8.1 Newton lock (0.658
    +/- 0.06; the conforming Picard's own stall sits at 0.604), per-mesh
    locks only (fold discipline)."""
    mesh = _load(mesh_dir, "coarse")
    mvop = _build(mesh)
    r = solve_multivalued_transonic(mvop, mesh, 0.80, alpha_deg=1.25,
                                    farfield="neumann", n_outer_level=3000)
    L = r["levels"][-1]
    assert L["n_limited"] == 0 and L["n_floored"] == 0
    # element_mach2 default flipped side->main 2026-07-14: re-read on the
    # cached M0.80 states is bit-identical (M0 1.3916 / M3 1.3851 both
    # readings -- quasi-2D has no beyond-tip cells). Band unchanged.
    assert L["mach_max"] < 1.6
    # Gamma within the Picard-quality band of the Newton solution (~0.2295);
    # the -4% B5 outlet truncation is inside this band by construction.
    assert abs(r["gamma"] / 0.2295 - 1.0) < 0.10
    cp = wall_cp_levelset(mesh, mvop, r["phi_ext"], m_inf=0.80)
    xs, cps = surface_curve_levelset(cp, "upper")
    sh = shock_metrics(xs, cps, 0.80)
    assert sh["has_shock"]
    assert abs(sh["x_shock"] - 0.658) < 0.06
