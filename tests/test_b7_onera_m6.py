"""
Track B / B7 gate: the ONERA M6 3D gate for the level-set solver
(docs/roadmap.md Track B B7; docs/design_track_b.md).

B7 is the first Track B phase that exercises the genuinely 3D machinery. B1-B6
ran on quasi-2D meshes, where the wake sheet is a flat strip with no tip and no
sweep, so three pieces of the level-set path were never tested at all:

  * the TE-**polyline** ruled level set (D9) -- the sheet is ruled between a
    swept root->tip TE polyline, so its per-segment (v, d_hat, n_hat) frame is
    OBLIQUE (the span axis is not perpendicular to the wake direction). B1 found
    and fixed a real defect here (an orthogonal projection wrongly clipped ~60%
    of the M6 cut set).
  * the **spanwise clip** 0 <= q <= span_length, which is what makes Gamma go to
    zero at the tip discretely: elements the wake PLANE crosses outboard of the
    sheet are rejected (`CutElementMap.beyond_tip_elems`), so no jump is carried
    there. This is the level-set analogue of the conforming path's free-edge
    rule (mesh/wake_cut.py, M1).
  * the g2 wake-BC component with a FREE spanwise jump gradient -- the trailing
    vortex DOF (D1). On a quasi-2D mesh there is no spanwise direction to be
    free in.

Far field = "neumann" (roadmap B7 / design_track_b.md): the B-path vortex
(`picard_ls._farfield_main`) is a SPAN-UNIFORM 2D point vortex whose branch cut
is the ray y=0, x>0 at every z. On a 3D wing that is wrong in two independent
ways, both measured on M6 coarse at M0.5 (demo cases/demo/b7_onera_m6, and
`test_farfield_vortex_is_contraindicated_in_3d` below):
  (a) the alpha-aimed sheet is NOT coplanar with that cut (by the outlet the
      sheet has climbed to y ~ x*tan(alpha)), so the outlet carries a prescribed
      Gamma jump no cut supports -- B3's recorded coplanarity rule, in 3D;
  (b) even re-aimed coplanar, a span-uniform Gamma cannot match Gamma(z) (which
      decays to 0 at the tip) -- P5's branch-ray artifact, whose conforming fix
      was the Gamma(z) taper (farfield_dirichlet(spanwise_gamma=True)).
Neumann carries no vortex, so it has nothing to misapply: neither defect exists.
Cost: the B5-measured O(Gamma/R) outlet truncation (a few % of lift).

Runtime policy follows P5: only fast, solve-free checks of the new 3D level-set
post-processing run always-on; the M6 transonic solves are behind
PYFP3D_TRANSONIC_GATES=1 (each is minutes), with the committed evidence in
cases/demo/b7_onera_m6/. Both M6 mesh families are gitignored, so these tests
skip until the generators are run:
    python cases/meshes/onera_m6/generate_onera_m6.py --level coarse
    python cases/meshes/onera_m6_wakefree/generate_onera_m6_wakefree.py
"""

import os
from pathlib import Path

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.meshgen.wing3d import B_SEMI, x_te
from pyfp3d.post.section_cut import section_cp_curve
from pyfp3d.post.shock import shock_report
from pyfp3d.post.surface import cl_kj_3d, planform_area
from pyfp3d.post.surface_ls import (
    cl_pressure_3d_levelset,
    section_cp_curve_levelset,
    wall_cp_levelset,
)
from pyfp3d.solve.picard_ls import (
    solve_multivalued_lifting,
    solve_multivalued_transonic,
)
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

REPO_ROOT = Path(__file__).parent.parent
MESHES = REPO_ROOT / "cases" / "meshes"
M1_DIR = MESHES / "onera_m6"            # wake-embedded (A/B vs conforming)
M4_DIR = MESHES / "onera_m6_wakefree"   # wake-free (the workflow target)

M_INF = 0.84
ALPHA = 3.06
ETAS = (0.44, 0.65, 0.90)               # AGARD AR-138 gate stations

GATES = os.environ.get("PYFP3D_TRANSONIC_GATES", "0") == "1"
run_gates = pytest.mark.skipif(
    not GATES,
    reason="M6 3D level-set solves take minutes; set PYFP3D_TRANSONIC_GATES=1",
)


def _require(directory: Path, level: str = "coarse") -> Path:
    p = directory / f"{level}.msh"
    if not p.exists():
        pytest.skip(f"{directory.name}/{level}.msh not generated (gitignored); "
                    "see this module's header")
    return p


def _m6_levelset(alpha_deg: float = ALPHA) -> WakeLevelSet:
    """The swept ruled sheet: TE polyline root->tip, aimed at the incidence."""
    a = np.radians(alpha_deg)
    te = np.array([[x_te(0.0), 0.0, 0.0], [x_te(B_SEMI), 0.0, B_SEMI]])
    return WakeLevelSet(te, direction=(np.cos(a), np.sin(a), 0.0))


def _setup(path: Path, alpha_deg: float = ALPHA):
    mesh = read_mesh(path)
    wls = _m6_levelset(alpha_deg)
    cm = CutElementMap(mesh.nodes, mesh.elements, wls,
                       wall_nodes=np.unique(mesh.boundary_faces["wall"]))
    mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)
    return mesh, cm, mvop


def _gamma_of_z(mesh, cm, phi_ext, mvop):
    """(z, Gamma) per TE station, sorted root -> tip."""
    z = mesh.nodes[cm.te_nodes, 2]
    g = mvop.te_jump(phi_ext)
    o = np.argsort(z)
    return z[o], g[o]


# ---------------------------------------------------------------------------
# Always-on: the 3D post-processing the B path needs (solve-free).
# ---------------------------------------------------------------------------

def test_cl_kj_3d_elliptic_sanity():
    """cl_kj_3d on an analytic elliptic Gamma(z) reproduces the closed form
    CL = 2 * (pi/4) * G0 * b_semi / s_ref. Pins the B7 lift reduction to the
    same helper the conforming P5 path uses (no Track-B copy)."""
    b, g0, s_ref = 1.1963, 0.1, 0.7532
    z = np.linspace(0.0, b, 60)
    g = g0 * np.sqrt(np.clip(1.0 - (z / b) ** 2, 0.0, None))
    cl = cl_kj_3d(g, z, s_ref=s_ref, b_semi=b)
    exact = 2.0 * (np.pi / 4.0) * g0 * b / s_ref
    assert abs(cl - exact) / exact < 5e-3, (cl, exact)


def test_te_stations_are_index_aligned_with_te_nodes():
    """te_jump is one value per TE node, index-aligned with cm.te_nodes -- the
    contract Gamma(z) and cl_kj_3d rely on. Guards a silent mis-pairing that
    would corrupt the spanwise integral without ever raising."""
    mesh, cm, mvop = _setup(_require(M1_DIR))
    phi = np.arange(mvop.n_total, dtype=np.float64)   # arbitrary state
    g = mvop.te_jump(phi)
    assert g.shape == (len(cm.te_nodes),)
    assert len(cm.te_nodes) > 20, "M6 must have many spanwise TE stations"
    z = mesh.nodes[cm.te_nodes, 2]
    assert z.min() >= -1e-9 and z.max() <= B_SEMI + 1e-9


@pytest.mark.parametrize("directory", [M1_DIR, M4_DIR], ids=["M1", "M4"])
def test_section_extractor_tracks_the_swept_planform(directory):
    """section_cp_curve_levelset derives the LOCAL chord/x_le from its own cut,
    so on the swept, tapered M6 both must follow the planform: chord shrinks and
    the LE moves aft toward the tip. Both surfaces must be populated."""
    mesh, cm, mvop = _setup(_require(directory))
    phi = np.zeros(mvop.n_total)
    phi[: mvop.n_main] = mesh.nodes[:, 0]             # phi = x (uniform flow)
    chords, x_les = [], []
    for eta in ETAS:
        c = section_cp_curve_levelset(mesh, mvop, phi, eta=eta, b_semi=B_SEMI,
                                      m_inf=M_INF)
        assert len(c["x_upper"]) >= 5 and len(c["x_lower"]) >= 5
        chords.append(c["chord"])
        x_les.append(c["x_le"])
    assert chords[0] > chords[1] > chords[2], chords     # taper
    assert x_les[0] < x_les[1] < x_les[2], x_les         # sweep


def test_section_extractor_raises_beyond_the_tip():
    """A cut plane outboard of the tip has no wing to cut -- the guard must
    raise rather than silently return an empty/garbage curve."""
    mesh, cm, mvop = _setup(_require(M1_DIR))
    phi = np.zeros(mvop.n_total)
    with pytest.raises(ValueError, match="too sparse"):
        section_cp_curve_levelset(mesh, mvop, phi, eta=1.05, b_semi=B_SEMI,
                                  m_inf=M_INF)


def test_d11_upper_surface_equals_the_main_based_section():
    """D11 in 3D: the UPPER surface reads the TE node's MAIN dof (= its own
    side), so the level-set upper section is BIT-IDENTICAL to the conforming
    extractor fed the main potential -- which is why the gate's shock metrics
    (all upper-surface) need no new machinery. The LOWER surface is where the
    per-side mapping bites: it must differ, at the TE triangles, by O(1) in Cp
    (reading `main` there is the junk that made cl_pressure = -3.35 in B3).
    """
    mesh, cm, mvop = _setup(_require(M1_DIR))
    # A physically-scaled multivalued state: freestream on the main dofs, plus a
    # constant jump on the aux copies (the discrete shape of a circulation). A
    # random field would NOT do -- its q^2 is so large that the isentropic Cp
    # saturates at the vacuum limit on every triangle, and both sides come out
    # at the same constant (measured), which would pass the test vacuously.
    phi = np.zeros(mvop.n_total)
    phi[: mvop.n_main] = mesh.nodes[:, 0]
    cut = np.flatnonzero(mvop.cm.ext_dof_of_node >= 0)
    phi[mvop.cm.ext_dof_of_node[cut]] = phi[cut] + 0.12
    for eta in ETAS:
        ls = section_cp_curve_levelset(mesh, mvop, phi, eta=eta, b_semi=B_SEMI)
        mn = section_cp_curve(mesh, mvop.main_potential(phi), eta=eta,
                              b_semi=B_SEMI)
        assert np.array_equal(ls["cp_upper"], mn["cp_upper"])
        assert ls["cp_lower"].shape == mn["cp_lower"].shape
        assert np.max(np.abs(ls["cp_lower"] - mn["cp_lower"])) > 1e-6


# ---------------------------------------------------------------------------
# Gated: the M6 3D solves (minutes each).
# ---------------------------------------------------------------------------

# Committed baselines: roadmap P5 ledger (conforming PICARD) and P8 ledger
# (conforming NEWTON -- the true discrete solution of that discretization).
#
# ★ The lift is gated against the NEWTON value, per the B6 user arbitration
# (2026-07-12): the conforming Picard UNDER-circulates (P8 measured its M6
# medium lift +7.9% below the Newton truth; the P4-erratum bias), so it was
# never a valid A/B target. B7 reproduces that inversion in 3D -- the level-set
# Picard sits within a few % of the Newton truth while the conforming Picard
# (P5) sits ~9% below it -- so gating on P5 would penalise the B path for being
# CLOSER to the truth. P5 remains the anchor for the SHOCK positions (P8's M6
# run left them essentially unchanged: 0.596/0.541/0.362 vs 0.596/0.570/0.425).
P5_CL_KJ = 0.24788
P5_SHOCKS = {0.44: 0.596, 0.65: 0.570, 0.90: 0.425}
P8_CL_KJ = 0.2692          # the Newton truth -- the lift reference


@run_gates
@pytest.mark.parametrize("directory", [M1_DIR, M4_DIR], ids=["M1", "M4"])
def test_b7_m6_coarse_subsonic_3d_machinery(directory):
    """The 3D-only machinery, isolated from the transonic fold (M0.5): the swept
    ruled sheet carries a circulation that decays monotonically to ~0 at the tip
    (the spanwise clip), and the pressure and circulation lifts agree (V6)."""
    mesh, cm, mvop = _setup(_require(directory))
    r = solve_multivalued_lifting(mvop, mesh, 0.5, alpha_deg=ALPHA,
                                  farfield="neumann", n_outer_max=60,
                                  tol_residual=1e-7)
    assert r["converged"], r["residual_norm"]
    z, g = _gamma_of_z(mesh, cm, r["phi_ext"], mvop)

    assert g.min() > -1e-3, "circulation must not go negative"
    assert abs(g[-1]) < 0.02, f"tip Gamma must vanish, got {g[-1]:.4f}"
    half = len(g) // 2
    assert int((np.diff(g[half:]) > 1e-3).sum()) == 0, "outer half must decay"

    s_ref = planform_area(mesh.nodes, mesh.boundary_faces["wall"])
    clkj = cl_kj_3d(g, z, s_ref=s_ref, b_semi=B_SEMI)
    cp = wall_cp_levelset(mesh, mvop, r["phi_ext"], m_inf=0.5)
    clp = cl_pressure_3d_levelset(mesh, cp["cp"], cp["area"], cp["n_out"],
                                  ALPHA, s_ref)
    v6 = abs(clp - clkj) / abs(clkj)
    assert v6 < 0.05, f"V6 consistency {v6*100:.2f}%"


@run_gates
def test_farfield_vortex_is_contraindicated_in_3d():
    """★ The B7 far-field decision, measured (demo cases/demo/b7_onera_m6).

    The B-path vortex is a SPAN-UNIFORM 2D point vortex whose branch cut is the
    ray y=0, x>0 at every z. On M6 that prescribes a Dirichlet jump no cut
    supports, and the outlet -- where the sheet leaves the domain -- goes
    spuriously near-sonic. Two independent causes, separated here:
      (a) the alpha-aimed sheet is not coplanar with the y=0 cut (B3's rule);
      (b) a scalar Gamma cannot match Gamma(z) -> 0 (P5's branch-ray artifact).
    Re-aiming the sheet coplanar removes (a) but not (b), so the artifact
    shrinks without vanishing. Neumann carries no vortex => neither exists.
    This is the evidence for `farfield="neumann"` being B7's recipe, and it
    fails loudly if someone flips the default back.
    """
    path = _require(M1_DIR)
    mesh = read_mesh(path)
    c = mesh.nodes[mesh.elements].mean(axis=1)
    outlet = (c[:, 0] > 6.0) & (np.abs(c[:, 1]) < 0.3) & (c[:, 2] < B_SEMI)

    m_out = {}
    for aim, aim_tag in ((ALPHA, "alpha_aimed"), (0.0, "coplanar")):
        _, _, mvop = _setup(path, alpha_deg=aim)
        for ff in ("neumann", "vortex"):
            r = solve_multivalued_lifting(mvop, mesh, 0.5, alpha_deg=ALPHA,
                                          farfield=ff, n_outer_max=60,
                                          tol_residual=1e-7)
            m = np.sqrt(mvop.element_mach2(r["phi_ext"], 0.5))
            m_out[(aim_tag, ff)] = float(m[outlet].max())

    neu = m_out[("alpha_aimed", "neumann")]
    vor = m_out[("alpha_aimed", "vortex")]
    cop = m_out[("coplanar", "vortex")]
    assert vor > neu + 0.2, f"(a) not reproduced: vortex {vor:.3f} vs neu {neu:.3f}"
    assert neu < cop < vor, (
        f"(b) not reproduced: neumann {neu:.3f} < coplanar-vortex {cop:.3f} "
        f"< alpha-aimed-vortex {vor:.3f} expected")
    # neumann itself must stay well subsonic at the outlet
    assert neu < 0.7, neu


@run_gates
@pytest.mark.parametrize("directory", [M1_DIR, M4_DIR], ids=["M1", "M4"])
def test_b7_m6_coarse_transonic_gate(directory):
    """B7 gate: M6 coarse M0.84 / alpha 3.06 on the level-set path, A/B against
    the committed P5 (Picard) / P8 (Newton) conforming baseline.

    Per-mesh locks, NOT a grid-convergence claim (the P8 fold-zone discipline).

    Convergence semantics = the recorded transonic Picard regime (P4/B6): the
    top Mach levels park on a bounded residual tail (~1e-6) rather than reaching
    tol_residual, so the assertion is on a BOUNDED, PHYSICAL field (0 limited /
    0 floored, physical M_max) plus the gate metrics -- not on `converged`. The
    LS Newton is the cure for the tail and is deferred on M6 (it uses a plain
    splu, and P8's N6 measured true-3D LU fill at ~100x the 2.5D cost -- it needs
    the lagged-LU treatment first). See roadmap B7.
    """
    mesh, cm, mvop = _setup(_require(directory))
    r = solve_multivalued_transonic(
        mvop, mesh, M_INF, alpha_deg=ALPHA, farfield="neumann",
        m_start=0.60, dm=0.04,
        n_outer_seed=120, n_outer_level=600, tol_residual=1e-7,
    )
    phi = r["phi_ext"]

    # bounded + physical field (NOT `converged` -- see the docstring)
    assert r["n_limited"] == 0 and r["n_floored"] == 0, (
        r["n_limited"], r["n_floored"])
    assert r["residual_norm"] < 1e-4, r["residual_norm"]   # tail, not divergence
    m_max = float(np.sqrt(r["mach2_max"]))
    assert 1.0 < m_max < 2.5, m_max          # P5 1.398 / P8 Newton 2.13

    # Gamma(z): the 3D-only machinery
    z, g = _gamma_of_z(mesh, cm, phi, mvop)
    assert g.min() > -1e-3
    assert abs(g[-1]) < 0.02, f"tip Gamma {g[-1]:.4f}"
    half = len(g) // 2
    assert int((np.diff(g[half:]) > 1e-3).sum()) == 0

    # lift against the conforming NEWTON truth (see the P8_CL_KJ note above)
    s_ref = planform_area(mesh.nodes, mesh.boundary_faces["wall"])
    clkj = cl_kj_3d(g, z, s_ref=s_ref, b_semi=B_SEMI)
    assert abs(clkj - P8_CL_KJ) / P8_CL_KJ < 0.10, (clkj, P8_CL_KJ)

    # upper-surface shocks: present, monotone, no expansion shock, near P5,
    # and migrating FORWARD toward the tip (the P5 G5.1 clause)
    x_shock = {}
    for eta in ETAS:
        cur = section_cp_curve_levelset(mesh, mvop, phi, eta=eta, b_semi=B_SEMI,
                                        m_inf=M_INF)
        s = shock_report(cur, M_INF)["upper"]
        assert s["has_shock"] and s["monotone"] and not s["expansion_shock"], (eta, s)
        assert abs(s["x_shock"] - P5_SHOCKS[eta]) <= 0.06, (eta, s["x_shock"])
        x_shock[eta] = s["x_shock"]
    assert x_shock[0.90] <= x_shock[0.65] + 0.05, x_shock
