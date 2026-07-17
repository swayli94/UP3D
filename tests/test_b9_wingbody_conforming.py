"""
Track B / B9 (re-spec 2026-07-17): the CONFORMING wing-body capability.

B9 runs BOTH wake models on the M2 wing-body geometry (M0.5, coarse+medium)
and compares them. The conforming half is a NEW capability -- until this
phase there was no conforming wing-body mesh and cut_wake ValueError'd on
the wake-free family. The mesh generator now grows an EMBEDDED wake variant
(pyfp3d/meshgen/wingbody.py::onera_m6_wingbody_mesh(embed_wake=True)); the
solver-side plumbing (cut_wake, the wake Gamma constraint, the P14 pressure
Kutta) is UNCHANGED -- the fuselage waterline duplicates under the same
boundary-edge rule as the wing-alone symmetry root edge.

These are the GB9.1/GB9.3 checks. The coarse solves stay UN-GATED (they are
minutes, like the wing-alone M6 coarse); the M0.5 Newton pair is behind
PYFP3D_TRANSONIC_GATES=1 and lives in the B9 demo, since it is the expensive
run. The mesh is gitignored, so everything here skips until:

    python cases/meshes/onera_m6_wingbody_conforming/generate_onera_m6_wingbody_conforming.py --levels coarse
"""

import os
from pathlib import Path

import numpy as np
import pytest

from pyfp3d.constraints.te_pressure import TEControlVolumes
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.meshgen.fuselage import FuselageParams, radius_at
from pyfp3d.meshgen.wing3d import B_SEMI, x_te
from pyfp3d.meshgen.wingbody import junction_z
from pyfp3d.post.surface import planform_area, wall_force_coefficients
from pyfp3d.solve.picard import solve_laplace_lifting

REPO_ROOT = Path(__file__).parent.parent
MESH_DIR = REPO_ROOT / "cases" / "meshes" / "onera_m6_wingbody_conforming"

ALPHA = 3.06
FUSELAGE = FuselageParams()
Z_JUNC = junction_z(FUSELAGE)

GATES = os.environ.get("PYFP3D_TRANSONIC_GATES", "0") == "1"


def _require(level: str = "coarse") -> Path:
    p = MESH_DIR / f"{level}.msh"
    if not p.exists():
        pytest.skip(f"{MESH_DIR.name}/{level}.msh not generated (gitignored); "
                    "see this module's header")
    return p


@pytest.fixture(scope="module")
def coarse_cut():
    path = _require("coarse")
    mesh = read_mesh(str(path))
    mc, wc = cut_wake(mesh)
    return mesh, mc, wc


# ---------------------------------------------------------------------------
# GB9.1 (topology) -- the wake sheet is stitched to the wing TE, rides the
# fuselage waterline, and terminates at the tip as the only free edge.
# ---------------------------------------------------------------------------

def test_group_set(coarse_cut):
    mesh, mc, wc = coarse_cut
    assert set(mesh.boundary_faces) == {"wall", "fuselage", "farfield",
                                        "symmetry", "wake"}


def test_free_nodes_only_at_the_tip(coarse_cut):
    """The crack detector: a sheet-body stitch failure turns waterline or
    aft-symmetry edges into interior free edges at z < B_SEMI. Every free
    node must sit at the tip edge (z ~ B_SEMI)."""
    mesh, mc, wc = coarse_cut
    if len(wc.free_nodes) == 0:
        pytest.fail("no free nodes at all -- even the tip edge should be free")
    zf = mesh.nodes[wc.free_nodes, 2]
    assert zf.min() > B_SEMI - 1e-6, (
        f"free (single-valued) sheet nodes down to z={zf.min():.4f} < B_SEMI "
        f"({B_SEMI}) -- the sheet is not stitched somewhere inboard"
    )


def test_junction_te_is_a_station(coarse_cut):
    mesh, mc, wc = coarse_cut
    st = np.sort(wc.station_z)
    assert abs(st[0] - Z_JUNC) < 0.05, (
        f"innermost Kutta station z={st[0]:.4f}, expected the junction "
        f"{Z_JUNC:.4f}"
    )
    # one station per TE node (3D swept TE)
    assert wc.n_stations == len(wc.te_nodes)


def test_te_nodes_are_wing_only(coarse_cut):
    """TE (Kutta) stations must be wing-only: none inboard of the junction,
    which is what keeps wall_tag='wall' from minting waterline stations."""
    mesh, mc, wc = coarse_cut
    te_z = mesh.nodes[wc.te_nodes, 2]
    assert te_z.min() > Z_JUNC - 1e-6


def test_waterline_nodes_all_duplicated(coarse_cut):
    """The wake nodes shared with the fuselage skin (the waterline) carry
    the jump onto the body, so they must all be duplicated (masters), by the
    same boundary-edge rule as the wing-alone symmetry root edge."""
    mesh, mc, wc = coarse_cut
    wake_nodes = np.unique(mesh.boundary_faces["wake"])
    fus_nodes = set(np.unique(mesh.boundary_faces["fuselage"]).tolist())
    waterline = np.array([n for n in wake_nodes.tolist() if n in fus_nodes],
                         dtype=np.int64)
    assert len(waterline) > 0, "no waterline (wake n fuselage) nodes"
    # on the revolution surface
    wl = mesh.nodes[waterline]
    R = np.array([radius_at(FUSELAGE, float(x)) for x in wl[:, 0]])
    assert np.abs(np.abs(wl[:, 2]) - R).max() < 0.01 * FUSELAGE.r_f + 3e-3
    master = set(wc.master_nodes.tolist())
    n_dup = sum(1 for n in waterline.tolist() if n in master)
    assert n_dup == len(waterline), \
        f"only {n_dup}/{len(waterline)} waterline nodes duplicated"


def test_inboard_strip_maps_to_junction_station(coarse_cut):
    """Wake nodes inboard of the junction (the below-symmetry strip, z<0)
    take the innermost station's Gamma (constant carry-across the body -- no
    shed vorticity between the symmetry plane and the junction)."""
    mesh, mc, wc = coarse_cut
    innermost = int(np.argmin(wc.station_z))
    # masters at z < junction all belong to the innermost station
    mz = mesh.nodes[wc.master_nodes, 2]
    strip = mz < Z_JUNC - 1e-6
    if strip.any():
        assert np.all(wc.node_station[strip] == innermost), (
            "an inboard-strip wake node maps to a station other than the "
            "junction"
        )


# ---------------------------------------------------------------------------
# GB9.1 (freestream + lifting sanity) -- the cut mesh solves.
# ---------------------------------------------------------------------------

def test_fixed_gamma_jump_reaches_the_waterline(coarse_cut):
    """A prescribed-Gamma Laplace solve must realize the jump slave-minus-
    master == Gamma EXACTLY, INCLUDING the waterline and aft-symmetry slaves
    (not just the wing TE) -- the FE-space check that the body genuinely
    carries the branch cut."""
    mesh, mc, wc = coarse_cut
    g0 = 0.3
    r = solve_laplace_lifting(mc, wc, alpha_deg=ALPHA, gamma_fixed=g0)
    phi = r["phi"]
    jump = phi[wc.slave_nodes] - phi[wc.master_nodes]
    # each master's station Gamma
    g_of = g0 * np.ones(len(wc.master_nodes))
    assert np.abs(jump - g_of).max() < 1e-9, (
        f"prescribed-Gamma jump off by {np.abs(jump - g_of).max():.2e} "
        "(some slave does not carry Gamma -- a broken duplication)"
    )


@pytest.mark.skipif(not GATES, reason="coarse lifting solve is ~minutes; "
                    "set PYFP3D_TRANSONIC_GATES=1")
def test_laplace_lifting_loads_the_junction(coarse_cut):
    """Kutta loop converges and Gamma > 0 at every station INCLUDING the
    junction (the B8 lift-loss analogue detector -- a mis-terminated sheet
    would unload the innermost stations)."""
    mesh, mc, wc = coarse_cut
    r = solve_laplace_lifting(mc, wc, alpha_deg=ALPHA, max_kutta_updates=20)
    assert r["kutta_converged"], "Kutta loop did not converge in 20 updates"
    g = np.asarray(r["gamma"])
    o = np.argsort(wc.station_z)
    assert np.all(g[o][:-1] > 0.0), "a non-tip station carries Gamma <= 0"
    # fuselage carries ~no lift vs the wing (GB9.4, coarse RECORDED band)
    s_ref = planform_area(mesh.nodes, mesh.boundary_faces["wall"])
    cl_w = wall_force_coefficients(
        mesh.nodes, mesh.elements, mesh.boundary_faces["wall"], r["phi"],
        alpha_deg=ALPHA, s_ref=s_ref, m_inf=0.0)["cl"]
    cl_f = wall_force_coefficients(
        mesh.nodes, mesh.elements, mesh.boundary_faces["fuselage"], r["phi"],
        alpha_deg=ALPHA, s_ref=s_ref, m_inf=0.0)["cl"]
    assert abs(cl_f) < 0.15 * abs(cl_w), (
        f"fuselage cl {cl_f:.4f} not small vs wing {cl_w:.4f}"
    )


# ---------------------------------------------------------------------------
# GB9.3 -- junction TE control volumes take only WING-side elements (the
# M2 open verification item, conforming side). TEControlVolumes builds
# wall-adjacency from boundary_faces["wall"] only, so a tet touching the
# fuselage enters a junction fan only through its WING face.
# ---------------------------------------------------------------------------

def test_te_control_volumes_construct_and_are_wing_side(coarse_cut):
    mesh, mc, wc = coarse_cut
    cvs = TEControlVolumes(mc, wc)          # raises on any empty/mis-sided fan
    stats = cvs.fan_stats()
    assert stats["fan_u_min"] >= 1 and stats["fan_l_min"] >= 1, (
        f"a TE control volume is empty: {stats}"
    )

    # the junction (innermost) TE node's fans, from the packed storage
    j = int(np.argmin(mesh.nodes[wc.te_nodes, 2]))
    up = cvs._u["elems"][cvs._u["off"][j]:cvs._u["off"][j + 1]]
    lo = cvs._l["elems"][cvs._l["off"][j]:cvs._l["off"][j + 1]]
    assert len(up) >= 1 and len(lo) >= 1

    # No junction-fan element is fuselage-only: TEControlVolumes builds
    # wall-adjacency from EXACT wing-'wall' face ownership, so every fan
    # element must own a wing wall face (a >=3-node subset in the wing face
    # set). This is the conforming half of GB9.3.
    wing_faces = {frozenset(f.tolist())
                  for f in np.asarray(mc.boundary_faces["wall"], np.int64)}
    el = np.asarray(mc.elements, np.int64)

    def owns_wing_face(e):
        n = el[e]
        return any(frozenset(n[list(c)].tolist()) in wing_faces
                   for c in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)))

    for e in np.concatenate([up, lo]):
        assert owns_wing_face(int(e)), (
            f"junction-fan element {e} owns no wing wall face -- a fuselage "
            "element polluted the TE control volume"
        )
