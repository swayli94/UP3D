"""
Track M / M2 — ONERA M6 wing-body (wing + simplified axisymmetric fuselage),
wake-free, for the level-set (Track B) solver path.

The mesh itself is GITIGNORED (M1/M5 policy: minutes to regenerate). Tests that
need it skip unless it has been generated locally:

    python cases/meshes/onera_m6_wingbody/generate_onera_m6_wingbody.py

Pure-geometry tests (FuselageParams / radius_at / te_polyline / junction_z)
have no mesh or gmsh dependency and always run.
"""

import math
from pathlib import Path

import numpy as np
import pytest

from pyfp3d.mesh.metrics import compute_aspect_ratios, compute_min_dihedral_angles
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.meshgen.fuselage import (
    FuselageParams,
    profile_points,
    radius_at,
)
from pyfp3d.meshgen.wing3d import B_SEMI, TIP_CAP_RADIUS, x_te
from pyfp3d.meshgen.wingbody import junction_z, te_polyline
from pyfp3d.post.surface import wall_crease_angles
from pyfp3d.wake import CutElementMap, WakeLevelSet

MESH_DIR = Path(__file__).resolve().parents[1] / "cases" / "meshes" / "onera_m6_wingbody"
COARSE = MESH_DIR / "coarse.msh"
MEDIUM = MESH_DIR / "medium.msh"

QUALITY_BOUNDS = {"min_dihedral_deg": 2.0, "max_aspect_ratio": 60.0}


def _load(path: Path):
    if not path.exists():
        pytest.skip(
            f"onera_m6_wingbody/{path.name} not generated "
            "(run cases/meshes/onera_m6_wingbody/generate_onera_m6_wingbody.py)"
        )
    return read_mesh(str(path))


@pytest.fixture(scope="module")
def coarse():
    return _load(COARSE)


@pytest.fixture(scope="module")
def medium():
    return _load(MEDIUM)


# --------------------------------------------------------------------------
# Fuselage geometry (no mesh, no gmsh)
# --------------------------------------------------------------------------

def test_fuselage_profile_closes_on_the_axis():
    """R = 0 at both tips, so the revolved solid pinches to a point there."""
    p = FuselageParams()
    assert radius_at(p, p.x_nose_tip) == 0.0
    assert radius_at(p, p.x_tail_tip) == 0.0
    xs, rs = profile_points(p, 80)
    assert rs[0] == 0.0 and rs[-1] == 0.0
    assert (rs[1:-1] > 0.0).all(), "interior profile must have positive radius"
    assert math.isclose(rs.max(), p.r_f, rel_tol=1e-12)


def test_fuselage_sections_match_the_piecewise_law():
    p = FuselageParams()
    # cylinder section is exactly r_f
    for x in (0.0, 0.5 * p.x_body_end, p.x_body_end):
        assert math.isclose(radius_at(p, x), p.r_f, rel_tol=1e-12)
    # nose is an ellipse: R(-l_nose/2) = r_f * sqrt(3)/2
    assert math.isclose(radius_at(p, -0.5 * p.l_nose),
                        p.r_f * math.sqrt(3.0) / 2.0, rel_tol=1e-12)
    # cone tapers linearly r_f -> r_tail
    mid = 0.5 * (p.x_body_end + p.x_tail_start)
    assert math.isclose(radius_at(p, mid), 0.5 * (p.r_f + p.r_tail), rel_tol=1e-12)
    # tail cap is a sphere of radius r_tail
    assert math.isclose(radius_at(p, p.x_tail_start), p.r_tail, rel_tol=1e-12)


def test_fuselage_params_reject_nonsense():
    with pytest.raises(ValueError):
        FuselageParams(r_f=-1.0)
    with pytest.raises(ValueError):
        FuselageParams(r_tail=0.5, r_f=0.15)   # afterbody must taper


def test_wing_root_chord_lies_in_the_constant_radius_section():
    """This is what makes junction_z = r_f EXACT: the M6 root chord
    (x in [0, C_ROOT]) never leaves the cylinder, where R(x) = r_f."""
    from pyfp3d.meshgen.wing3d import C_ROOT, x_le
    p = FuselageParams()
    assert x_le(0.0) >= 0.0
    assert x_te(0.0) <= p.x_body_end
    assert math.isclose(radius_at(p, x_le(0.0)), p.r_f, rel_tol=1e-12)
    assert math.isclose(radius_at(p, x_te(0.0)), p.r_f, rel_tol=1e-12)


def test_te_polyline_runs_from_the_junction_to_the_tip():
    """The level-set wake's inboard end is the junction (z = r_f), NOT the
    symmetry plane -- inboard of the junction there is no wing TE, only body."""
    p = FuselageParams()
    te = te_polyline(p)
    assert te.shape == (2, 3)
    assert math.isclose(junction_z(p), p.r_f, rel_tol=1e-12)
    # inboard end: on the chord plane, at the junction, on the wing TE line
    assert math.isclose(te[0, 2], p.r_f, rel_tol=1e-12)
    assert te[0, 1] == 0.0
    assert math.isclose(te[0, 0], x_te(p.r_f), rel_tol=1e-12)
    # outboard end: the tip TE corner (unmoved by the fuse -- the round cap
    # radius vanishes at the TE)
    assert math.isclose(te[1, 2], B_SEMI, rel_tol=1e-12)
    assert math.isclose(te[1, 0], x_te(B_SEMI), rel_tol=1e-12)


# --------------------------------------------------------------------------
# Mesh gates
# --------------------------------------------------------------------------

def test_mesh_is_wake_free(coarse):
    """The LS path never sees a wake tag -- that is the whole point of doing
    M2 on the Track B path (roadmap M2 note)."""
    assert "wake" not in coarse.boundary_faces
    assert set(coarse.boundary_faces) == {"wall", "fuselage", "farfield",
                                          "symmetry"}


def test_mesh_quality_within_bounds(coarse):
    dihedral = compute_min_dihedral_angles(coarse.nodes, coarse.elements)
    aspect = compute_aspect_ratios(coarse.nodes, coarse.elements)
    assert dihedral.min() >= QUALITY_BOUNDS["min_dihedral_deg"]
    assert aspect.max() <= QUALITY_BOUNDS["max_aspect_ratio"]


def test_fuselage_group_lies_on_the_surface_of_revolution(coarse):
    """The wing/fuselage split is geometric, so assert it actually holds:
    every 'fuselage' node satisfies sqrt(y^2+z^2) = R(x)."""
    p = FuselageParams()
    fus = np.unique(coarse.boundary_faces["fuselage"])
    rr = np.hypot(coarse.nodes[fus, 1], coarse.nodes[fus, 2])
    R = np.array([radius_at(p, float(x)) for x in coarse.nodes[fus, 0]])
    assert np.abs(rr - R).max() < 0.01 * p.r_f + 3e-3


def test_wing_group_is_not_on_the_fuselage(coarse):
    """The converse: the wing skin must NOT be mistaken for body. Most wall
    nodes sit clear of the revolution surface (they coincide only at the seam)."""
    p = FuselageParams()
    wall = np.unique(coarse.boundary_faces["wall"])
    rr = np.hypot(coarse.nodes[wall, 1], coarse.nodes[wall, 2])
    R = np.array([radius_at(p, float(x)) for x in coarse.nodes[wall, 0]])
    off = np.abs(rr - R) > 0.01
    assert off.mean() > 0.9, "wall group looks like it absorbed fuselage faces"


def test_wake_polyline_endpoints_are_exact_wall_nodes(coarse):
    """BOTH ends of the level-set TE polyline land on real wall nodes, exactly.

    Outboard: the tip TE corner (the round cap's radius vanishes there, so the
    M5 fuse never moved it -- measured 0.0). Inboard: the wing TE (y = 0) meets
    the fuselage skin (z = r_f in the constant-radius section) at an exact OCC
    vertex -- measured 1.5e-9 at coarse and medium alike. This is what makes
    the analytic te_polyline() a faithful description of the discrete TE.
    """
    p = FuselageParams()
    wall = coarse.nodes[np.unique(coarse.boundary_faces["wall"])]
    corner = np.array([x_te(B_SEMI), 0.0, B_SEMI])
    assert np.linalg.norm(wall - corner, axis=1).min() < 1e-9
    te_in = np.array([x_te(junction_z(p)), 0.0, junction_z(p)])
    assert np.linalg.norm(wall - te_in, axis=1).min() < 1e-6


def test_round_tip_cap_survives_the_fuse(coarse):
    """The fuselage fuse must not disturb the tip: the cap still reaches its
    apex, and the tip TE corner (where the cap radius vanishes) is unmoved."""
    wall = coarse.nodes[np.unique(coarse.boundary_faces["wall"])]
    z_max = wall[:, 2].max()
    apex = B_SEMI + TIP_CAP_RADIUS
    assert B_SEMI < z_max <= apex + 1e-9
    assert z_max > B_SEMI + 0.8 * TIP_CAP_RADIUS


def test_fuselage_skin_has_no_seam_crease(coarse):
    """The body is ONE splined surface of revolution, not four fused
    primitives, so its crease angles are O(h) faceting of a smooth surface --
    not a fixed seam angle. (The M5 tip-cap argument, applied to the body.)"""
    ang, _ = wall_crease_angles(coarse.nodes, coarse.elements,
                                coarse.boundary_faces["fuselage"])
    assert np.percentile(ang, 99) <= 25.0
    # a real seam would show as a ridge of near-constant large angle; a smooth
    # faceted surface has a small median
    assert np.median(ang) < 8.0


# --------------------------------------------------------------------------
# Level-set ingest census (locks the f3c7989 / roadmap-M2 claim as tests;
# committed evidence artifact = cases/meshes/onera_m6_wingbody/
# ls_ingest_census.csv, regenerated by census_ls_ingest.py in that directory)
# --------------------------------------------------------------------------

def _wingbody_levelset(alpha_deg: float = 0.0) -> WakeLevelSet:
    """The ingest B9 will use: WakeLevelSet on the analytic TE polyline
    (junction -> tip, chord plane), same alpha=0 convention as the B1 M6
    census tests."""
    a = np.radians(alpha_deg)
    return WakeLevelSet(te_polyline(FuselageParams()),
                        direction=(np.cos(a), np.sin(a), 0.0))


def _cut_map(mesh, alpha_deg: float = 0.0, *, include_fuselage: bool = False):
    wall = mesh.boundary_faces["wall"].ravel()
    if include_fuselage:
        wall = np.concatenate([wall, mesh.boundary_faces["fuselage"].ravel()])
    return CutElementMap(mesh.nodes, mesh.elements,
                         _wingbody_levelset(alpha_deg),
                         wall_nodes=np.unique(wall))


# junction TE node is an OCC vertex 1.5e-9 off exact z = r_f (measured, both
# levels), so "inboard of the junction" uses the same 1e-6 tolerance as
# test_wake_polyline_endpoints_are_exact_wall_nodes.
INBOARD_TOL = 1e-6


def test_ls_ingest_census_coarse_exact(coarse):
    """Locks the numbers quoted by commit f3c7989 / roadmap M2 (measured,
    alpha=0): 1,415 cut elements, 76 TE nodes."""
    cm = _cut_map(coarse)
    assert len(cm.cut_elems) == 1415
    assert len(cm.te_nodes) == 76


def test_ls_ingest_census_medium_counts(medium):
    """Medium-level census (measured 2026-07-14; NOT in the f3c7989 prose,
    which quoted the coarse level only)."""
    cm = _cut_map(medium)
    assert len(cm.cut_elems) == 29108
    assert len(cm.te_nodes) == 150


@pytest.mark.parametrize("level_fixture", ["coarse", "medium"])
def test_ls_te_nodes_span_junction_to_tip(level_fixture, request):
    """TE nodes span z in [junction_z = r_f, B_SEMI] = junction -> tip, with
    ZERO TE nodes inboard of the junction: B1's spanwise clip
    (0 <= q <= span_length) excludes the fuselage region with no new code."""
    mesh = request.getfixturevalue(level_fixture)
    cm = _cut_map(mesh)
    te_z = mesh.nodes[cm.te_nodes, 2]
    z_junc = junction_z(FuselageParams())
    assert abs(te_z.min() - z_junc) < INBOARD_TOL       # inboard end = junction
    assert math.isclose(te_z.max(), B_SEMI, rel_tol=1e-12)  # outboard end = tip
    assert int(np.sum(te_z < z_junc - INBOARD_TOL)) == 0


@pytest.mark.parametrize("level_fixture", ["coarse", "medium"])
def test_ls_te_nodes_immune_to_fuselage_wall_nodes(level_fixture, request):
    """TE-node identification is identical whether wall_nodes is wing-only
    or wing+fuselage: the fuselage does not contaminate the TE set."""
    mesh = request.getfixturevalue(level_fixture)
    cm_wing = _cut_map(mesh, include_fuselage=False)
    cm_full = _cut_map(mesh, include_fuselage=True)
    assert np.array_equal(np.sort(cm_wing.te_nodes),
                          np.sort(cm_full.te_nodes))


def test_ls_te_nodes_independent_of_wake_aim(coarse):
    """TE nodes live on the TE line itself (s ~ 0 AND d ~ 0 there for any
    aim), so re-aiming the sheet to incidence must not change the TE set
    (the cut census legitimately does change: 1,415 at alpha=0 vs 1,364 at
    alpha=3.06, measured)."""
    cm0 = _cut_map(coarse, alpha_deg=0.0)
    cm3 = _cut_map(coarse, alpha_deg=3.06)
    assert np.array_equal(np.sort(cm0.te_nodes), np.sort(cm3.te_nodes))


def test_symmetry_plane_and_farfield(coarse):
    """Symmetry at z = 0; far field on the sphere. NOTE the center is the
    geometric one, NOT the far-field nodes' centroid -- they are only the
    z >= 0 HALF sphere, so their centroid sits well off the axis."""
    from pyfp3d.meshgen.wing3d import MAC

    sym = coarse.nodes[np.unique(coarse.boundary_faces["symmetry"])]
    assert np.abs(sym[:, 2]).max() < 1e-6

    p = FuselageParams()
    r_far = 15.0 * MAC
    xc = 0.5 * (p.x_nose_tip + max(x_te(B_SEMI), p.x_tail_tip))
    far = coarse.nodes[np.unique(coarse.boundary_faces["farfield"])]
    r = np.linalg.norm(far - np.array([xc, 0.0, 0.0]), axis=1)
    assert np.abs(r - r_far).max() < 1e-3 * r_far
