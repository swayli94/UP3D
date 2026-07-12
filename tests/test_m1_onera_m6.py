"""
M1 gate checks for the ONERA M6 swept/tapered half-wing mesh family
(cases/meshes/onera_m6, generate_onera_m6.py; docs/roadmap.md Track M).

Mesh-side M1 gate items covered here:
  - solver preprocessor ingests the family; tags complete (wall /
    farfield / symmetry / wake); spherical far field at ~15 MAC
  - wake is one planar interior sheet in the chord plane y = 0, swept
    from the TE line to the far field, ending at the tip (wake-tip
    closure: tip edge is one open chain starting at the exact tip TE
    corner)
  - element quality within bounds (min dihedral / max aspect -- the same
    bounds generate_onera_m6.py enforces at generation time)
  - the P2 topology asserts pass through cut_wake with the M1-specific
    semantics: per-node TE stations on the swept TE, single-valued
    free-edge nodes on the tip edge (Gamma(tip) = 0), Kutta probes found
    on the unstructured wing surface
  - G2.1-style freestream preservation on the CUT coarse mesh (the
    "solver preprocessor ingests it" acceptance at V0 level)

Heavier per-mesh topology asserts on every wake-tagged mesh (hard rule 7)
run in test_p2_wake_cut.py::test_topology_asserts_all_wake_meshes, which
picks up the M6 family automatically.
"""

from pathlib import Path

import numpy as np
import pytest

from pyfp3d.kernels.residual import assemble_stiffness_matrix
from pyfp3d.mesh.metrics import (
    compute_aspect_ratios,
    compute_min_dihedral_angles,
)
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.constraints.wake import WakeConstraint
from pyfp3d.meshgen.wing3d import B_SEMI, MAC, x_te

LEVELS = ["coarse", "medium"]
R_FAR = 15.0 * MAC
# Mirror generate_onera_m6.QUALITY_BOUNDS (the case script is not an
# importable package module).
QUALITY_BOUNDS = {"min_dihedral_deg": 2.0, "max_aspect_ratio": 60.0}


@pytest.fixture(scope="module")
def meshes():
    """The M6 .msh files are large and deliberately NOT committed
    (cases/meshes/onera_m6/*.msh is gitignored); regenerate them with
    `python cases/meshes/onera_m6/generate_onera_m6.py` (~30 s)."""
    mesh_dir = Path(__file__).parent.parent / "cases" / "meshes" / "onera_m6"
    missing = [lv for lv in LEVELS if not (mesh_dir / f"{lv}.msh").exists()]
    if missing:
        pytest.skip(
            f"onera_m6 meshes not generated ({', '.join(missing)}); run "
            "cases/meshes/onera_m6/generate_onera_m6.py"
        )
    return {level: read_mesh(mesh_dir / f"{level}.msh") for level in LEVELS}


@pytest.fixture(scope="module")
def coarse_cut(meshes):
    return cut_wake(meshes["coarse"])


class TestM6MeshFamily:
    @pytest.mark.parametrize("level", LEVELS)
    def test_ingestion_tags(self, meshes, level):
        mesh = meshes[level]
        assert set(mesh.boundary_faces) == {"wall", "farfield", "symmetry", "wake"}

    @pytest.mark.parametrize("level", LEVELS)
    def test_domain_geometry(self, meshes, level):
        """Spherical far field ~15 MAC; symmetry plane at z = 0; wake sheet
        planar in the chord plane y = 0, from the TE to the far field."""
        mesh = meshes[level]
        far = np.unique(mesh.boundary_faces["farfield"])
        xc = 0.5 * x_te(B_SEMI)  # sphere center x used by the builder
        r = np.linalg.norm(mesh.nodes[far] - np.array([xc, 0.0, 0.0]), axis=1)
        assert np.all(np.abs(r - R_FAR) < 1e-3 * R_FAR)

        sym = np.unique(mesh.boundary_faces["symmetry"])
        assert np.abs(mesh.nodes[sym, 2]).max() < 1e-9

        wake = np.unique(mesh.boundary_faces["wake"])
        assert np.abs(mesh.nodes[wake, 1]).max() < 1e-9
        assert mesh.nodes[wake, 2].min() > -1e-9
        assert mesh.nodes[wake, 2].max() < B_SEMI + 1e-9
        # sheet reaches the far-field sphere downstream
        assert mesh.nodes[wake, 0].max() > xc + R_FAR - 1e-3

    @pytest.mark.parametrize("level", LEVELS)
    def test_wake_attaches_to_swept_te(self, meshes, level):
        """Wall-and-wake nodes lie on the swept TE line x = x_te(z)."""
        mesh = meshes[level]
        te = np.intersect1d(np.unique(mesh.boundary_faces["wake"]),
                            np.unique(mesh.boundary_faces["wall"]))
        assert len(te) > 10
        xz = mesh.nodes[te][:, [0, 2]]
        x_expected = np.array([x_te(z) for z in xz[:, 1]])
        np.testing.assert_allclose(xz[:, 0], x_expected, rtol=0, atol=1e-6)

    @pytest.mark.parametrize("level", LEVELS)
    def test_wake_tip_closure(self, meshes, level):
        """M1 gate: the sheet's tip edge is ONE open chain whose wall-side
        endpoint is the exact tip TE corner (no crack, no dangling piece)."""
        mesh = meshes[level]
        wake = np.asarray(mesh.boundary_faces["wake"], dtype=np.int64)
        edges = np.sort(np.concatenate([wake[:, [0, 1]], wake[:, [1, 2]],
                                        wake[:, [2, 0]]]), axis=1)
        uniq, counts = np.unique(edges, axis=0, return_counts=True)
        boundary_edges = uniq[counts == 1]
        at_tip = np.abs(mesh.nodes[:, 2] - B_SEMI) < 1e-6
        tip_edges = boundary_edges[at_tip[boundary_edges].all(axis=1)]
        assert len(tip_edges) > 0

        nodes, deg = np.unique(tip_edges, return_counts=True)
        assert deg.max() <= 2, "tip edge self-intersects"
        endpoints = nodes[deg == 1]
        assert len(endpoints) == 2, "tip edge is not one open chain"
        wall_nodes = np.unique(mesh.boundary_faces["wall"])
        corner = endpoints[np.isin(endpoints, wall_nodes)]
        assert len(corner) == 1
        assert abs(mesh.nodes[corner[0], 0] - x_te(B_SEMI)) < 1e-9

    @pytest.mark.parametrize("level", LEVELS)
    def test_element_quality_within_bounds(self, meshes, level):
        mesh = meshes[level]
        dihedral = compute_min_dihedral_angles(mesh.nodes, mesh.elements)
        aspect = compute_aspect_ratios(mesh.nodes, mesh.elements)
        assert dihedral.min() >= QUALITY_BOUNDS["min_dihedral_deg"]
        assert aspect.max() <= QUALITY_BOUNDS["max_aspect_ratio"]

    def test_refinement_is_monotone(self, meshes):
        n = {level: len(meshes[level].elements) for level in LEVELS}
        assert n["medium"] > 3 * n["coarse"]


class TestM6WakeCut:
    """M1-specific cut_wake semantics (topology asserts run inside)."""

    def test_stations_and_free_edge(self, meshes, coarse_cut):
        mesh = meshes["coarse"]
        _, wc = coarse_cut

        # Swept TE: every TE node is its own spanwise station.
        assert wc.n_stations == len(wc.te_nodes)
        assert wc.n_stations > 20

        # Tip edge nodes are single-valued; the tip TE corner is one of
        # them and is therefore NOT a Kutta station (Gamma(tip) = 0).
        assert len(wc.free_nodes) > 0
        assert np.abs(mesh.nodes[wc.free_nodes, 2] - B_SEMI).max() < 1e-6
        assert wc.station_z.max() < B_SEMI - 1e-6
        assert not np.isin(wc.free_nodes, wc.te_nodes).any()
        assert not np.isin(wc.free_nodes, wc.master_nodes).any()

        # Kutta probes exist on both sides for every station and sit off
        # the chord plane on the correct side.
        assert np.all(wc.kutta_upper >= 0) and np.all(wc.kutta_lower >= 0)
        assert np.all(mesh.nodes[wc.kutta_upper, 1] > 0)
        assert np.all(mesh.nodes[wc.kutta_lower, 1] < 0)

    def test_g21_freestream_on_cut_mesh(self, coarse_cut):
        """G2.1 analogue on the cut M6 coarse mesh: phi = x gives machine-
        zero reduced residual away from Dirichlet/wall rows."""
        mesh_cut, wc = coarse_cut
        A = assemble_stiffness_matrix(mesh_cut.nodes, mesh_cut.elements)
        con = WakeConstraint(A, wc)
        R = con.A_reduced @ mesh_cut.nodes[: con.n_reduced, 0]

        to_master = np.arange(len(mesh_cut.nodes))
        to_master[wc.slave_nodes] = wc.master_nodes
        check = np.ones(con.n_reduced, dtype=bool)
        for tag in ("wall", "farfield", "symmetry"):
            faces = mesh_cut.boundary_faces[tag]
            check[np.unique(to_master[np.unique(faces)])] = False
        # Free-edge nodes end the sheet inside the domain; their rows sum
        # fluxes from both sides like any interior node and must vanish too.
        assert np.max(np.abs(R[check])) < 1e-10


# ---------------------------------------------------------------------------
# Mesh-family SELF-SIMILARITY (P13/G13.3, 2026-07-13)
#
# The M6 family derives every length from h_wall -- except that the far field
# was CLAMPED (h_far = min(2.5, 120*h_wall)), and the clamp bit ONLY at
# `coarse` (120*0.030 = 3.6 > 2.5). So coarse->medium refined the far field by
# 1.39x while the wall refined by 2x: `coarse` was NOT on the same refinement
# ray, and ANY three-point Richardson over (coarse, medium, fine) was invalid
# -- including the one P9/G9.1 attempted. These tests lock the fix:
#   * the SHIPPED levels keep the clamp, so they stay bit-identical and the
#     P5 / P8-G8.2 / B7 / M1 regression locks are untouched;
#   * RICHARDSON_LADDER is self-similar by exactly 2 in EVERY length scale.
# ---------------------------------------------------------------------------
def _generator():
    import importlib.util
    from .conftest import REPO_ROOT
    p = REPO_ROOT / "cases" / "meshes" / "onera_m6" / "generate_onera_m6.py"
    spec = importlib.util.spec_from_file_location("_gen_m6", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_shipped_levels_keep_the_clamp_and_stay_bit_identical():
    """coarse/medium/fine are the dev/gate meshes and carry regression locks;
    the self-similarity fix must NOT move them."""
    g = _generator()
    for level, h_far in (("coarse", 2.5), ("medium", 1.8), ("fine", 0.9)):
        assert g.LEVELS[level]["h_far"] == pytest.approx(h_far, rel=1e-12), level


def test_richardson_ladder_is_self_similar():
    """EVERY length scale must refine by exactly 2 at EVERY step, or a
    three-point Richardson over the ladder is not a grid-convergence study."""
    g = _generator()
    ladder = [g.LEVELS[k] for k in g.RICHARDSON_LADDER]
    assert len(ladder) == 3
    for key in ("h_wall", "h_edge", "h_wake", "h_far"):
        for a, b in zip(ladder[:-1], ladder[1:]):
            assert a[key] / b[key] == pytest.approx(2.0, rel=1e-12), (
                f"{key} refines by {a[key] / b[key]:.3f}, not 2 -- the ladder "
                f"is not self-similar and no Richardson over it is valid")


def test_the_clamp_is_what_broke_the_old_family():
    """Regression-documents the defect itself: with the clamp ON, the coarse
    end is off the refinement ray (far field refines 1.39x, wall 2x)."""
    g = _generator()
    clamped = [g._level_params(h) for h in (0.030, 0.015)]
    assert clamped[0]["h_wall"] / clamped[1]["h_wall"] == pytest.approx(2.0)
    ratio = clamped[0]["h_far"] / clamped[1]["h_far"]
    assert ratio == pytest.approx(2.5 / 1.8, rel=1e-9)
    assert ratio < 1.5, "the clamp must still be the thing we are guarding against"
