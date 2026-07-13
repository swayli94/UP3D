"""
M5 gate checks for the ROUNDED-TIP ONERA M6 family (cases/meshes/
onera_m6_roundtip, generate_onera_m6_roundtip.py; docs/roadmap.md Track M M5).

WHAT M5 IS FOR. P13/G13.3 localized the last thing blocking 3D grid convergence
to the FLAT tip cap of the M1 family: a flat cap meets the upper/lower surfaces
at a sharp convex edge, which in potential flow is a singularity, and its box
peak Mach DIVERGES under uniform refinement (p = +0.321) while the wake free
edge and the wing interior stay bounded. M5 replaces the cap with the half body
of revolution swept by the tip section about its own chord line.

The gate has to establish two separate things, and the tests below are grouped
that way:

  (1) THE CAP IS ROUND, DISCRETELY. Not "the CAD is round" -- the solver only
      ever sees the triangulation. The discriminator is the crease angle on the
      tip-section seam (the locus that IS the sharp edge in M1): it must decay
      like O(h) -- facets approximating a smooth surface -- rather than sitting
      at the edge's turning angle. M1 measures ~92 deg at every level, which is
      what a real edge does; M5 measures 47 -> 25 -> 12 deg. That contrast IS
      the gate, so both families are measured by the same code here.

  (2) NOTHING ELSE MOVED. The cap radius vanishes at the LE and the TE, so the
      TE line, the wake sheet, the tip TE corner, the Kutta stations and B_SEMI
      must all be exactly what they were -- otherwise the A/B against M1 (and
      every P13 conclusion drawn from it) is confounded.

The heavy per-mesh topology asserts on every wake-tagged mesh (hard rule 7) run
in test_p2_wake_cut.py::test_topology_asserts_all_wake_meshes, which picks this
family up automatically from cases/meshes/*/*.msh.
"""

from pathlib import Path

import numpy as np
import pytest

from pyfp3d.constraints.wake import WakeConstraint
from pyfp3d.kernels.residual import assemble_stiffness_matrix
from pyfp3d.mesh.metrics import compute_aspect_ratios, compute_min_dihedral_angles
from pyfp3d.mesh.reader import read_mesh
from pyfp3d.mesh.wake_cut import cut_wake
from pyfp3d.meshgen.wing3d import (
    B_SEMI, C_TIP, MAC, TIP_CAP_RADIUS, x_le, x_te,
)
from pyfp3d.post.surface import wall_crease_angles

LEVELS = ["coarse", "medium"]
R_FAR = 15.0 * MAC
QUALITY_BOUNDS = {"min_dihedral_deg": 2.0, "max_aspect_ratio": 60.0}

#: Mirrors generate_onera_m6_roundtip.SEAM_XI / SEAM_MAX_DEG (the case scripts
#: are not importable package modules).
SEAM_XI = (0.05, 0.95)
SEAM_MAX_DEG = {"coarse": 60.0, "medium": 35.0, "fine": 20.0}

#: What a genuinely sharp edge measures, for contrast: the M1 flat cap, at any
#: refinement level. The turning angle of a flat cap against the section surface
#: is a property of the GEOMETRY, so refinement resolves it and never shrinks it.
FLAT_CAP_CREASE_DEG = 90.0


def _mesh_dir(name):
    return Path(__file__).parent.parent / "cases" / "meshes" / name


@pytest.fixture(scope="module")
def meshes():
    """The .msh files are large and deliberately NOT committed (gitignored, the
    M1 policy); regenerate with
    `python cases/meshes/onera_m6_roundtip/generate_onera_m6_roundtip.py`
    (~1 min for coarse + medium)."""
    d = _mesh_dir("onera_m6_roundtip")
    missing = [lv for lv in LEVELS if not (d / f"{lv}.msh").exists()]
    if missing:
        pytest.skip(
            f"onera_m6_roundtip meshes not generated ({', '.join(missing)}); "
            "run cases/meshes/onera_m6_roundtip/generate_onera_m6_roundtip.py"
        )
    return {level: read_mesh(d / f"{level}.msh") for level in LEVELS}


@pytest.fixture(scope="module")
def coarse_cut(meshes):
    return cut_wake(meshes["coarse"])


def seam_crease(mesh) -> np.ndarray:
    """Crease angles on the tip-section seam at z = B_SEMI, away from the LE and
    the TE (sharp by design in BOTH families -- the TE carries the Kutta
    condition, so it is not the cap's business)."""
    ang, mid = wall_crease_angles(mesh.nodes, mesh.elements,
                                  mesh.boundary_faces["wall"])
    xi = (mid[:, 0] - x_le(B_SEMI)) / C_TIP
    on_seam = ((np.abs(mid[:, 2] - B_SEMI) < 1e-9)
               & (xi > SEAM_XI[0]) & (xi < SEAM_XI[1]))
    return ang[on_seam]


# ---------------------------------------------------------------------------
# (1) the cap is round -- discretely, which is the only way the solver sees it
# ---------------------------------------------------------------------------
class TestTheCapIsRound:
    @pytest.mark.parametrize("level", LEVELS)
    def test_seam_is_not_an_edge(self, meshes, level):
        seam = seam_crease(meshes[level])
        assert len(seam) > 50, "seam locus not found -- did the tip geometry move?"
        assert seam.max() <= SEAM_MAX_DEG[level], (
            f"{level}: seam crease {seam.max():.1f} deg -- a flat cap measures "
            f"~{FLAT_CAP_CREASE_DEG:.0f} deg, so this cap is not round (or "
            "h_tip stopped resolving it)"
        )

    def test_seam_crease_decays_like_h(self, meshes):
        """★ THE GATE. A smooth surface approximated by facets creases by
        O(h * curvature) and halves when h halves. A sharp EDGE creases by its
        own turning angle and does not move. This is the difference between
        removing the singularity and merely re-meshing it."""
        coarse, medium = (seam_crease(meshes[lv]).max() for lv in ("coarse", "medium"))
        ratio = coarse / medium
        assert 1.5 < ratio < 2.6, (
            f"seam crease went {coarse:.1f} -> {medium:.1f} deg (ratio "
            f"{ratio:.2f}); O(h) decay would halve it. A ratio near 1 means the "
            "cap still has an edge"
        )

    def test_the_flat_cap_is_the_thing_we_replaced(self):
        """Regression-documents the defect: on the SAME metric and the SAME
        code, M1's flat cap creases by ~90 deg and does NOT improve with
        refinement. Without this the M5 numbers have nothing to be better than.
        """
        d = _mesh_dir("onera_m6")
        levels = [lv for lv in LEVELS if (d / f"{lv}.msh").exists()]
        if not levels:
            pytest.skip("onera_m6 (flat) meshes not generated")
        for level in levels:
            seam = seam_crease(read_mesh(d / f"{level}.msh"))
            assert seam.max() > 80.0, (
                f"flat {level}: seam crease {seam.max():.1f} deg -- expected the "
                "~90 deg sharp cap edge P13/G13.3 blamed"
            )

    @pytest.mark.parametrize("level", LEVELS)
    def test_cap_reaches_its_apex(self, meshes, level):
        """The cap is a body of revolution of radius t(x), so the wall reaches
        z = B_SEMI + TIP_CAP_RADIUS at the max-thickness station -- and the mesh
        must actually resolve that, not clip it to a couple of facets."""
        wall = np.unique(meshes[level].boundary_faces["wall"])
        z_max = meshes[level].nodes[wall, 2].max()
        apex = B_SEMI + TIP_CAP_RADIUS
        assert z_max <= apex + 1e-9
        assert z_max > apex - 0.01 * TIP_CAP_RADIUS, (
            f"{level}: wall reaches z = {z_max:.6f}, apex is {apex:.6f}"
        )


# ---------------------------------------------------------------------------
# (2) nothing else moved -- or the A/B against M1 is confounded
# ---------------------------------------------------------------------------
class TestNothingElseMoved:
    @pytest.mark.parametrize("level", LEVELS)
    def test_tip_te_corner_is_unmoved(self, meshes, level):
        """The cap radius vanishes at the TE, so the corner the wake attaches to
        must be the exact same point it is in M1."""
        mesh = meshes[level]
        wall = np.unique(mesh.boundary_faces["wall"])
        corner = np.array([x_te(B_SEMI), 0.0, B_SEMI])
        assert np.linalg.norm(mesh.nodes[wall] - corner, axis=1).min() < 1e-9

    @pytest.mark.parametrize("level", LEVELS)
    def test_wake_sheet_is_untouched(self, meshes, level):
        """Chord plane, attached to the TE, ending at the tip: the cap must not
        have perturbed the sheet (it never reaches aft of the local TE)."""
        mesh = meshes[level]
        assert "wake" in mesh.boundary_faces
        wake = np.unique(mesh.boundary_faces["wake"])
        assert np.abs(mesh.nodes[wake, 1]).max() < 1e-9, "wake left the chord plane"
        assert mesh.nodes[wake, 2].max() <= B_SEMI + 1e-9, "wake past the tip"
        te = np.intersect1d(wake, np.unique(mesh.boundary_faces["wall"]))
        assert len(te) >= 2, "wake sheet does not attach to the wing TE"

    @pytest.mark.parametrize("level", LEVELS)
    def test_ingestion_tags_and_far_field(self, meshes, level):
        mesh = meshes[level]
        assert set(mesh.boundary_faces) == {"wall", "farfield", "symmetry", "wake"}
        sym = np.unique(mesh.boundary_faces["symmetry"])
        assert np.abs(mesh.nodes[sym, 2]).max() < 1e-9
        far = np.unique(mesh.boundary_faces["farfield"])
        xc = 0.5 * (x_le(0.0) + x_te(B_SEMI))
        r = np.linalg.norm(mesh.nodes[far] - np.array([xc, 0.0, 0.0]), axis=1)
        assert np.all(np.abs(r - R_FAR) < 1e-3 * R_FAR)

    @pytest.mark.parametrize("level", LEVELS)
    def test_element_quality_within_bounds(self, meshes, level):
        """The cap is a small-radius feature and it does cost some quality --
        the gate is that it stays inside the bounds the solver has always run
        on (the same bounds M1 enforces)."""
        mesh = meshes[level]
        assert (compute_min_dihedral_angles(mesh.nodes, mesh.elements).min()
                >= QUALITY_BOUNDS["min_dihedral_deg"])
        assert (compute_aspect_ratios(mesh.nodes, mesh.elements).max()
                <= QUALITY_BOUNDS["max_aspect_ratio"])

    def test_cut_wake_keeps_the_m1_semantics(self, coarse_cut, meshes):
        """Per-node TE stations on the swept TE; the tip edge stays single-valued
        (Gamma(tip) = 0 discretely); Kutta probes found on both sides."""
        mesh_cut, wc = coarse_cut
        assert len(wc.station_z) > 50
        assert wc.station_z.min() >= -1e-9
        assert wc.station_z.max() < B_SEMI, (
            "the tip TE corner must stay OUT of the Kutta stations (free edge)"
        )
        assert np.all(wc.kutta_upper >= 0) and np.all(wc.kutta_lower >= 0)
        assert np.all(mesh_cut.nodes[wc.kutta_upper, 1] > 0)
        assert np.all(mesh_cut.nodes[wc.kutta_lower, 1] < 0)

    def test_g21_freestream_on_cut_mesh(self, coarse_cut):
        """G2.1 analogue on the cut round-tip mesh: phi = x gives a machine-zero
        reduced residual away from the Dirichlet/wall rows. This is the "solver
        ingests it" acceptance -- and it exercises the cap, whose wall rows are
        excluded but whose neighbouring fluid rows are not."""
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
        assert np.max(np.abs(R[check])) < 1e-10


# ---------------------------------------------------------------------------
# the ladder (no gmsh needed -- pure parameter algebra)
# ---------------------------------------------------------------------------
def _generator():
    import importlib.util
    from .conftest import REPO_ROOT
    p = (REPO_ROOT / "cases" / "meshes" / "onera_m6_roundtip"
         / "generate_onera_m6_roundtip.py")
    spec = importlib.util.spec_from_file_location("_gen_m5", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_ladder_is_self_similar_including_the_cap():
    """M1b had to bolt a `coarse_ss` level on because h_far was clamped. This
    family is self-similar from the start -- and h_tip has to refine with
    everything else, or the cap would be better resolved at fine than at coarse
    and the ladder would be measuring the cap's discretization, not the flow's.
    """
    g = _generator()
    ladder = [g.LEVELS[k] for k in g.RICHARDSON_LADDER]
    assert len(ladder) == 3
    for key in ("h_wall", "h_edge", "h_wake", "h_far", "h_tip"):
        for a, b in zip(ladder[:-1], ladder[1:]):
            assert a[key] / b[key] == pytest.approx(2.0, rel=1e-12), (
                f"{key} refines by {a[key] / b[key]:.3f}, not 2"
            )


def test_h_far_is_not_clamped():
    """The M1b defect, explicitly guarded against in the new family."""
    g = _generator()
    for level, p in g.LEVELS.items():
        assert p["h_far"] == pytest.approx(120.0 * p["h_wall"], rel=1e-12), level


def test_the_cap_is_off_by_default():
    """`tip_cap` defaults to "flat": every existing mesh family, and every P5 /
    P8-G8.2 / B7 / M1 lock anchored to one, must be bit-identical."""
    import inspect

    from pyfp3d.meshgen.wing3d import onera_m6_wing_mesh

    sig = inspect.signature(onera_m6_wing_mesh)
    assert sig.parameters["tip_cap"].default == "flat"
    assert sig.parameters["h_tip"].default is None
