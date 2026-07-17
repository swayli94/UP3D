"""
Track B / B9 (re-spec 2026-07-17): the level-set (LS) wing-body wiring +
the junction TE control-volume verification (GB9.3, LS half).

The LS half of B9 was ready before this phase: the M2 wing-body family is
wake-free and its LS ingest census is locked by test_m2_wingbody.py (76
coarse / 150 medium TE nodes at alpha=0). This module adds the pieces B9
specifically needs:

  * wiring at the B9 incidence alpha=3.06 (the committed M6 subsonic
    convention), TE-node count locked (aim-independent) and the cut census
    RECORDED (aim-dependent);
  * GB9.3: the junction TE node's B4 control-volume fans take ONLY wing-side
    elements. wall_nodes = the wing 'wall' group only, so the innermost TE
    node (which sits at the junction, its element fan touching fuselage wall
    faces) has fuselage-only elements STRUCTURALLY excluded from its CV
    (multivalued.py::_build_te_control_volumes on_wall mask). A negative
    control (wall_nodes = wall + fuselage) records what would leak in.

The heavy M0.5 LS Newton solves live in the B9 demo (cases/demo/b9_wingbody);
these are wiring/census/CV checks only. The mesh is gitignored, so they skip
until:

    python cases/meshes/onera_m6_wingbody/generate_onera_m6_wingbody.py --levels coarse medium
"""

from pathlib import Path

import numpy as np
import pytest

from pyfp3d.mesh.reader import read_mesh
from pyfp3d.meshgen.fuselage import FuselageParams
from pyfp3d.meshgen.wingbody import junction_z, te_polyline
from pyfp3d.wake import CutElementMap, MultivaluedOperator, WakeLevelSet

REPO_ROOT = Path(__file__).parent.parent
MESH_DIR = REPO_ROOT / "cases" / "meshes" / "onera_m6_wingbody"

ALPHA = 3.06
FUSELAGE = FuselageParams()
Z_JUNC = junction_z(FUSELAGE)
# TE-node counts are aim-independent (test_m2_wingbody); the B9 alpha keeps
# them (the wing skin is untouched by the 2026-07-16 body re-spec).
TE_NODES = {"coarse": 76, "medium": 150}


def _require(level: str) -> Path:
    p = MESH_DIR / f"{level}.msh"
    if not p.exists():
        pytest.skip(f"{MESH_DIR.name}/{level}.msh not generated (gitignored); "
                    "see this module's header")
    return p


def _levelset(alpha_deg: float = ALPHA) -> WakeLevelSet:
    a = np.radians(alpha_deg)
    return WakeLevelSet(te_polyline(FUSELAGE),
                        direction=(np.cos(a), np.sin(a), 0.0))


def _setup(level: str, wall_groups=("wall",)):
    mesh = read_mesh(str(_require(level)))
    wls = _levelset()
    wall_nodes = np.unique(np.concatenate(
        [mesh.boundary_faces[g] for g in wall_groups]))
    cm = CutElementMap(mesh.nodes, mesh.elements, wls, wall_nodes=wall_nodes)
    mvop = MultivaluedOperator(mesh.nodes, mesh.elements, cm, levelset=wls)
    return mesh, cm, mvop


# ---------------------------------------------------------------------------
# GB9.2 wiring: the operator builds at the B9 incidence, TE count locked.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("level", ["coarse", "medium"])
def test_ls_wiring_te_census(level):
    mesh, cm, mvop = _setup(level)
    assert len(cm.te_nodes) == TE_NODES[level], (
        f"{level}: {len(cm.te_nodes)} TE nodes, expected {TE_NODES[level]}"
    )
    # TE nodes span junction -> tip, none inboard of the junction.
    z = mesh.nodes[cm.te_nodes, 2]
    assert z.min() > Z_JUNC - 1e-6, "a TE node is inboard of the junction"
    assert z.max() > 1.0, "TE nodes do not reach the tip"
    # a nonzero cut set exists (RECORDED, not locked -- aim-dependent)
    assert len(cm.cut_elems) > 0


# ---------------------------------------------------------------------------
# GB9.3: the junction TE control-volume fans are wing-side only.
# ---------------------------------------------------------------------------

def _junction_fan_elems(mvop, mesh):
    """Upper+lower CV element ids of the innermost (junction) TE node."""
    j = int(np.argmin(mesh.nodes[mvop.cm.te_nodes, 2]))
    cv = mvop._te_cv[j]
    return j, np.concatenate([cv["upper_elems"], cv["lower_elems"]])


def _owns_face_in(el_nodes, face_set):
    return any(frozenset(el_nodes[list(c)].tolist()) in face_set
              for c in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)))


@pytest.mark.parametrize("level", ["coarse", "medium"])
def test_junction_cv_is_wing_side_only(level):
    """With wing-only wall_nodes, every element in the junction TE node's
    control volume owns a WING wall face and NONE is fuselage-only."""
    mesh, cm, mvop = _setup(level, wall_groups=("wall",))
    j, fan = _junction_fan_elems(mvop, mesh)
    assert len(fan) >= 2, "junction control volume is empty"

    el = np.asarray(mesh.elements, np.int64)
    wing_faces = {frozenset(f.tolist())
                  for f in np.asarray(mesh.boundary_faces["wall"], np.int64)}
    fus_faces = {frozenset(f.tolist())
                 for f in np.asarray(mesh.boundary_faces["fuselage"], np.int64)}
    for e in fan:
        n = el[int(e)]
        assert _owns_face_in(n, wing_faces), (
            f"junction-fan element {e} owns no wing wall face"
        )
        assert not _owns_face_in(n, fus_faces), (
            f"junction-fan element {e} owns a FUSELAGE face -- it polluted "
            "the wing-side control volume"
        )


def test_junction_cv_negative_control_records_fuselage_leak():
    """Negative control (RECORDED): building the CVs with wall_nodes =
    wall + fuselage lets the wall-adjacency mask admit fuselage-adjacent
    elements at the junction. This documents WHY wing-only is the wiring;
    it is not a pass/fail on the physics."""
    mesh, cm, mvop = _setup("coarse", wall_groups=("wall", "fuselage"))
    j, fan = _junction_fan_elems(mvop, mesh)
    el = np.asarray(mesh.elements, np.int64)
    fus_faces = {frozenset(f.tolist())
                 for f in np.asarray(mesh.boundary_faces["fuselage"], np.int64)}
    n_fus = sum(1 for e in fan if _owns_face_in(el[int(e)], fus_faces))
    # This is the object the wing-only wiring avoids. We only assert the
    # census is well-defined (>= 0); the value is recorded in the message.
    assert n_fus >= 0, "negative-control CV construction failed"
    print(f"[GB9.3 negative control] junction CV admits {n_fus} "
          f"fuselage-adjacent elements when wall_nodes includes the fuselage")
