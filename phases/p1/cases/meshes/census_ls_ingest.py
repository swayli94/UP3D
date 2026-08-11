"""
Track M / M2 — level-set ingest census of the wing-body meshes (evidence
artifact for the f3c7989 claim; CLAUDE.md: "A claim without a committed
artifact is not evidence").

Builds the SAME level-set ingest B9 will use — WakeLevelSet on the analytic
TE polyline (wingbody.te_polyline, junction -> tip, chord plane) +
CutElementMap with wall_nodes — on the local gitignored coarse.msh /
medium.msh, and writes ls_ingest_census.csv with, per level:

  level, alpha_deg, n_cut_elems, n_te_nodes, te_z_min, te_z_max,
  n_te_inboard_of_junction, te_set_identical_wing_only_vs_full

Structural claims being evidenced (roadmap M2 / commit f3c7989):
  * TE nodes span z in [junction_z = r_f = 0.15, B_SEMI = 1.1963] exactly,
  * ZERO TE nodes inboard of the junction (B1's spanwise clip excludes the
    fuselage region with no new code),
  * the TE-node set is identical whether wall_nodes is wing-only
    ("wall") or wing+fuselage ("wall" + "fuselage").

Run (meshes must exist locally; see generate_onera_m6_wingbody.py):
    python cases/meshes/onera_m6_wingbody/census_ls_ingest.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT))

from pyfp3d.mesh.reader import read_mesh                       # noqa: E402
from pyfp3d.meshgen.fuselage import FuselageParams             # noqa: E402
from pyfp3d.meshgen.wingbody import junction_z, te_polyline    # noqa: E402
from pyfp3d.wake import CutElementMap, WakeLevelSet            # noqa: E402

MESH_DIR = Path(__file__).resolve().parent
LEVELS = ("coarse", "medium")
# alpha = 0: the chord-plane ruled sheet, the same convention as the B1 M6
# census tests (tests/test_b1_cut_elements.py::_m6_levelset default).
ALPHA_DEG = 0.0


def wingbody_levelset(alpha_deg: float = ALPHA_DEG) -> WakeLevelSet:
    a = np.radians(alpha_deg)
    te = te_polyline(FuselageParams())
    return WakeLevelSet(te, direction=(np.cos(a), np.sin(a), 0.0))


def census(mesh, alpha_deg: float = ALPHA_DEG) -> dict:
    wls = wingbody_levelset(alpha_deg)
    wall_wing = np.unique(mesh.boundary_faces["wall"])
    wall_full = np.unique(np.concatenate([
        mesh.boundary_faces["wall"].ravel(),
        mesh.boundary_faces["fuselage"].ravel(),
    ]))

    cm = CutElementMap(mesh.nodes, mesh.elements, wls, wall_nodes=wall_wing)
    cm_full = CutElementMap(mesh.nodes, mesh.elements, wls,
                            wall_nodes=wall_full)

    te_z = mesh.nodes[cm.te_nodes, 2]
    z_junc = junction_z(FuselageParams())
    # Inboard tolerance 1e-6: the junction TE node is an OCC vertex sitting
    # 1.5e-9 off exact z = r_f (commit f3c7989; same tolerance as
    # tests/test_m2_wingbody.py::test_wake_polyline_endpoints_are_exact_wall_nodes),
    # so it must count as ON the junction, not inboard of it.
    return {
        "alpha_deg": alpha_deg,
        "n_cut_elems": int(len(cm.cut_elems)),
        "n_te_nodes": int(len(cm.te_nodes)),
        "te_z_min": float(te_z.min()),
        "te_z_max": float(te_z.max()),
        "n_te_inboard_of_junction": int(np.sum(te_z < z_junc - 1e-6)),
        "te_set_identical_wing_only_vs_full":
            bool(np.array_equal(np.sort(cm.te_nodes),
                                np.sort(cm_full.te_nodes))),
    }


def main() -> int:
    rows = []
    for level in LEVELS:
        path = MESH_DIR / f"{level}.msh"
        if not path.exists():
            print(f"SKIP {level}: {path} not generated "
                  "(gitignored; run generate_onera_m6_wingbody.py)")
            continue
        mesh = read_mesh(str(path))
        row = {"level": level, **census(mesh)}
        rows.append(row)
        print(f"{level}: " + ", ".join(f"{k}={v}" for k, v in row.items()
                                       if k != "level"))

    if not rows:
        print("No meshes found; nothing written.")
        return 1

    out = MESH_DIR / "ls_ingest_census.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
